import os.path
import sys
from pandas import read_table
from tempfile import NamedTemporaryFile
from pandas import read_table,read_csv,concat
from tcbf.run_command import run_command
def parse_block(workdir,query,target,out):
    query_boundries = read_csv(os.path.join(workdir,"Step1",f"{query}.boundary.bed"))
    query_gene = read_table(os.path.join(workdir,"Step2",f"{query}.bed")).iloc[:,:4]
    query_gene.columns = "chromosome start end gene_name".split()

    target_gene = read_table(os.path.join(workdir,"Step2",f"{target}.bed")).iloc[:,:4]
    target_gene.columns = "chromosome start end target".split()
    def get_gene(info):
        chrom, start, end = info[1:]
        f = query_gene.query(f"chromosome == '{chrom}' & end >= {start} & start <= {end}")
        return f
    result = []
    for value in query_boundries.itertuples():
        res = get_gene(value[1:])
        if res.shape[0] > 0:
            res["tad_name"] = value[1]
            result.append(res)

    query_boundary_gene = concat(result)

    block = read_table(os.path.join(workdir,"Step2",f"{query}.{target}.lifted.anchors"),header = None,comment = "#")
    block.columns = "gene_name target score".split()
    tmp = query_boundary_gene.merge(block).loc[:,["tad_name","target"]]
    final = tmp.merge(target_gene).loc[:,["tad_name","chromosome","start","end"]]
    final.columns = "seq_id seq_id2 start2 end2".split()
    final.to_csv(out,index=False)
    # return final


def extract_pep(workdir,species):
    genome = os.path.join(workdir,"Step1",f"{species}.genome.fa")
    gff3 = os.path.join(workdir,"Step1",f"{species}.gff3")
    out = os.path.join(workdir,"Step2",f"{species}.pep")
    command = f"gffread -g {genome} -y {out} {gff3};" \
              f"cd {os.path.join(workdir,'Step2')};" \
              f" lastdb -p {species} {species}.pep"
    run_command(command)
def extract_mRNA_bed(workdir,species):

    if not os.path.exists(os.path.join(workdir,"Step2")):
        os.mkdir(os.path.join(workdir,"Step2"))
    gff_file= os.path.join(workdir,"Step1",f"{species}.gff3")
    out =  os.path.join(workdir,"Step2",f"{species}.bed")
    def parse_key(fileds):
        if "ID" in fileds:
            return "ID"
        elif "Name" in fileds:
            return "Name"
        elif "transcript_id" in fileds:
            return "transcript_id"
        else:
            sys.exit("gff file mRNA field error")
    with open(gff_file)as f:
        for line in f:
            if line.startswith("#") or not line:
                continue
            info = line.split()
            if info[2] == "mRNA":
                type = "mRNA"
                fields = parse_key(info[8])
                break
            elif info[2] == "transcript":
                type = "transcript"
                fields = parse_key(info[8])
                break
    try:
        if type == "mRNA":
            command = f"python -m jcvi.formats.gff bed --primary_only --type=mRNA  --key={fields} {gff_file} -o {out}"
        elif type == "transcript":
            command = f"python -m jcvi.formats.gff bed  --type=transcript --key={fields} {gff_file} -o {out} "
        else:
            raise Exception("missing valid mRNA feature ")
    except:
        sys.exit("missing valid mRNA feature")
    run_command(command)
    pep_file =  os.path.join(workdir,"Step2",f"{species}.pep")
    peps_id = {i.split()[0][1:] for i in open(pep_file) if i.startswith(">")}
    bed = read_table(out,header=None)
    bed = bed[bed[3].isin(peps_id)]
    bed.to_csv(out,header=None,index=False,sep="\t")

def get_pep_pair(workdir,query,target):
    step2 = os.path.join(workdir,"Step2")
    command = f"cd {step2};" \
              f"python -m jcvi.compara.catalog ortholog --dbtype=prot {query} {target} --no_strip_names --no_dotplot"
    run_command(command)


# def run_pep(workdir):
#     from tcbf.network_construct import get_species
#     species = get_species(workdir)
#     for i in species:
#         extract_mRNA_bed(workdir,i)
#         extract_pep(workdir,i)
#     from itertools import combinations
#
#     for a,b in combinations(species,2):
#         get_pep_pair(workdir,a,b)
#         with NamedTemporaryFile("w+t") as Collinearity:
#             parse_block(workdir, a, b, Collinearity.name)
#             bound_bed = os.path.join(workdir,"Step1",f"{b}.bound.bed")
#             networkout = os.path.join(workdir,"Step2",f"{a}-{b}.gene.network")
#             command = f"tcbf_syn_process  {Collinearity.name}  {bound_bed} {networkout} 40000"
#             run_command(command)


def check_pep_bed(workdir,species):
    if not os.path.exists(os.path.join(workdir,"Step2")):
        os.mkdir(os.path.join(workdir,"Step2"))
    file2 = os.path.join(workdir,"Step2",f"{species}.pep")
    if not os.path.exists(file2):
        extract_pep(workdir,species)
    file1 = os.path.join(workdir,"Step2",f"{species}.bed")
    if not os.path.exists(file1):
        extract_mRNA_bed(workdir,species)



def align_gene(workdir,query,target,maxgap):

    get_pep_pair(workdir, query, target)

    with NamedTemporaryFile("w+t") as Collinearity:
        parse_block(workdir, query, target, Collinearity.name)
        boundary = os.path.join(workdir, "Step1", f"{target}.boundary.bed")
        networkout = os.path.join(workdir, "Step2", f"{query}-{target}.gene.pair")
        command = f"tcbf_syn_process  {Collinearity.name}  {boundary} {networkout} {maxgap}"
        run_command(command)