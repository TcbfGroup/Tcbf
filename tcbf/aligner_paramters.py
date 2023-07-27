import os.path
from tempfile import NamedTemporaryFile

import pysam

from tcbf.run_command import run_command
from pandas import read_table,concat



def mash_distance(genome1, genome2):
    if not os.path.exists(genome1 + ".msh"):
        run_command(f"mash sketch {genome1}")
    if not os.path.exists(genome2 + ".msh"):
        run_command(f"mash sketch {genome1}")
    distance = float(run_command(f"mash dist {genome1}.msh {genome2}.msh").split()[2])
    return distance


def minimap2_align(workdir,bound_query, target, output_file, query,threads,map_length = 50, parameter=None):
    def process_minimap_result(paf_file, out,map_length ):
        col_names = "seq_id length start end strand seq_id2 " \
                    "length2 start2 end2 map_match map_length map_quality".split()

        data = read_table(paf_file, header=None, names=col_names, usecols=range(12))
        data = data.query(f"map_length >= {map_length}")
        need = "seq_id start end seq_id2 start2 end2".split()

        data.loc[:, need].sort_values(["seq_id", "seq_id2", "start2", "end2"]).to_csv(out,index=False)

    genome1 = os.path.join(workdir,"Step1",f"{query}.genome.fa")
    distance = mash_distance(genome1, target)
    if not parameter:

        if distance <= 0.001:
            parameter = " -x asm5 "
        elif distance <= 0.1:
            parameter = " -x asm10 "
        elif distance <= 0.2:
            parameter = " -x asm20 "
        elif distance <= 0.5:
            parameter = ""
        else:
            parameter = "-x sr"
        parameter += " -p 0.3 "
    with NamedTemporaryFile("w+t") as paf:
        command = f"minimap2  {parameter} -t {threads}  {target} {bound_query}  > {paf.name}"
        run_command(command)
        process_minimap_result(paf.name, output_file,map_length)



def parallel_lastz(query, target, parameters, threads):
    def run_lastz(seq_id):
        from tcbf.extract_TAD_boundary import format_seq
        seq = sequences[seq_id]
        with NamedTemporaryFile("w+t") as result_tmp:
            with NamedTemporaryFile("w+t") as tmp_file:
                tmp_file.write(f">{seq_id}\n{format_seq(seq)}")
                command = f"lastz {tmp_file.name} {query} {parameters} > {result_tmp.name} "
                run_command(command)
                result = read_table(result_tmp.name)
                return result

    from concurrent.futures import ThreadPoolExecutor
    target_file = pysam.FastaFile(target)
    sequences = {i:target_file[i] for i in target_file.references}
    with ThreadPoolExecutor(max_workers=threads)as executor:
        results = list(executor.map(run_lastz, sequences.keys()))
        results = concat(results)
    return results

def lastz_align(workdir,bound_query, target, output_file, query,threads,map_length = 50, parameter=None):

    if not parameter:
#         parameter = f"E=30 H=3000 K=5000 L=5000 M=10 O=400 T=1 Q={os.path.join(workdir,'Step2','general.q')} --notransition --step=20"
#     parameter += " --ambiguous=iupac    --format=general:name1,start1,end1,name2,start2,end2  "
#     general_q = """A C G T
# 91 -114 -31 -123
# -114 100 -125 -31
# -31 -125 100 -114
# -123 -31 -114 91"""

        # with open(os.path.join(workdir,"Step2","general.q"),"w")as f:
        #     f.write(general_q)
        parameter = " --notransition --step=20  --format=general:name1,start1,end1,name2,start2,end2  --ambiguous=iupac  --nogapped "
    results = parallel_lastz(bound_query,target,parameter,threads)
    results.index =  "seq_id start end seq_id2 start2 end2".split()

    results.to_csv(output_file,index=False)

def last_align(workdir,bound_query, target, output_file,threads,map_length = 50):
    names = "queryid subjectid  identity alignment_length mismatches gapopens " \
            "qstart qend sstart send evalue bitscore querylength subjectlength rawscore".split()
    target_species = os.path.basename(target).split("_")[0]
    db_dir = os.path.join(workdir, "Step1", f"{target_species}_lastdb")
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
    def run_last(seq_id):
        from tcbf.extract_TAD_boundary import format_seq

        with NamedTemporaryFile("w+t") as result_tmp:
            with NamedTemporaryFile("w+t") as fasta_tmp_file:
                if not os.path.exists(os.path.join(db_dir,f"{seq_id}.prj")):
                    seq = sequences[seq_id]
                    fasta_tmp_file.write(f">{seq_id}\n{format_seq(seq)}")
                    command1 = f"lastdb -uMAM4 {os.path.join(db_dir,seq_id)} {fasta_tmp_file.name};"
                    run_command(command1)


                command = f"lastal -f BlastTab+ {os.path.join(db_dir,seq_id)} {bound_query}  > {result_tmp.name} "
                run_command(command)

                result = read_table(result_tmp.name,comment="#", names=names).query(f"alignment_length >= {map_length}").loc[:,["queryid","qstart","qend","subjectid","sstart","send"]]
        return result
    from concurrent.futures import ThreadPoolExecutor
    target_file = pysam.FastaFile(target)
    sequences = {i:target_file[i] for i in target_file.references}
    with ThreadPoolExecutor(max_workers=threads)as executor:
        results = concat(list(executor.map(run_last, sequences)))
    results.index =  "seq_id start end seq_id2 start2 end2".split()
    results.to_csv(output_file,index=False)

def blat_align(workdir,bound_query, target, output_file, query,threads,map_length = 50, parameter=None):
    genome1 = os.path.join(workdir,"Step1",f"{query}.genome.fa")
    distance = mash_distance(genome1, target)
    if not parameter:

        if distance <= 0.1:
            parameter = " -tileSize=11 -minScore=100 -minIdentity=98"
        elif distance <= 0.2:
            parameter = " -tileSize=11 -stepSize=11 -oneOff=0 -minMatch=2 -minScore=30 -minIdentity=90 -maxGap=2 -maxIntron=75000 "
        else:
            parameter = " -tileSize=12 -oneOff=1 -minMatch=1 -minScore=30 -minIdentity=80 -maxGap=3 -maxIntron=75000 "
        parameter += " -t=dna -q=dna -fastMap -noHead  "

