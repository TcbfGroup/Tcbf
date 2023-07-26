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



def parallel_lastz(query,target,paramters,threads):
    def run_lastz(seq_id):
        from tcbf.extract_TAD_boundary import format_seq
        seq = sequences[seq_id]
        with NamedTemporaryFile("w+t") as result_tmp:
            with NamedTemporaryFile("w+t") as tmp_file:
                tmp_file.write(f">{seq_id}\n{format_seq(seq)}")
                command = f"lastz {tmp_file.name} {query} {paramters} > {result_tmp} "
                run_command(command)
                result = read_table(result_tmp.name)
                return result

    from concurrent.futures import ProcessPoolExecutor
    target_file = pysam.FastaFile(target)
    sequences = {i:target_file for i in target_file.references}
    with ProcessPoolExecutor(max_workers=threads)as executor:
        results = [result for result in executor.map(run_lastz, sequences)]
        results = concat([read_table(i) for i in results])
    return results






def lastz_align(workdir,bound_query, target, output_file, query,threads,map_length = 50, parameter=None):

    if not parameter:
        parameter = "E=30 H=3000 K=5000 L=5000 M=10 O=400 T=1 Q=general.q --notransition "
    parameter += " --ambiguous=iupac    --format=general:name1,start1,end1,name2,start2,end2  "
    result = parallel_lastz(bound_query,target,parameter,threads)
    result.to_csv(output_file,index=False)


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

