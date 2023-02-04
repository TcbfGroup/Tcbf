import os.path
from tempfile import NamedTemporaryFile
from tcbf.run_command import run_command
from pandas import read_table



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
