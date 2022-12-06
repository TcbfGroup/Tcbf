import os.path
from tempfile import NamedTemporaryFile
from tcbf.run_command import run_command
from pandas import read_table
import pyarrow.feather as feather
aligner_paramter = {

    "minimap2": {"near": "-cx asm5",
                 "medium": "-cx asm10",
                 "far": "-cx asm20"},
}


def mash_distance(genome1, genome2):
    if not os.path.exists(genome1 + ".msh"):
        run_command(f"mash sketch {genome1}")
    if not os.path.exists(genome2 + ".msh"):
        run_command(f"mash sketch {genome1}")
    distance = float(run_command(f"mash dist {genome1}.msh {genome2}.msh").split()[2])
    return distance


def minimap2_align(bound_query, target, output_file,threads, parameter=None):
    def process_minimap_result(paf_file, out,distance):
        col_names = "seq_id length start end strand seq_id2 length2 start2 end2 map_match map_length map_quality".split()

        data = read_table(paf_file, header=None, names=col_names, usecols=range(12))
        if distance <= 0.02:
            map_length = 3000
        elif distance <= 0.08:
            map_length = 1000
        elif distance <= 0.15:
            map_length = 200
        else:
            map_length = 100
        data = data.query(f"map_length >= {map_length}")
        need = "seq_id start end seq_id2 start2 end2".split()

        feather.write_feather(data.loc[:, need].sort_values(["seq_id", "seq_id2", "start2", "end2"]), out)
    genome1 = bound_query.replace("bound.fasta","genome.fa")

    if not parameter:
        distance = mash_distance(genome1, target)
        if distance <= 0.02:
            parameter = " -x asm5 "
        elif distance <= 0.08:
            parameter = " -x asm10 "
        elif distance <= 0.15:
            parameter = " -x asm20 "
        else:
            parameter = ""
    with NamedTemporaryFile("w+t") as paf:
        command = f"minimap2  {parameter} -t {threads}  {target} {bound_query}  > {paf.name}"
        run_command(command)
        process_minimap_result(paf.name, output_file,distance)