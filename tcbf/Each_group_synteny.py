import os.path
from tcbf.extract_TAD_bound import format_seq

from tempfile import NamedTemporaryFile
from tcbf.run_command import run_command
import pysam
def pysam_get_sequence(genome, name):
    fa = pysam.FastaFile(genome)


    sequence = fa.fetch(name)

    return sequence

def plot_synteny(workdir, TAD_names: list, output_file):
    workdir = os.path.abspath(workdir)
    TAD_file_dir = os.path.join(workdir, "Step1")
    seqs = []

    for line in TAD_names:


        genome = line.split("_")[0]
        file = os.path.join(TAD_file_dir, f"{genome}.bound.fasta")

        seq = pysam_get_sequence(file,line)
        seqs.append((line,seq))
    with NamedTemporaryFile("wt") as temp_paf:
        with NamedTemporaryFile("wt") as temp_file:
            for line in seqs:
                temp_file.write(f">{line[0]}\n{format_seq(line[1])}\n")
            command = f"minimap2 -X -N 50  -c  -p 0.1 {temp_file.name}  {temp_file.name}  > {temp_paf.name}"
            run_command(command)
            print(command)
            command2 = f"plot_TAD_bound_synteny {temp_paf.name} {output_file}"
            run_command(command2)
