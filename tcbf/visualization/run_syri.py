import os.path
import unittest
import pysam
from tcbf.run_command import run_command
from tempfile import NamedTemporaryFile
from tcbf.aligner_paramters import mash_distance
from tcbf.extract_TAD_bound import format_seq,bed2position
from tempfile import TemporaryDirectory
from pandas import read_table


def reverse_complement(seq):
    from Bio.Seq import Seq
    s = str(Seq(seq).reverse_complement())
    return s

def get_sequence(genome,chrom,start,end,orientation =1):
    start,end = int(start),int(end)
    position = bed2position(chrom, start, end)



    fasta = pysam.Fastafile(genome)

    seq = fasta.fetch(region=position)
    if orientation == -1:
        seq = reverse_complement(seq)
    return format_seq(seq.strip("N"))

def parse_syri_out(outfile):


    VARS = ['SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
    data = read_table(outfile,header=None)
    data = data[data[10].isin(VARS)].loc[:,[0, 1, 2, 5, 6, 7, 10]]
    data[[0, 5, 10]] = data[[0, 5, 10]].astype(str)
    data[[1, 2, 6, 7]] = data[[1, 2, 6, 7]].astype(int)
    data.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend', 'type']
    return data





def SV_identify(name1,position1,
                name2,position2,
                workdir,
                ):
    genome1 = os.path.join(workdir,"Step1",f"{name1}.genome.fa")
    genome2 = os.path.join(workdir, "Step1", f"{name2}.genome.fa")
    distance = mash_distance(genome1, genome2)
    print(distance)
    if distance <= 0.02:
        parameter = " -ax asm5 "
    elif distance <= 0.15:
        parameter = " -ax asm10"

    elif distance <= 0.8:
        parameter = " -ax asm20 "
    else:
        parameter = "-a"

    name1 = os.path.basename(genome1).split(".",1)[0]
    name2 = os.path.basename(genome2).split(".",1)[0]
    with NamedTemporaryFile("w+t") as f1:
        with NamedTemporaryFile("w+t") as f2:
            with NamedTemporaryFile("w+t")as sam:
                with TemporaryDirectory() as tmpdirname:

                    f1.write(f">{name1}\n{get_sequence(genome1,*position1)}\n")
                    f2.write(f">{name2}\n{get_sequence(genome2,*position2)}\n")

                    command = f"minimap2  --eqx {parameter}  {f1.name} {f2.name}  > {sam.name}"
                    run_command(command)

                    command2 = f"syri -c {sam.name} -r {f1.name} -q {f2.name} -k -F S   --dir  {tmpdirname} --nosnp "
                    run_command(command2)
                    outfile = os.path.join(tmpdirname,"syri.out")
                    if os.path.exists(outfile):
                        return parse_syri_out(outfile)




class TestSyri(unittest.TestCase):
    """
    测试syri function
    """
    def test_syri(self):


        SV_identify("A2",("A2_Chr06",103080000,110720000,1),
                    "F1",("F1_Chr06",24400000,31000000,1),
                    "/root/run"
                    )
