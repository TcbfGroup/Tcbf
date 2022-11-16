import os.path
import pyfastx
from pandas import DataFrame
from pandas import read_table
import pandas as pd
import click
# from run_command import run_command


import subprocess
import sys
def run_command(command):
    """
    此段代码借鉴了  109-130行  https://github.com/davidemms/OrthoFinder/blob/master/scripts_of/__main__.py
    command : 需要执行的命令
    """

    capture = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    stdout, stderr = capture.communicate()
    try:
        stdout = stdout.decode()
        stderr = stderr.decode()
    except(UnicodeError, AttributeError):
        stdout = stdout.encode()
        stderr = stderr.encode()
    n_stderr_lines = stderr.count('\n')
    if capture.returncode != 0:
        print(f"Returned error :{capture.returncode}")
        print(f"Command: {command}")
        print("stdout:\n-------")
        print(stdout)
        if n_stderr_lines > 0:
            print("stderr:\n-------")
            print(stderr)
        sys.exit()
    return stdout

"""
python extract_TAD_bound.py  -t /data/A_genome_tad.txt \
-g /data/Garboreum_genome_HAU_v1.0/Lachesis_assembly_changed.fa -d 150000 -p A -o ~/data

"""


def check_dependency():
    try:
        import pysam
        from pysam import Fastafile
        return "pysam"
    except ImportError:
        try:
            import pyfastx
            from pyfastx import Fasta
            return "pyfastx"
        except ImportError:
            from Bio import SeqIO
            return


class TAD:
    def __init__(self,
                 chromosome: str,
                 start: int,
                 end: int,
                 bound_percent: float,
                 position=None,
                 distance: int = None,
                 name=None):
        self.name = name
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.bound_percent = bound_percent
        self.distance = distance
        self.length = self.end - self.start + 1
        self._position = position
        if self.bound_percent:
            self._left_coord = (max(0, self.start - int(self.bound_percent * self.length)),
                                self.start + (self.bound_percent * self.length))
            self._right_coord = (max(0, self.end - int(self.bound_percent * self.length)),
                                 self.end + (self.bound_percent * self.length))

        else:
            self._left_coord = (max(0, self.start - self.distance),
                                self.start + self.distance)
            self._right_coord = (max(0, self.end - self.distance),
                                 self.end + self.distance)
        self._left_coord = tuple(int(i) for i in self._left_coord)
        self._right_coord = tuple(int(i) for i in self._right_coord)

    @property
    def position(self):
        return self._position

    def set_position(self, position=None):
        self._position = position

    @property
    def left_bound(self):
        return self.chromosome, *self._left_coord

    @property
    def right_bound(self):
        return self.chromosome, *self._right_coord

    def __repr__(self):
        return f"{self.name}\t {self.chromosome}\t{self.start}\t{self.end} Length:{self.length}\t\n" \
               f"Left bound:\t{self.left_bound}Bound length:\t{self.right_bound[2] - self.left_bound[1] + 1}\n" \
               f"Right bound:\t{self.right_bound}\n"


def bed2position(chrom, start, end):
    return f"{chrom}:{start + 1}-{end}"


def format_seq(seq):
    return "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))


def add_prefix(fa, prefix, out):
    seq = pyfastx.Fasta(fa)
    with open(out, "w") as f:
        for line in seq:
            f.write(f">{prefix}_{line.name}\n{format_seq(line.seq)}")


def mash_genome(genome_fasta):
    command = f"mash sketch {genome_fasta}"
    run_command(command)


class TADs:
    def __init__(self, tad_file: str,
                 fasta_file: str,
                 bound_percent: float = None,
                 distance: int = None,
                 prefix: str = "TAD",
                 ):
        if (not bound_percent) and (not distance):
            raise TypeError
        self.TAD_bound_table = None
        self.prefix = prefix

        tad_table = read_table(tad_file, header=None)
        tad_table.columns = "chromosome start end".split()
        tad_table["chromosome"] = prefix + "_" + tad_table["chromosome"]
        self.bound_percent = bound_percent
        self.distance = distance
        # 列表存储每一个TAD信息
        self.TADs = []
        # 构建一个dataframe
        tad_infos = []
        num = 0
        for line in tad_table.groupby("chromosome"):
            data = []
            tad_info = []
            line[1].index = range(line[1].shape[0])
            for info in line[1].itertuples():
                tad = TAD(*list(info[1:]),
                          bound_percent=bound_percent,
                          position=None,
                          distance=distance,
                          name=self.prefix + f"_{num}")
                num += 1
                data.append(tad)
                tad_info.append((tad.name, *tad.left_bound[1:], *tad.right_bound[1:]))
            tad_bound = DataFrame(tad_info)
            tad_bound.columns = "tad_name left_start left_end right_start right_end".split()
            ss = pd.concat([line[1], tad_bound], axis=1)
            tad_infos.append(ss)

            data[0].set_position("first")
            data[-1].set_position("last")
            self.TADs.extend(data)

        self.TAD_table = pd.concat(tad_infos)

        self.imported_module = check_dependency()
        if self.imported_module == "pysam":
            from pysam import Fastafile
            self.fasta = Fastafile(fasta_file)
        elif self.imported_module == "pyfastx":
            from pyfastx import Fasta
            self.fasta = Fasta(fasta_file)
        else:
            from Bio import SeqIO
            self.fasta = {i.id: i for i in SeqIO.parse(fasta_file, "fasta")}

    def export_table(self, FileName):
        self.TAD_table.to_csv(FileName + ".csv", index=False)

    def _extract_seq(self, chrom, start, end):
        if self.imported_module == "pysam":
            position = bed2position(chrom, start, end)
            seq = self.fasta.fetch(region=position)

        elif self.imported_module == "pyfastx":
            seq = self.fasta[chrom][start:end].seq
        else:
            seq = self.fasta[chrom][start:end].seq

        return format_seq(seq)

    def extract_bound_bed(self, handle):
        s1 = self.TAD_table.loc[:, ["chromosome", "left_start", "left_end"]]
        s2 = self.TAD_table.loc[:, ["chromosome", "right_start", "right_end"]]
        s1.columns = ("chromosome", "start", "end")
        s2.columns = s1.columns
        self.TAD_bound_table = pd.concat([s1, s2], axis=0).drop_duplicates().sort_values(["chromosome", "start"]). \
            reset_index(drop=True)
        self.TAD_bound_table.index = self.TAD_bound_table.index.set_names(["tad_name"])
        self.TAD_bound_table.reset_index(inplace=True)
        self.TAD_bound_table["tad_name"] = self.prefix + "_bound_" + self.TAD_bound_table["tad_name"].astype(str)

        self.TAD_bound_table.to_csv(handle, index=False)

    def extract_bound_seq(self, out_fasta):
        for tad_bound in self.TAD_bound_table.itertuples():
            seq = self._extract_seq(*tad_bound[2:])
            out_fasta.write(f">{tad_bound[1]} {tad_bound[2:]}\n{seq}\n")

    # def extract_bound_seq(self, out_fasta: NamedTemporaryFile):
    #     """
    #     提取TAD的边界序列，对于每一条染色体的TAD ，如果是按照TAD边界上下游固定距离
    #     提取每一个TAD的左边界和最后一个TAD的右边界
    #     如果是
    #     :param out_fasta: 文件句柄，可以是一个临时文件，也可以是一个打开的文件
    #     :return:
    #     """
    #     for index, tad in enumerate(self.TADs):
    #         left_seq = self._extract_seq(*tad.left_bound)
    #         if self.distance:
    #             if tad.position == "first":
    #                 out_fasta.write(f">{tad.name}_left_first    {tad.left_bound}\n{left_seq}\n")
    #             elif tad.position is None:
    #                 out_fasta.write(f">{tad.name}_left_middle   {tad.left_bound}\n{left_seq}\n")
    #             elif tad.position == "last":
    #                 right_seq = self._extract_seq(*tad.right_bound)
    #                 out_fasta.write(f">{tad.name}_left_middle   {tad.left_bound}\n{left_seq}\n")
    #                 out_fasta.write(f">{tad.name}_right_last     {tad.right_bound}\n{right_seq}\n")
    #         else:
    #             right_seq = self._extract_seq(*tad.right_bound)
    #
    #             out_fasta.write(f">{tad.name}_left    {tad.left_bound}\n{left_seq}\n")
    #             out_fasta.write(f">{tad.name}_right    {tad.right_bound}\n{right_seq}\n")
    #
    # def export_bound_bed(self, handle):
    #     handle.write("tad_name\tchromosome\tstart\tend\n")
    #     for index, tad in enumerate(self.TADs):
    #         data = f"{tad.left_bound[0]}\t{tad.left_bound[1]}\t{tad.left_bound[2]}"
    #         if tad.position == "first":
    #             handle.write(f"{tad.name}_left_first\t{data}\n")
    #         elif tad.position is None:
    #             handle.write(f"{tad.name}_left_middle\t{data}\n")
    #         elif tad.position == "last":
    #             right_data = f"{tad.right_bound[0]}\t{tad.right_bound[1]}\t{tad.right_bound[2]}"
    #             handle.write(f"{tad.name}_left_middle\t{data}\n")
    #             handle.write(f"{tad.name}_right_last\t{right_data}\n")


def extract_TAD_boundary(tad: str,
                         genome: str,
                         distance: int,
                         prefix: str,
                         output: str):
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(os.path.join(output, "Step1")):
        os.mkdir(os.path.join(output, "Step1"))

    output = os.path.join(output, "Step1")
    genome_file = os.path.join(output, f"{prefix}.genome.fa")

    add_prefix(genome, prefix, genome_file)
    mash_genome(genome_file)
    G1_TAD = TADs(tad,
                  genome_file,
                  distance=distance,
                  prefix=prefix,
                  )

    G1_TAD.export_table(os.path.join(output, f"{prefix}.TAD"))

    with open(os.path.join(output, f"{prefix}.bound.bed", ), "w") as f:
        G1_TAD.extract_bound_bed(f)
    with open(os.path.join(output, f"{prefix}.bound.fasta", ), "w") as f:
        G1_TAD.extract_bound_seq(f)


@click.command()
@click.option("-t", '--TAD', type=click.Path(exists=True), help="TAD file ", required=True)
@click.option("-g", '--genome', type=click.Path(exists=True), help="genome sequence file  ", required=True)
@click.option('-d', "--distance", type=int, default=150000, required=False)
@click.option('-p', "--prefix", type=str, default=None, help="输出文件的前缀", required=True)
@click.option('-o', "--output", type=str, default='.', help="输出文件的目录", required=False)
def main(tad: str,
         genome: str,
         distance: int,
         prefix: str,
         output: str):
    extract_TAD_boundary(tad,genome,distance,prefix,output)


if __name__ == '__main__':
    pd.options.mode.chained_assignment = None
    main()
