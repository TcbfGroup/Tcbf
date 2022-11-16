<h1 align="center">Tcbf</h1>
<p align="center">
    <em> ❤️ Tcbf = Topological association structure(TAD) Conservative Boundary Finder</em>
</p>
<p>
    <a href="https://opensource.org/licenses/MIT">
        <img src="https://img.shields.io/badge/License-MIT-brightgreen.svg" alt="License">
    </a>
</p>

[中文 README](README.ch.md)

## 📣 Introduction
___
TADs are fundamental regulatory chromatin structures and are
largely conserved across tissues and species. We developed 
a python pipeline Tcbf to identify the conservative TAD boundary between
multiple genome.
___
## ✨ Pre-requisite:
### [MCL](https://github.com/micans/mcl)
The mcl clustering algorithm is available in the
repositories of some Linux distributions and so can be
installed in the same way as any other package. 
For example, on Ubuntu, Debian, Linux Mint:

`sudo apt-get install mcl`


Alternatively, it can be built from source which 
will likely require the 'build-essential'
or equivalent package on the Linux distribution being 
used. Instructions are provided on the MCL webpage, 
http://micans.org/mcl/.

### [Minimap2](https://github.com/lh3/minimap2)
The default aligner is minimap2. In the future, we will consider supporting more aligner.

### [Mash](https://github.com/marbl/Mash)
Mash is used to measure the genetics distance between different genome sequences and adjust the alignment parameters.


### R language and the following packages

```
install.packages("BiocManager")
install.packages("feather")
BiocManager::install("GenomicRanges")
BioManager::install("plyranges")
install.packages("tidyverse")
```

### Python >= 3.6
This package is heavily dependent on f-string. As of Python 3.6, f-strings are a great new way to format strings. Not only are they more readable, more concise, and less prone to error than other ways of
formatting, but they are also faster!
___
## 🔰 Installation

**pip install**
```shell
pip install Tcbf
```

**Install from source**
```shell
git clone https://github.com/hexin010101/Tcbf
cd Tcbf
pip install -r requirements.txt
python setup.py install
```

### Quick installation using docker

`docker pull Tcbf`

Singularity container source

`singularity pull Tcbf.sif docker://Tcbf`
___
## 📝 Usage

You can test the Tcbf pipeline with a toy project. which takes about five mins.
To run Tcbf on the example Data:
```
Tcbf.py run -c config.txt  -o test
```

### Inputs
The config.txt is a table file with three columns. 
The first columns is path of genome sequence path.

The second columns is TAD annotaed file.


The third column is the individual name.


Of note, pay attention to the contigs or scaffolds. For organisms have well assembly,
we usually recommand to remove all scaffolds.
```
genome1.fasta   genome1_tad.txt g1
genome2.fasta   genome2_tad.txt g2
genome3.fasta   genome3_tad.txt g3
genome4.fasta   genome4_tad.txt g4
```

we use the results from [HiTAD](https://academic.oup.com/nar/article/45/19/e163/4093166) 
as the TAD annotate file 
```
#cat genome1_tad.txt
Chr01   0       2300000
Chr01   2300000 3950000
Chr01   3950000 4600000
Chr01   4600000 5900000
Chr01   5900000 6350000
```
if your TAD annotated file from other software, we also provide a [method]() 
to convert.


### Output
For conserved TAD boundaries among species, we provide a table file to show the 
the clustering results. This table has one TAD boundary group per line and one spcies per column and is ordered
from largest orthogroup to smallest.



| g1                                                  | g2                                   |
|-----------------------------------------------------|--------------------------------------|
| g1_0_left_first,g1_12_left_middle,g1_25_left_middle | g2_12_left_middle                    |
| g1_34_left_first,g1_34_left_middlee                 | g2_124_right_last,g2_14_left_middle, |


For each pair of genome comparisons, we provided a JSON file 
to explore the conserved TAD structure. For More detail, see [Here]().
```
{"A_0":["B_10"], ## highly conserved
"A_1":["B_2","B_3"] ## TAD splited              
}
```
### Advanced
For HPC users have multiple nodes, we also provide a step by step to accelerate work.
```
tcbf run -c config.txt  -o test --only_print_command
```
For more details, see [Here]()



## 😉 Author
Tbcf are maintained by: * [@HeXin](https://github.com/hexin010101)


For more Help, Please leave a message in the issue, 
I will reply as soon as possible.


## 📃 License

MIT 
