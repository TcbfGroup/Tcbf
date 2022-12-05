import collections
import os.path
import pycircos
import matplotlib.pyplot as plt
from pandas import read_table
from pandas import concat
from pandas import read_csv
from tcbf.network_construct import get_species


def get_chromosome(workdir, genome):
    file = os.path.join(workdir, "Step1", f"{genome}.genome.fa.fai")
    result = []
    with open(file) as f:
        for line in f:
            line = line.split()
            name = line[0]
            length = int(line[1])
            result.append((name, length))
    return result


def parse_TAD_count(workdir, reference):
    count_file = os.path.join(workdir, "Result", "TAD_groups_count.tsv")
    count = read_table(count_file, index_col=0)
    count[count >= 1] = 1

    group_file = os.path.join(workdir, "Result", "TAD_groups.tsv")
    group_data = read_table(group_file, index_col=0)
    group_data["present_genome"] = count.sum(axis=1)

    unassign_file = os.path.join(workdir, "Result", "Unassignd_TAD.tsv")
    unassign = read_table(unassign_file)
    unassign["present_genome"] = 1

    merge = concat([group_data, unassign])

    merge.dropna(subset=[reference], inplace=True)

    dic = {}
    for k, v in merge.groupby("present_genome"):
        d = ";".join(v[reference]).split(";")
        dic[k] = d
    return dic


def get_TAD_position(workdir, reference, TAD_list):
    arcdata_dict = collections.defaultdict(dict)
    file = os.path.join(workdir, "Step1", f"{reference}.bound.bed")
    data = read_csv(file)
    data = data[data["tad_name"].isin(TAD_list)]
    for line in data.itertuples():
        # print(line[0])
        name = line[2]
        start = line[3]
        end = line[4]
        width = (end - start) *3
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"] = []
            arcdata_dict[name]["values"] = []
        arcdata_dict[name]["positions"].append(start)
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(100)

    return arcdata_dict


color = """
#a3dbd2
#ffe058
#cdcae3
#fc9f95
#fc427b
#fdb86a
#b3de69
#6060f7
""".split()


def plot_circos(workdir, reference,output):
    workdir = os.path.abspath(workdir)

    species_num = len(get_species(workdir))

    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle

    circle = Gcircle(figsize=(30, 30))
    chrom_info = get_chromosome(workdir, reference)

    min_circos = 350
    max_circos = 800
    every = int((max_circos - min_circos) / species_num)

    for line in chrom_info:
        if line[1] <= 2e+6:
            continue
        arc = Garc(arc_id=line[0],
                   size=line[1],
                   interspace=2,
                   raxis_range=(max_circos, max_circos + every * 0.4),
                   labelposition=70, label_visible=True,
                   labelsize=30,
                   facecolor="#ffbbb8",
                   edgecolor=None,
                   linewidth=0,
                   )
        circle.add_garc(arc)
    circle.set_garcs(0, 360)
    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id, raxis_range=(max_circos + every * 0.4, max_circos + every * 0.6),

                        tickinterval=20000000,
                        ticklabels=None)

    TAD_info = parse_TAD_count(workdir, reference)


    for i in range(species_num ,0,-1 ):
        data1 = get_TAD_position(workdir, reference, TAD_info[species_num -i + 1])

        min_position = min_circos + every * (i - 1)
        max_position = min_position + every * 0.8
        for key in data1:
            circle.barplot(key, data=data1[key]["values"],
                           positions=data1[key]["positions"],
                           width=data1[key]["widths"],
                           spine=False,
                           edgecolor="white",
                           raxis_range=[min_position, max_position],
                           rlim=[0, 100],

                           facecolor=color[i-1])





    circle.figure.savefig(output)

