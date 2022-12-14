import os.path
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import pandas as pd
import matplotlib.patches as mpatches
from pandas import read_csv, read_table
import matplotlib.path as mpath
from numpy import full_like, float64
import streamlit as st


def plot_syn(ax, point1, point2, style="curve", color="#e8effb", zorder=1, alpha=0.9):
    Path = mpath.Path
    A = point1
    B = point2
    A1, A2 = A
    B1, B2 = B
    ax1, ay1 = A1
    ax2, ay2 = A2
    bx1, by1 = B1
    bx2, by2 = B2
    M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
    ymid = (ay1 + by1) / 2
    if style == "curve":
        path_data = [
            (M, A1),
            (C4, (ax1, ymid)),
            (C4, (bx1, ymid)),
            (C4, B1),
            (L, B2),
            (C4, (bx2, ymid)),
            (C4, (ax2, ymid)),
            (C4, A2),
            (CP, A1),
        ]
    else:
        path_data = [(M, A1), (L, B1), (L, B2), (L, A2), (CP, A1)]

    codes, verts = zip(*path_data)
    path = Path(verts, codes)

    pp1 = mpatches.PathPatch(
        path,
        fc=color, transform=ax.transData,
        lw=0.2, zorder=zorder,
        alpha=alpha,
        ec="k")
    ax.add_patch(pp1)

    return


def coord_compression(X, x_range,
                      x_min=None,
                      x_max=None):
    if x_min is None:
        x_min = X.min()
    if x_max is None:
        x_max = X.max()

    x_size = x_max - x_min
    ratio = (X - x_min) / x_size
    ratio[ratio <= 0] = 0
    ratio[ratio >= 1] = 1
    x_range_size = x_range[1] - x_range[0]
    new_value = full_like(X, x_range[0])
    new_value =  new_value + x_range_size * ratio
    return new_value


# def position_to_boundary(workdir, reference, chrom, start, end, tads=None):
#     try:
#         if tads is not None:
#             start = tads.iloc[0, 4]
#             end = tads.iloc[-1, 7]
#     except:
#         st.write(reference, start, end)
#     file = os.path.join(workdir, "Step1", f"{reference}.bound.bed")
#     bounds = read_csv(file).query(f"chromosome == '{chrom}'")
#     mask = (bounds["start"] >= start) & (bounds["end"] <= end)
#     bounds = bounds[mask]
#
#     return bounds

def position_to_boundary_TAD(workdir, reference, chrom, start, end,
                             xrange=(0, 1)):
    def parse_data(df):
        df = df.query(f"chromosome == '{chrom}' & start <= {end} & end >= {start}")
        d = df.iloc[:, -2:]
        df.iloc[:, -2:] = coord_compression(d.astype(float64).to_numpy(), x_range=xrange,
                                            x_max=end,
                                            x_min=start)
        return df
    bound_file = os.path.join(workdir, "Step1", f"{reference}.bound.bed")
    bounds = parse_data(read_csv(bound_file))

    tad_file = os.path.join(workdir, "Step1", f"{reference}.TAD.csv")
    tads = parse_data(read_csv(tad_file).iloc[:, [3, 0, 1, 2]])
    return bounds, tads


def print_coordinate(chrom, start, end):
    return f"{chrom.split('_')[0]}  {chrom.split('_')[1]}\n{start / 1e+6}-{end / 1e+6} Mb"


s = (print_coordinate("G1_Chr02", 11.2e+6, 22.345e+6))



def plot_TAD(workdir, ax, reference, chrom, start, end, y_start, y_end, x_range, orientation=1, ):


    height = y_end - y_start
    bounds, tads = position_to_boundary_TAD(workdir, reference, chrom, start, end, x_range)

    if orientation == -1:
        columns = bounds.columns
        bounds.iloc[:, -2:] = x_range[1] + x_range[0] - bounds.iloc[:, -2:]
        bounds = bounds.iloc[:, [0, 1, 3, 2]]
        bounds.columns = columns

        tads.iloc[:, -2:] = x_range[1] + x_range[0] - tads.iloc[:, -2:]
        tads = tads.iloc[:, [0, 1, 3, 2]]
        tads.columns = columns

    coord = print_coordinate(chrom,start,end) if orientation == 1 else print_coordinate(chrom,end,start)
    ax.text(-0.15,y_start,coord,alpha = 0.8,color = "#798a9c",horizontalalignment  = "center")

    tads["width"] = tads["end"] - tads["start"]
    bounds["width"] = bounds["end"] - bounds["start"]

    tad_collection = []
    for item in tads.itertuples():
        rect = Rectangle(
            (item[3], y_start),
            item[-1], height)
        tad_collection.append(rect)
    tad_collection = PatchCollection(tad_collection, zorder=2, color="#f19f9c")
    ax.add_collection(tad_collection)

    bound_coord_dict = {}

    bound_collection = []

    for index, item in enumerate(bounds.itertuples()):
        rect = Rectangle((item[3], y_start), item[-1], height,
                         )
        bound_collection.append(rect)

        bound_coord_dict[item[1]] = (item[3], y_start, item[-1], height)

    bound_collection = PatchCollection(bound_collection, zorder=3, edgecolors="black", facecolors="yellow")
    ax.add_collection(bound_collection)

    rect = Rectangle((x_range[0],y_start),x_range[1] ,height,color = "grey",alpha = 0.7)
    ax.add_patch(rect)
    return bound_coord_dict


#
# def plot_TAD_position(workdir, reference, chrom, start, end,
#                       ax, y_start, y_end, orientation=1,
#                       xrange=(0, 1)):
#     # height = y_end - y_start
#     # file = os.path.join(workdir, "Step1", f"{reference}.TAD.csv")
#     # tads = read_csv(file).query(f"chromosome == '{chrom}'")
#     # mask = (tads["left_start"] >= start) & (tads["right_end"] <= end)
#     # tads = tads[mask]
#     # boundary_name = list(position_to_boundary(workdir, reference, chrom, start, end, tads)["tad_name"])
#     #
#     # number = tads.iloc[:, [1, 2, 4, 5, 6, 7]]
#     # new = coord_compression(number.astype(float64).to_numpy(), x_range=xrange)
#     # new_coords = pd.DataFrame(new, columns=number.columns)
#     # new_coords["width"] = new_coords["end"] - new_coords["start"]
#     #
#     # s1 = new_coords.loc[:, ["left_start", "left_end"]]
#     # s2 = new_coords.loc[:, ["right_start", "right_end"]]
#     # s2.columns = s1.columns
#     # bound_coord = pd.concat([s1, s2], axis=0).drop_duplicates().sort_values(["left_start"]). \
#     #     reset_index(drop=True)
#     # bound_coord["width"] = bound_coord["left_end"] - bound_coord["left_start"]
#     if orientation == -1:
#         columns = bound_coord.columns
#
#         bound_coord.iloc[:, :2] = xrange[1] + xrange[0] - bound_coord.iloc[:, :2]
#         bound_coord = bound_coord.iloc[:, [1, 0, 2]]
#         bound_coord.columns = columns
#
#         columns = new_coords.columns
#         new_coords.iloc[:, :2] = xrange[1] + xrange[0] - new_coords.iloc[:, :2]
#         new_coords.iloc[:, [0, 1]] = new_coords.iloc[:, [1, 0]]
#         new_coords.columns = columns
#
#     tad_collection = []
#     for item in new_coords.itertuples():
#         rect = Rectangle(
#             (item[1], y_start),
#             item[-1], height)
#         tad_collection.append(rect)
#     tad_collection = PatchCollection(tad_collection, zorder=2, color="#f19f9c")
#     ax.add_collection(tad_collection)
#
#     bound_collection = []
#
#     bound_coord_collection = []
#     bound_coord_dict = {}
#
#     assert bound_coord.shape[0] == len(boundary_name)
#     for index, item in enumerate(bound_coord.itertuples()):
#         rect = Rectangle((item[1], y_start), item[-1], height)
#         bound_collection.append(rect)
#         bound_coord_collection.append(((item[1], y_start), (item[2], y_start)))
#         bound_coord_dict[boundary_name[index]] = (item[1], y_start, item[-1], height)
#     bound_collection = PatchCollection(bound_collection, zorder=3, color="green")
#     ax.add_collection(bound_collection)
#     return bound_coord_dict


def get_orientation(ia: [int], ib: [int]) -> int:
    """
    Infer the orientation of a pairwise block.
    Args:
        ia (List[int]): List a
        ib (List[int]): List b
    Returns:
        str: plus (+) or minus (-)
    """
    if len(ia) != len(ib) or len(ia) < 2:
        return 1  # Just return a default orientation
    from numpy import polyfit
    slope, _ = polyfit(ia, ib, deg=1)
    return 1 if slope >= 0 else -1


def get_pair_TAD_boundary_region(workdir, pair_df, only_primary):
    if pair_df.shape[1] != 2:
        sys.exit()
    paired_boundary = [j for i in pair_df.iloc[:, 1].dropna() for j in i.split(";")]

    reference = pair_df.columns[1]

    file = os.path.join(workdir, "Step1", f"{reference}.bound.bed")
    boundary = read_csv(file)
    boundary.index = boundary["tad_name"]
    need = boundary.loc[paired_boundary]

    need.reset_index(drop=True, inplace=True)
    result = []

    ress = None
    num = 0
    for index, res in need.groupby("chromosome"):
        bcol = res["tad_name"].str.split("_").str.get(-1).astype(int)

        orientation = get_orientation(range(len(bcol)), bcol)
        res = res.sort_values("start")
        start = res.iloc[0, 2]
        end = res.iloc[-1, 3]
        result.append([index, start, end, orientation])
        if len(bcol) >= num:
            ress = index, start, end, orientation
            num = len(bcol)
    if only_primary:
        return [ress]
    return result


def get_pairwise_genome_region(workdir, reference,
                               chrom,
                               start,
                               end,
                               only_primary=True):
    # file = os.path.join(workdir, "Step1", f"{reference}.bound.bed")
    # boundary = read_csv(file)

    chrom = reference + "_" + chrom
    # boundary = boundary.query(f"chromosome == '{chrom}'")
    # mask = (boundary["end"] >= start) & (boundary["start"] <= end)
    # boundary = boundary[mask]["tad_name"]
    boundary = position_to_boundary_TAD(workdir, reference, chrom, start, end, xrange=(0, 1))[0]["tad_name"]

    one_to_many = os.path.join(workdir, "Result", "one_to_many", f"{reference}.txt")
    pair = read_table(one_to_many)
    pair = pair[pair[reference].isin(boundary)]

    res = {}
    for index in range(1, pair.shape[1]):
        pair_genome_region = get_pair_TAD_boundary_region(workdir, pair.iloc[:, [0, index]], only_primary)

        res[pair.columns[index]] = pair_genome_region
    return res





def plot_SV(workdir, name1, name2, ax, x1_ranges, x2_ranges,
            y1_ranges, y2_ranges, genome1_pos, genome2_pos):
    from tcbf.visualization.run_syri import SV_identify
    SV = SV_identify(name1, genome1_pos,
                     name2, genome2_pos, workdir)

    SV.iloc[:, [1, 2]] = coord_compression(SV.iloc[:, [1, 2]].astype(float).to_numpy(),
                                           x1_ranges, x_min=1, x_max=genome1_pos[2] - genome1_pos[1])
    SV.iloc[:, [4, 5]] = coord_compression(SV.iloc[:, [4, 5]].astype(float).to_numpy(),
                                           x2_ranges, x_min=1, x_max=genome2_pos[2] - genome2_pos[1])

    SYN = SV.query("type == 'SYN'")
    for line in SYN.itertuples():
        point1 = (line[2], y1_ranges[1]), (line[3], y1_ranges[1])
        point2 = (line[5], y2_ranges[0]), (line[6], y2_ranges[0])

        plot_syn(ax, point1, point2, color="gray", alpha=0.4)

    SYN = SV.query("type == 'INV'")
    for line in SYN.itertuples():
        point1 = (line[3], y1_ranges[1]), (line[2], y1_ranges[1])
        point2 = (line[5], y2_ranges[0]), (line[6], y2_ranges[0])

        plot_syn(ax, point1, point2, color="yellow")
    # def plot_SYN(dataframe):
    #     pass
    # st.stop()


def plot_boundary_pair(workdir, ax, reference, boundary_position, species_order, style="curve"):
    from itertools import product
    # file = os.path.join(workdir, "Result", "one_to_many", f"{reference}.txt")
    # data = read_table(file)
    # reference_boundary = boundary_position.get(reference)
    # data = data[data[reference].isin(reference_boundary.keys())].loc[:, species_order]
    # st.write(data)
    all_need_boundary = set(j for i in boundary_position.values() for j in list(i.keys()))
    file = os.path.join(workdir,"Result","TAD_groups.tsv")
    import re
    datas = []
    with open(file)as f:
        head = (f.readline().split())
        for item in f:
            d = set(re.split("\t|;",item.strip()))
            if all_need_boundary & d :
                datas.append(item.split("\t")[1:])
    df = pd.DataFrame(datas)
    df.columns = head
    data = df.loc[:,species_order]

    data[data.applymap(lambda x:"bound" not in x)] = None

    for s in species_order[1:]:
        tmp = data.loc[:, [reference, s]].dropna().apply(lambda x:x.str.strip()).applymap(lambda x: x.split(";")).apply(
            lambda x: list(product(*x)), axis=1).to_list()
        # st.write(data.loc[:, [reference, s]].dropna())

        pair = [j for i in tmp for j in i]

        for s1, s2 in pair:
            d1 = boundary_position.get(reference)
            if not d1:
                continue
            d2 = boundary_position.get(s)
            if not d2:
                continue

            x1 = d1.get(s1, False)
            x2 = d2.get(s2, False)
            if x1 and x2:

                x, y, width, height = x1
                point1 = (x, y + height), (x + width, y + height)
                x, y, width, height = x2
                point2 = (x, y), (x + width, y)

                plot_syn(ax, point1, point2, style=style, zorder=2)

        reference = s






def plot_main(workdir,genome,chrom,start,end,species_order,):
    fig, ax = plt.subplots()
    pair_genome_region = get_pairwise_genome_region(workdir=workdir,
                                                    reference=genome,
                                                    chrom=chrom, start=start, end=end)

    # reference_span = end - start
    # var = max([i[0][2] - i[0][1] for i in pair_genome_region.values() if i[0]])
    # max_span = max(reference_span, var)
    species_num = len(species_order)
    one_space = 1 / species_num
    bound_position = {}

    bound_position[genome] = plot_TAD(workdir=workdir,
                                      reference=genome,
                                      chrom=f"{genome}_{chrom}",
                                      start=start,
                                      end=end,
                                      ax=ax, y_start=0,
                                      y_end=0.01,
                                      x_range=(0, 1)
                                      # x_range=(0.5 - ratio, 0.5 + ratio)
                                      )

    # reference = genome
    # position = (f"{genome}_{chrom}", start, end, 1)
    # x_ranges = (0.5 - ratio, 0.5 + ratio)
    # y_ranges = (0, 0.01)
    for index in range(1, len(species_order)):
        value = pair_genome_region[species_order[index]]
        species, v = species_order[index], value[0]
        if not v:
            continue
        bound_position[species] = plot_TAD(workdir=workdir,
                                           reference=species,
                                           chrom=v[0],
                                           start=v[1],
                                           end=v[2],
                                           ax=ax,
                                           y_start=index * one_space,
                                           y_end=(index * one_space) + 0.01,
                                           orientation=v[-1],
                                           x_range=(0, 1)
                                           # x_range=(0.5 - ratio, 0.5 + ratio)
                                           )





    plot_boundary_pair(workdir,
                       ax, genome, bound_position, species_order, style="fdg")
    ax.set_xlim(-0.2, 1)
    ax.set_ylim(-0.2, 1)
    ax.set_axis_off()
    return fig

def plot_synteny(workdir,reference,chrom,start,end,species_order,output):
    fig = plot_main(workdir,reference,chrom,start,end,species_order)
    fig.savefig(output,dpi = 300)


if __name__ == "__main__":
    with st.sidebar:
        with st.form("my form"):

            genome_list = " B1 C1 D5 E1 F1 G1 K2 A2".split()
            genome_list2 = st.multiselect("Select you want deplayed genome",genome_list,
                                         default=genome_list[:5])

            workdir = "/root/run"
            genome = st.radio("Pick one",genome_list)
            genome_list2.remove(genome)

            chrom_list = [f"Chr{str(i).zfill(2)}" for i in range(1,14)]
            chrom = st.radio("Select a Chrom",chrom_list)


            start,end  = st.slider("double ended slider",
                                   min_value=0.1,
                                   max_value=100.,
                                   value = (12.0,23.4))
            end = end * 1000000
            start = start * 1000000
            genome_list2.insert(0,genome)

            submitted = st.form_submit_button("Submit")
            if not submitted:
                st.stop()
    fig = plot_main("/root/run",genome,chrom,start,end,genome_list2)
    st.pyplot(fig)