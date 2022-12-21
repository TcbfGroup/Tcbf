from tcbf.visualization.plot_tad_structure import *
import matplotlib.pyplot as plt
import streamlit as st
import cooler
import numpy as np
import itertools


def plot_heatmap(ax, reference, cool, chrom, start, end, xrange, yrange, up_down, orientation=1):
    clr = cooler.Cooler(cool)
    M = clr.matrix(balance=False, sparse=False).fetch((chrom, start, end))
    M[np.isnan(M)] = 0
    n = M.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    product = itertools.product(range(n, -1, -1), range(0, n + 1, 1))
    A = np.dot(np.array([(ii[1], ii[0]) for ii in product]), t)
    cmap = 'RdYlBu_r'
    x = A[:, 1].reshape(n + 1, n + 1)
    y = A[:, 0].reshape(n + 1, n + 1)

    if up_down == -1:
        y[y > 0] = -y[y > 0]
    else:
        y[y < 0] = -y[y < 0]
    vmax = np.percentile(M[M.nonzero()], 95)
    vmin = M.min()
    vmin, vmax = float(vmin), float(vmax)
    vmin, vmax = st.slider("Select value", min_value=vmin, max_value=vmax, value=(vmin, vmax))
    x = coord_compression(x, x_range=xrange)
    y = coord_compression(y, x_range=yrange)
    if orientation == -1:
        x = 1 - x

    def plot_TAD(reference):
        tad_file = os.path.join(workdir, "Step1", f"{reference}.TAD.csv")

        tads = read_csv(tad_file).query(f"chromosome == '{reference}_{chrom}'")
        mask = (tads["end"] > start) & (tads['start'] <= end)
        tads = tads[mask]
        st.write(tads)
        for s, e in zip(tads['start'], tads['end']):
            si = int(s // res - start // res)
            ei = int(e // res - start // res)
            if si < 0:
                si = 0
            if ei > n - 1:
                ei = n - 1

            if ei - si < 2:
                continue

            x_data = [x[:-1, :-1][n - 1 - si, si],
                      x[:-1, :-1][n - 1 - si, ei],
                      x[:-1, :-1][n - 1 - ei, ei]]

            y_data = [y[:-1, :-1][n - 1 - si, si],
                      y[:-1, :-1][n - 1 - si, ei],
                      y[:-1, :-1][n - 1 - ei, ei]]

            ax.plot(x_data, y_data, color="gray", linestyle='-',
                    linewidth=2)

    ax.pcolormesh(x, y, np.flipud(M), cmap=cmap, vmin=vmin, vmax=vmax)
    # plot_TAD(reference)


def plot_pair_heat_map(workdir, query, target, chrom, start, end, query_cool, target_cool):
    fig, ax = plt.subplots(figsize=(9, 9))
    x_range = (0, 1)
    bound_position = {}
    bound_position[query] = plot_TAD(workdir, ax, query, f"{query}_{chrom}", start, end, 0.2, 0.21, x_range,
                                     orientation=1, )
    plot_heatmap(ax, query, query_cool, chrom, start, end, xrange=(0, 1), yrange=(0.22, 1), up_down=1)
    pair_genome_region = get_pairwise_genome_region(workdir=workdir,
                                                    reference=query,
                                                    chrom=chrom, start=start, end=end)
    chrom2, start2, end2, orientation2 = pair_genome_region[target][0]
    bound_position[target] = plot_TAD(workdir, ax, target, chrom2, start2, end2, -0.21, -0.2, (0, 1),
                                      orientation=orientation2, )
    plot_heatmap(ax, target, target_cool, chrom2.split("_")[1], start2, end2, xrange=(0, 1), yrange=(-1, -0.22),
                 up_down=-1,
                 orientation=orientation2)
    plot_boundary_pair(workdir, ax, query, bound_position, [query, target], style="d", color="red")
    # plot_SV(workdir,query,target,ax,(0,1),(0,1),(0.2 ,0.21),(-0.21 ,-0.2),
    #         genome1_pos = (f"{query}_{chrom}",start,end,1),
    #         genome2_pos = (chrom2,start2,end2,orientation2))

    ax.set_xlim(-0.2, 1)
    ax.set_ylim(-1, 1)
    ax.set_axis_off()

    return fig


res = 40e+3
workdir = "/root/run/"
query = "D5"
target = "C1"
chrom = "Chr03"
start = 42.96e+6
end = 45.74e+6

query_cool = f"/data/{query}.cool"
target_cool = f"/data/{target}.cool"

st.pyplot(plot_pair_heat_map(workdir, query, target, chrom, start, end, query_cool, target_cool))


# fig, ax = plt.subplots(figsize=(9, 9))
#
# x_range = (0, 1)
# bound_position = {}
# bound_position[query] = plot_TAD(workdir, ax, query, f"{query}_{chrom}", start, end, 0.2, 0.21, x_range,
#                                  orientation=1, )
#
# fig,ax = plt.subplots(figsize=(9, 9))
# plot_heatmap(ax, "","/data/matrix_10000.cool", "3", 154450000, 157850000, xrange=(0, 1), yrange=(0.22, 1), up_down=1)
# ax.set_xlim(-0.2, 1)
# ax.set_ylim(-1, 1)
# ax.set_axis_off()
# fig.saveplot("/data/")
# st.pyplot(fig)
