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
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('interaction',
                                             ['#FFFFFF', '#FFDFDF', '#FF7575', '#FF2626', '#F70000'])
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
        max_value = 0
        tad_file = os.path.join(workdir, "Step1", f"{reference}.TAD.csv")

        tads = read_csv(tad_file).query(f"chromosome == '{reference}_{chrom}'")
        mask = (tads["end"] > start) & (tads['start'] <= end)
        tads = tads[mask]

        for s, e in zip(tads['start'], tads['end']):
            si = int(s // res - start // res)
            ei = int(e // res - start // res)
            if si <= 0:
                si = 0
            if ei > n - 1:
                ei = n - 1

            if ei - si < 2:
                continue

            x_data = [x[:-1, :-1][n - 1-si, si],
                      x[:-1, :-1][n - 1- si, ei],
                      x[:-1, :-1][n - 1- ei, ei]]

            y_data = [y[:-1, :-1][n -1- si , si],
                      y[:-1, :-1][n -1- si, ei],
                      y[:-1, :-1][n -1- ei, ei]]

            ax.plot(x_data, y_data, color="gray", linestyle='-',
                    linewidth=2)
            m = abs(y_data[1])
            if m > max_value:
                max_value = m
        return max_value
    #
    # max_value = plot_TAD(reference)
    # y = y[:-1, :-1]
    # x = x[:-1, :-1]
    #
    # plot_M = np.flipud(M)[:x.shape[0],:y.shape[0]]
    # if up_down == -1:
    #     max_value = -max_value
    #     plot_M[y < max_value] = 0
    # else:
    #     plot_M[y>max_value] = 0
    plot_TAD(reference)
    plot_M = np.flipud(M)



    ax.pcolormesh(x, y,plot_M, cmap=cmap, vmin=vmin, vmax=vmax)


def plot_pair_heat_map(workdir, query, target, chrom, start, end, query_cool, target_cool):
    fig, ax = plt.subplots(figsize=(9, 9))
    x_range = (0, 1)
    bound_position = {query: plot_TAD(workdir, ax, query, f"{query}_{chrom}", start, end, 0.2, 0.21, x_range,
                                      orientation=1, )}
    pair_genome_region = get_pairwise_genome_region(workdir=workdir,
                                                    reference=query,
                                                    chrom=chrom, start=start, end=end)
    chrom2, start2, end2, orientation2 = pair_genome_region[target][0]
    bound_position[target] = plot_TAD(workdir, ax, target, chrom2, start2, end2, 0.0, 0.01, (0, 1),
                                      orientation=orientation2)
    plot_heatmap(ax, target, target_cool, chrom2.split("_")[1], start2, end2, xrange=(0, 1), yrange=(-0.01, -0.8),
                 up_down=1,
                 orientation=orientation2)
    plot_heatmap(ax, query, query_cool, chrom, start, end, xrange=(0, 1), yrange=(0.21, 1), up_down=1)


    plot_boundary_pair(workdir, ax, query, bound_position, [query, target], style="curve")
    ax.set_xlim(-0.2, 1)
    ax.set_ylim(-1, 1)
    ax.set_axis_off()

    return fig


res = 20e+3
workdir = "/data/cool/o2"
query = "Ahypogaea"
target = "Pvalgaris"
chrom = "11"
start =5280000
end = 12200000

query_cool = f"/data/cool/peanut.cool::20000"
target_cool = f"/data/cool/common_bean.cool::20000"

st.pyplot(plot_pair_heat_map(workdir, query, target, chrom, start, end, query_cool, target_cool))
