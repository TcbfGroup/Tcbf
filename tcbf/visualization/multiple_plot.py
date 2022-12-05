import streamlit as st
import itertools, cooler
import numpy as np
import pandas as pd
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt




def coord_compression(X, x_range):
    x_min = X.min()
    x_max = X.max()
    x_size = x_max - x_min
    ratio = (X - x_min) / x_size

    x_range_size = x_range[1] - x_range[0]

    new_value = np.full_like(X, x_range[0])
    new_value += x_range_size * ratio
    return new_value


def plot_heatmap(ax, matrix, region,chrom,start,end,bound,TAD_pos,genome):
    n = matrix.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    A = np.dot(np.array([(i[1], i[0]) for i in itertools.product(range(n, -1, -1), range(0, n + 1, 1))]), t)
    cmap = "RdYlBu_r"
    x = A[:, 1].reshape(n + 1, n + 1)

    y = A[:, 0].reshape(n + 1, n + 1)
    y[y < 0] = -y[y < 0]

    x = coord_compression(x, region[0])
    y = coord_compression(y, region[1])

    vmax = int(np.percentile(matrix[matrix.nonzero()], 90))
    vmin = int(matrix.min())

    ax.pcolormesh(x, y, np.flipud(matrix),

                  cmap=cmap, vmin=vmin, vmax=vmax,
                  edgecolor='none', snap=False, linewidth=.001, rasterized=True,
                  )

    def plot_TAD():
        tadtype = np.dtype({'names': ['chr', 'start', 'end'],
                            'formats': ['U5', np.int, np.int]})
        tads = np.loadtxt(TAD_pos, dtype=tadtype, usecols=[0, 1, 2])

        tads = tads[(tads['chr'] == chrom)]
        mask = (tads['end'] > start) & (tads['start'] < end)
        tads = tads[mask]

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

            ax.plot(x_data, y_data, color="red", linestyle='-',
                    linewidth=5)
    syn = []
    def get_boundary_location():
        boundary = pd.read_csv(bound)
        boundary = boundary[(boundary['chromosome'] == f"{genome}_" + chrom)]
        mask = (boundary['end'] > start) & (boundary['start'] < end)
        boundary = boundary[mask]
        for s, e in zip(boundary['start'], boundary['end']):

            si = int(s // res - start // res)
            ei = int(e // res - start // res)
            if si < 0:
                si = 0
            if ei > n - 1:
                ei = n - 1

            if ei - si < 2:
                continue
            syn.append([[x[:-1, :-1][n - 1 - si, si],region[1][0]],
                        [x[:-1, :-1][n - 1 - ei, ei], region[1][0]],
                        ])


    plot_TAD()
    get_boundary_location()
    return syn

def plot_syn(ax, point1, point2):
    Path = mpath.Path

    B = point2
    A = point1
    A1, A2 = A
    B1, B2 = B
    ax1, ay1 = A1
    ax2, ay2 = A2
    bx1, by1 = B1
    bx2, by2 = B2
    M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
    ymid = (ay1 + by1) / 2

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

    codes, verts = zip(*path_data)
    path = Path(verts, codes)

    pp1 = mpatches.PathPatch(
        path,
        fc="gainsboro", transform=ax.transData,
        lw=0.2, zorder=1,
        alpha=0.5,
        ec="k")

    ax.add_patch(pp1)
    ax.plot([2], [2], "")
    return




clr = cooler.Cooler("/root/data/D5.cool")
chrom1 = "Chr05"
start1 = 20.68e+6
end1 = 25.62e+6

M = clr.matrix(balance=False, sparse=False).fetch((chrom1, start1, end1))
M[np.isnan(M)] = 0

clr2 = cooler.Cooler("/data/A2.cool")
chrom2 = "Chr02"
start2 = 30.64e+6
end2 = 35.99e+6
M2 = clr2.matrix(balance=False, sparse=False).fetch((chrom2, start2, end2))
M2[np.isnan(M2)] = 0

res = clr.binsize
heatmap_pos1 = [0, 0, 1, 1]

fig = plt.figure(figsize=(7.5, 8))

h_ax = fig.add_axes(heatmap_pos1)

syn1 = plot_heatmap(h_ax, M, ([0.1, 0.9], [0, 0.4]),chrom1,start1,end1,"/root/run/Step1/D5.bound.bed",
                    "/root/data/D5.txt","D5")
syn2 = plot_heatmap(h_ax, M2, ([0.1, 0.9], [0.45, 0.85]),chrom2,start2,end2,"/root/run/Step1/A2.bound.bed",
                    "/root/data/A2.txt","A2")

button = st.button("rerun")
if button:
    st.experimental_rerun()


for line in zip(syn1,syn2):
    plot_syn(h_ax,line[0],line[1])


    # st.write(line)
h_ax.set_xlim(0, 1)
h_ax.set_ylim(0, 1)
h_ax.set_axis_off()

st.pyplot(fig)
