from tcbf.visualization.plot_tad_structure import *
import matplotlib.pyplot as plt
import streamlit as st
import cooler
import numpy as np
import itertools

def plot_heatmap(ax, cool, chrom, start, end,xrange,yrange,up_down,orientation = 1):
    clr = cooler.Cooler(cool)
    M = clr.matrix(balance=False, sparse=False).fetch((chrom, start, end))
    M[np.isnan(M)] = 0
    n = M.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    product = itertools.product(range(n, -1, -1), range(0, n + 1, 1))
    A = np.dot(np.array([(ii[1], ii[0]) for ii in product]), t)
    cmap = 'RdYlBu_r'
    x =  A[:, 1].reshape(n + 1, n + 1)
    y = A[:, 0].reshape(n + 1, n + 1)

    if up_down == -1:
        y[y > 0] = -y[y > 0]
    else:
        y[y < 0] = -y[y < 0]
    vmax = np.percentile(M[M.nonzero()], 95)
    vmin = M.min()
    vmin,vmax = float(vmin),float(vmax)
    vmin,vmax = st.slider("Select value",min_value=vmin,max_value=vmax,value=(vmin,vmax))
    x = coord_compression(x,x_range=xrange)
    y = coord_compression(y,x_range=yrange)
    if orientation == -1:
        x = 1-x
    ax.pcolormesh(x, y, np.flipud(M), cmap=cmap, vmin=vmin, vmax=vmax)


def plot_pair_heat_map(workdir,query,target,chrom,start,end,query_cool,target_cool):



    fig, ax = plt.subplots(figsize=(9, 9))


    x_range = (0,1)
    bound_position = {}
    bound_position[query] = plot_TAD(workdir, ax, query, f"{query}_{chrom}", start, end, 0.2, 0.21, x_range, orientation=1,)



    plot_heatmap(ax, query_cool, chrom, start, end,xrange=(0,1),yrange=(0.22 ,1),up_down=1)


    pair_genome_region = get_pairwise_genome_region(workdir=workdir,
                                                        reference=query,
                                                        chrom=chrom, start=start, end=end)



    chrom2,start2,end2,orientation2 = pair_genome_region[target][0]

    bound_position[target] = plot_TAD(workdir, ax, target, chrom2, start2, end2, -0.21, -0.2, (0,1), orientation=orientation2,)
    plot_heatmap(ax,target_cool, chrom2.split("_")[1], start2, end2,xrange=(0,1),yrange=(-1 ,-0.22),up_down = -1,
                 orientation=orientation2)


    plot_boundary_pair(workdir, ax, query, bound_position, [query,target], style="curve")
    # plot_SV(workdir,query,target,ax,(0,1),(0,1),(0.2 ,0.21),(-0.21 ,-0.2),
    #         genome1_pos = (f"{query}_{chrom}",start,end,1),
    #         genome2_pos = (chrom2,start2,end2,orientation2))
    ax.set_xlim(-0.2, 1)
    ax.set_ylim(-1, 1)
    ax.set_axis_off()
    # st.pyplot(fig)
    return fig


workdir = "/root/run/"
query = "A2"
target = "D5"
chrom = "Chr02"
start = 3.44e+6
end = 6e+6

query_cool = f"/data/{query}.cool"
target_cool = f"/data/{target}.cool"

st.pyplot(plot_pair_heat_map(workdir,query,target,chrom,start,end,query_cool,target_cool))

fig, ax = plt.subplots(figsize=(9, 9))

x_range = (0, 1)
bound_position = {}
bound_position[query] = plot_TAD(workdir, ax, query, f"{query}_{chrom}", start, end, 0.2, 0.21, x_range,
                                 orientation=1, )

fig,ax = plt.subplots(figsize=(9, 9))
plot_heatmap(ax, "/data/matrix_10000.cool", "14", 8200000, 10200000, xrange=(0, 1), yrange=(0.22, 1), up_down=1)
ax.set_xlim(-0.2, 1)
ax.set_ylim(-1, 1)
ax.set_axis_off()
st.pyplot(fig)