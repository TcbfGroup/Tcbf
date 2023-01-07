import os.path
import sys
from tempfile import NamedTemporaryFile
import pandas as pd
from tcbf.run_command import run_command
from tcbf.aligner_paramters import minimap2_align


def process_tad_paris(file):
    pd.options.mode.chained_assignment = None
    df = pd.read_feather(file)
    df["bound"] = df["bound"].astype(int)
    s1 = df["seq_id"].str.split("_|-", expand=True)
    s1 = s1.replace({"middle": None})
    s1[1] = s1[1].astype(int)
    direction_convert = {"left": 0,
                         "right": 1}

    def determine_status(tad_index1, direction, islast_first):
        if islast_first is not None:
            result = ((tad_index1, direction_convert[direction]),)
        else:
            result = (tad_index1, direction_convert[direction]), (tad_index1 - 1, 1)
        return result

    def flatten(lst):
        result = []
        for i in lst:
            for j in i:
                result.append(j)
        return result

    df["align2"] = s1.iloc[:, 1:].apply(lambda x: determine_status(*x), axis=1)

    align_pool = {}
    for line in df.drop_duplicates().groupby(["tad_name", "bound"]):
        status1 = line[0]
        align_pool.setdefault(status1[0], [[] for _ in range(2)])
        status2_lst = flatten(line[1]["align2"])
        align_pool[status1[0]][status1[1]].extend(status2_lst)

    final_lst = {}
    for k, v in align_pool.items():
        left_align = v[0]
        right_align = v[1]
        # 必须左右两边的TAD边界都要有比对结果
        if not (left_align and right_align):
            continue
        d1 = pd.DataFrame(left_align, columns="left_hit left_hit_direction".split())
        d2 = pd.DataFrame(right_align, columns="right_hit right_hit_direction".split())

        # conserved
        # 获取左右各自比对到的其他所有TAD
        left_hit = set(i[0] for i in left_align)
        right_hit = set(i[0] for i in right_align)
        # 检查是否有重合的
        if len(left_hit & right_hit) != 0:

            merge = d1.merge(d2, left_on="left_hit", right_on="right_hit")
            # 如果两边都有能比对上的，说明是保守的
            final = merge.query("left_hit_direction != right_hit_direction")
            if final.shape[0] != 0:
                for aligned in final["left_hit"].unique():
                    final_lst[(k, aligned)] = 1
        # splited
        # left-left right-right
        l1 = d1.query("left_hit_direction == 0")
        r1 = d2.query("right_hit_direction == 1")
        l1["left_hit"] += 1
        merge = l1.merge(r1, left_on="left_hit", right_on="right_hit")
        if merge.shape[0] >= 1:
            for aligned in merge.itertuples():
                final_lst[(k, aligned[1] - 1, aligned[3])] = 0.5
        else:
            l1["left_hit"] += 1
            merge = l1.merge(r1, left_on="left_hit", right_on="right_hit")
            if merge.shape[0] >= 1:
                for aligned in merge.itertuples():
                    final_lst[(k, aligned[1] - 1, aligned[3])] = 0.5

        # split
        # right-left left-right
        l1 = d1.query("left_hit_direction == 1")
        r1 = d2.query("right_hit_direction == 0")
        l1["left_hit"] -= 1
        merge = l1.merge(r1, left_on="left_hit", right_on="right_hit")
        if merge.shape[0] >= 1:
            for aligned in merge.itertuples():
                final_lst[(k, aligned[1] + 1, aligned[3])] = 0.5
        else:
            l1["left_hit"] -= 1
            merge = l1.merge(r1, left_on="left_hit", right_on="right_hit")
            if merge.shape[0] >= 1:
                for aligned in merge.itertuples():
                    final_lst[(k, aligned[1] + 1, aligned[3])] = 0.5

    datas = []

    query_genome = df["seq_id"][0].split("_")[0]

    for k, v in final_lst.items():
        for line in k[1:]:
            row = (k[0], line, v)

            datas.append(row)

    d = pd.DataFrame(datas)

    d[1] = query_genome + "_" + d[1].astype(str)
    return d


# def process_nucmer_coord(coord_file, out):
#     from pandas import read_table
#     need = "seq_id start end seq_id2 start2 end2".split()
#     df = read_table(coord_file, header=None,
#                     usecols=[0, 1, 2, 3, 7, 8]).loc[:, [8, 2, 3, 7, 0, 1]]
#     df.columns = need
#     feather.write_feather(df, out)
#
#
# def nucmer_align(bound_fasta, reference, output,
#                  threads=4):
#     with NamedTemporaryFile("w+t") as delta:
#         with NamedTemporaryFile("w+t") as filter_delta:
#             with NamedTemporaryFile("w+t") as coord:
#                 command = f"nucmer {reference} {bound_fasta} --delta={delta.name} -t {threads} &&" \
#                           f"delta-filter -m -i 90 -l 1000 {delta.name} > {filter_delta.name} &&" \
#                           f"show-coords -THq {filter_delta.name} > {coord.name} "
#                 run_command(command)
#                 process_nucmer_coord(coord.name, output)


def align_genome(query, target, workdir, threads, aligner, minioverlap,
                 xdist=10,ydist=10,N = 3):
    if threads <= 0:
        from multiprocessing import cpu_count
        threads = cpu_count()

    result_dir = os.path.join(workdir, "Step2")
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)


    query_bound = os.path.join(workdir, "Step1", f"{query}.bound.fasta")
    target_genome = os.path.join(workdir, "Step1", f"{target}.genome.fa")

    bound_bed = os.path.join(workdir, "Step1", f"{target}.bound.bed")


    network_out = os.path.join(result_dir, f"{query}_{target}.network.bed")


    with NamedTemporaryFile("w+t") as Collinearity:
        if aligner == "nucmer":
            sys.exit("Nucmer is not supported temporarily, and will be added in the future")
            # nucmer_align(query_bound, target_genome, Collinearity.name, threads=threads)
        elif aligner == "minimap2":
            minimap2_align(query_bound,
                           target_genome,
                           Collinearity.name,
                           threads=threads
                           )
        else:
            raise TypeError(f"Wrong aligner for {aligner}!!!")

        command = f"tcbf_syn_process {Collinearity.name}  " \
                  f" {bound_bed}  {network_out} {minioverlap}"
        run_command(command)
    from tcbf.construct_synteny_block import construct_block
    construct_block(workdir,query,target,xdist,ydist,N)


