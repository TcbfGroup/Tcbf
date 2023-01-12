import glob
import os.path
from itertools import combinations
from tcbf.run_command import run_command
from collections import defaultdict
from pandas import read_csv, concat, DataFrame, read_table, options


def get_species(workdir):
    step1 = os.path.join(workdir, "Step1")
    species = [i.split(".")[0] for i in os.listdir(step1) if i.endswith(".fa")]
    return species


def sub_network_construct(workdir, need_syn):
    step2 = os.path.join(workdir, "Step2")
    species = get_species(workdir)
    merge_file = os.path.join(workdir, "Step3", "merge.network.txt")
    if os.path.exists(merge_file):
        os.unlink(merge_file)
    network_file = os.path.join(workdir, "Step3", "out.clean.network.txt")
    file_template = "{}_{}.block.txt" if need_syn else "{}_{}.network.bed"
    for s1, s2 in combinations(species, 2):
        file1 = os.path.join(step2, file_template.format(s1, s2))
        file2 = os.path.join(step2, file_template.format(s2, s1))
        max_score = get_max_score(file1, file2)
        if max_score is not None:
            max_score.to_csv(merge_file, header=False, index=False, sep="\t", mode="a")

    command = f"mcl {merge_file} --abc -o {network_file}"
    run_command(command)


def get_total_TAD_boundary(workdir):
    rundir = os.path.join(workdir, "Step1")
    lst = glob.glob(f"{rundir}/*bed")
    total = concat([read_csv(i) for i in lst])["tad_name"]
    return set(total)


def extract_ortho_group(abspath):
    total_boundary = get_total_TAD_boundary(abspath)

    network_file = os.path.join(abspath, "Step3", "out.clean.network.txt")
    Un_Assign = os.path.join(abspath, "Result", "Unassignd_TAD.tsv")
    Group_file = os.path.join(abspath, "Result", "TAD_groups.tsv")
    TAD_count = os.path.join(abspath, "Result", "TAD_groups_count.tsv")

    dic = {}
    unassign = defaultdict(lambda: {})
    unassign_num = 1

    in_ortho_boundary = set()

    with open(network_file) as f:
        for index, line in enumerate(f, 1):
            if len(line.split()) > 1:
                group = f"Group_{index}"
                for k in line.split():
                    species = k.split("_")[0]
                    dic.setdefault(group, {}).setdefault(species, []).append(k)
                    in_ortho_boundary.add(k)

            else:
                species = line.strip().split("_")[0]
                unassign[species][unassign_num] = line.strip()
                in_ortho_boundary.add(line.strip())
                unassign_num += 1

    unmapped = total_boundary - in_ortho_boundary
    for i in unmapped:
        species = i.split("_")[0]
        unassign[species][unassign_num] = i
        unassign_num += 1
    data = DataFrame(dic).T.fillna("")
    count = data.apply(lambda x: x.str.len())
    data = data.applymap(lambda x: ";".join(x) if x else "")
    data.to_csv(Group_file, sep="\t")
    unassign_df = DataFrame(unassign)
    unassign_df.to_csv(Un_Assign, sep="\t", index=False)
    count.to_csv(TAD_count, sep="\t")

    #
    # one_to_many = os.path.join(abspath, "Result", "one_to_many")
    # if not os.path.exists(one_to_many):
    #     os.mkdir(one_to_many)
    # all_species = get_species(abspath)
    #
    # merge = concat([data, unassign_df])
    # for species in all_species:
    #     boundary_file = os.path.join(abspath, "Step1", f"{species}.bound.bed")
    #     reference_TAD_boundary_order = read_csv(boundary_file)["tad_name"]
    #     tmp = merge.dropna(subset=[species])
    #     rr = []
    #     for i in reference_TAD_boundary_order:
    #         rr.append(tmp[tmp[species].str.contains(f"{i}(;|$)")])
    #
    #     result_file = os.path.join(one_to_many, f"{species}.txt")
    #     res = concat(rr)
    #     res.pop(species)
    #     res.insert(0, species, list(reference_TAD_boundary_order))
    #     res.to_csv(result_file, index=False, sep="\t")


def get_max_score(network1, network2):
    n1 = read_table(network1)
    n2 = read_table(network2)
    if (n1.shape[0] == 0) or (n2.shape[0] == 0):
        return
    n1["combine"] = n1["seq_id"].str.cat(n1["tad_name"], sep="-")
    n2["combine"] = n2["tad_name"].str.cat(n2["seq_id"], sep="-")
    n1 = n1.iloc[:, 2:]
    n2 = n2.iloc[:, 2:]
    n = n1.merge(n2, how="outer", on="combine")
    max_score = n.iloc[:, [0, 2]].max(axis=1)
    data = n["combine"].str.split("-", expand=True)
    data["max_score"] = max_score
    data.columns = "genome1 genome2 score".split()
    return data

def construct_one_to_many(workdir,need_syn):
    step2 = os.path.join(workdir, "Step2")
    species = get_species(workdir)
    file_template = "{}_{}.block.txt" if need_syn else "{}_{}.network.bed"
    step3 = os.path.join(workdir,"Step3")
    for s1 in species:

        for s2 in species:
            if s1 == s2:
                continue
            bound_file = os.path.join(workdir,"Step1",f"{s1}.bound.bed")
            file1 = os.path.join(step2, file_template.format(s1, s2))
            file2 = os.path.join(step2, file_template.format(s2, s1))
            max_score = get_max_score(file1, file2)
            boundaries = read_csv(bound_file)[["tad_name"]]
            boundaries.columns = ["genome1"]
            if max_score is not None:
                max_score = max_score.iloc[:,:2]
                m = boundaries.merge(max_score, how="left").fillna(0)
                s = defaultdict(list)
                for k in m.itertuples():
                    key,value = k[1:]
                    s[key].append(value)
                with open(os.path.join(step3,f"{s1}-{s2}.txt"),"w")as f:
                    for k,v in s.items():
                        if not v[0]:
                            f.write(f"{k}\t\n")
                        else:
                            f.write(f"{k}\t{';'.join(v)}\n")
        anchor_file = [os.path.join(step3,i) for i in os.listdir(step3) if i.startswith(f"{s1}-")]
        first = read_table(anchor_file[0],names=os.path.basename(anchor_file[0])[:-4].split("-"))
        if len(anchor_file) >1:
            result = concat([read_table(i,names = i[:-4].split("-"))[[i[:-4].split("-")[1]]] for i in anchor_file[1:]],axis = 1)
            final = concat([first, result], axis=1)
        else:
            final = first

        final.to_csv(os.path.join(step3,f"{s1}.join.txt"),index=False,sep = "\t")







def network_construct(workdir, need_syn):
    options.mode.chained_assignment = None
    abs_path = os.path.abspath(workdir)
    run_dir = os.path.join(abs_path, "Step3")

    if not os.path.exists(run_dir):
        os.mkdir(run_dir)
    Result = os.path.join(abs_path, "Result")
    if not os.path.exists(Result):
        os.mkdir(Result)
    sub_network_construct(abs_path, need_syn)

    extract_ortho_group(abs_path)
    construct_one_to_many(workdir,need_syn)
