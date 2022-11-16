import os.path
from itertools import combinations
from tcbf.run_command import run_command
import pandas as pd
from collections import defaultdict




def get_species(workdir):
    step1 = os.path.join(workdir, "Step1")
    species = [i.split(".")[0] for i in os.listdir(step1) if i.endswith(".fa")]
    return species


def sub_network_construct(workdir):
    step2 = os.path.join(workdir, "Step2")
    species = get_species(workdir)
    merge_file = os.path.join(workdir,"Step3","merge.network.txt")
    if os.path.exists(merge_file):
        os.unlink(merge_file)
    network_file = os.path.join(workdir, "Step3", "out.clean.network.txt")
    for s1,s2 in combinations(species,2):
        print(f"processing {s1} {s2}")
        get_max_score(os.path.join(step2,f"{s1}_{s2}.network.bed"),
                      os.path.join(step2, f"{s2}_{s1}.network.bed"),
                      ).to_csv(merge_file,header = False,index=False,sep = "\t",mode = "a")

    command = f"mcl {merge_file} --abc -o {network_file}"
    run_command(command)


def extract_ortho_group(abspath):
    network_file = os.path.join(abspath, "Step3", "out.clean.network.txt")
    Un_Assign = os.path.join(abspath,"Result","Unassignd_TAD.tsv")
    Group_file = os.path.join(abspath,"Result","TAD_groups.tsv")
    TAD_count = os.path.join(abspath,"Result","TAD_groups_count.tsv")



    dic = {}
    unassign = defaultdict(lambda :{})
    unassign_num = 1
    with open(network_file)as f:
        for index, line in enumerate(f,1):
            if len(line.split()) > 1:
                group = f"Group_{index}"
                for k in line.split():
                    species = k.split("_")[0]
                    dic.setdefault(group,{}).setdefault(species,[]).append(k)

            else:
                species = line.strip().split("_")[0]
                unassign[species][unassign_num] = line.strip()
                unassign_num += 1
    data = pd.DataFrame(dic).T.fillna("")
    count = data.apply(lambda x:x.str.len())
    data = data.applymap(lambda x:";".join(x) if x else "")
    data.to_csv(Group_file,sep = "\t")
    pd.DataFrame(unassign).to_csv(Un_Assign,sep="\t")
    count.to_csv(TAD_count,sep = "\t")






def get_max_score(network1, network2):
    n1 = pd.read_table(network1)
    n2 = pd.read_table(network2)
    n1["combine"] = n1["seq_id"].str.cat(n1["tad_name"], sep="-")
    n2["combine"] = n2["tad_name"].str.cat(n2["seq_id"], sep="-")
    n1 = n1.iloc[:, 2:]
    n2 = n2.iloc[:, 2:]
    n = n1.merge(n2, how="outer", on="combine")
    max_score = n.max(axis=1)
    data = n["combine"].str.split("-", expand=True)
    data["max_score"] = max_score
    data.columns = "genome1 genome2 score".split()
    return data



def network_construct(workdir):
    pd.options.mode.chained_assignment = None
    abs_path = os.path.abspath(workdir)
    run_dir = os.path.join(abs_path, "Step3")

    if not os.path.exists(run_dir):
        os.mkdir(run_dir)

    sub_network_construct(abs_path)

    Result = os.path.join(abs_path,"Result")
    if not os.path.join(Result):
        os.mkdir(Result)

    extract_ortho_group(abs_path)




