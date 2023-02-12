import glob
import os.path

import pandas as pd
files = glob.glob("*join.txt")
species = [i[:-9] for i in files]
if not os.path.exists("result"):
    os.mkdir("result")

def shift(data):
    # data["end"] = data["start"].shift(-1)
    # data = data.dropna().iloc[:,1:]
    data = data.iloc[:,1:]
    data.columns = "Chr Start End Value".split()
    data["End"] = data["End"].astype(int)
    return data

max_number = len(species)

for index,s in enumerate(species):

    join = pd.read_table(files[index]).apply(lambda x:len(x.dropna()),axis = 1)
    bound = pd.read_csv(f"../Step1/{s}.bound.bed")
    bound["chromosome"] = bound["chromosome"].str.split("_").str.get(-1)
    bound["present_number"] = join
    karyotype = bound.sort_values(["chromosome","end"]).groupby("chromosome").tail(1).loc[:, "chromosome start end".split()]
    karyotype.columns = "Chr Start End".split()
    karyotype.to_csv(f"result/{s}.karyotype.csv",index=False)

    density = bound.groupby("chromosome").apply(shift)
    density.to_csv(f"result/{s}.density.csv",index=False)

    t = bound[bound["present_number"].isin([1, max_number])]
    t["Type"] = ""
    t["Shape"] = ""
    t["color"] = ""
    t.loc[t["present_number"] == 1,"Type"] = "unique"
    t.loc[t["present_number"] != 1, "Type"] = "conserved"
    t.loc[t["present_number"] == 1, "Shape"] = "box"
    t.loc[t["present_number"] != 1, "Shape"] = "triangle"

    t.loc[t["present_number"] == 1, "color"] = "6a3d9a"
    t.loc[t["present_number"] != 1, "color"] = "ff7f00"

    t = t.loc[:,"Type Shape chromosome start end color".split()]
    t.columns = "Type Shape Chr Start End color".split()
    t.to_csv(f"result/{s}.marker.csv",index=False)


