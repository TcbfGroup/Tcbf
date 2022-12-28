import os.path
from pandas import read_csv
from collections import defaultdict
from tcbf.Grouper import Grouper


class AlignmentLine(object):
    def __init__(self, line):
        arg = line.split()
        self.query = arg[0]
        self.target = arg[1]
        self.hitLen = int(arg[2])
        self.query_index = None
        self.target_index = None

    def __repr__(self):
        return f"AlignmentLine {self.query} to {self.target}, Length = {self.hitLen}"

    def qi(self):
        return self.query_index

    def ti(self):
        return self.target_index


class Alignment(object):
    def __init__(self, filename):
        self.file = filename

    def __iter__(self):
        with open(self.file) as fp:
            fp.readline()
            for row in fp:
                yield AlignmentLine(row)


def get_boundary_order(workdir, genome):
    file = os.path.join(workdir, "Step1", f"{genome}.bound.bed")
    data = read_csv(file)
    data["index_number"] = range(data.shape[0])

    result = dict(zip(data["tad_name"], list(zip(data["index_number"], data["chromosome"]))))

    return result


def read_alignments(workdir, genome1, genome2):
    result = []
    seen = set()
    file = os.path.join(workdir, "Step2", f"{genome1}_{genome2}.network.bed")

    qorder = get_boundary_order(workdir, genome1)
    sorder = get_boundary_order(workdir, genome2)
    b1 = Alignment(file)
    for b in b1:
        query, target = b.query, b.target
        query_index, q = qorder.get(query)
        target_index, s = sorder.get(target)
        key = query, target
        if key in seen:
            continue
        seen.add(key)

        b.qseqid = q
        b.tseqid = s

        b.qi, b.ti = query_index, target_index

        result.append(b)
    return result


def group_hits(alignments):
    all_hits = defaultdict(list)
    for b in alignments:
        all_hits[(b.qseqid, b.tseqid)].append((b.qi, b.ti, b.hitLen))
    return all_hits


def _score(cluster):
    """
    score of the cluster, in this case, is the number of non-repetitive matches
    """
    x, y = list(zip(*cluster))[:2]
    return min(len(set(x)), len(set(y)))


def synteny_scan(points, xdist, ydist, N):
    cluster = Grouper()
    n = len(points)
    points.sort()
    for i in range(n):
        for j in range(i - 1, -1, -1):
            del_x = points[i][0] - points[j][0]
            if del_x > xdist:
                break
            del_y = points[i][1] - points[j][1]
            if abs(del_y) > ydist:
                continue
            cluster.join(points[i], points[j])
    clusters = [sorted(cluster) for cluster in list(cluster) if _score(cluster) >= N]
    return clusters


def batch_scan(alignment, xdist, ydist, N):
    chr_pair_points = group_hits(alignment)
    clusters = []
    for chr_pair in sorted(chr_pair_points.keys()):
        points = chr_pair_points.get(chr_pair)
        clusters.extend(
            synteny_scan(points, xdist, ydist, N))
    return clusters


def construct_block(workdir,query,target,xdist,ydist,N):
    query_number_boundary = {v[0]:k for k,v in get_boundary_order(workdir,query).items()}
    target_number_boundary = {v[0]:k for k,v in get_boundary_order(workdir,target).items()}

    alignments = read_alignments(workdir,query,target)
    cluster = batch_scan(alignments,xdist=xdist,ydist = ydist,N = N)

    cluster = [i for i in cluster if max(j[2] for j in i) >= 3000]
    result_file = os.path.join(workdir,"Step2",f"{query}_{target}.block.txt")
    with open(result_file,"w")as f:
        f.write("seq_id\ttad_name\tscore\n")
        for c in cluster:

            for query_index,target_index,hitLen in c:
                query_boundary = query_number_boundary.get(query_index)
                target_boundary = target_number_boundary.get(target_index)

                f.write("\t".join((query_boundary,target_boundary,str(hitLen))) + "\n")

