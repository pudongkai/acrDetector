# coding=utf-8
import os


def extractHTH(faa):
    hthFile = faa.replace(".faa", ".hth")
    hmmscan = "feature/lib/hmmscan"
    hmmdb = 'feature/utils/pfarmHTHdatabase/HTH.hmm'
    if not os.path.isfile(hthFile):
        print(hmmdb)
        os.system('%s --tblout %s %s %s' % (hmmscan, hthFile, hmmdb, faa))

    hth_proteins = set()
    with open(hthFile, 'r') as f:
        for line in f:
            if line[0] != '#':
                temp = line.split()
                hth_proteins.add(temp[2].split('|')[0])
    return hth_proteins


def caculateHTHgap(protein_id_list, hth_list, target):

    # 获取目标蛋白下游的蛋白
    downstream = getOrderList(protein_id_list, target)

    # 计算目标蛋白与Aca之间的间隔
    gap = 1
    for protein in downstream:
        if protein in hth_list:
            return gap
        gap += 1

    return -999


# 获得目标蛋白上下游的n个蛋白
# 结果list顺序从目标蛋白开始
def getOrderList(protein_id_list, target):
    n = 300
    total = len(protein_id_list)
    index = protein_id_list.index(target)
    start = index + 1
    end = start + n
    if end > total:
        downstream = protein_id_list[start:] + protein_id_list[: n-(total-start)+1]
    else:
        downstream = protein_id_list[start: end]

    return downstream
