from ete3 import Tree
from sys import argv, stderr

tree_file = argv[1]
lineages_file = argv[2]

lineages = {}
with open(lineages_file, 'r') as f:
    next(f) # With header
    for line in f:
        ID, Host = line.strip().split('\t')
        if Host not in lineages.keys():
            lineages[Host] = [ID]
        else:
            lineages[Host] = lineages[Host] + [ID]

monophyl_out = {l:[] for l in lineages.keys()}
with open(tree_file, 'r') as f:
    for line in f:
        tre = Tree(line)
        mid = tre.get_midpoint_outgroup()
        tre.set_outgroup(mid)
        lineages_tree = {}
        for ID in tre.get_leaf_names():
            Host = list(lineages.keys())[[ID in lineages[l] for l in lineages].index(True)]
            if Host not in lineages_tree.keys():
                lineages_tree[Host] = [ID]
            else:
                lineages_tree[Host] = lineages_tree[Host] + [ID]
        for lineage in lineages_tree.keys():
            common_anc = len(tre.get_common_ancestor(lineages_tree[lineage]).get_leaves())
            if common_anc == len(lineages_tree[lineage]):
                monophyl_out[lineage].append(1)
            else:
                monophyl_out[lineage].append(0)

for l in monophyl_out:
    P = sum(monophyl_out[l]) / len(monophyl_out[l])
    print(l, P, sep = '\t', file = stderr)


for l in monophyl_out:
    print(l, ''.join([str(i) for i in monophyl_out[l]]), sep = '\t')

out = []
for i in range(len(monophyl_out[list(monophyl_out.keys())[0]])):
    shape = [monophyl_out[l][i] for l in list(monophyl_out.keys())]
    out.append(''.join([str(x) for x in shape]))

print(list(lineages.keys()))
unique_shapes = set(out)
print(unique_shapes)
print([out.count(i) for i in unique_shapes])
