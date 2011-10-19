import os
import re
import itertools

data = [["Held-Karp bound", "HELD_KARP_BOUND"], ["Nearest Neighbor", "NEAREST_NEIGHBOR_TOUR"], ["Greedy", "GREEDY_TOUR"], ["Greedy + 2-opt", "GREEDY_2OPT_TOUR"], ["Greedy + best-2-opt", "GREEDY_BEST2OPT_TOUR"]]

print "\\documentclass{article}"
print "\\begin{document}"
print "\\begin{tabular}{|" + "".join(itertools.repeat('r', len(data) + 1)) + "|}"
print "\\hline"
print "$n$",
for item in data:
    print "&", item[0],
print "\\\\"
print "\\hline"
for n in range(10, 210, 10):
    file = "../logs/log" + str(n) + ".txt"
    current_item = dict()
    entries = []
    for line in open(file):
        if "----------" in line:
            if len(current_item):
                entries.append(current_item) 
            current_item = dict()
        else:
            tokens = line.split()
            current_item[tokens[0]] = tokens[1]
    if len(current_item):
        entries.append(current_item)
    sums = dict()
    num_entries = 0
    for row in entries:
        num_entries += 1
        for key in row.keys():
            if not (key in sums):
                sums[key] = 0.0
            sums[key] += float(row[key])
    for stat in sums:
        sums[stat] /= num_entries
    print "%d" % n,
    for item in data:
        print "& %.3f" % sums[item[1]],
    print "\\\\"
print "\\hline"
print "\\end{tabular}"
print "\\end{document}"
