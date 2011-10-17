import os
import re

print "\\documentclass{article}"
print "\\begin{document}"
print "\\begin{tabular}{|rrrrr|}"
print "\\hline"
print "$n$ & 2-matching bound & Held-Karp bound & cutting planes & fractional variables \\\\"
print "\\hline"
for n in range(10, 510, 10):
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
    for row in entries:
        for key in row.keys():
            if not (key in sums):
                sums[key] = 0.0
            sums[key] += float(row[key])
    for stat in sums:
        sums[stat] /= len(entries)
    print "%d & %.3f & %.3f & %.3f & %.3f \\\\" % (n, sums["BIMATCHING_BOUND"], sums["HELD_KARP_BOUND"], sums["CUTTING_PLANES"], sums["FRACTIONAL_VARIABLES"])
print "\\hline"
print "\\end{tabular}"
print "\\end{document}"
