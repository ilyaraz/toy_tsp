import math
import random
import sys

def sqr(x):
    return x * x

n = int(sys.argv[1])
print n

x = []
y = []
for i in range(n):
    x.append(random.random())
    y.append(random.random())
for i in range(n):
    for j in range(n):
        print "%.20lf" % math.sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j])),
    print
