#!/usr/bin/python3

import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vertices', dest='vertices', help='number of vertices', required=True, type=int)
parser.add_argument('--max_weight', dest='max_weight', help='maximum weight, defaults to 10', default=10, type=int)
parser.add_argument('--sparse', dest='sparse', help='sparseness of the graph, defaults to 0.5', type=float, default=0.5)

args = parser.parse_args()
N = args.vertices
MAXWEIGHT = args.max_weight
STEP = int(args.sparse * N)
if STEP == 0:
	STEP = 1

n_lines = 0
lines = []
for i in range(1, N+1):
	for j in range(1, N+1, random.randint(1, STEP)):
		if i == j:
			continue
		weight = random.randrange(1, MAXWEIGHT)
		lines.append("{}-{}-{}".format(i, weight, j))
		n_lines += 1

print(N)
print(n_lines)

for line in lines:
	print(line)
