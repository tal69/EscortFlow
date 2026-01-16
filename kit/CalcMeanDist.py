"""
Calculate the mean distance to the closest output cell - VERY SIMPLE
"""

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
parser.add_argument("-O", "--output_cells", nargs='+', type=int, help="List of output locations", required=True)

args = parser.parse_args()

Lx = args.Lx
Ly = args.Ly

ez = [int(a) for a in args.output_cells]
if len(ez) % 2:
    print(f"Warning: output cell coordinate list must be of even length - {ez}")
    exit(1)

O = []
for i in range(len(ez) // 2):
    O.append((ez[i * 2], ez[i * 2 + 1]))

total = 0

for x in range(Lx):
    for y in range(Ly):
        min_val = Lx + Ly
        for o in O:
            min_val = min(min_val, abs(o[0] - x) + abs(o[1] - y))
        total += min_val

print(f"Mean distance to output cell: {total / (Lx * Ly)}")