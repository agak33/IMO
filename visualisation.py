import argparse
import glob
import os
from math import ceil

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file_dir", help="Directory with vertices coords", required=True)
parser.add_argument("-i", "--input-dir", help="Input directory", required=True)
parser.add_argument("-o", "--output-dir", help="Output directory", required=True)
parser.add_argument(
    "-n", "--cycles-number", type=int, help="Number of separable cycles.", default=2
)

INSTANCES = ("kroa100.tsp", "krob100.tsp")

args = parser.parse_args()

if __name__ == "__main__":
    for input_file in glob.glob(f"{args.input_dir}/*"):
        output_filename = os.path.split(input_file)[-1]
        instance_filename = INSTANCES[0] if INSTANCES[0] in output_filename else INSTANCES[1]

        x, y = [], []
        cycles = []
        with open(f"{args.file_dir}/{instance_filename}") as file:
            for line in file:
                line = line.strip().split()
                if line[0].isnumeric():
                    x.append(int(line[1]))
                    y.append(int(line[2]))

        with open(input_file) as file:
            vertexes = [int(vertex) - 1 for vertex in file]
            vertexes_in_one_cycle = len(vertexes) // args.cycles_number

            cycles = [
                vertexes[
                    i * vertexes_in_one_cycle : i * vertexes_in_one_cycle + vertexes_in_one_cycle
                ]
                for i in range(args.cycles_number)
            ]

        for cycle in cycles:
            plt.plot(
                [x[c] for c in cycle] + [x[cycle[0]]],
                [y[c] for c in cycle] + [y[cycle[0]]],
            )
            plt.scatter([x[c] for c in cycle], [y[c] for c in cycle])

        plt.title(f"{instance_filename}")
        plt.savefig(f"{args.output_dir}/{output_filename}.png")
        plt.close()
