import argparse

import matplotlib.pyplot as plt

INPUT_DIR = "input"
RESULTS_DIR = "output"
OUTPUT_DIR = "plots"

RESULTS_FILE_PREFIXES = ("min", "max", "median")
ALGORITHMS = ("cycle", "regret", "neighbour")
INSTANCES = ("kroa100", "krob100")

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--instance",
    help="Graph instance filename (kroa100, krob100). "
    "If not set, the script will show results for both.",
)
parser.add_argument(
    "-p",
    "--prefix",
    help="File prefix (min, max, median). If not set, the script will show results for all.",
)
parser.add_argument(
    "-a",
    "--algo",
    help="Algorithm to show (cycle, regret, neighbour). "
    "If not set, the script will show results for all.",
)
parser.add_argument("-n", "--cycles-number", type=int, help="Number of separable cycles.")

args = parser.parse_args()

if __name__ == "__main__":
    algos_to_show = ALGORITHMS if args.algo is None else args.algo
    prefixes_to_show = RESULTS_FILE_PREFIXES if args.prefix is None else args.prefix

    file_prefixes = [f"{prefix}_{algo}" for prefix in prefixes_to_show for algo in algos_to_show]

    for instance in INSTANCES if args.instance is None else args.instance:
        x, y = [], []
        cycles = []
        with open(f"{INPUT_DIR}/{instance}.tsp") as file:
            for line in file:
                line = line.strip().split()
                if line[0].isnumeric():
                    x.append(int(line[1]))
                    y.append(int(line[2]))

        for file_prefix in file_prefixes:
            with open(f"{RESULTS_DIR}/{file_prefix}_{instance}.tsp") as file:
                vertexes = [int(vertex) - 1 for vertex in file]
                vertexes_in_one_cycle = len(vertexes) // args.cycles_number

                cycles = [
                    vertexes[
                        i * vertexes_in_one_cycle : i * vertexes_in_one_cycle
                        + vertexes_in_one_cycle
                    ]
                    for i in range(args.cycles_number)
                ]

            for cycle in cycles:
                plt.plot(
                    [x[c] for c in cycle] + [x[cycle[0]]],
                    [y[c] for c in cycle] + [y[cycle[0]]],
                )
                plt.scatter([x[c] for c in cycle], [y[c] for c in cycle])

            plt.title(f"{file_prefix}_{instance}")
            plt.savefig(f"{OUTPUT_DIR}/{file_prefix}_{instance}")
            plt.close()
