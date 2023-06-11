import os

import matplotlib.pyplot as plt
import numpy as np

DIR = "output_tests"
INSTANCES = ["kroa200.tsp", "krob200.tsp", "kroa100.tsp", "krob100.tsp"]
BEST_SCORES = [34259, 35539, 24480, 25268]
BEST_SCORES_NRS = [300, 137, 929, 84]


def similarity_vertexes(cycles1: list[set[int]], cycles2: list[set[int]]) -> int:
    common_vertex_count = [
        # pierwsza możliwość (1szy z 1szym, 2gi z 2gim)
        len(cycles1[0].intersection(cycles2[0])),
        len(cycles1[1].intersection(cycles2[1])),
        # druga możliwość (1szy z 2gim, 2gi z 1szym)
        len(cycles1[0].intersection(cycles2[1])),
        len(cycles1[1].intersection(cycles2[0])),
    ]

    argmax_count = np.argmax(common_vertex_count)
    if argmax_count % 2:
        return common_vertex_count[argmax_count] + common_vertex_count[argmax_count - 1]

    return common_vertex_count[argmax_count] + common_vertex_count[argmax_count + 1]


def similarity_edges(solution1: np.ndarray, solution2: np.ndarray) -> int:
    similarity = (solution1.astype(int) & solution2.astype(int)).sum()

    if similarity % 2:
        raise ValueError("XD")

    return similarity // 2


def get_neighbourhood_matrix(vertexes: list[int]) -> np.ndarray:
    neighbourhood_matrix = np.zeros(shape=(len(vertexes), len(vertexes)))

    for i in range(len(vertexes)):
        v1, v2 = vertexes[i], vertexes[(i + 1) % len(vertexes)]

        neighbourhood_matrix[v1][v2] = 1
        neighbourhood_matrix[v2][v1] = 1

    return neighbourhood_matrix


if __name__ == "__main__":
    for instance, best_score, best_score_nr in zip(INSTANCES, BEST_SCORES, BEST_SCORES_NRS):
        # wczytanie wierzchołków najlepszej instancji
        with open(os.path.join(DIR, f"{instance}_{best_score_nr}"), "r") as file:
            best_vertexes = [int(vertex) - 1 for vertex in file]
            best_cycles = [
                set(best_vertexes[: len(best_vertexes) // 2]),
                set(best_vertexes[len(best_vertexes) // 2 :]),
            ]

        best_nm = get_neighbourhood_matrix(best_vertexes)

        # wczytanie wyników wszystkich instancji
        with open(os.path.join(DIR, f"{instance}_results"), "r") as file:
            all_scores = [int(score) for x, score in enumerate(file)]

        s_v = []
        s_e = []

        all_cycles = []
        all_nm = []
        for i in range(1, 1001):
            with open(os.path.join(DIR, f"{instance}_{i}"), "r") as file:
                vertexes = [int(vertex) - 1 for vertex in file]
                cycles = [
                    set(vertexes[: len(vertexes) // 2]),
                    set(vertexes[len(vertexes) // 2 :]),
                ]

                nm = get_neighbourhood_matrix(vertexes)

                all_cycles.append(cycles)
                all_nm.append(nm)

                sv_val = similarity_vertexes(best_cycles, cycles)
                se_val = similarity_edges(best_nm, nm)

                if i == best_score_nr:
                    try:
                        assert sv_val == len(best_vertexes)
                        assert se_val == len(best_vertexes)
                    except AssertionError:
                        print(sv_val, se_val)
                        raise
                else:
                    s_v.append(sv_val)
                    s_e.append(se_val)

        all_scores_reduced = all_scores[: best_score_nr - 1] + all_scores[best_score_nr:]
        plt.title(
            f"{instance}, liczba wspólnych krawędzi\nwspółczynnik korelacji: {round(np.corrcoef(all_scores_reduced, s_e)[0, 1], 3)}"
        )
        plt.scatter(all_scores_reduced, s_e)
        plt.xlabel("Wartość funkcji celu")
        plt.ylabel(f"Podobieństwo do najlepszego rozwiązania ({best_score})")
        plt.savefig(f"{instance}, edges.png")
        plt.close()

        plt.title(
            f"{instance}, liczba wspólnych wierzchołków\nwspółczynnik korelacji: {round(np.corrcoef(all_scores_reduced, s_v)[0, 1], 3)}"
        )
        plt.scatter(all_scores_reduced, s_v)
        plt.xlabel("Wartość funkcji celu")
        plt.ylabel(f"Podobieństwo do najlepszego rozwiązania ({best_score})")
        plt.savefig(f"{instance}, vertexes.png")
        plt.close()

        s_v_mean = [
            np.mean(
                [
                    similarity_vertexes(cycles1, cycles2)
                    for x, cycles1 in enumerate(all_cycles)
                    if x != y
                ]
            )
            for y, cycles2 in enumerate(all_cycles)
        ]
        s_e_mean = [
            np.mean([similarity_edges(nm1, nm2) for x, nm1 in enumerate(all_nm) if x != y])
            for y, nm2 in enumerate(all_nm)
        ]

        plt.title(
            f"{instance}, liczba wspólnych krawędzi\nwspółczynnik korelacji: {round(np.corrcoef(all_scores, s_e_mean)[0, 1], 3)}"
        )
        plt.scatter(all_scores, s_e_mean)
        plt.xlabel("Wartość funkcji celu")
        plt.ylabel("Średnie podobieństwo")
        plt.savefig(f"{instance}, edges, mean.png")
        plt.close()

        plt.title(
            f"{instance}, liczba wspólnych wierzchołków\nwspółczynnik korelacji: {round(np.corrcoef(all_scores, s_v_mean)[0, 1], 3)}"
        )
        plt.xlabel("Wartość funkcji celu")
        plt.ylabel("Średnie podobieństwo")
        plt.scatter(all_scores, s_v_mean)
        plt.savefig(f"{instance}, vertexes, mean.png")
        plt.close()
