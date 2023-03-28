#include "utils.hpp"

string INPUT_DIR = "../greedy_tsp/input";
string INSTANCES[] = { "kroa100.tsp", "krob100.tsp" };

string INPUT_DIR_CYCLES = "../greedy_tsp/output";
string RESULT_PREFIX = "min";
string ALGORITHMS[] = { "random", "regret_no_weight" };


void swap_vertexes(vector<vector<int>> cycles) {
    /*
    A function to swap vertices between cycles.
    */

}

void swap_vertexes(vector<int> cycle) {
    /*
    A function to swap vertices inside the cycle.
    */
}

void find_best_move(vector<vector<int>> cycles) {

    for(int vertex1 : cycles[0]) {
        for (int vertex2 : cycles[1]) {

        }
    }
}


int steepest_local_search(
    vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    const bool swap_edges
) {
    int solution_delta = 0;
    int curr_delta = 0;

    while(true) {
        // zamiana wierzchołków między cyklami

        if (swap_edges) {
            // zamiana krawędzi wewnątrz cykli
        } else {
            // zamiana wierzchołków wewnątrz cykli
        }


        if (curr_delta > solution_delta) {
            solution_delta = curr_delta;
        } else {
            break;
        }
    }

    return solution_delta;
}


void evaluation_algorithm(
    const string& algorithm,
    const string& instance_name,
    const vector<vector<int>>& matrix,
    const vector<vector<int>>& cycles,
    const bool swap_edges=false, // if set to ``true``, then the algorithm should swap the edges instead of vertices inside the cycle
    int n=100
)
{
	for (int i = 0; i < n; i++)
	{
		vector<vector<int>> cycles_to_modify(cycles);

        if (algorithm == "greedy")
        {
            //TODO
        }
        else if (algorithm == "steepest")
        {
            steepest_local_search(cycles_to_modify, matrix, swap_edges); //TODO
        }
	}
    // TODO
}


int main() {
    srand(time(nullptr));

    for (string instance : INSTANCES)
	{
		int n;
		vector<vector<int>> v;
		read_file(INPUT_DIR + "/" + instance, &n, v);
		vector<vector<int>> matrix(n);
		for (int i = 0; i < matrix.size(); i++)
		{
			matrix[i] = vector<int>(n);
		}
		make_distance_matrix(v, matrix);

        for(string cycle_algo : ALGORITHMS)
        {
            vector<vector<int>> cycles(2);

            read_file(INPUT_DIR_CYCLES + "/" + RESULT_PREFIX + "_" + cycle_algo + "_" + instance, cycles);

            evaluation_algorithm("greedy", instance, matrix, cycles);
            evaluation_algorithm("greedy", instance, matrix, cycles, true);

            evaluation_algorithm("steepest", instance, matrix, cycles);
            evaluation_algorithm("steepest", instance, matrix, cycles, true);
            cout << "-----------------------------------" << "\n";
        }
	}
	return 0;
}
