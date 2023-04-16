#include <chrono>

#include "algorithms.hpp"

string INPUT_DIR = "../greedy_tsp/input";
string INSTANCES[] = { "kroa100.tsp", "krob100.tsp" };

string INPUT_DIR_CYCLES = "../greedy_tsp/output";
string RESULT_PREFIX = "min";
string ALGORITHMS[] = { "random", "regret" };

string OUTPUT_DIR = "output";

void swap_edges(vector<int>& cycle, int i, int j) // `i` i `j` to indeksy wierzchołków w cyklu które rozpoczynają krawędź
{
    int n = cycle.size();
    if(min(i,j) == 0 && max(i,j) == n-1)
    {
        swap(cycle[i], cycle[j]);
    }
    reverse(cycle.begin() + i + 1, cycle.begin() + j);
}

void swap_vertexes(vector<vector<int>>& cycles, int i, int j) {
    /*
    A function to swap vertices between cycles.
    `i` - vertex index in the first cycle
    `j` - vertex index in the second cycle
    */
    vector<int>* first_cycle = &(cycles[0]);
    vector<int>* second_cycle = &(cycles[1]);

    swap((*first_cycle)[i], (*second_cycle)[j]);
}

void swap_vertexes(vector<int>& cycle, int i, int j) {
    /*
    A function to swap vertices inside the cycle.
    `i`, `j` - vertices to swap
    */
    swap(cycle[i], cycle[j]);
}

int delta_swap_edges(const vector<vector<int>> & matrix,vector<int> cycle, int i, int j)// `i` i `j` to indeksy wierzchołków w cyklu które rozpoczynają krawędź
{
    int n = cycle.size();
    int a,b,c,d;
    if(min(i,j) == 0 && max(i,j) == n-1)
    {
        a = cycle[i];
		b = cycle[(i+1)%n];
		c = cycle[(j-1+n)%n];
		d = cycle[j];
    }
    else
    {
        a = cycle[(i-1+n)%n];
		b = cycle[i];
		c = cycle[j];
		d = cycle[(j+1)%n];
    }
    return matrix[a][c] + matrix[b][d] - matrix[a][b]-matrix[c][d];
}

int delta_swap_vertexes(const vector<vector<int>>& matrix, vector<vector<int>> cycles, int i, int j) {
    /*
    A function to swap vertices between cycles.
    `i` - vertex index in the first cycle
    `j` - vertex index in the second cycle
    */
    vector<int>* first_cycle = &(cycles[0]);
    vector<int>* second_cycle = &(cycles[1]);

    int n1 = (*first_cycle).size();
    int n2 = (*second_cycle).size();

    int i_prev = (*first_cycle)[(i - 1 + n1) % n1];
    int i_curr = (*first_cycle)[i];
    int i_next = (*first_cycle)[(i + 1) % n1];

    int j_prev = (*second_cycle)[(j - 1 + n2) % n2];
    int j_curr = (*second_cycle)[j];
    int j_next = (*second_cycle)[(j + 1) % n2];

    return (
        matrix[i_prev][j_curr] + matrix[j_curr][i_next] + matrix[j_prev][i_curr] + matrix[i_curr][j_next]
    ) - (
        matrix[i_prev][i_curr] + matrix[i_curr][i_next] + matrix[j_prev][j_curr] + matrix[j_curr][j_next]
    );
}

int delta_swap_vertexes(const vector<vector<int>>& matrix, vector<int> cycle, int i, int j) {
    /*
    A function to swap vertices inside the cycle.
    `i`, `j` - vertices to swap
    */
    int n = cycle.size();

    int i_prev = cycle[(i - 1 + n) % n];
    int i_curr = cycle[i];
    int i_next = cycle[(i + 1) % n];

    int j_prev = cycle[(j - 1 + n) % n];
    int j_curr = cycle[j];
    int j_next = cycle[(j + 1) % n];

    if (abs(i - j) == 1) { // wierzchołki sąsiadujące ze sobą
        return (
            matrix[i_prev][j_curr] + matrix[i_curr][j_next]
        ) - (
            matrix[i_prev][i_curr] + matrix[j_curr][j_next]
        );
    }

    if (i == 0 && j == cycle.size() - 1) { // wierzchołki na końcach cyklu
        return (
            matrix[j_curr][i_next] + matrix[j_prev][i_curr]
        ) - (
            matrix[i_curr][i_next] + matrix[j_prev][j_curr]
        );
    }

    return (
        matrix[i_prev][j_curr] + matrix[j_curr][i_next] + matrix[j_prev][i_curr] + matrix[i_curr][j_next]
    ) - (
        matrix[i_prev][i_curr] + matrix[i_curr][i_next] + matrix[j_prev][j_curr] + matrix[j_curr][j_next]
    );
}


int move_swap_edges(const vector<vector<int>> & matrix,vector<int>& cycle, bool steepest, int* v1=NULL, int* v2=NULL)
{
    /*
    steepest == true -> we cannot make a move ;)
    */
    std::default_random_engine rng(std::random_device{}());
    int n = cycle.size();
    vector<vector<int>>moves_edge_swap;
    for(int i = 0; i < n; i++)
    {
        for(int j = i + 1; j < n; j++)
        {
            moves_edge_swap.push_back({i,j});
        }
    }
    std::shuffle(moves_edge_swap.begin(), moves_edge_swap.end(), rng); // posorotwane ruchy
    int best_value = 0;
    for(vector<int> v : moves_edge_swap)
    {
        int delta = delta_swap_edges(matrix, cycle,v[0], v[1]);

        if( delta < 0)
        {
            if (steepest == false)
            {
                swap_edges(cycle, v[0], v[1]);
                return delta;
            }
            else
            {
                if(best_value > delta)
                {
                    best_value = delta;
                    *(v1) = v[0];
                    *(v2) = v[1];
                }
            }
        }
    }
    return best_value;
}

int move_swap_vertexes(const vector<vector<int>>& matrix, vector<vector<int>>& cycles, bool steepest=false, int* v1=NULL, int* v2=NULL)
{
    /*
    Look for a best move to swap vertices between two cycles.

    If a swap move will be performed, the function will return ``true``.
    */
    vector<int>* first_cycle = &(cycles[0]);
    vector<int>* second_cycle = &(cycles[1]);

    if ((*first_cycle).size() > (*second_cycle).size()) {
        first_cycle = &(cycles[1]);
        second_cycle = &(cycles[0]);
    }

    // generate all possible moves
    vector< vector<int> > moves;

    for (int i = 0; i < (*first_cycle).size(); ++i) {
        for (int j = i; j < (*second_cycle).size(); ++j) {
            moves.push_back({i, j});
        }
    }

    // shuffle
    default_random_engine rng = default_random_engine { random_device {}() };
    shuffle(moves.begin(), moves.end(), rng);

    // check moves
    int best_delta = 0;
    int delta;

    for (const vector<int>& move : moves) {
        delta = delta_swap_vertexes(matrix, cycles, move[0], move[1]);

        if (delta < 0) {
            if (!steepest) {
                swap_vertexes(cycles, move[0], move[1]);
                return delta;
            }

            if (delta < best_delta) {
                *(v1) = move[0];
                *(v2) = move[1];
                best_delta = delta;
            }
        }
    }
    return best_delta;
}

int move_swap_vertexes(const vector<vector<int>>& matrix, vector<int>& cycle, bool steepest=false, int* v1=NULL, int* v2=NULL)
{
    /*
    Look for a best move to swap vertices in cycle.

    If a swap move will be performed, the function will return ``true``.
    */

    // generate all possible moves
    vector< vector<int> > moves;

    for (int i = 0; i < cycle.size(); ++i) {
        for (int j = i + 1; j < cycle.size(); ++j) {
            moves.push_back({i, j});
        }
    }

    // shuffle
    default_random_engine rng = default_random_engine { random_device {}() };
    shuffle(moves.begin(), moves.end(), rng);

    // check moves
    int best_delta = 0;
    int delta;

    for (const vector<int>& move : moves) {
        delta = delta_swap_vertexes(matrix, cycle, move[0], move[1]);

        if (delta < 0) {
            if (!steepest) {
                swap_vertexes(cycle, move[0], move[1]);
                return delta;
            }

            if (delta < best_delta) {
                *(v1) = move[0];
                *(v2) = move[1];
                best_delta = delta;
            }
        }
    }

    return best_delta;
}


int steepest_local_search(
    vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    const bool perform_swap_edges
) {
    int solution_delta = 0;

    int min_curr_delta;
    int curr_delta_type1 = 0;

    int curr_delta1_type2 = 0;
    int curr_delta2_type2 = 0;

    int v1_type1 = 0, v2_type1 = 0;

    int v1_type2 = 0, v2_type2 = 0;
    int v3_type2 = 0, v4_type2 = 0;

    while(true) {
        // zamiana wierzchołków między cyklami
        curr_delta_type1 = move_swap_vertexes(matrix, solution, true, &v1_type1, &v2_type1);

        if (perform_swap_edges) {
            // zamiana krawędzi wewnątrz cykli
            curr_delta1_type2 = move_swap_edges(matrix, solution[0], true, &v1_type2, &v2_type2);
            curr_delta2_type2 = move_swap_edges(matrix, solution[1], true, &v3_type2, &v4_type2);
        } else {
            // zamiana wierzchołków wewnątrz cykli
            curr_delta1_type2 = move_swap_vertexes(matrix, solution[0], true, &v1_type2, &v2_type2);
            curr_delta2_type2 = move_swap_vertexes(matrix, solution[1], true, &v3_type2, &v4_type2);
        }

        min_curr_delta = min({
            curr_delta_type1,
            curr_delta1_type2,
            curr_delta2_type2
        });

        if (min_curr_delta >= 0) {
            break;
        }

        solution_delta += min_curr_delta;
        if (min_curr_delta == curr_delta_type1) {
            swap_vertexes(solution, v1_type1, v2_type1);
            // cout << "Swap: " << min_curr_delta << " " << v1_type1 << " " <<  v2_type1 << endl;
        }
        else if (min_curr_delta == curr_delta1_type2) {
            perform_swap_edges ? swap_edges(solution[0], v1_type2, v2_type2) : swap_vertexes(solution[0], v1_type2, v2_type2);
            // cout << "Swap in 1st cycle: " << min_curr_delta << " " << v1_type2 << " " <<  v2_type2 << " " << perform_swap_edges << endl;
        }
        else {
            perform_swap_edges ? swap_edges(solution[1], v3_type2, v4_type2) : swap_vertexes(solution[1], v3_type2, v4_type2);
            // cout << "Swap in 2nd cycle: " << min_curr_delta << " " << v3_type2 << " " <<  v4_type2 << " " << perform_swap_edges << endl;
        }
    }

    return solution_delta;
}

int greedy_local_search(
    vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    const bool flag_swap_edges
)
{
    int solution_delta = 0;

    while(true)
    {
        bool flag_end = true;
        int cycle_index = rand()%2;     // od którego cylu zaczniemy
        int edges_or_vertex = rand()%2; // zaczynbamy od zamiany krawędzi albo wierzchołków
        int outside_insde = rand()%2;

        if(flag_swap_edges == true)    // mamy zamianę krawędzi
        {
            if(edges_or_vertex == true)
            {
                // 1. zamiana krawedzi w obu cyklach
                // 2. zamiana wierzcholkow miedzy cyklami
                int delta = move_swap_edges(matrix, solution[cycle_index],false);

                if(!delta)
                {
                    delta = move_swap_edges(matrix, solution[1-cycle_index],false);
                }

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
            else
            {
                // 1. zamiana wierzcholkow miedzy cyklami
                // 2. zamiana krawedzi w obu cyklach
                int delta = move_swap_vertexes(matrix, solution);

                if(!delta)
                {
                    delta = move_swap_edges(matrix, solution[cycle_index],false);
                }

                if(!delta)
                {
                    delta = move_swap_edges(matrix, solution[1-cycle_index],false);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
        }
        else
        {
            if(edges_or_vertex == true)
            {
                // 1. zamiana wierzcholkow w obu cyklach
                // 2. zamiana wierzcholkow miedzy cyklami
                int delta = move_swap_vertexes(matrix, solution[cycle_index]);

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution[1-cycle_index]);
                }

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
            else
            {
                // 1. zamiana wierzcholkow miedzy cyklami
                // 2. zamiana wierzcholkow w obu cyklach
                int delta = move_swap_vertexes(matrix, solution);

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution[cycle_index]);
                }

                if(!delta)
                {
                   delta = move_swap_vertexes(matrix, solution[1-cycle_index]);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
        }
    }
}

int random_walk(vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    int time_ms
) {
    int result = 0;
    int v1, v2;
    auto start_time = chrono::high_resolution_clock::now();

    while(chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start_time).count() < time_ms) {
        int outside = rand() % 2;

        if (outside) { // wykonujemy ruch między cyklami (zamiana wierzchołków)
            v1 = rand() % solution[0].size();
            v2 = rand() % solution[1].size();

            result += delta_swap_vertexes(matrix, solution, v1, v2);
            swap_vertexes(solution, v1, v2);
        } else {
            int edges = rand() % 2;
            int cycle = rand() % 2;

            v1 = rand() % solution[cycle].size();
            do {
                v2 = rand() % solution[cycle].size();
            } while(v1 == v2);

            if (edges) { // zamiana krawędzi
                result += delta_swap_edges(matrix, solution[cycle], v1, v2);
                swap_edges(solution[cycle], v1, v2);
            } else {
                result += delta_swap_vertexes(matrix, solution[cycle], v1, v2);
                swap_vertexes(solution[cycle], v1, v2);
            }
        }
    }

    return result;
}


class Result {
    public:
        vector<vector<int>> solution;
        float time;
        int score;
        int delta;

        Result() {
            this->solution = vector< vector<int> >(2);
            this->time = .0;
            this->score = 0;
            this->delta = 0;
        }

        Result(vector<vector<int>> solution, float time, int score, int delta) {
            this->solution = solution;
            this->time = time;
            this->score = score;
            this->delta = delta;
        }

        void compute_score(const vector<vector<int>>& matrix) {
            this->score = score_cycle(matrix, this->solution[0]) + score_cycle(matrix, this->solution[1]);
        }
};


void evaluation_algorithm(
    const string& algorithm,
    const string& instance_name,
    const string& cycle_algo,
    const vector<vector<int>>& matrix,
    const bool swap_edges=false, // if set to ``true``, then the algorithm should swap the edges instead of vertices inside the cycle
    int n=100,
    int max_time_ms=184 // for random walk
)
{
    vector<Result> results(n); // wszystkie wyniki

    for (int i = 0; i < n; ++i) {

        // krok 1 - wygeneruj rozwiązanie
        if (cycle_algo == "random") {
            random_solution(matrix, results[i].solution);
        } else {
            regrest_heuristics(matrix, n, results[i].solution, i);
        }

        // krok 2 - zapisz wynik
        // results[i].compute_score(matrix);

        // krok 3 - zmodyfikuj rozwiązanie, zapisz deltę i czas
        int d = 0;
        float time;

        auto start = chrono::high_resolution_clock::now();
        if (algorithm == "greedy")
        {
            d = greedy_local_search(results[i].solution, matrix, swap_edges);
        }
        else if (algorithm == "steepest")
        {
            d = steepest_local_search(results[i].solution, matrix, swap_edges);
        }
        else if(algorithm == "random walk")
        {
            d = random_walk(results[i].solution, matrix, max_time_ms);
        }
        auto stop = chrono::high_resolution_clock::now();

        results[i].time = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
        results[i].delta = d;
        results[i].compute_score(matrix);

        // cout << i << " INSTANCE: " << instance_name << "; ALGO: " << algorithm << "; CYCLE ALGO: " << cycle_algo << endl;
        // cout << "DELTA: " << d << "; SCORE: " << results[i].score << "; TIME: " << results[i].time << endl;
    }
    string neighborhood = swap_edges ? "edges" : "vertices";
    cout << "INSTANCE: " << instance_name << "; ALGO: " << algorithm << "; CYCLE ALGO: " << cycle_algo << "; NEIGHBOURHOOD: " << neighborhood << endl;

    // score
    auto max_score = max_element(results.begin(), results.end(), [](Result a, Result b){ return a.score < b.score; });
    auto min_score = min_element(results.begin(), results.end(), [](Result a, Result b){ return a.score < b.score; });
    int mean_score = accumulate(results.begin(), results.end(), 0, [](const int& a, Result b){ return a + b.score; }) / n;
    cout << "SCORE (mean (min - max))" << endl;
    cout << mean_score << " (" << (*min_score).score << " - " << (*max_score).score << ")" << endl;

    // time
    auto max_time = max_element(results.begin(), results.end(), [](Result a, Result b){ return a.time < b.time; });
    auto min_time = min_element(results.begin(), results.end(), [](Result a, Result b){ return a.time < b.time; });
    float mean_time = accumulate(results.begin(), results.end(), 0, [](const float& a, Result b){ return a + b.time; }) / n;
    cout << "TIME (mean (min - max)) [ms]" << endl;
    cout << round(mean_time) << " (" << (*min_time).time << " - " << (*max_time).time << ")" << endl;

    // delta
    auto max_delta = max_element(results.begin(), results.end(), [](Result a, Result b){ return a.delta < b.delta; });
    auto min_delta = min_element(results.begin(), results.end(), [](Result a, Result b){ return a.delta < b.delta; });
    int mean_delta = accumulate(results.begin(), results.end(), 0, [](const int& a, Result b){ return a + b.delta; }) / n;
    cout << "DELTA (mean (min - max))" << endl;
    cout << mean_delta << " (" << (*max_delta).delta << " - " << (*min_delta).delta << ")" << endl;

    // best solution
    save_cycles_to_file(
        OUTPUT_DIR + "/" + cycle_algo + "_" + algorithm + "_" + neighborhood + "_" + instance_name,
        (*min_score).solution
    );
    cout << "----------------------------" << endl;
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
            evaluation_algorithm("greedy", instance, cycle_algo, matrix);
            evaluation_algorithm("greedy", instance, cycle_algo, matrix, true);

            evaluation_algorithm("steepest", instance, cycle_algo, matrix);
            evaluation_algorithm("steepest", instance, cycle_algo, matrix, true);
            cout << "-----------------------------------" << "\n";
        }
        evaluation_algorithm("random walk", "kroa100.tsp", "random", matrix, NULL, 100, 201);
        evaluation_algorithm("random walk", "krob100.tsp", "random", matrix, NULL, 100, 184);

        evaluation_algorithm("random walk", "kroa100.tsp", "regret", matrix, NULL, 100, 17);
        evaluation_algorithm("random walk", "krob100.tsp", "regret", matrix, NULL, 100, 19);
	}

}
