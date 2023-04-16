#include <chrono>

#include "algorithms.hpp"

using namespace std;

string INPUT_DIR = "./input";
string INSTANCES[] = { "kroA200.tsp", "kroB200.tsp" };
// string ALGORITHMS[] = { "regret", "steepest", "evaluation memory", "candidates" };
string ALGORITHMS[] = { "regret", "steepest", "evaluation memory" };
string OUTPUT_DIR = "output";

enum applicable {
    // gdy brakuje krawędzi, wierzchołków, etc.
    No,

    // wszystko ok, ale różny kierunek krawędzi (dla ruchu zamiany krawędzi)
    Maybe,

    // można użyć tu i teraz
    Yes
};

class Move {
    public:
        // wierzchołek pierwszy: indeks, numer poprzednika, numer, numer następnika
        int i_index, i_prev, i, i_succ;

        // wierzchołek drugi: indeks, numer poprzednika, numer, numer następnika
        int j_index, j_prev, j, j_succ;

        // 0 lub 1 dla ruchu wewnątrztrasowego (zamiany krawędzi), -1 dla ruchu międzytrasowego (zamiana wierzchołków)
        int cycle_number;

        // ocena ruchu
        int score;

        Move(
            const vector<vector<int>>& matrix,
            const vector<vector<int>>& cycles,
            int i, int j,
            int cycle_number
        ) {
            int n = cycles[cycle_number].size();

            this->i_index = i;
            this->j_index = j;

            this->i_prev = cycles[cycle_number][(i - 1 + n) % n];
            this->j_prev = cycles[cycle_number][(j - 1 + n) % n];

            this->i = cycles[cycle_number][i];
            this->j = cycles[cycle_number][j];

            this->i_succ = cycles[cycle_number][(i + 1) % n];
            this->j_succ = cycles[cycle_number][(j + 1) % n];

            this->score = delta_swap_edges(matrix, cycles[cycle_number], i_index, j_index);
            this->cycle_number = cycle_number;
        }

        Move(
            const vector<vector<int>>& matrix,
            const vector<vector<int>>& cycles,
            int i, int j
        )
        {
            int n1 = cycles[0].size();
            int n2 = cycles[1].size();

            this->i_index = i;
            this->j_index = j;

            this->i_prev = cycles[0][(i - 1 + n1) % n1];
            this->j_prev = cycles[1][(j - 1 + n2) % n2];

            this->i = cycles[0][i];
            this->j = cycles[1][j];

            this->i_succ = cycles[0][(i + 1) % n1];
            this->j_succ = cycles[1][(j + 1) % n2];

            this->score = delta_swap_vertexes(matrix, cycles, i_index, j_index);
            this->cycle_number = -1;
        }

        applicable is_applicable(const vector<vector<int>>& cycles){
            /*
                Checks if a move is applicable.

                For edges applicable, when edges have the same direction
                (i_prev -> i -> i_succ, j_prev -> j -> j_succ)
                OR
                (i_succ -> i -> i_prev, j_succ -> j -> j_prev).
                If edges are in different directions, the `Maybe` value is
                returned - cause we can use this move in the future.

                For vertexes applicable, when vertexes are in separable cycles
                and egdes i_prev -- i -- i_succ, j_prev -- j -- j_succ
                exist, in any direction.
            */
            vector<int>::const_iterator v_i, v_j;
            if (cycle_number == -1) { // zamiana wierzchołków
                const vector<int>* cycle_i;
                const vector<int>* cycle_j;

                if( (v_i = find(cycles[0].begin(), cycles[0].end(), i)) == cycles[0].end() )
                {
                    v_i = find(cycles[1].begin(), cycles[1].end(), i);
                    v_j = find(cycles[0].begin(), cycles[0].end(), j);

                    cycle_i = &(cycles[1]);
                    cycle_j = &(cycles[0]);
                }
                else
                {
                    v_j = find(cycles[1].begin(), cycles[1].end(), j);

                    cycle_i = &(cycles[0]);
                    cycle_j = &(cycles[1]);
                }

                // jeżeli wierzchołki są w 1 cyklu
                if (v_i == cycle_i->end() || v_j == cycle_j->end()) {
                    return No;
                }

                // jeżeli istnieją krawędzie i_prev - i - i_succ lub j_prev - j - j_succ
                this->i_index = v_i - cycle_i->begin();
                this->j_index = v_j - cycle_j->begin();
                int n_i = cycle_i->size(), n_j = cycle_j->size();

                int i_prev = (*cycle_i)[(i_index - 1 + n_i) % n_i];
                int i_succ = (*cycle_i)[(i_index + 1) % n_i];
                int j_prev = (*cycle_j)[(j_index - 1 + n_j) % n_j];
                int j_succ = (*cycle_j)[(j_index + 1) % n_j];
                if ((
                        (i_prev == this->i_prev && i_succ == this->i_succ) ||
                        (i_prev == this->i_succ && i_succ == this->i_prev)
                    ) && (
                        (j_prev == this->j_prev && j_succ == this->j_succ) ||
                        (j_prev == this->j_succ && j_succ == this->j_prev)
                    )
                ) {
                    return Yes;
                }

                // krawędzie nie istnieją, zwróć No
                return No;
            } else
            {
                const vector<int>* cycle = &(cycles[cycle_number]);
                int n = cycle->size();

                v_i = find(cycle->begin(), cycle->end(), i);
                v_j = find(cycle->begin(), cycle->end(), j);

                if (v_i == cycle->end() || v_j == cycle->end()) {
                    return No;
                }

                this->i_index = v_i - cycle->begin();
                this->j_index = v_j - cycle->begin();

                int i_prev = (*cycle)[(i_index - 1 + n) % n];
                int i_succ = (*cycle)[(i_index + 1) % n];
                int j_prev = (*cycle)[(j_index - 1 + n) % n];
                int j_succ = (*cycle)[(j_index + 1) % n];

                if (
                    ((i_prev == this->i_prev && i_succ == this->i_succ) && (j_prev == this->j_prev && j_succ == this->j_succ)) ||
                    ((i_prev == this->i_succ && i_succ == this->i_prev) && (j_prev == this->j_succ && j_succ == this->j_prev))
                )
                {
                    return Yes;
                }

                if ((
                        (i_prev == this->i_prev && i_succ == this->i_succ) ||
                        (i_prev == this->i_succ && i_succ == this->i_prev)
                    ) && (
                        (j_prev == this->j_prev && j_succ == this->j_succ) ||
                        (j_prev == this->j_succ && j_succ == this->j_prev)
                    )
                ) {
                    return Maybe;
                }

                return No;
            }

        }

        void apply(vector<vector<int>>& cycles) {
            if (cycle_number == -1) {
                swap_vertexes(cycles, this->i_index, this->j_index);
            } else {
                swap_edges(cycles[cycle_number], this->i_index, this->j_index);
            }
        }
};


void init_improving_moves(
    vector<Move>& moves,
    const vector<vector<int>>& cycles,
    const vector<vector<int>>& matrix
) {
    // 1. ruchy wewnątrztrasowe
    for (int cycle_index = 0; cycle_index < cycles.size(); ++cycle_index) {
        for (int i = 0; i < cycles[cycle_index].size(); ++i) {
            for (int j = i + 1; j < cycles[cycle_index].size(); ++j) {
                Move new_move(
                    matrix, cycles, i, j, cycle_index
                );

                // dodajemy do listy tylko jeżeli przynosi poprawę
                if (new_move.score < 0) {
                    moves.push_back(new_move);
                }
            }
        }
    }

    // 2. ruchy międzytrasowe
    const vector<int>* first_cycle = &(cycles[0]);
    const vector<int>* second_cycle = &(cycles[1]);

    if ((*first_cycle).size() > (*second_cycle).size()) {
        first_cycle = &(cycles[1]);
        second_cycle = &(cycles[0]);
    }

    for (int i = 0; i < (*first_cycle).size(); ++i) {
        for (int j = i; j < (*second_cycle).size(); ++j) {
            Move new_move(matrix, cycles, i, j);

            // dodajemy do listy tylko jeżeli przynosi poprawę
            if (new_move.score < 0) {
                moves.push_back(new_move);
            }
        }
    }
}


void add_new_moves(
    Move applied_move,
    vector<Move>& moves,
    const vector<vector<int>>& cycles,
    const vector<vector<int>>& matrix
) {
    if(applied_move.cycle_number == -1) {

        int i_prev = (applied_move.i_index - 1 + cycles[0].size()) % cycles[0].size();
        int i_index = applied_move.i_index;

        int j_prev = (applied_move.j_index - 1 + cycles[1].size()) % cycles[1].size();
        int j_index = applied_move.j_index;

        // 2 nowe krawędzie wychodzące z i_prev, i (zamiana krawędzy)
        for(int i = 0; i < cycles[0].size(); ++i) {
            if (i != i_prev) {
                Move new_move = i < i_prev ? Move(matrix, cycles, i, i_prev, 0) : Move(matrix, cycles, i_prev, i, 0);
                if(new_move.score < 0) moves.push_back(new_move);
            }
            if (i != i_index) {
                Move new_move = i < i_index ? Move(matrix, cycles, i, i_index, 0) : Move(matrix, cycles, i_index, i, 0);
                if(new_move.score < 0) moves.push_back(new_move);
            }

            // zamiana wierzchołków
            Move new_move(matrix, cycles, i, j_prev);
            if(new_move.score < 0) moves.push_back(new_move);

            new_move = Move(matrix, cycles, i, j_index);
            if(new_move.score < 0) moves.push_back(new_move);
        }

        // 2 nowe krawędzie wychodzące z j_prev, j (zamiana krawędzi)
        for(int j = 0; j < cycles[1].size(); ++j) {
            if (j != j_prev) {
                Move new_move = j < j_prev ? Move(matrix, cycles, j, j_prev, 1) : Move(matrix, cycles, j_prev, j, 1);
                if(new_move.score < 0) moves.push_back(new_move);
            }
            if (j != j_index) {
                Move new_move = j < j_index ? Move(matrix, cycles, j, j_index, 1) : Move(matrix, cycles, j_index, j, 1);
                if(new_move.score < 0) moves.push_back(new_move);
            }

            // zamiana wierzchołków
            Move new_move(matrix, cycles, i_prev, j);
            if(new_move.score < 0) moves.push_back(new_move);

            new_move = Move(matrix, cycles, i_index, j);
            if(new_move.score < 0) moves.push_back(new_move);
        }
    } else {
        int i_succ = (applied_move.i_index + 1) % cycles[applied_move.cycle_number].size();
        int i_index = applied_move.i_index;
        int j_index = applied_move.j_index;

        // 2 nowe krawędzie wychodzące z i_prev, i (zamiana krawędzy)
        for(int i = 0; i < cycles[applied_move.cycle_number].size(); ++i) {
            if (i != i_succ) {
                Move new_move = i < i_succ ? Move(matrix, cycles, i, i_succ, applied_move.cycle_number) : Move(matrix, cycles, i_succ, i, applied_move.cycle_number);
                if(new_move.score < 0) moves.push_back(new_move);
            }
            if (i != i_index) {
                Move new_move = i < i_index ? Move(matrix, cycles, i, i_index, applied_move.cycle_number) : Move(matrix, cycles, i_index, i, applied_move.cycle_number);
                if(new_move.score < 0) moves.push_back(new_move);
            }
            if (i != j_index) {
                Move new_move = i < j_index ? Move(matrix, cycles, i, j_index, applied_move.cycle_number) : Move(matrix, cycles, j_index, i, applied_move.cycle_number);
                if(new_move.score < 0) moves.push_back(new_move);
            }
        }

        // zamiana wierzchołków
        int other_cycle = (applied_move.cycle_number + 1) % 2;
        for(int i = 0; i < cycles[other_cycle].size(); ++i) {
            if(other_cycle == 0) {
                Move new_move = Move(matrix, cycles, i, i_index);
                if(new_move.score < 0) moves.push_back(new_move);

                new_move = Move(matrix, cycles, i, i_succ);
                if(new_move.score < 0) moves.push_back(new_move);

                new_move = Move(matrix, cycles, i, j_index);
                if(new_move.score < 0) moves.push_back(new_move);
            } else {
                Move new_move = Move(matrix, cycles, i_index, i);
                if(new_move.score < 0) moves.push_back(new_move);

                new_move = Move(matrix, cycles, i_succ, i);
                if(new_move.score < 0) moves.push_back(new_move);

                new_move = Move(matrix, cycles, j_index, i);
                if(new_move.score < 0) moves.push_back(new_move);
            }
        }
    }
}


void local_search_evaluation_memory(
    vector<vector<int>>& cycles,
    const vector<vector<int>>& matrix
) {
    vector<Move> moves;
    init_improving_moves(moves, cycles, matrix);

    bool move_was_found = true;
    applicable is_applicable;
    int i;

    while (move_was_found) {
        move_was_found = false;

        // sort moves by scores, from best to worst
        sort(moves.begin(), moves.end(), [](Move a, Move b){ return a.score < b.score; });

        i = 0;
        while (i < moves.size()) {
            is_applicable = moves[i].is_applicable(cycles);

            if (is_applicable == Yes) {
                // cout << moves[i].score << " " << moves[i].cycle_number << " " << moves[i].i << " " << moves[i].i_succ << "   " << moves[i].j << " " << moves[i].j_succ << endl;
                // cout << moves[i].i_index << " " << moves[i].j_index << endl;

                move_was_found = true;
                moves[i].apply(cycles);
                add_new_moves(moves[i], moves, cycles, matrix);

                moves.erase(moves.begin() + i);
                break;
            } else if (is_applicable == No) {
                moves.erase(moves.begin() + i);
                continue;
            }
            ++i;
        }
    }
}


class Result {
    public:
        vector<vector<int>> solution;
        float time;
        int score;

        Result() {
            this->solution = vector< vector<int> >(2);
            this->time = .0;
            this->score = 0;
        }

        Result(vector<vector<int>> solution, float time, int score, int delta) {
            this->solution = solution;
            this->time = time;
            this->score = score;
        }

        void compute_score(const vector<vector<int>>& matrix) {
            this->score = score_cycle(matrix, this->solution[0]) + score_cycle(matrix, this->solution[1]);
        }
};


void evaluation_algorithm(
    const string& algorithm,
    const string& instance_name,
    const vector<vector<int>>& matrix,
    int n=100
)
{
    vector<Result> results(n); // wszystkie wyniki
    for (int i = 0; i < n; ++i) {
        // krok 1 - wygeneruj rozwiązanie
        if (algorithm != "regret"){
            random_solution(matrix, results[i].solution);
        }

        // krok 2 - uruchom algorytm dla rozwiązania losowego
        float time;

        auto start = chrono::high_resolution_clock::now();
        if (algorithm == "regret")
        {
            regrest_heuristics(matrix, results[i].solution);
        }
        else if (algorithm == "steepest")
        {
            steepest_local_search(results[i].solution, matrix, true);
        }
        else if(algorithm == "evaluation memory")
        {
            local_search_evaluation_memory(results[i].solution, matrix);
        }
        else if(algorithm == "candidates")
        {
            candidates_algorithm(results[i].solution, matrix);
        }
        auto stop = chrono::high_resolution_clock::now();

        results[i].time = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
        results[i].compute_score(matrix);

        // cout << i << " INSTANCE: " << instance_name << "; ALGO: " << algorithm << "; CYCLE ALGO: " << cycle_algo << endl;
        // cout << "DELTA: " << d << "; SCORE: " << results[i].score << "; TIME: " << results[i].time << endl;
    }
    cout << "INSTANCE: " << instance_name << "; ALGO: " << algorithm << endl;

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
    // auto max_delta = max_element(results.begin(), results.end(), [](Result a, Result b){ return a.delta < b.delta; });
    // auto min_delta = min_element(results.begin(), results.end(), [](Result a, Result b){ return a.delta < b.delta; });
    // int mean_delta = accumulate(results.begin(), results.end(), 0, [](const int& a, Result b){ return a + b.delta; }) / n;
    // cout << "DELTA (mean (min - max))" << endl;
    // cout << mean_delta << " (" << (*max_delta).delta << " - " << (*min_delta).delta << ")" << endl;

    // best solution
    save_cycles_to_file(
        OUTPUT_DIR + "/" + algorithm + "_" + instance_name,
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
    
        for(string algorithm : ALGORITHMS)
        {
            evaluation_algorithm(algorithm, instance, matrix, 5);
            cout << "-----------------------------------" << "\n";
        }
    }

}