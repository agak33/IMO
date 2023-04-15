#include <chrono>

#include "algorithms.hpp"

using namespace std;

string INPUT_DIR = "./input";
string INSTANCES[] = { "kroA200.tsp", "kroB200.tsp" };
string ALGORITHMS[] = { "regret", "steepest", "evaluation memory", "candidates" };
string OUTPUT_DIR = "output";


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

        }
        else if(algorithm == "candidates")
        {
            // TODO
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