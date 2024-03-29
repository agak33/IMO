#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <tuple>
#include <map>

#include "utils.hpp"
#include <random>

using namespace std;


void random_solution(const vector<vector<int>>& matrix, vector<vector<int>>& cycles)
{
	vector<int> remaining;
	for(int i = 0; i < matrix.size(); i++) {
		remaining.push_back(i);
	}

	auto random_number_generator = std::default_random_engine { std::random_device {}() };
	shuffle(remaining.begin(), remaining.end(), random_number_generator);

	int middle_index = ceil(remaining.size() / 2.0);
    copy(remaining.begin(), remaining.begin() + middle_index, back_inserter(cycles[0]));
    copy(remaining.begin() + middle_index, remaining.end(), back_inserter(cycles[1]));
}


void nearest_neighbour(const vector<vector<int>>& matrix, vector<vector<int>>& cycles, int start_cycle = -1)
{
	set<int> remaining;
	for(int i = 0; i < matrix.size(); i++) {
		remaining.insert(i);
	}
	int v1 = start_cycle == -1 ? rand() % matrix.size() : start_cycle;
	int v2 = find_farthest_vertex(matrix, v1, remaining);

	add_vertex_to_cycle(v1, cycles[0], remaining);
	add_vertex_to_cycle(v2, cycles[1], remaining);

	while(!remaining.empty())
	{
		for(int i = 0; i < 2; i++){
			int vertex_begin = cycles[i][0];
			int vertex_end = cycles[i][cycles[i].size() - 1];

			int nearest_vertex_begin = find_nearest_vertex(matrix, vertex_begin, remaining);
			int nearest_vertex_end = find_nearest_vertex(matrix, vertex_end, remaining);

			if(matrix[vertex_begin][nearest_vertex_begin] <= matrix[vertex_end][nearest_vertex_end]){
				add_vertex_to_cycle(nearest_vertex_begin, cycles[i], remaining, 0);
			} else {
				add_vertex_to_cycle(nearest_vertex_end, cycles[i], remaining);
			}
		}
	}
}


pair<int, int> score_extention_cycle(const vector<vector<int>>& matrix, vector<int>& cycle, set<int>& remaining)
{
	int min_score = 1e9;
	int min_index = -1;
	int min_vertex = -1;
	int n = cycle.size();
	for (int i = 0; i < cycle.size(); i++)
	{
		for (auto it = remaining.begin(); it != remaining.end(); ++it)
		{
			int score_diff = matrix[cycle[i]][*it] + matrix[*it][cycle[(i + 1) % n]] - matrix[cycle[i]][cycle[(i + 1) % n]];
			if (score_diff < min_score)
			{
				min_score = score_diff;
				min_index = (i + 1) % n;
				min_vertex = *it;
			}
		}
	}
	return { min_index, min_vertex };
}

void greedy_cycle(const vector<vector<int>>& matrix, int  n, vector<vector<int>>& cycles, int start_cycle = -1)
{
	int v1;
	if (start_cycle == -1)
	{
		v1 = rand() % 100;
	}
	else
	{
		v1 = start_cycle;
	}
	set<int> remaining;
	for (int i = 0; i < n; i++)
	{
		remaining.insert(i);
	}

	remaining.erase(v1);
	cycles[0].push_back(v1);
	int v2 = find_farthest_vertex(matrix, v1, remaining);
	cycles[1].push_back(v2);
	remaining.erase(v2);
	v1 = find_nearest_vertex(matrix, cycles[0][0], remaining);
	cycles[0].push_back(v1);
	remaining.erase(v1);
	v2 = find_nearest_vertex(matrix, cycles[1][0],remaining);
	cycles[1].push_back(v2);
	remaining.erase(v2);
	while (!remaining.empty())
	{
		for (int i = 0; i < 2; i++)
		{
			int best_index, best_vertex;
			tie(best_index, best_vertex) = score_extention_cycle(matrix, cycles[i], remaining);
			cycles[i].insert(cycles[i].begin() + best_index, best_vertex);
			remaining.erase(best_vertex);
		}
	}
	return;
}
pair<int, int> regret2(const vector<vector<int>>& matrix, set<int>& remaining, vector<int>& cycle, float weight=1.37)
{
	vector<vector<float>> regrets(remaining.size());
	int n = cycle.size();
	for (int i = 0; i < remaining.size(); i++)
	{
		regrets[i].resize(n);
	}
	map<int, int>m;
	int cc = 0;
	for (int i = 0; i < n; i++)
	{
		for (auto it = remaining.begin(); it != remaining.end(); ++it)
		{
			if (m.find(*it) == m.end())
			{
				m[*it] = cc;
				cc++;
			}
			regrets[m[*it]][i] = matrix[cycle[i]][*it] + matrix[*it][cycle[(i + 1) % n]] - matrix[cycle[i]][cycle[(i + 1) % n]];

		}
	}
	vector<vector<float>> tmp_regrets;
	copy(regrets.begin(), regrets.end(), back_inserter(tmp_regrets));
	float max_regret = -1e9;
	int best_index = -1;
	float tmp_regret_min = 0;
	for (int j = 0; j < tmp_regrets.size(); j++)
	{
		sort(tmp_regrets[j].begin(), tmp_regrets[j].end());
		float tmp_regret;
		if (tmp_regrets[j].size() <= 1)
		{
			 tmp_regret =  tmp_regrets[j][0];
		}
		else
		{
			 tmp_regret = tmp_regrets[j][1] - weight *  tmp_regrets[j][0];
		}

		if (max_regret <= tmp_regret)
		{
			max_regret = tmp_regret;
			best_index = j;
			tmp_regret_min = tmp_regrets[j][0];
		}
	}
	for(int i =0;i<n; i++)
	{
		if (regrets[best_index][i] == tmp_regret_min)
		{
			for (const auto& pair : m)
			{
				if (pair.second == best_index)
				{
					return { (i+1)%n, (pair.first) };
				}
			}
		}
	}

}
void regrest_heuristics(const vector<vector<int>>& matrix, int  n, vector<vector<int>>& cycles, int start_vertex = -1, float weight=1.37)
{
	set<int> remaining;
	for (int i = 0; i < n; i++)
	{
		remaining.insert(i);
	}
	int v1;
	if (start_vertex == -1)
	{
		int v1 = rand() % 100;
	}
	else
	{
		v1 = start_vertex;
	}
	cycles[0].push_back(v1);
	remaining.erase(v1);
	int v2 = find_farthest_vertex(matrix, v1, remaining);

	cycles[1].push_back(v2);
	remaining.erase(v2);

	int licznik = 0;

	v1 = find_nearest_vertex(matrix, v1, remaining);
	cycles[0].push_back(v1);
	remaining.erase(v1);

	v2 = find_nearest_vertex(matrix, v2, remaining);
	cycles[1].push_back(v2);
	remaining.erase(v2);

	while (!remaining.empty())
	{
		for (int i = 0; i < 2; i++)
		{
			int best_index, best_vertex;
			tie(best_index, best_vertex) = regret2(matrix, remaining, cycles[i], weight);
			cycles[i].insert(cycles[i].begin() + best_index, best_vertex);
			remaining.erase(best_vertex);
		}
	}
	return;

}


void evaluation_algorithm(string algorithm, string instance,  const vector<vector<int>>& matrix, int n)
{
	vector < pair<int, vector<vector<int>>>> results; // vector<pair<ocena, cykle>>
	for (int i = 0; i < n; i++)
	{
		vector<vector<int>>cycles(2);
		results.push_back({0, cycles});
		if (algorithm == "regret")
		{
			regrest_heuristics(matrix, n, results[i].second, i);

		}
		if (algorithm == "regret_no_weight")
		{
			regrest_heuristics(matrix, n, results[i].second, i, 1.0);
		}
		if (algorithm == "cycle")
		{
			greedy_cycle(matrix, n, results[i].second, i);
		}
		if (algorithm == "neighbour")
		{
			nearest_neighbour(matrix, results[i].second, i);
		}
		if (algorithm == "random")
		{
			random_solution(matrix, results[i].second);
		}
		results[i].first = score_cycle(matrix, results[i].second[0]) + score_cycle(matrix, results[i].second[1]);
	}
	sort(results.begin(), results.end());
	int mean = 0;
	for (int i = 0; i < results.size(); i++)
	{
		mean += results[i].first;
	}
	cout << "-----------------------------------"<<"\n";
	cout << algorithm +"  "<< instance << "\n";
	cout << "najlepszy wynik: " << results[0].first << "\n";
	cout << "najgorszy wynik: " << results[n - 1].first << "\n";
	cout << "średni wynik: " << mean/n << "\n";
	save_cycles_to_file("output/min_" + algorithm+"_" + instance, results[0].second);
	save_cycles_to_file("output/max_" + algorithm + "_" + instance, results[n - 1].second);
	save_cycles_to_file("output/median_" + algorithm + "_" + instance, results[n / 2].second);
}

int main()
{
	srand(time(nullptr));
	string instances[] = { "kroa100.tsp", "krob100.tsp" };
	string dir = "input";
	for (string instance : instances)
	{
		int n;
		vector<vector<int>> v;
		read_file(dir + "/" + instance, &n, v);
		vector<vector<int>> matrix(n);
		for (int i = 0; i < matrix.size(); i++)
		{
			matrix[i] = vector<int>(n);
		}
		make_distance_matrix(v, matrix);

		evaluation_algorithm("regret", instance, matrix, n);
		evaluation_algorithm("regret_no_weight", instance, matrix, n);
		evaluation_algorithm("cycle", instance, matrix, n);
		evaluation_algorithm("neighbour", instance, matrix, n);
		evaluation_algorithm("random", instance, matrix, n);
	}
	cout << "-----------------------------------" << "\n";
	return 0;
}