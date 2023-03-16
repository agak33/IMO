#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <tuple>
#include<map>
using namespace std;

void read_file(string file_name, int* n, vector<vector<int>>& cords)
{
	ifstream file(file_name);
	string line;
	int line_number = 0;
	bool tmp = file.good();
	while (getline(file, line))
	{
		line_number++;
		if (line_number == 4)
		{
			line = line.substr(10);
			*n = stoi(line);
		}
		if (line_number == 6)
		{
			break;
		}
	}
	cords.resize(*n);
	for (int i = 0; i < *n; i++)
	{
		int a, b, c;
		file >> a >> b >> c;
		cords[i].push_back(b);
		cords[i].push_back(c);
	}
	file.close();
}

int find_farthest_vertex(vector<vector<int>>& matrix, int  n, int v, set<int>& remaining)
{
	int max_distance = -1;
	int max_vertex = -1;
	for (auto it = remaining.begin(); it != remaining.end(); ++it)
	{
		if (matrix[*it - 1][v - 1] > max_distance)
		{
			max_distance = matrix[*it - 1][v - 1];
			max_vertex = *it;
		}
	}
	return max_vertex;
}

int find_nearest_vertex(vector<vector<int>>& matrix, int v, set<int>& remaining)
{
	int min_distance = 1e9;
	int min_vertex = -1;
	for (auto it = remaining.begin(); it != remaining.end(); ++it)
	{
		if (matrix[*it - 1][v - 1] < min_distance && matrix[*it - 1][v - 1] > 0)
		{
			min_distance = matrix[*it-1][v - 1];
			min_vertex = *it;
		}
	}
	return min_vertex;
}


void make_distance_matrix(vector<vector<int>>& coords, vector<vector<int>>& matrix, int  n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrix[i][j] = round(sqrt(pow(coords[i][0] - coords[j][0], 2) + pow(coords[i][1] - coords[j][1], 2)));
		}
	}
}
pair<int, int> score_extention_cycle(vector<vector<int>>& matrix, vector<int>& cycle, set<int>& remaining)
{
	int min_score = 1e9;
	int min_index = -1;
	int min_vertex = -1;
	int n = cycle.size();
	for (int i = 0; i < cycle.size(); i++)
	{
		for (auto it = remaining.begin(); it != remaining.end(); ++it)
		{
			int score_diff = matrix[cycle[i] - 1][*it - 1] + matrix[*it - 1][cycle[(i + 1) % n] - 1] - matrix[cycle[i] - 1][cycle[(i + 1) % n] - 1];
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
void greedy_cycle(vector<vector<int>>& matrix, int  n, vector<vector<int>>& cycles, int start_cycle = -1)
{
	int v1;
	if (start_cycle == -1)
	{
		v1 = rand() % 100 + 1;
	}
	else
	{
		v1 = start_cycle;
	}
	set<int> remaining;
	for (int i = 1; i <= n; i++)
	{
		remaining.insert(i);
	}
	
	remaining.erase(v1);
	cycles[0].push_back(v1);
	int v2 = find_farthest_vertex(matrix, n, v1, remaining);
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
pair<int, int> regret2(vector<vector<int>>& matrix, set<int>& remaining, vector<int>& cycle)
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
			if (m.find(*it - 1) == m.end())
			{
				m[*it - 1] = cc;
				cc++;
			}
			regrets[m[*it - 1]][i] = matrix[cycle[i] - 1][*it - 1] + matrix[*it - 1][cycle[(i + 1) % n] - 1] - matrix[cycle[i] - 1][cycle[(i + 1) % n] - 1];
			
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
			 tmp_regret = tmp_regrets[j][1] - 1.37 *  tmp_regrets[j][0];
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
					return { (i+1)%n, (pair.first +1) };
				}
			}
		}
	}

}
void regrest_heuristics(vector<vector<int>>& matrix, int  n, vector<vector<int>>& cycles, int start_vertex = -1)
{
	set<int> remaining;
	for (int i = 1; i <= n; i++)
	{
		remaining.insert(i);
	}
	int v1;
	if (start_vertex == -1)
	{
		int v1 = rand() % 100 + 1;
	}
	else
	{
		v1 = start_vertex;
	}
	cycles[0].push_back(v1);
	remaining.erase(v1);
	int v2 = find_farthest_vertex(matrix, n, v1, remaining);
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
			tie(best_index, best_vertex) = regret2(matrix, remaining, cycles[i]);
			cycles[i].insert(cycles[i].begin() + best_index, best_vertex);
			remaining.erase(best_vertex);
		}
	}
	return;

}
int score_cycle(vector<vector<int>>& matrix, vector<int>& cycle)
{
	int score = 0;
	int n = cycle.size();
	for (int i = 0; i < cycle.size(); i++)
	{
		score += matrix[cycle[i] - 1][cycle[(i + 1) % n] - 1];
	}
	return score;
}
void save_cycles_to_file(string file_name, vector<vector<int>>& cycles)
{
	ofstream File(file_name);
	for (int i = 0; i < cycles.size(); i++)
	{
		for (int j = 0; j < cycles[i].size(); j++)
		{
			File << cycles[i][j] << "\n";
		}
	}
}

void evaluation_algorithm(string algorithm, string instance,  vector<vector<int>>& matrix, int n)
{
	int min_vertex=1;
	int max_vertex=1;
	vector < pair<int, vector<vector<int>>>> results; // vector<pair<ocena, cykle>>
	for (int i = 1; i <= n; i++)
	{
		vector<vector<int>>cycles(2);
		results.push_back({0, cycles});
		if (algorithm == "regret")
		{
			regrest_heuristics(matrix, n, results[i - 1].second, i);
			
		}
		if (algorithm == "cycle")
		{
			greedy_cycle(matrix, n, results[i - 1].second, i);
		}
		results[i - 1].first = score_cycle(matrix, results[i - 1].second[0]) + score_cycle(matrix, results[i - 1].second[1]);
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
	cout << "œredni wynik: " << mean/n << "\n";
	save_cycles_to_file("min_" + algorithm+"_" + instance, results[0].second);
	save_cycles_to_file("max_" + algorithm + "_" + instance, results[n - 1].second);
	save_cycles_to_file("median" + algorithm + "_" + instance, results[n / 2].second);
}

int main()
{
	srand(time(nullptr));
	string instances[] = { "kroa100.tsp", "krob100.tsp" };
	for (string instance : instances)
	{
		int n;
		vector<vector<int>> v;
		read_file(instance, &n, v);
		vector<vector<int>> matrix(n);
		for (int i = 0; i < matrix.size(); i++)
		{
			matrix[i] = vector<int>(n);
		}
		make_distance_matrix(v, matrix, n);
		evaluation_algorithm("regret", instance, matrix, n);
		evaluation_algorithm("cycle", instance, matrix, n);
	}
	cout << "-----------------------------------" << "\n";
	return 0;
}