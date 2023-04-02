#include <vector>
#include <algorithm>
#include <tuple>
#include <set>
#include <map>
#include <random>

#include "utils.hpp"

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
