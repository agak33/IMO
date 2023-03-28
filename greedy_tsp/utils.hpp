#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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


void save_cycles_to_file(string file_name, vector<vector<int>>& cycles)
{
	ofstream File(file_name);
	for (vector<int> cycle : cycles)
	{
		for (int vertex_index : cycle)
		{
			File << vertex_index + 1 << "\n";
		}
	}
}


int find_farthest_vertex(const vector<vector<int>>& matrix, int v, set<int>& remaining)
{
	int max_distance = -1;
	int max_vertex = -1;
	for (auto it = remaining.begin(); it != remaining.end(); ++it)
	{
		if (matrix[*it][v] > max_distance)
		{
			max_distance = matrix[*it][v];
			max_vertex = *it;
		}
	}
	return max_vertex;
}


int find_nearest_vertex(const vector<vector<int>>& matrix, int v, set<int>& remaining)
{
	int min_distance = 1e9;
	int min_vertex = -1;
	for (auto it = remaining.begin(); it != remaining.end(); ++it)
	{
		if (matrix[*it][v] < min_distance && matrix[*it][v] > 0)
		{
			min_distance = matrix[*it][v];
			min_vertex = *it;
		}
	}
	return min_vertex;
}


void add_vertex_to_cycle(int vertex, vector<int>& cycle, set<int>& remaining, int insert_index=-1)
{
	remaining.erase(vertex);
	if (insert_index == -1){
		cycle.push_back(vertex);
	} else {
		cycle.insert(cycle.begin() + insert_index, vertex);
	}
}


void make_distance_matrix(const vector<vector<int>>& coords, vector<vector<int>>& matrix)
{
	for (int i = 0; i < coords.size(); i++)
	{
		for (int j = 0; j < coords.size(); j++)
		{
			matrix[i][j] = round(sqrt(pow(coords[i][0] - coords[j][0], 2) + pow(coords[i][1] - coords[j][1], 2)));
		}
	}
}

int score_cycle(const vector<vector<int>>& matrix, vector<int>& cycle)
{
	int score = 0;
	int n = cycle.size();
	for (int i = 0; i < cycle.size(); i++)
	{
		score += matrix[cycle[i]][cycle[(i + 1) % n]];
	}
	return score;
}