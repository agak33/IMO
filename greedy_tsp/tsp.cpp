#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
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

int main()
{
	int n;
	vector<vector<int>> v;
	read_file("kroa100.tsp", &n, v);
	vector<vector<int>> matrix(n);
	for (int i = 0; i < matrix.size(); i++)
	{
		matrix[i] = vector<int>(n);
	}
	make_distance_matrix(v, matrix, n);
	cout << matrix[1][2];
	return 0;
}