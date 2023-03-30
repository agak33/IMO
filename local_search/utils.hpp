#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <random> 
using namespace std;

void read_file(string path, vector<vector<int>>& cycles) {
    fstream input_file;
    input_file.open(path, ios::in);

    vector<int> data;
    int vertex;

    if(input_file.good()) {
        while(input_file >> vertex) {
            data.push_back(vertex);
        }
    } else {
        cout << "Error while reading the file: " << path << endl;
		input_file.close();
		return;
    }
    input_file.close();

    int middle_index = ceil(data.size() / 2.0);
    copy(data.begin(), data.begin() + middle_index, back_inserter(cycles[0]));
    copy(data.begin() + middle_index, data.end(), back_inserter(cycles[1]));
}

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

void make_distance_matrix(const vector<vector<int>>& coords, vector<vector<int>>& matrix)
{
	for (int i = 1; i <=coords.size(); i++)
	{
		for (int j = 1; j <=coords.size(); j++)
		{
			matrix[i][j] = round(sqrt(pow(coords[i-1][0] - coords[j-1][0], 2) + pow(coords[i-1][1] - coords[j-1][1], 2)));
		}
	}
}

int score_cycle(const vector<vector<int>>& matrix,const vector<int>& cycle)
{
	int score = 0;
	int n = cycle.size();
	for (int i = 0; i < cycle.size(); i++)
	{
		score += matrix[cycle[i]][cycle[(i + 1) % n]];
	}
	return score;
}
