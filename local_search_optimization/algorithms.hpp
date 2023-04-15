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
void regrest_heuristics(const vector<vector<int>>& matrix, vector<vector<int>>& cycles, int start_vertex = -1, float weight=1.37)
{
	set<int> remaining;
    int n = matrix.size();
	for (int i = 0; i < n; i++)
	{
		remaining.insert(i);
	}
	int v1;
	if (start_vertex == -1)
	{
		v1 = rand() % n;
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

void swap_edges(vector<int>& cycle, int i, int j) // `i` i `j` to indeksy wierzchołków w cyklu które rozpoczynają krawędź
{
    int n = cycle.size();
    if(min(i,j) == 0 && max(i,j) == n-1)
    {
        swap(cycle[i], cycle[j]);
    }
    reverse(cycle.begin() + i, cycle.begin() + (j+1)%n);
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
