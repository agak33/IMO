#include <vector>
#include <algorithm>
#include <tuple>
#include <set>
#include <map>
#include <random>
#include <numeric>
#include<cassert>

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

void swap_edges(vector<int> &cycle, int i, int j) // `i` i `j` to indeksy wierzchołków w cyklu które rozpoczynają krawędź
{
    int ii = min(i, j);
    int jj = max(i, j);
    int n = cycle.size();
    ii = (ii+1)%n;
    int d = (abs(jj-ii)+n)%n;
    for(int k=0;k<(d/2) +1; k++)
    {
        int a =(ii+k)%n;
        int b = (ii+d-k+n)%n;
        swap(cycle[a], cycle[b]);
    }
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

int delta_swap_edges(const vector<vector<int>> & matrix,const vector<int>& cycle, int i, int j)// `i` i `j` to indeksy wierzchołków w cyklu które rozpoczynają krawędź
{
    int n = cycle.size();
    int a, b, c, d;
    a = cycle[i];
    b = cycle[(i+1)%n];
    c = cycle[j%n];
    d = cycle[(j+1)%n];
    return matrix[a][c] + matrix[b][d] - matrix[a][b] - matrix[c][d];
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

void make_candidate_matrix(const vector<vector<int>> &matrix, vector<vector<int>> &candidate_matrix, int n = 10)
{
    for (int i = 0; i < candidate_matrix.size(); i++)
    {
        vector<pair<int, int>> indexed_vec(matrix[i].size());
        for (int j = 0; j < matrix[i].size(); j++)
        {
            indexed_vec[j] = make_pair(matrix[i][j], j);
        }
        indexed_vec[i].first = 1e9; // nie bierzemy pod uwagę krawędzi do siebie
        sort(indexed_vec.begin(), indexed_vec.end());
        for (int j = 0; j < n; j++)
        {
            candidate_matrix[i].push_back(indexed_vec[j].second);
        }
    }
}

int find_vertex(const vector<vector<int>> &solution, int v)
{
    if(find(solution[0].begin(), solution[0].end(),v) != solution[0].end())
    {
        return 0;
    }
    return 1;
}

int find_index(const vector<int> solution, int v)
{
    for(int i=0;i<solution.size(); i++)
    {
        if(solution[i] == v)
        {
            return i;
        }
    }
    return -1;
}

void candidates_algorithm( vector<vector<int>> &solution,const vector<vector<int>> &matrix, int k=10)
{
    vector<vector<int>> candidates(matrix.size());
    make_candidate_matrix(matrix, candidates, k);
    while(true)
    {
      //  cout<<score_cycle(matrix,solution[0])+score_cycle(matrix,solution[1])<<"\n";
        pair<int,int> best_move; //zawiera krawędź pomiędzy dwoma wierzchołkami którą chcemy dodać
        int delta_best_move = 1e9; //ocena najlepszego ruchu
        int type_best_move = -1;
        int n = solution[0].size();
        int best_cycle = -1;
        pair<int, int> tmp_v = { -1,-1 };
        for(int i=0;i < candidates.size(); i++)
        {
            for(int j=0; j<candidates[i].size(); j++)
            {
                // dla każdej krawędzi kandydackiej będziemy próbowali ją dodać
                int c1 = find_vertex(solution,i);
                int c2 = find_vertex(solution, candidates[i][j]);
                int v1 = i;
                int v2 = candidates[i][j];
                if (c1 == c2)
                {
                    // wierzchołki są w tym samym cyklu, robimy zamianę krawędzi
                    //znajdumy indeksy wierzhcołków w cyklu
                    int index_1 = find_index(solution[c1], v1);
                    int index_2 = find_index(solution[c2], v2);
                    int local_delta = delta_swap_edges(matrix, solution[c1], index_1, index_2);
                    if(local_delta < delta_best_move)
                    {
                        delta_best_move = local_delta;
                        best_move = make_pair(index_1, index_2);
                        type_best_move = 1;
                        best_cycle = c1;
                    }
                   local_delta = delta_swap_edges(matrix, solution[c1], (index_1 - 1 + n) % n, (index_2 - 1 + n) % n);
                    if(local_delta < delta_best_move)
                    {
                        delta_best_move = local_delta;
                        best_move = make_pair((index_1-1+n)%n, (index_2-1+n)%n);
                        type_best_move = 1;
                        best_cycle = c1;
                    }
                }
                else
                {
                    //wierzchołki są w różnych cyklach, robię zamianę wierzchołków wprowadzając
                    //naszą krawędź
                    int index_1 = find_index(solution[c1], v1);
                    int index_2 = find_index(solution[c2], v2);
                    vector<pair<int, pair<int, int>>> deltas;
                    if (c1 == 0)
                    {
                        // wektor par <ocena, ruch>
                        deltas.push_back({ delta_swap_vertexes(matrix,solution,(index_1 + 1) % n, index_2), {(index_1 + 1) % n, index_2} });
                        deltas.push_back({ delta_swap_vertexes(matrix,solution,index_1,(index_2 - 1 + n) % n), {index_1,(index_2 - 1 + n) % n} });
                    }
                    else
                    {
                        deltas.push_back({ delta_swap_vertexes(matrix,solution,(index_2 - 1+n) % n, index_1), {(index_2 - 1 + n) % n, index_1} });
                        deltas.push_back({ delta_swap_vertexes(matrix,solution,index_2,(index_1 + 1 + n) % n), {index_2,(index_1 + 1 + n) % n } });
                    }
                    sort(deltas.begin(), deltas.end());
                    if(deltas[0].first < delta_best_move)
                    {
                        delta_best_move = deltas[0].first;
                        best_move = deltas[0].second;
                        type_best_move = 2;
                        tmp_v = { v1,v2 };
                    }
                }
            }
        }
        if(type_best_move > 0 && delta_best_move < 0)
        {
           
            if(type_best_move == 1)
            {
                assert(delta_swap_edges(matrix, solution[best_cycle], best_move.first, best_move.second) < 0 && "błąd delty krawedzie");
                swap_edges(solution[best_cycle], best_move.first, best_move.second);
            }
            else
            {
                int i1 = best_move.first;
                int i2 = best_move.second;
                assert(delta_swap_vertexes(matrix,solution, i1, i2) < 0 && "błąd delty");
                swap_vertexes(solution, i1, i2);
                /*
                v1 = tmp_v.first;
                v2 = tmp_v.second;
                int c1 = find_vertex(solution, v1);
                int c2 = find_vertex(solution, v2);
                assert(c1 == c2 && "błąd cykle");
                i1 = find_index(solution[c1], v1);
                i2 = find_index(solution[c2], v2);
                assert((abs(i1 - i2)== 1 || abs(i1 - i2) == 99) && "błąd indeksy");
                assert(solution[c1][i1] == v1 && "błąd wierzchołek");
                assert(solution[c2][i2] == v2 && "błąd wierzchołek");
                */

            }
        }
        else
        {
            break;
        }
    }
    return;
}


