#include <vector>
#include <algorithm>
#include <tuple>
#include <set>
#include <map>
#include <random>
#include <numeric>
#include<cassert>
#include <climits>
#include "utils.hpp"

using namespace std;

void random_solution(const vector<vector<int>>& matrix, vector<vector<int>>& cycles)
{
    vector<int> remaining;
	for(int i = 0; i < matrix.size(); i++) {
        remaining.push_back(i);
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle (remaining.begin(), remaining.end(), std::default_random_engine(seed));

    int middle_index = ceil(remaining.size() / 2.0);
    if (cycles[0].size() == 0) {
        copy(remaining.begin(), remaining.begin() + middle_index, back_inserter(cycles[0]));
        copy(remaining.begin() + middle_index, remaining.end(), back_inserter(cycles[1]));
    } else {
        copy(remaining.begin(), remaining.begin() + middle_index, cycles[0].begin());
        copy(remaining.begin() + middle_index, remaining.end(), cycles[1].begin());
    }
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

int greedy_local_search(
    vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    const bool flag_swap_edges
)
{
    int solution_delta = 0;

    while(true)
    {
        bool flag_end = true;
        int cycle_index = rand()%2;     // od którego cylu zaczniemy
        int edges_or_vertex = rand()%2; // zaczynbamy od zamiany krawędzi albo wierzchołków
        int outside_insde = rand()%2;

        if(flag_swap_edges == true)    // mamy zamianę krawędzi
        {
            if(edges_or_vertex == true)
            {
                // 1. zamiana krawedzi w obu cyklach
                // 2. zamiana wierzcholkow miedzy cyklami
                int delta = move_swap_edges(matrix, solution[cycle_index],false);

                if(!delta)
                {
                    delta = move_swap_edges(matrix, solution[1-cycle_index],false);
                }

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
            else
            {
                // 1. zamiana wierzcholkow miedzy cyklami
                // 2. zamiana krawedzi w obu cyklach
                int delta = move_swap_vertexes(matrix, solution);

                if(!delta)
                {
                    delta = move_swap_edges(matrix, solution[cycle_index],false);
                }

                if(!delta)
                {
                    delta = move_swap_edges(matrix, solution[1-cycle_index],false);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
        }
        else
        {
            if(edges_or_vertex == true)
            {
                // 1. zamiana wierzcholkow w obu cyklach
                // 2. zamiana wierzcholkow miedzy cyklami
                int delta = move_swap_vertexes(matrix, solution[cycle_index]);

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution[1-cycle_index]);
                }

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
            else
            {
                // 1. zamiana wierzcholkow miedzy cyklami
                // 2. zamiana wierzcholkow w obu cyklach
                int delta = move_swap_vertexes(matrix, solution);

                if(!delta)
                {
                    delta = move_swap_vertexes(matrix, solution[cycle_index]);
                }

                if(!delta)
                {
                   delta = move_swap_vertexes(matrix, solution[1-cycle_index]);
                }

                if(!delta)
                {
                    return solution_delta;
                }
                solution_delta += delta;
            }
        }
    }
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
void make_random_moves(int n, vector<vector<int>> &solution )
{
     auto now = std::chrono::system_clock::now();
    auto seed = now.time_since_epoch().count();
    default_random_engine rng = default_random_engine(seed);
    std::uniform_int_distribution<int> dist1(0, 1);
    std::uniform_int_distribution<int> distn(0, solution[0].size()-1);
    for(int i=0; i<n; i++)
    {
        if(dist1(rng) ==0) // zamian wierzchołków między cyklami
        {
            int i = distn(rng);
            int j = distn(rng);
            swap_vertexes(solution, i, j);
        }
        else //zamian krawędzi
        {
            std::uniform_int_distribution<int> dist(0, solution[0].size()-1);
            int i = distn(rng);
            int j = distn(rng);
            int cycle_index = dist1(rng);
            if( i != j)
            {
                swap_edges(solution[cycle_index], i,j);
            }
        }
    }
}

void ILS1(vector<vector<int>> &solution,const vector<vector<int>> &matrix, int & iteration, int time_limit)
{
    // iteration - zwracam liczbę iteracji
    // zakładam że w solution jest rozwiązanie losowe
    auto start = std::chrono::high_resolution_clock::now();
    int score_solution = score_cycle(matrix, solution[0]) +  score_cycle(matrix, solution[1]);
    //cout<<"score: (przed) "<<score_solution <<"\n";
    candidates_algorithm(solution, matrix);

    while(true)
    {
        iteration++;
        vector<vector<int>> tmp_solution(solution.size());
        copy(solution.begin(), solution.end(), tmp_solution.begin());
        make_random_moves(10,tmp_solution);
        score_solution = score_cycle(matrix, solution[0]) +  score_cycle(matrix, solution[1]);
        //regrest_heuristics(matrix, tmp_solution);
        candidates_algorithm(tmp_solution, matrix);
        int score_tmp_solution = score_cycle(matrix, tmp_solution[0]) +  score_cycle(matrix, tmp_solution[1]);
        if (score_tmp_solution < score_solution)
        {
            copy(tmp_solution.begin(), tmp_solution.end(), solution.begin());
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(elapsed_time >= time_limit)
        {
            break;
        }
    }
}
void make_big_perturbation(vector<vector<int>> &solution, set<int> & deleted, float moves_fraction = 0.2)
{
    int vertex_to_delete = moves_fraction*solution[0].size();
    auto now = std::chrono::system_clock::now();
    auto seed = now.time_since_epoch().count();
    default_random_engine rng = default_random_engine(seed);
    std::uniform_int_distribution<int> dist1(0, 1);
    for(int i=0; i<vertex_to_delete; i=i+2) // dla łatwości implementacji usuwamy z każdego cyklu po tyle samo wierzchołków
    {
        int cycle_index = dist1(rng);
        std::uniform_int_distribution<int> dist(0, solution[cycle_index].size()-1);
        int vertex_index = dist(rng);
        deleted.insert(solution[cycle_index][vertex_index]);
        solution[cycle_index].erase(solution[cycle_index].begin() + vertex_index, solution[cycle_index].begin()+vertex_index+1);
        vertex_index = dist(rng);
        cycle_index = 1-cycle_index;
        deleted.insert(solution[cycle_index][vertex_index]);
        solution[cycle_index].erase(solution[cycle_index].begin() + vertex_index, solution[cycle_index].begin()+vertex_index+1);
    }
}
void init_regret2(vector<vector<int>> &cycles,const vector<vector<int>> &matrix, set<int> &remaining, float weight = 1.37, int n = 200)
{

    while (!remaining.empty())
    {
        for (int i = 0; i < 2; i++)
        {
            if(cycles[i].size() == n) //pełen cykl
            {
                continue;
            }
            int best_index, best_vertex;
            tie(best_index, best_vertex) = regret2(matrix, remaining, cycles[i], weight);
            cycles[i].insert(cycles[i].begin() + best_index, best_vertex);
            remaining.erase(best_vertex);
        }
    }
    assert(cycles[0].size() == cycles[1].size());
    return;
}
void ILS2(vector<vector<int>> &solution, const vector<vector<int>> &matrix, int & iteration, int time_limit, bool with_local = false)
{
// iteration - zwracam liczbę iteracji
    // zakładam że w solution jest rozwiązanie losowe
    auto start = std::chrono::high_resolution_clock::now();
    candidates_algorithm(solution, matrix);
    while(true)
    {
        iteration++;
        vector<vector<int>> tmp_solution(solution.size());
        copy(solution.begin(), solution.end(), tmp_solution.begin());
        set <int> deleted;
        make_big_perturbation(tmp_solution, deleted);
        init_regret2(tmp_solution, matrix, deleted);
        int score_solution = score_cycle(matrix, solution[0]) +  score_cycle(matrix, solution[1]);
        if(with_local)
        {
            candidates_algorithm(tmp_solution, matrix);
        }
        int score_tmp_solution = score_cycle(matrix, tmp_solution[0]) +  score_cycle(matrix, tmp_solution[1]);
        if (score_tmp_solution < score_solution)
        {
            copy(tmp_solution.begin(), tmp_solution.end(), solution.begin());
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(elapsed_time >= time_limit)
        {
            break;
        }
    }
}


void MSLS(vector<vector<int>> &solution, const vector<vector<int>> &matrix, int iterations=100)
{
    vector<vector<vector<int>>> solution_candidates(iterations);
    for(int i = 0; i < iterations; ++i) {
        solution_candidates[i] = vector<vector<int>>(solution);
        random_solution(matrix, solution_candidates[i]);
    }

    int best_solution_index = 0;
    int best_solution_score = INT_MAX;
    int curr_score = INT_MAX;

    for(int i = 0; i < iterations; ++i){
        // cout << (i + 1) << "... ";
        candidates_algorithm(solution_candidates[i], matrix);

        curr_score = score_cycle(matrix, solution_candidates[i][0]) +  score_cycle(matrix, solution_candidates[i][1]);
        if (curr_score < best_solution_score) {
            best_solution_index = i;
            best_solution_score = curr_score;
        }
    }
    solution = solution_candidates[best_solution_index];
    // cout << endl;
}


void crossover(const vector<vector<int>>& matrix, const vector<vector<int>> &parent_solution, const vector<vector<int>> &solution, vector<vector<int>> &new_solution)
{
    vector<vector<int>>parent_solution_copy(2);
    copy(parent_solution.begin(), parent_solution.end(), parent_solution_copy.begin());
    int solution_size = parent_solution_copy[0].size(); // 100 wierzchołków
    vector<pair<int,int>> deleted_edges;
    vector<int> deleted_vertex;
    vector<vector<int>> counter(2);

    vector<vector<bool>> has_egde(200);
    for (int i = 0; i < 200; ++i){
        has_egde[i] = vector<bool>(200, false);
    }

    // usuwanie kawędzi (?)
    for(int i=0;i<2; i++)
    {
        int n = solution[i].size();
        for(int j=0; j<n; j++)
        {
            counter[i].push_back(0);
            int e1 = min(solution[i][j], solution[i][(j+1)%n]);
            int e2 = max(solution[i][j], solution[i][(j+1)%n]);
            // deleted_edges.push_back({e1,e2});

            has_egde[e1][e2] = true;
            has_egde[e2][e1] = true;
        }
    }
    // dla każdego wierzchołka oznaczamy ile krawędzi z nim związanych zostaje usunięte
    for(int i=0;i<2; i++)
    {
        int n = parent_solution_copy[i].size();
        for(int j=0; j<n; j++)
        {
            int e1 = min(parent_solution_copy[i][j], parent_solution_copy[i][(j+1)%n]);
            int e2 = max(parent_solution_copy[i][j], parent_solution_copy[i][(j+1)%n]);
            // if(find(deleted_edges.begin(), deleted_edges.end(), make_pair(e1,e2)) == deleted_edges.end())
            if(!has_egde[e1][e2])
            {
                counter[i][j] -=1;
                counter[i][(j + 1) % n] -= 1;
                if(counter[i][j] == -2)
                {
                    deleted_vertex.push_back(parent_solution_copy[i][j]);// zbieramy wierzchołki które zostaną usunięte
                }
            }
        }
    }
    for(auto v : deleted_vertex)
    {
        for(int i=0;i<2; i++)
        {
            for(int j=0; j<parent_solution_copy[i].size(); j++)
            {
                if(parent_solution_copy[i][j] == v)
                {
                    parent_solution_copy[i].erase(parent_solution_copy[i].begin()+j, parent_solution_copy[i].begin()+j+1);
                }
            }
        }
    }
    set<int> remaining(deleted_vertex.begin(), deleted_vertex.end());

    init_regret2(parent_solution_copy, matrix,remaining, 1.37,solution_size);
    make_random_moves(3, parent_solution_copy);

    new_solution = vector<vector<int>>(parent_solution_copy);

    // copy(parent_solution_copy.begin(), parent_solution_copy.end(), new_solution.begin());
   // cout << score_cycle(matrix, parent_solution[0]) + score_cycle(matrix, parent_solution[1]) - score_cycle(matrix, new_solution[0]) - score_cycle(matrix, new_solution[0]);
}


class Solution {
    public:
        vector<vector<int>> solution;
        int score;

    Solution() {
        this->solution = vector<vector<int>>(2);
    }
    Solution(vector<vector<int>> solution, int score) {
        this->solution = solution;
        this->score = score;
    }

    bool operator==(const int score) {
        return this->score == score;
    }
};


void genetic(vector<vector<int>> &solution, const vector<vector<int>> &matrix, int& iteration, int time_limit, const bool use_local = false, const int population_size = 20)
{
    vector<Solution> population(population_size);
    for (int i = 0; i < population.size(); ++i) {
        do {
            random_solution(matrix, population[i].solution);
            candidates_algorithm(population[i].solution, matrix);
            population[i].score = score_cycles(matrix, population[i].solution);
        } while(find(population.begin(), population.begin() + i, population[i].score) < population.begin() + i);
    }

    int words_solution_index = max_element(population.begin(), population.end(), [](const Solution& a, const Solution& b){ return a.score < b.score; }) - population.begin();
    // cout << "population generated" << endl;
    int sol1, sol2;
    iteration = 0;

    auto start = std::chrono::high_resolution_clock::now();
    while(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() < time_limit)
    {
        ++iteration;
        // cout << "iteration: " << iteration << endl;
        sol1 = rand() % population_size;
        do {
            sol2 = rand() % population_size;
        } while(sol1 == sol2);

        // cout << "Chosen: " << sol1 << " " << sol2 <<endl;

        vector<vector<int>> new_solution;
        crossover(matrix, population[sol1].solution, population[sol2].solution, new_solution);

        if(use_local) {
            candidates_algorithm(new_solution, matrix);
        }
        int new_score = score_cycles(matrix, new_solution);

        if (new_score < population[words_solution_index].score && find(population.begin(), population.end(), new_score) == population.end()) {
            population[words_solution_index] = Solution(new_solution, new_score);
            // cout << "New score was added: " << new_score << endl;
            words_solution_index = max_element(population.begin(), population.end(), [](const Solution& a, const Solution& b){ return a.score < b.score; }) - population.begin();
        }
    }

    int best_solution_index = min_element(population.begin(), population.end(), [](const Solution& a, const Solution& b){ return a.score < b.score; }) - population.begin();
    solution = population[best_solution_index].solution;
    cout << population[best_solution_index].score << endl;
}


