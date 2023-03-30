#include "utils.hpp"

string INPUT_DIR = "../greedy_tsp/input";
string INSTANCES[] = { "kroa100.tsp", "krob100.tsp" };

string INPUT_DIR_CYCLES = "../greedy_tsp/output";
string RESULT_PREFIX = "min";
string ALGORITHMS[] = { "random", "regret_no_weight" };

void swap_edges(vector<int>& cycle, int i, int j) // `i` i `j` to indeksy wierzchołków w cyklu które rozpoczynają krawędź
{
    int n =cycle.size();
    if(min(i,j) == 0 && max(i,j) == n-1)
    {
        swap(cycle[i], cycle[j]);
    }
    reverse(cycle.begin() + i, cycle.begin() + (j+1)%n);
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


bool move_swap_edges(const vector<vector<int>> & matrix,vector<int>& cycle, bool steepest)
{
    std::default_random_engine rng(std::random_device{}());
    int n = cycle.size();
    vector<vector<int>>moves_edge_swap;
    for(int i=0; i<n-1; i++)
    {
        for(int j=i+1; j<n-1; j++)
        {
            moves_edge_swap.push_back({i,j});
        }
    }
    std::shuffle(moves_edge_swap.begin(), moves_edge_swap.end(), rng); // posorotwane ruchy
    vector<int> best_edges = {-1,-1};
    int best_value = 1e9;
    for(vector<int> v : moves_edge_swap)
    {
        int delta = delta_swap_edges(matrix, cycle,v[0], v[1]);
        if( delta < 0)
        {
            if (steepest == false)
            {
                swap_edges(cycle, v[0], v[1]);
                return true;
            }
            else
            {
                if(best_value > delta)
                {
                    best_value = delta;
                    best_edges[0] = v[0];
                    best_edges[1] = v[1];
                }
            }
        }
    }
    if(best_value < 0)
    {
        swap_edges(cycle, best_edges[0], best_edges[1]);
        return true;
    }
    return false;
}
bool swap_vertexes(vector<vector<int>> cycles, bool steepest ) {
    /*
    A function to swap vertices between cycles.
    */
   return false;

}

bool swap_vertexes(vector<int> cycle, bool steepest) {
    /*
    A function to swap vertices inside the cycle.
    */
   return false;
}

void find_best_move(vector<vector<int>> cycles) {

    for(int vertex1 : cycles[0]) {
        for (int vertex2 : cycles[1]) {

        }
    }
}


int steepest_local_search(
    vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    const bool swap_edges
) {
    int solution_delta = 0;
    int curr_delta = 0;

    while(true) {
        // zamiana wierzchołków między cyklami

        if (swap_edges) {
            // zamiana krawędzi wewnątrz cykli
        } else {
            // zamiana wierzchołków wewnątrz cykli
        }


        if (curr_delta > solution_delta) {
            solution_delta = curr_delta;
        } else {
            break;
        }
    }

    return solution_delta;
}

void greedy_local_search(
    vector<vector<int>>& solution,
    const vector<vector<int>>& matrix,
    const bool flag_swap_edges
)
{
    
    while(true)
    {
        bool flag_end = true;
        int cycle_index = rand()%2; // od którego cylu zaczniemy
        int edges_or_vertex = rand()%2; //zaczynbamy od zamiany krawędzi albo wierzchołków
        int outside_insde = rand()%2;
        if(flag_swap_edges == true) // mamy zamianę krawędzi
        {
            if(edges_or_vertex == true)
            {
                bool flag_change = move_swap_edges(matrix, solution[cycle_index],false);
                if(flag_change == false)
                {
                    flag_change = move_swap_edges(matrix, solution[1-cycle_index],false);
                }
                if(flag_change == false)
                {
                    bool flag_change = swap_vertexes(solution, false);
                }
                if(flag_change == false)
                {
                    return;
                }
            }
            else
            {
                bool flag_change = swap_vertexes(solution, false);
                if(flag_change == false)
                {
                    flag_change = move_swap_edges(matrix, solution[cycle_index],false);
                }
                if(flag_change == false)
                {
                    bool flag_change = move_swap_edges(matrix, solution[1-cycle_index],false);
                }
                if(flag_change == false)
                {
                    return;
                }
            }
        }
        else
        {
            if(edges_or_vertex == true)
            {
                bool flag_change = swap_vertexes( solution[cycle_index],false);
                if(flag_change == false)
                {
                    flag_change = swap_vertexes(solution[1-cycle_index],false);
                }
                if(flag_change == false)
                {
                    bool flag_change = swap_vertexes(solution, false);
                }
                if(flag_change == false) // nie można wykonać żadnego ruchu poprawiającego wynik
                {
                    return;
                }
            }
            else
            {
                bool flag_change = swap_vertexes(solution, false);
                if(flag_change == false)
                {
                    flag_change = swap_vertexes(solution[cycle_index],false);
                }
                if(flag_change == false)
                {
                   flag_change = swap_vertexes(solution[1-cycle_index],false);
                }
                if(flag_change == false) // nie można wykonać żadnego ruchu poprawiającego wynik
                {
                    return;
                }
            }
        }
    }
}


void evaluation_algorithm(
    const string& algorithm,
    const string& instance_name,
    const vector<vector<int>>& matrix,
    const vector<vector<int>>& cycles,
    const bool swap_edges=false, // if set to ``true``, then the algorithm should swap the edges instead of vertices inside the cycle
    int n=100
)
{
	for (int i = 0; i < n; i++)
	{
		vector<vector<int>> cycles_to_modify(cycles);

        if (algorithm == "greedy")
        {
            //cout<<score_cycle(matrix,cycles_to_modify[0])+score_cycle(matrix,cycles_to_modify[1])<<"<-  wejście \n";
            greedy_local_search(cycles_to_modify, matrix, true);
            //cout<<score_cycle(matrix,cycles_to_modify[0])+score_cycle(matrix,cycles_to_modify[1])<<"<-wyjście \n";

        }
        else if (algorithm == "steepest")
        {
            steepest_local_search(cycles_to_modify, matrix, swap_edges); //TODO
        }
	}
    // TODO
}


int main() {
    srand(time(nullptr));
    for (string instance : INSTANCES)
	{
		int n;
		vector<vector<int>> v;
		read_file(INPUT_DIR + "/" + instance, &n, v);
		vector<vector<int>> matrix(n+1);
		for (int i = 0; i < matrix.size(); i++)
		{
			matrix[i] = vector<int>(n+1);
		}
		make_distance_matrix(v, matrix);

        for(string cycle_algo : ALGORITHMS)
        {
            vector<vector<int>> cycles(2);

            read_file(INPUT_DIR_CYCLES + "/" + RESULT_PREFIX + "_" + cycle_algo + "_" + instance, cycles);

            evaluation_algorithm("greedy", instance, matrix, cycles);
            evaluation_algorithm("greedy", instance, matrix, cycles, true);

            evaluation_algorithm("steepest", instance, matrix, cycles);
            evaluation_algorithm("steepest", instance, matrix, cycles, true);
            cout << "-----------------------------------" << "\n";
        }
	}
	return 0;
}
