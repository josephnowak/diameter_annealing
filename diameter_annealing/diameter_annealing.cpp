#include <bits/stdc++.h>
#include <random>
#include <math.h>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm> 
#include <chrono> 
#include <map>
#include <iostream> 
#include <assert.h>
using namespace std;
using namespace std::chrono;

#define um unordered_map
#define us unordered_set
//default_random_engine generator;
random_device rd;
mt19937 generator(rd());
uniform_real_distribution<double> uniform_real(0, 1);


/*
TODO:
1. Document all the different attributes and the general logic of every method
2. Add an option to select random paths instead of use the dijkstra
*/



struct Comp
{
  bool operator()(const tuple<int, int, int> &x, const tuple<int, int, int> &y) const
  {
  	if(get<0>(x) == get<0>(y)){
  		return get<1>(x) < get<1>(y);
  	}
    return get<0>(x) > get<0>(y);
  }
};


class BaseRouting {
	/**
	 *
	 * @brief Abstract class used on the simulated annealing class/algorithm to control how to find neighborhood_solution
	 */
	public:
		int sample_size;
		/**
		 * This is equivalent to calculate the R[node].size() * charge, is only useful to simplify the use of symmetric paths
		 *
		 * @brief Charge of everynode on the routing or how many paths use this node as internal 
		 */
		vector<int> charge_per_node;

		/**
		 * The data structure is equivalent to adjacency list representation but using arrays/vectors
		 * @brief Graph that represent the connection between every pair of nodes
		 */
		vector<vector<int>> g;

		/**
		 * This is equivalent to the sum of the charge of every internal node on the routing
		 *
		 * @brief Charge of every path
		 */
		vector<int> total_charge_per_path;

		/**
		 * This can be implemented in multiple ways
		 * @brief Find a neighborhood solution
		 */
		virtual void neighborhood_solution() = 0;
		
		/**
		 * It is useful for expensive computations like update the distances using a dijkstra algorithm
		 * @brief Method executed before reduce the temperature on the simulated annealing algorithm
		 */
		virtual void update_routing() = 0;

		/**
		 * This method is called everytime the simulated annealing algorithm decide to not use the new neighborhood solution
		 *
		 * @brief Method used to forget the new solution. 
		 */
		virtual void forget_neighborhood_solution() = 0;

		/**
		 * This method is called everytime the simulated annealing algorithm decide to not use the new neighborhood solution
		 *
		 * @brief Method used to forget the new solution. 
		 */
		virtual vector<vector<vector<int>>> get_standard_format_routing() = 0;


		/**
		 *
		 * @brief Return a pair with the forwarding diameter (1) and any other value that is going to break a tie (2)
		 */
		virtual pair<int, int> evaluate() = 0;

		bool validate_routing(vector<vector<vector<int>>> &routing, int expected_forwarding_diameter){
			
			int n = this->g.size();
			vector<int> charge_per_node = vector<int>(n, 0);
			vector<vector<int>> charge_per_path = vector<vector<int>>(n, vector<int>(n, 0));
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					int prev_node = i;
					bool incorrect = false;
					for(int node: routing[i][j]){
						if(std::find(g[prev_node].begin(), g[prev_node].end(), node) == g[prev_node].end()) {
							incorrect = true;
						}
						charge_per_node[node] += 1;
						prev_node = node;
					}
					if(incorrect){
						cout << "The path between the nodes: " << i << " and " << j << " is incorrect" << endl;
						for(int node: routing[i][j])
							cout << node << " ";
						cout << endl;
						return false;
					}
				}
			}
			int forwarding_diameter = 0;
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					for(int node: routing[i][j]){
						charge_per_path[i][j] += charge_per_node[node];
					}
					forwarding_diameter = max(forwarding_diameter, charge_per_path[i][j]);
				}
				
			}

			if(forwarding_diameter != expected_forwarding_diameter){
				cout << "Routing forwarding diameter " << forwarding_diameter << ", expected: " << expected_forwarding_diameter << endl;
				return false;
			}
			return true;
		}

};


class Routing: virtual public BaseRouting {
	/**
	 * Implementation of the BaseRounting using symmetric paths and dijkstra to find neighborhood solutions.
	 *
	 * @brief Routing
	 */
	private:
		/**
		 * @brief Ids of all the path, they are generated as autoincrement indexes.
		 */
		unordered_set<int> paths_ids;

		/**
		 * This is useful for speed up the recalculation of the forwarding diameter and is also useful
		 * on the forget_neighborhood_solution method.
		 * @brief Paths added during the execution of the neighborhood_solution method.
		 */
		vector<unordered_map<int, int>> added_paths;
		/**
		 * This is useful for speed up the recalculation of the forwarding diameter and is also useful
		 * on the forget_neighborhood_solution method.
		 * @brief Paths deleted during the execution of the neighborhood_solution method.
		 */
		vector<unordered_map<int, int>> erased_paths;

		/**
		 * @brief Vector that contains all the paths that are going to be use to modify the actual routing 
		 */
		vector<int> selected_paths;

		/**
		 * @brief Dijkstra distance of every path
		 */
		vector<int> dijkstra_path_distances;

		/**
		 * Useful only for symmetric paths
		 * @brief Represent the charge generated by every use of a node as a internal path
		 */
		int charge;

		/**
		 * @brief Number of node on the graph
		 */
		int num_nodes;

		void generate_initial_routing()
		{
			for(int i = 0; i < this->num_nodes; i++)
				get_random_paths(i);
			
			for(int i = 0; i < this->num_nodes; i++){
				this->charge_per_node[i] = this->R[i].size() * this->charge;
				for(const pair<int, int> path_map: this->R[i]){
					this->total_charge_per_path[path_map.first] += this->charge_per_node[i];	
					this->paths_ids.insert(path_map.first);
				}
			}
			this->update_routing();
			
		}
		vector<int> build_path(vector<int> &prev, int x, int y)
		{

			vector<int> path;
			int act = y;
			while(prev[act] != act)
			{
				if(prev[act] != x)
					path.push_back(prev[act]);
				act = prev[act];
			}
			reverse(path.begin(), path.end());
			return path;
		}

		pair<vector<int>, vector<int>> dijkstra(int x, int y)
		{
			set<tuple<int,float, int>> q;
			q.insert(make_tuple(0, uniform_real(generator), x));
			vector<int> dist_per_node(this->num_nodes, INT_MAX);
			vector<int> prev(this->num_nodes);
			prev[x] = x;
			dist_per_node[x] = 0;
			while(not q.empty())
			{
				tuple<int, int, int> temp = *q.begin();
				int dist = get<0>(temp), random_value = get<1>(temp), node = get<2>(temp);
				if(dist < 0)
					break;
				dist_per_node[node] = dist;
				q.erase(q.begin());
				if(node == y) 
					break;

				for(int &adj: this->g[node])
				{
					if(dist_per_node[adj] > dist + this->charge_per_node[adj] + 2)
					{
						prev[adj] = node;
						q.erase(make_tuple(dist_per_node[adj], random_value, adj));
						dist_per_node[adj] = dist + this->charge_per_node[adj] + 2;
						q.insert(make_tuple(dist_per_node[adj], uniform_real(generator), adj));
					}
				}
			}
			return {dist_per_node, prev};
		}
		void modify_routing(vector<unordered_map<int, int>> &new_paths, vector<unordered_map<int, int>> &old_paths)
		{
			unordered_map<int, int> new_paths_charge, changed_charge_paths;

			for(int i = 0; i < this->num_nodes; i++)
			{
				int new_charge = this->charge * (this->R[i].size() + new_paths[i].size() - old_paths[i].size());
				int change_value = this->charge * (new_paths[i].size() - old_paths[i].size());
				// calculate the new charge for every node based on the amount of new 
				this->charge_per_node[i] = new_charge;

				for(const pair<int, int> &id_path_map: old_paths[i])
					this->R[i].erase(id_path_map.first);

				if(change_value != 0)
					for(const pair<int, int> &path_map: this->R[i])
						changed_charge_paths[path_map.first] += change_value;
				
				for(const pair<int, int> &id_path_map: new_paths[i]){
					this->R[i].insert(id_path_map);
					new_paths_charge[id_path_map.first] += new_charge;
				}
			}
			for(const pair<int, int> &path_charge: new_paths_charge)
				this->total_charge_per_path[path_charge.first] = path_charge.second;

			for(const pair<int, int> &path_charge: changed_charge_paths)
				this->total_charge_per_path[path_charge.first] += path_charge.second;
			
		}
		void get_random_paths(int source)
		{
			queue<pair<int, int>> q;
			q.push({source, 0});
			vector<int> distances(this->num_nodes, -1);
			while(!q.empty())
			{
				pair<int, int> temp = q.front();
				int node = temp.first;
				int len = temp.second;
				q.pop();
				if(distances[node] != -1) continue;
				distances[node] = len;
				for(int v: this->g[node])
				{
					if(distances[v] != -1) continue;
					q.push({v, len + 1});
				}
			}
			for(int i = source + 1; i < this->num_nodes; i++)
			{
				int best_distance = distances[i];
				int act_node = i;
				int path_id = source * this->num_nodes + i;
				for(int dist = best_distance; dist > 1; dist--)
				{
					vector<int> valid;
					for(int &node: this->g[act_node])
						if(distances[node] == dist - 1)
							valid.push_back(node);
					uniform_int_distribution<int> uniform_int_path(0, valid.size() - 1);
					int pos = uniform_int_path(generator);
					act_node = valid[pos];
					this->R[act_node].insert({path_id, dist - 2});
				}
			}
			assert (1 == 1);
		}

	public:
		/**
		 * Every cell/node of the vector is a set with all the paths (IDs) that use that node as an internal one.
		 * @brief Routing of the network
		 */
		vector<unordered_map<int, int>> R;

		
		/**
		 * This attribute is useful to restrict the randomness of the algorithm, an smaller integer (call it N)
		 * means that the algorithm is only going to take into account the top N paths based on the three metrics
		 * described on the neighborhood_solution method, this parameter is combined with the sample_size
		 * 
		 * @brief Determine the number of best possible paths to select randonmly
		 */
		int universe_size;

		/**
		 * Once the numbre of paths is restricted based on the universe_size, they are sorted randomly and then M of them
		 * are selected (this m is the sample_size)
		 * @brief Number of paths that are going to be modified every time the neighborhood_solution method is called
		 */
		int sample_size;
	
		Routing(){}

		Routing(vector<vector<int>> &g, int sample_size, int universe_size)
		{
			this->g = g;
			this->num_nodes = g.size();
			this->R = vector<unordered_map<int, int>>(this->num_nodes);
			this->charge_per_node = vector<int>(this->num_nodes);
			this->added_paths = vector<unordered_map<int, int>>(this->num_nodes, unordered_map<int, int>());
			this->erased_paths = vector<unordered_map<int, int>>(this->num_nodes, unordered_map<int, int>());
			this->total_charge_per_path = vector<int>(this->num_nodes * this->num_nodes, 0);
			this->dijkstra_path_distances = vector<int>(this->num_nodes * this->num_nodes, 0);
			this->charge = 2;
			this->universe_size = universe_size;
			this->sample_size = sample_size;
			this->selected_paths = vector<int>(this->universe_size);

			this->generate_initial_routing();
		}
		Routing(const Routing &a)
		{
			this->num_nodes = a.num_nodes;
			this->R = a.R;
			this->charge_per_node = a.charge_per_node;
			this->g = a.g;
			this->total_charge_per_path = a.total_charge_per_path;
			this->paths_ids = a.paths_ids;
		}
		void neighborhood_solution()
		{
			for(int i = 0; i < this->num_nodes; i++)
			{
				this->erased_paths[i].clear();
				this->added_paths[i].clear();
			}

			set<tuple<int, int, int>, Comp> ordered_charge_per_path;
			
			for(int id_path: this->paths_ids)
			{
				int x = id_path / this->num_nodes, y = id_path % this->num_nodes;
				/*
					Every path is sorted based on three elements:
					1. The charge, which is the sum of the charge of the nodes used as internal nodes on the path
					2. The dijkstra distance, explained in the update_routing method
					3. The id of the path, used to break ties

				*/
				ordered_charge_per_path.insert(
					make_tuple(this->total_charge_per_path[id_path], 
					this->dijkstra_path_distances[id_path], 
					id_path
				));

				// preserve only universe_size elements
				if(ordered_charge_per_path.size() > this->universe_size){
					ordered_charge_per_path.erase(prev(ordered_charge_per_path.end()));
				}
			}
			
			// Select and sort randomly the paths on the ordered_charge_per_path 
			int i = 0;
			for(const tuple<int, int, int> &path: ordered_charge_per_path){
				this->selected_paths[i++] = get<2>(path);
			}
			shuffle(this->selected_paths.begin(), this->selected_paths.end(), generator);

			for(int i = 0; i < min(this->sample_size, (int)ordered_charge_per_path.size()); i++)
			{
				int id_path = this->selected_paths[i];

				// remove the previous path first to improve the results of the dijkstra
				for(int j = 0; j < this->num_nodes; j++)
				{
					// verify if the path was using the jth node
					auto r_element = this->R[j].find(id_path);
					if(r_element == this->R[j].end())
						continue;
					// if it was using the node then reduce the charge of that node (only useful for the dijkstra) 
					// and add it to the erased_paths set
					this->charge_per_node[j] -= this->charge;
					this->erased_paths[j].insert(*r_element);
				}

				int x = id_path / this->num_nodes, y = id_path % this->num_nodes;
				// get the dijkstra path from the node X to Y
				vector<int> prev_path_node = this->dijkstra(x, y).second;
				vector<int> path = this->build_path(prev_path_node, x, y);
				for(int i = 0; i < path.size(); i++)
				{
					// for every node on the path increment the charge (only useful for the dijkstra) 
					// and add it to the added_paths set
					int v = path[i];
					this->charge_per_node[v] += this->charge;
					this->added_paths[v][id_path] = i;
				}	
				
			}
			// modify the routing based on the added and erased paths
			this->modify_routing(this->added_paths, this->erased_paths);

		}
		void update_routing()
		{
			for(int i = 0; i < this->num_nodes; i++){
				vector<int> distances_dijkstra = this->dijkstra(i, -1).first;
				for(int j = 0; j < this->num_nodes; j++)
					this->dijkstra_path_distances[i * this->num_nodes + j] = distances_dijkstra[j];
			}
		}
		void forget_neighborhood_solution()
		{			
			this->modify_routing(this->erased_paths, this->added_paths);
		}
		pair<int, int> evaluate()
		{
			return {
				*max_element(this->total_charge_per_path.begin(), this->total_charge_per_path.end()), 
				*max_element(this->charge_per_node.begin(), this->charge_per_node.end())
			};
		}

		vector<vector<vector<int>>> get_standard_format_routing(){
			
			vector<vector<vector<int>>> routing = vector<vector<vector<int>>>(
				this->num_nodes, 
				vector<vector<int>>(
					this->num_nodes, 
					vector<int>()
				)
			);

			for(int id_path: this->paths_ids)
			{
				int x = id_path / this->num_nodes, y = id_path % this->num_nodes;
				vector<pair<int, int>> elements = vector<pair<int, int>>();

				for(int i = 0; i < this->num_nodes; i++){
					auto r_element = this->R[i].find(id_path);
					if(r_element == this->R[i].end())
						continue;
					elements.push_back({r_element->second, i});
				}
				sort(elements.begin(), elements.end());

				vector<int> path = vector<int>(elements.size());
				for(int i = 0; i < path.size(); i++)
					path[i] = elements[i].second;
				
				routing[x][y] = path;
				vector<int> inverse_path (path.rbegin(), path.rend());
				routing[y][x] = inverse_path;
			}
			return routing;

		}
};

template<typename T>
class SimulatedAnnealing
{
	double temperature;
	double beta; 
	double min_temperature;
	T routing;
	T best_routing;

public:
	SimulatedAnnealing(T &routing, double temperature, double min_temperature, double beta)
	{
		this->temperature = temperature;
		this->routing = routing;
		this->best_routing = routing;
		this->min_temperature = min_temperature;
		this->beta = beta;
	}
	T optimize()
	{
		double temperature = this->temperature;
		pair<int,int> act_cost = this->routing.evaluate();
		pair<int,int> best_cost = act_cost;
		while(temperature > this->min_temperature)
		{
			cout << "temperature: " << temperature << endl;
			auto start_total = chrono::steady_clock::now();
			for(int i = 0; i < 100; i++)
			{
				auto start = chrono::steady_clock::now();
				this->routing.neighborhood_solution();
				auto end = chrono::steady_clock::now();
				// cout << "time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
				pair<int, int> new_cost = this->routing.evaluate();

				// cout << new_cost.first << endl;
				float dif_cost = ((act_cost.first - new_cost.first) / 2) / (float)this->routing.sample_size;
				//  cout << dif_cost << endl;
				double x = uniform_real(generator);
				if(dif_cost > 0)
					act_cost = new_cost;
				else if(dif_cost == 0 and (act_cost.second >= new_cost.second))
					act_cost = new_cost;
				else if(dif_cost < 0 and x < exp(dif_cost/temperature))
					act_cost = new_cost;
				else{
					cout << "act_sol: " << act_cost.first << " new_cost: " << new_cost.first << endl;
 					this->routing.forget_neighborhood_solution();
				}

				if(best_cost.first > new_cost.first  or 
					(new_cost.first == best_cost.first and (best_cost.second > new_cost.second))){
					this->best_routing = this->routing;
					best_cost = new_cost;
				}
				// cout << new_cost.first << endl;
			}
			this->routing.update_routing();
			auto end = chrono::steady_clock::now();
			cout << "total time: " << chrono::duration_cast<chrono::microseconds>(end - start_total).count() << endl;
			cout << this->routing.evaluate().first << endl;
			temperature = this->reduce_temperature(temperature);
		}
		/*for(int i = 0; i < this->best_routing.num_nodes; i++)
			for(int j = i + 1; j < this->best_routing.num_nodes; j++){
				this->best_routing.R[j][i] = this->best_routing.R[i][j];
				reverse(this->best_routing.R[j][i].begin(), this->best_routing.R[j][i].end());
			}*/
		return this->best_routing;
	}
	double reduce_temperature(double temperature)
	{
		return temperature - this->beta;
	}
};


int main()
{
	
	int n, m;
	cin >> n >> m;
	vector<vector<int>> g(n, vector<int>());
	for(int i = 0; i < m; i++)
	{
		int x, y, w;
		cin >> x >> y;
		g[x].push_back(y);
		g[y].push_back(x);
	}

	
	Routing routing(g, 180, 240);

	SimulatedAnnealing<Routing> simulated_annealing(routing, 1, 1e-1, 1e-1);

	auto start = chrono::steady_clock::now();
	Routing sol = simulated_annealing.optimize();
	// Routing sol = optimize(routing, 1, 1e-1, 1e-1);
	auto end = chrono::steady_clock::now();
	cout << "time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
	
	cout << "annealing: " << sol.evaluate().first << endl;
	int sum = 0;
	for(int i = 0; i < n; i++){
		sum += sol.charge_per_node[i];
		cout << sol.charge_per_node[i] << " ";
	}
	cout << endl;

	vector<vector<vector<int>>> result_routing = sol.get_standard_format_routing();

	cout << sol.validate_routing(result_routing, sol.evaluate().first) << endl;
	return 0;
}