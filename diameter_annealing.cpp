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
3. Change the name of total_charge by total_change_per_path
4. Change the name of charge_per_path by ordered_paths.
5. Change the name of prev by prev_path_node, or something more clear.
6. Change the name of the universe_size and sample_size for something more congruent to what they do, which is basically
	 define how much elements of the total paths are going to be taken for the modification of paths (sample size right now)
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


class Routing {
	public:
		vector<unordered_set<int>> R, added_paths, erased_paths;
		vector<int> charge_per_node, path_per_node, selection_paths;
		unordered_set<int> total_paths;
		vector<vector<int>> g;
		vector<int> total_charge, total_distances;
		int charge, neighborhood_len, num_nodes, universe_size, sample_size;

		Routing(){}

		Routing(vector<vector<int>> &g, int sample_size, int universe_size)
		{
			this->g = g;
			this->num_nodes = g.size();
			this->R = vector<unordered_set<int>>(this->num_nodes);
			this->charge_per_node = vector<int>(this->num_nodes);
			this->added_paths = vector<unordered_set<int>>(this->num_nodes, unordered_set<int>());
			this->erased_paths = vector<unordered_set<int>>(this->num_nodes, unordered_set<int>());
			this->total_charge = vector<int>(this->num_nodes * this->num_nodes, 0);
			this->total_distances = vector<int>(this->num_nodes * this->num_nodes, 0);
			this->charge = 2;
			this->generate_initial_routing();
			this->universe_size = universe_size;
			this->sample_size = sample_size;
			this->selection_paths = vector<int>(this->universe_size);
		}
		Routing(const Routing &a)
		{
			this->num_nodes = a.num_nodes;
			this->R = a.R;
			this->charge_per_node = a.charge_per_node;
			this->g = a.g;
			this->total_charge = a.total_charge;
			this->total_paths = a.total_paths;
		}
		void generate_initial_routing()
		{
			for(int i = 0; i < this->num_nodes; i++)
				get_random_paths(i);
			
			for(int i = 0; i < this->num_nodes; i++){
				this->charge_per_node[i] = this->R[i].size() * this->charge;
				for(int path: this->R[i]){
					this->total_charge[path] += this->charge_per_node[i];	
					this->total_paths.insert(path);
				}
			}
			this->update_distances();
			
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
					this->R[act_node].insert(path_id);
				}
			}
		}
		void modify_routing(vector<unordered_set<int>> &new_paths, vector<unordered_set<int>> &old_paths)
		{
			unordered_map<int, int> new_paths_charge, changed_charge_paths;

			for(int i = 0; i < this->num_nodes; i++)
			{
				int new_charge = this->charge * (this->R[i].size() + new_paths[i].size() - old_paths[i].size());
				int change_value = this->charge * (new_paths[i].size() - old_paths[i].size());
				this->charge_per_node[i] = new_charge;

				for(const int &id_path: old_paths[i])
					this->R[i].erase(id_path);

				if(change_value != 0)
					for(int path: this->R[i])
						changed_charge_paths[path] += change_value;
				
				for(const int &id_path: new_paths[i]){
					this->R[i].insert(id_path);
					new_paths_charge[id_path] += new_charge;
				}
			}
			for(const pair<int, int> &path_charge: new_paths_charge)
				this->total_charge[path_charge.first] = path_charge.second;

			for(const pair<int, int> &path_charge: changed_charge_paths)
				this->total_charge[path_charge.first] += path_charge.second;
			
		}
		void neighborhood_solution()
		{
			// cout << "here " << endl;
			auto start = chrono::steady_clock::now();

			for(int i = 0; i < this->num_nodes; i++)
			{
				this->erased_paths[i].clear();
				this->added_paths[i].clear();
			}

			set<tuple<int, int, int>, Comp> charge_per_path;
			
			// cout << endl;
			for(int id_path: this->total_paths)
			{
				//if(this->used_paths.find(id_path) != this->used_paths.end())
				//	continue;
				int x = id_path / this->num_nodes, y = id_path % this->num_nodes;
				// cout << "(" << this->total_charge[id_path] << ", " << this->total_distances[id_path] << ", " << x << ", " << y << ")" << endl;
				
				charge_per_path.insert(make_tuple(this->total_charge[id_path], this->total_distances[id_path], id_path));
				// cout << charge_per_path.size() << endl;
				if(charge_per_path.size() > this->universe_size){
					// cout << "here " << endl;
					charge_per_path.erase(prev(charge_per_path.end()));
				}
			}
			
			// cout << "b" << endl;
			int i = 0;
			for(const tuple<int, int, int> &path: charge_per_path){
				// int x = get<2>(path) / this->num_nodes, y = get<2>(path) % this->num_nodes;
				// cout << "(" << x << ", " << y << ", " << get<0>(path) << ", " << get<1>(path) << ")" << endl;
				this->selection_paths[i++] = get<2>(path);
			}
			// cout << endl;
			shuffle(this->selection_paths.begin(), this->selection_paths.end(), generator);
			auto end = chrono::steady_clock::now();
			// cout << "first step: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
			// cout << charge_per_path.size() << endl;
			// cout << this->total_paths.size() << endl;
			start = chrono::steady_clock::now();
			// cout << this->sample_size << " " << (int)charge_per_path.size() << endl;
			for(int i = 0; i < min(this->sample_size, (int)charge_per_path.size()); i++)
			{
				/*for(const unordered_set<int> &paths: this->R){
					for(int path: paths)
						cout << "(" << path / this->num_nodes << ", " << path % this->num_nodes << ")" << " ";
					cout << endl;
				}
				cout << endl;
				cout << "i: " << i << endl;
				cout << "charge per node: ";
				for(int i = 0; i < this->num_nodes; i++)
					cout << this->charge_per_node[i] << " ";
				cout << endl;
				cout << "charge per path: ";

				cout << endl;*/
				int id_path = this->selection_paths[i];
				for(int j = 0; j < this->num_nodes; j++)
				{
					if(this->R[j].find(id_path) == this->R[j].end())
						continue;
					this->charge_per_node[j] -= this->charge;
					this->erased_paths[j].insert(id_path);
				}

				// cout << "before disjktra " << endl;
				// cout << "x: " << id_path / this->num_nodes << " y: " << id_path % this->num_nodes << endl;
				int x = id_path / this->num_nodes, y = id_path % this->num_nodes;
				vector<int> prev = this->dijkstra(x, y).second;
				vector<int> path = this->build_path(prev, x, y);
				// cout << "after disjktra " << endl;
				for(int v: path)
				{
					// cout << v << " ";
					this->charge_per_node[v] += this->charge;
					this->added_paths[v].insert(id_path);
				}	
				
			}
			end = chrono::steady_clock::now();
			// cout << "second step: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
			start = chrono::steady_clock::now();
			this->modify_routing(this->added_paths, this->erased_paths);			// cout << "fin " << endl;
			end = chrono::steady_clock::now();
			// cout << "third step: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;

		}
		void update_distances()
		{
			auto start = chrono::steady_clock::now();
			for(int i = 0; i < this->num_nodes; i++){
				vector<int> distances_dijkstra = this->dijkstra(i, -1).first;
				for(int j = 0; j < this->num_nodes; j++)
					this->total_distances[i * this->num_nodes + j] = distances_dijkstra[j];
			}
			auto end = chrono::steady_clock::now();
			// cout << "update step: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
		}
		void forget_neighborhood_solution()
		{			
			this->modify_routing(this->erased_paths, this->added_paths);
		}
		pair<int, int> evaluate()
		{
			return {
				*max_element(this->total_charge.begin(), this->total_charge.end()), 
				*max_element(this->charge_per_node.begin(), this->charge_per_node.end())
			};
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
};

class SimulatedAnnealing
{
	double temperature;
	double beta; 
	double min_temperature;
	Routing routing;
	Routing best_routing;

public:
	SimulatedAnnealing(Routing &routing, double temperature, double min_temperature, double beta)
	{
		this->temperature = temperature;
		this->routing = routing;
		this->best_routing = routing;
		this->min_temperature = min_temperature;
		this->beta = beta;
	}
	Routing optimize()
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
			this->routing.update_distances();
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
	/*vector<unordered_set<int>> R;
		vector<int> charge_per_node;
		vector<int> path_per_node;
		vector<int> selection_paths;
		vector<unordered_set<int>> added_paths;
		vector<unordered_set<int>> erased_paths;
		vector<vector<int>> g;
		set<pair<int, int>, greater<pair<int, int>>> charge_per_path;
		vector<int> total_charge;
		int charge;
		int neighborhood_len;
		int num_nodes;
		int universe_size;
		int sample_size;*/

	
	Routing routing(g, 180, 240);

	/*for(const unordered_set<int> &paths: routing.R){
		for(int path: paths)
			cout << "(" << path / n << ", " << path % n << ")" << " ";

		cout << endl;
	}
	cout << endl;
	for(int charge: routing.charge_per_node)
		cout << charge << " ";
	cout << endl;*/

	SimulatedAnnealing simulated_annealing(routing, 1, 1e-1, 1e-1);

	auto start = chrono::steady_clock::now();
	Routing sol = simulated_annealing.optimize();
	auto end = chrono::steady_clock::now();
	cout << "time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;

	/*for(const unordered_set<int> &paths: sol.R){
		for(int path: paths)
			cout << "(" << path / n << ", " << path % n << ")" << " ";

		cout << endl;
	}*/

	cout << "annealing: " << sol.evaluate().first << endl;
	int sum = 0;
	for(int i = 0; i < n; i++){
		sum += sol.charge_per_node[i];
		cout << sol.charge_per_node[i] << " ";
	}
	cout << endl;

	return 0;
}