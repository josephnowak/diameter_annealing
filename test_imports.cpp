
#include "diameter_annealing/simulated_annealing.hpp"
#include "diameter_annealing/routing.hpp"
#include <chrono>
#include <random>
#include <assert.h>

using namespace std;


int main()
{
	
	int n, m;
	cin >> n >> m;
	vector<vector<int>> g(n, vector<int>());
	for(int i = 0; i < m; i++)
	{
		int x, y;
		cin >> x >> y;
		g[x].push_back(y);
		g[y].push_back(x);
	}

	
	SymmetricRouting routing(g, 5, 80, 0.7);

	SimulatedAnnealing simulated_annealing(routing, 15, 1e-1, 1, 50);

	auto start = chrono::steady_clock::now();
	simulated_annealing.optimize();
	auto end = chrono::steady_clock::now();
	cout << "time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
	
	cout << "annealing: " << routing.evaluate().first << endl;
	// int sum = 0;
	// for(int i = 0; i < n; i++){
	// 	sum += sol.charge_per_node[i];
	// 	cout << sol.charge_per_node[i] << " ";
	// }
	// cout << endl;

	vector<vector<vector<int>>> result_routing = routing.get_standard_format_routing();

	assert(routing.validate_routing(result_routing, routing.evaluate().first));
	return 0;
}