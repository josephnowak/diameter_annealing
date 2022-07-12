#include "diameter_annealing/simulated_annealing.hpp"
#include "diameter_annealing/routing.hpp"

#include <iostream>
#include <chrono>

using namespace std;


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

	
	SymmetricRouting routing(g, 5, 50, 0.5);

	// SimulatedAnnealing<SymmetricRouting> simulated_annealing(routing, 10, 1e-1, 1);

	// auto start = chrono::steady_clock::now();
	// SymmetricRouting sol = simulated_annealing.optimize();
	// auto end = chrono::steady_clock::now();
	// cout << "time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
	
	// cout << "annealing: " << sol.evaluate().first << endl;
	// // int sum = 0;
	// // for(int i = 0; i < n; i++){
	// // 	sum += sol.charge_per_node[i];
	// // 	cout << sol.charge_per_node[i] << " ";
	// // }
	// // cout << endl;

	// vector<vector<vector<int>>> result_routing = sol.get_standard_format_routing();

	// cout << sol.validate_routing(result_routing, sol.evaluate().first) << endl;
	return 0;
}