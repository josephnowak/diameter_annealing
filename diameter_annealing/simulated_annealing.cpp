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


class SimulatedAnnealing
{
	double temperature;
	double beta; 
	double min_temperature;
	BaseRouting routing;
	BaseRouting best_routing;

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

