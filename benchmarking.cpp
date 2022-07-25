
#include "diameter_annealing/simulated_annealing.hpp"
#include "diameter_annealing/routing.hpp"
#include "graph_generators/cycle.hpp"
#include "graph_generators/de_bruijn.hpp"
#include "graph_generators/hypercube.hpp"
#include "graph_generators/wheel.hpp"

#include <chrono>
#include <random>
#include <assert.h>
#include <fstream>

using namespace std;


vector<vector<int>> read_graph(string file_path){
	ifstream myfile;
  	myfile.open(file_path);

	int n, m;
	myfile >> n >> m;
	vector<vector<int>> g(n, vector<int>());
	
	for(int i = 0; i < m; i++)
	{
		int x, y;
		myfile >> x >> y;
		g[x].push_back(y);
		g[y].push_back(x);
	}
	return g;
}


pair<int, int> run(
	vector<vector<int>> &g,
	int sample_size, 
	int universe_size,
	float regular_path_probability,
	double temperature, 
    double min_temperature, 
    double beta, 
    int n_iterations
){
	
	SymmetricRouting routing(g, sample_size, universe_size, regular_path_probability);

	SimulatedAnnealing simulated_annealing(
		routing, 
		temperature,
		min_temperature, 
		beta, 
		n_iterations,
		false
	);

	simulated_annealing.optimize();
	// useful to verify the solution
	// vector<vector<vector<int>>> result_routing = routing.get_standard_format_routing();

	// assert(routing.validate_routing(result_routing, routing.evaluate().first));
	// cout << "diameter: " << routing.evaluate().first << " initial " << simulated_annealing.initial_routing_diameter << endl;

	return {routing.evaluate().first, simulated_annealing.initial_routing_diameter};

}


int count_edges(vector<vector<int>> &g){
	int edge_sizes = 0;
	for(vector<int> &edges: g)
		edge_sizes += (int)edges.size();
	edge_sizes /= 2;
	return edge_sizes;
}


void run_wheel_graph(ofstream &myfile){
	// Graph Size, Sample Size, Universe Size, N iterations, path prob, tempe, min tempe, beta, test_repeats
	vector<vector<double>> simulated_inputs = {
		// small cases, trying to find the minimum with stable results and small time
		{10, 1, 3, 5, 0.5, 10.0, 0, 0.25, 50},
		{30, 2, 10, 15, 0.5, 10.0, 0, 0.25, 50},
		{51, 2, 15, 30, 0.5, 10.0, 0, 0.25, 50},
		{70, 3, 20, 35, 0.5, 10.0, 0, 0.25, 50},

		// Medium size cases, trying to find a good result / time  ratio
		{100, 5, 20, 35, 0.7, 10.0, 0, 0.25, 10},
		{300, 10, 70, 50, 0.7, 10.0, 0, 0.25, 10},
		{501, 30, 200, 45, 0.7, 10.0, 0, 0.25, 10},
		{700, 70, 350, 30, 0.8, 10.0, 0, 0.25, 10},

		// Big size cases, restricted by time
		{3000, 5000, 6000, 5, 0.99, 10.0, 0, 0.5, 3},
		{7000, 10000, 11000, 3, 0.99, 10.0, 0, 0.5, 3},

		// Biggest case tested



	};
	
	for(vector<double> &inputs: simulated_inputs){
		int graph_size = inputs[0];
		int sample_size = inputs[1];
		int universe_size = inputs[2];
		int n_iterations = inputs[3];
		float regular_path_probability = inputs[4];
		double temperature = inputs[5];
		double min_temperature = inputs[6];
		double beta = inputs[7];
		int test_repeats = inputs[8];

		write_wheel_graph(graph_size, "graph_files/wheel.txt");
		vector<vector<int>> g = read_graph("graph_files/wheel.txt");

		int edge_sizes = count_edges(g);

		for(int i = 0; i < test_repeats; i++){
			

			auto start = chrono::steady_clock::now();
			pair<int, int> diameter = run(
				g, 
				sample_size, 
				universe_size, 
				regular_path_probability, 
				temperature, 
				min_temperature, 
				beta, 
				n_iterations
			);
			auto end = chrono::steady_clock::now();

			chrono::duration<double> elapsed_time = end - start;
			myfile << "W(" << graph_size - 1 << ");";
			myfile << graph_size << ";";
			myfile << edge_sizes << ";";
			myfile << diameter.second << ";";
			myfile << diameter.first << ";";
			myfile << elapsed_time.count() << ";";	
			myfile << sample_size << ";";
			myfile << universe_size << ";";
			myfile << regular_path_probability << ";";
			myfile << temperature << ";";
			myfile << min_temperature << ";";
			myfile << beta << ";";
			myfile << n_iterations << "\n";
		}
	}
}

void run_cycle_graph(ofstream &myfile){
	// Graph Size, Sample Size, Universe Size, N iterations, path prob, tempe, min tempe, beta, test_repeats
	vector<vector<double>> simulated_inputs = {
		// small cases, trying to find the minimum with stable results and small time
		{10, 1, 3, 5, 0.5, 10.0, 0, 0.25, 50},
		{30, 2, 10, 10, 0.5, 10.0, 0, 0.25, 50},
		{50, 2, 15, 10, 0.5, 10.0, 0, 0.25, 50},
		{70, 3, 30, 20, 0.5, 10.0, 0, 0.25, 50},

		// Medium size cases, trying to find a good result / time  ratio
		{100, 5, 20, 5, 0.7, 10.0, 0, 0.25, 10},
		{300, 10, 70, 5, 0.7, 10.0, 0, 0.25, 10},
		{500, 10, 300, 5, 0.8, 10.0, 0, 0.5, 5},
		{700, 20, 500, 3, 0.8, 10.0, 0, 0.5, 5},

	};
	
	for(vector<double> &inputs: simulated_inputs){
		int graph_size = inputs[0];
		int sample_size = inputs[1];
		int universe_size = inputs[2];
		int n_iterations = inputs[3];
		float regular_path_probability = inputs[4];
		double temperature = inputs[5];
		double min_temperature = inputs[6];
		double beta = inputs[7];
		int test_repeats = inputs[8];

		write_cycle_graph(graph_size, "graph_files/cycle.txt");
		vector<vector<int>> g = read_graph("graph_files/cycle.txt");

		int edge_sizes = count_edges(g);

		for(int i = 0; i < test_repeats; i++){
			auto start = chrono::steady_clock::now();
			pair<int, int> diameter = run(
				g, 
				sample_size, 
				universe_size, 
				regular_path_probability, 
				temperature, 
				min_temperature, 
				beta, 
				n_iterations
			);
			auto end = chrono::steady_clock::now();

			chrono::duration<double> elapsed_time = end - start;
			myfile << "C(" << graph_size << ");";
			myfile << graph_size << ";";
			myfile << edge_sizes << ";";
			myfile << diameter.second << ";";
			myfile << diameter.first << ";";
			myfile << elapsed_time.count() << ";";	
			myfile << sample_size << ";";
			myfile << universe_size << ";";
			myfile << regular_path_probability << ";";
			myfile << temperature << ";";
			myfile << min_temperature << ";";
			myfile << beta << ";";
			myfile << n_iterations << "\n";
		}
	}
}


void run_hypercube_graph(ofstream &myfile){
	// Graph Size, Sample Size, Universe Size, N iterations, path prob, tempe, min tempe, beta, test_repeats
	vector<vector<double>> simulated_inputs = {
		// small cases, trying to find the minimum with stable results and small time
		{4, 1, 5, 15, 0.6, 10.0, 0, 0.25, 50},
		{5, 2, 10, 25, 0.6, 10.0, 0, 0.25, 50},
		{6, 3, 20, 35, 0.6, 10.0, 0, 0.25, 50},
		{7, 5, 20, 45, 0.6, 10.0, 0, 0.25, 50},

		// Medium size cases, trying to find a good result / time  ratio
		{8, 8, 70, 50, 0.8, 10.0, 0, 0.20, 10},
		{9, 15, 100, 35, 0.8, 10.0, 0, 0.20, 10},

		// Big size cases, restricted by time
		{10, 50, 130, 10, 0.99, 10.0, 0, 0.5, 5},
		{11, 100, 500, 5, 0.99, 10.0, 0, 1, 5},
		{12, 1000, 1500, 1, 0.99, 10.0, 0, 1, 2},
		{13, 3000, 3500, 1, 0.99, 10.0, 0, 1, 2},

		// Biggest case tested

	};
	
	for(vector<double> &inputs: simulated_inputs){
		int hypercube_size = inputs[0];
		int sample_size = inputs[1];
		int universe_size = inputs[2];
		int n_iterations = inputs[3];
		float regular_path_probability = inputs[4];
		double temperature = inputs[5];
		double min_temperature = inputs[6];
		double beta = inputs[7];
		int test_repeats = inputs[8];

		write_hypercube_graph(hypercube_size, "graph_files/hypercube.txt");
		vector<vector<int>> g = read_graph("graph_files/hypercube.txt");
		int graph_size = g.size();

		int edge_sizes = count_edges(g);

		for(int i = 0; i < test_repeats; i++){
			auto start = chrono::steady_clock::now();
			pair<int, int> diameter = run(
				g, 
				sample_size, 
				universe_size, 
				regular_path_probability, 
				temperature, 
				min_temperature, 
				beta, 
				n_iterations
			);
			cout << "result: " << diameter.first << endl;
			auto end = chrono::steady_clock::now();

			chrono::duration<double> elapsed_time = end - start;
			myfile << "Q(" << hypercube_size << ");";
			myfile << graph_size << ";";
			myfile << edge_sizes << ";";
			myfile << diameter.second << ";";
			myfile << diameter.first << ";";
			myfile << elapsed_time.count() << ";";	
			myfile << sample_size << ";";
			myfile << universe_size << ";";
			myfile << regular_path_probability << ";";
			myfile << temperature << ";";
			myfile << min_temperature << ";";
			myfile << beta << ";";
			myfile << n_iterations << "\n";
		}
	}
}



int main()
{
	ofstream myfile;
  	myfile.open("executions.csv");
	myfile << "Notation;" << "V;" << "E;" << "ID;" << "FD;";
	myfile << "Time;" << "SS;" << "US;" << "RPP;";
	myfile << "T;" << "MinT;" << "Beta;" << "IterN" << "\n";
	run_hypercube_graph(myfile);
	// run_hypercube_graph(myfile);
	myfile.close();
	// useful to verify the solution
	// vector<vector<vector<int>>> result_routing = routing.get_standard_format_routing();

	// assert(routing.validate_routing(result_routing, routing.evaluate().first));
	return 0;
}