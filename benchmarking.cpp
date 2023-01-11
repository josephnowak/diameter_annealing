
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
#include <math.h>

using namespace std;


void convert_to_geometry_reduction(
	vector<vector<double>> &simulated_inputs,
	float min_temp,
	int extra_index=0
){
	// Graph Size, Sample Size, Universe Size, N iterations, path prob, tempe, min tempe, beta, test_repeats
	for(vector<double> &v: simulated_inputs){
		double initial_temp = v[5 + extra_index];
		double beta = v[7 + extra_index];
		double num_iters = initial_temp / beta;
		double new_beta = pow(min_temp / initial_temp, 1.0/num_iters);
		v[6 + extra_index] = min_temp;
		v[7 + extra_index] = new_beta;
	}		
}


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
		{13500, 5, 5, 1, 0.99, 10.0, 0, 10, 1},

	};
	convert_to_geometry_reduction(simulated_inputs, 0.01);
	
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

		cout << "Processing wheel: " << graph_size << endl;

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
		{10, 1, 5, 5, 0.5, 10.0, 0, 0.25, 50},
		{30, 2, 10, 10, 0.5, 10.0, 0, 0.25, 50},
		{50, 2, 15, 10, 0.5, 10.0, 0, 0.25, 50},
		{70, 3, 30, 20, 0.5, 10.0, 0, 0.25, 50},

		// Medium size cases, trying to find a good result / time  ratio
		{100, 5, 20, 5, 0.7, 10.0, 0, 0.25, 10},
		{300, 10, 70, 5, 0.7, 10.0, 0, 0.25, 10},
		{500, 10, 300, 5, 0.8, 10.0, 0, 0.5, 5},
		{700, 20, 500, 3, 0.8, 10.0, 0, 0.5, 5},

	};
	convert_to_geometry_reduction(simulated_inputs, 0.01);

	
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

		cout << "Processing cycle: " << graph_size << endl;

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
		{4, 1, 5, 20, 0.6, 10.0, 0, 0.25, 50},
		{5, 2, 15, 30, 0.6, 10.0, 0, 0.25, 50},
		{6, 3, 25, 40, 0.6, 10.0, 0, 0.25, 50},
		{7, 5, 60, 50, 0.6, 10.0, 0, 0.25, 50},

		// Medium size cases, trying to find a good result / time  ratio
		{8, 8, 100, 50, 0.8, 10.0, 0, 0.20, 10},
		{9, 15, 200, 35, 0.8, 10.0, 0, 0.20, 10},

		// Big size cases, restricted by time
		{10, 50, 130, 10, 0.99, 10.0, 0, 0.5, 5},
		{11, 100, 500, 5, 0.99, 10.0, 0, 1, 5},
		{12, 1000, 1500, 1, 0.99, 10.0, 0, 1, 3},
		{13, 1, 2, 1, 0.99, 10.0, 0, 1, 1},

		// Biggest case tested

	};
	convert_to_geometry_reduction(simulated_inputs, 0.01);

	
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

		cout << "Processing hypercube: " << hypercube_size << endl;

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



void run_de_bruijn_graph(ofstream &myfile){
	// Graph Size, Sample Size, Universe Size, N iterations, path prob, tempe, min tempe, beta, test_repeats
	vector<vector<double>> simulated_inputs = {
		// small cases, trying to find the minimum with stable results and small time
		{4, 2, 1, 10, 20, 0.6, 10.0, 0, 0.25, 50},
		{5, 2, 2, 50, 50, 0.6, 10.0, 0, 0.25, 50},
		{4, 3, 3, 100, 70, 0.6, 10.0, 0, 0.25, 50},
		
		// Medium size cases, trying to find a good result / time  ratio
		{5, 3, 15, 130, 40, 0.8, 10.0, 0, 0.25, 30},
		{4, 4, 25, 300, 50, 0.8, 10.0, 0, 0.25, 30},
		{3, 7, 30, 350, 50, 0.9, 10.0, 0, 0.25, 30},
		{4, 5, 50, 500, 30, 0.95, 10.0, 0, 0.25, 10},
		{3, 9, 80, 700, 30, 0.99, 10.0, 0, 0.25, 5},

		// Big size cases, restricted by time
		{4, 8, 300, 500, 1, 0.99, 10.0, 0, 1, 3},
		{4, 9, 500, 700, 1, 0.99, 10.0, 0, 1, 3},
		
		// Biggest case tested
		{2, 100, 1, 2, 1, 0.99, 10.0, 0, 10, 1},
	};
	convert_to_geometry_reduction(simulated_inputs, 0.01, 1);
	
	for(vector<double> &inputs: simulated_inputs){
		int n = inputs[0];
		int m = inputs[1];
		int sample_size = inputs[2];
		int universe_size = inputs[3];
		int n_iterations = inputs[4];
		float regular_path_probability = inputs[5];
		double temperature = inputs[6];
		double min_temperature = inputs[7];
		double beta = inputs[8];
		int test_repeats = inputs[9];

		cout << "Processing de bruijn: " << n << ", " << m << endl;

		write_de_bruijn_graph(n, m, "graph_files/de_bruijn.txt");
		vector<vector<int>> g = read_graph("graph_files/de_bruijn.txt");
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
			auto end = chrono::steady_clock::now();

			chrono::duration<double> elapsed_time = end - start;
			myfile << "B(" << m << "," << n << ");";
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
	myfile << "Notation;" << "V;" << "E;" << "InitialDiameter;" << "FinalDiameter;";
	myfile << "Time;" << "SampleSize;" << "UniverseSize;" << "RegularPathProbability;";
	myfile << "T;" << "MinT;" << "Beta;" << "IterN" << "\n";

	run_wheel_graph(myfile);
	run_cycle_graph(myfile);
	run_hypercube_graph(myfile);
	run_de_bruijn_graph(myfile);
	myfile.close();
	// useful to verify the solution
	// vector<vector<vector<int>>> result_routing = routing.get_standard_format_routing();

	// assert(routing.validate_routing(result_routing, routing.evaluate().first));
	return 0;
}