#include <fstream>
#include <set>
#include <vector>

using namespace std;
#define ll long long int


void write_hypercube_graph(int n, string file_path){
	
	ofstream myfile;
  	myfile.open(file_path);

	myfile << (1 << n) << " " << (1 << (n - 1)) * n << "\n";

	// create the hypercube using a bottom up approach, this process takes two hypercube graph Qi and Pi with i being the grade
	// of the graphs and join them to create a new hypercube graph of grade i + 1 using the following logic:
	// Take Pi and insert it in the new graph.
	// Take the second graph and sum the len of the graph to the nodes and insert it in the new graph.
	// Create and edge from every node u of Qi to every node v of Pi, where u = v - len(Pi)
	vector<vector<int>> prev_graph = {{}};
	vector<vector<int>> graph;

	for(int i = 1; i <= n; i++){
		int prev_size = prev_graph.size();
		int size_g = prev_size * 2;
		graph = vector<vector<int>>(size_g, vector<int>());
		for(int j = 0; j < prev_size; j++){
			for(int node: prev_graph[j]){
				graph[j].push_back(node);
				graph[j + prev_size].push_back(node + prev_size);
			}
			graph[j].push_back(j + prev_size);
			graph[j + prev_size].push_back(j);
		}
		prev_graph = graph;
	}
	for(int i = 0; i < (int)graph.size(); i++){
		for(int node: graph[i]){
			if(i < node)
				myfile << i << " " << node << "\n";
		}
	}
	myfile.close();
}


