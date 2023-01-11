#include <fstream>
#include <set>
#include <vector>
#include <iostream>

using namespace std;
#define ll long long int


void write_cycle_graph(int n, string file_path){
	
	ofstream myfile;
  	myfile.open(file_path);

	myfile << n << " " << n << "\n";

    vector<vector<int>> g = vector<vector<int>>(n, vector<int>());
    for(int i = 0; i < n; i++){
        g[i].push_back((i + 1) % n);
    }
    for(int i = 0; i < n; i++){
        for(int node: g[i]){
            myfile << i << " " << node << "\n";
        }
    }
    
	myfile.close();
}

