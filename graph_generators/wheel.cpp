#include <fstream>
#include <set>
#include <vector>
#include <iostream>

using namespace std;
#define ll long long int


void write_wheel_graph(int n, string file_path){
	
	ofstream myfile;
  	myfile.open(file_path);

	myfile << n << " " << 2 * (n - 1) << "\n";

    vector<vector<int>> g = vector<vector<int>>(n, vector<int>());
    for(int i = 1; i < n; i++){
        g[0].push_back(i);
        if(i + 1 < n){
            g[i].push_back(i + 1);
        }
        else{
            g[1].push_back(i);
        }
            
    }
    for(int i = 0; i < n; i++){
        for(int node: g[i]){
            myfile << i << " " << node << "\n";
        }
    }
    
	myfile.close();
}


// int main(){
//     cout << "Here" << endl;
// 	write_wheel_graph(700, "graph_files/wheel.txt");
//     return 0;
// }
