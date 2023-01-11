#include <unordered_map>
#include <vector>

#include <fstream>
#include <set>
#include <map>
#include <math.h>
#include <assert.h>
#include <iostream>


using namespace std;

#define ll long long int
#define um unordered_map

void generate_de_bruijn(
	um<string, vector<string>> &g, 
	string &base, 
	string &set_char
)
{
	g[base] = vector<string>();
	for(char s: set_char)
	{
		string node = base.substr(1, base.size() - 1) + s;
		if(g.find(node) == g.end())
			generate_de_bruijn(g, node, set_char);
		
		if(node != base)
			g[base].push_back(node);
		// if(node != base)
		// 	g[node].push_back(base);
	}
}


void write_de_bruijn_graph(int n, int m, string file_path){
	// M is the the number of symbols (the alphabet) and N is the number of symbol in every word or node.
	// This graph has M ^ N nodes 
	
	ofstream myfile;
  	myfile.open(file_path);

	um<string, vector<string>> g;
	// m es el número de simbolos y n es la cantidad de ellos en cada nodo

	string set_char;
	for(int i = 0; i < m; i++)
		set_char += '0' + i;
	string base = "";
	for(int i = 0; i < n; i++)
		base += set_char[0];
	generate_de_bruijn(g, base, set_char);

	assert(g.size() == pow(m, n));

	set<string> keys;
	for(auto &p: g)
		keys.insert(p.first);

	map<string, int> conversion;
	for(string s: keys)
		conversion[s] = conversion.size();

	set<pair<int,int>> edges;
	for(auto &p: g){
		int x = conversion[p.first];
		for(string s: p.second){
			int y = conversion[s];
			edges.insert({min(x, y), max(x, y)});
		}
	}
	
	// cout << g.size() << " " << edges.size() << "\n";
	myfile << g.size() << " " << edges.size() << "\n";

	for(const pair<int, int> &p: edges){
		myfile << p.first << " " << p.second << endl;
	}
	myfile.close();
}

// int main(){
// 	write_de_bruijn_graph(2, 8, "graph_files/de_bruijn.txt");
// }