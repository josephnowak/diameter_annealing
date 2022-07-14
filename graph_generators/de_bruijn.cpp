#include <bits/stdc++.h>
#include <unordered_map>
#include <set>
using namespace std;
#define ll long long int
#define um unordered_map

void generate_de_bruijn(um<string, vector<string>> &g, 
						um<string, int> &conversion,
						string &base, string &set_char)
{
	g[base] = vector<string>();
	conversion[base] = conversion.size();
	for(char s: set_char)
	{
		string node = base.substr(1, base.size() - 1) + s;
		if(conversion.find(node) == conversion.end())
			generate_de_bruijn(g, conversion, node, set_char);
		g[base].push_back(node);
		if(node != base)
			g[node].push_back(base);
	}
}

int main()
{
	um<string, vector<string>> g;
	um <string, int> conversion;
	// m es el número de simbolos y n es la cantidad de ellos en un string
	int n, m;
	cin >> n >> m;
	string set_char;
	for(int i = 0; i < m; i++)
		set_char += 'a' + i;
	string base = "";
	for(int i = 0; i < n; i++)
		base += set_char[0];
	generate_de_bruijn(g, conversion, base, set_char);
	
	set<pair<int,int>> set_pairs;
	// cout << conversion.size() << endl;
	// cout << pow(m, n) << endl;
	for(auto p: g)
	{
		int node_a = conversion[p.first] - 1;
		for(string s: p.second){
			int node_b = conversion[s] - 1;
			if((set_pairs.find({node_a, node_b}) == set_pairs.end()) and 
				(set_pairs.find({node_b, node_a}) == set_pairs.end()))
			{
				set_pairs.insert({node_a, node_b});
				if (node_a != node_b)
					cout << node_a + 1 << " " << node_b + 1 << endl;
				// cout << p.first << " " << s << endl;
			}
		}
	}
	
}

