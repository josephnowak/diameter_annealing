#include <bits/stdc++.h>
#include <set>
using namespace std;
#define ll long long int


int main()
{
	int n;
	cin >> n;
	vector<vector<vector<int>>> graphs = vector<vector<vector<int>>>(n + 1);
	graphs[2] = {{2, 4}, {1, 3}, {2, 4}, {1, 3}};
	
	for(int i = 3; i <= n; i++)
	{
		graphs[i] = graphs[i - 1];
		int m = graphs[i].size();
		int size_g = (1 << i);
		for(int j = m; j < size_g; j++)
		{
			graphs[i].push_back(graphs[i - 1][j - m]);
			
			for(int p = 0; p < graphs[i - 1][j - m].size(); p++)
				graphs[i][graphs[i].size() - 1][p] += m;
				graphs[i][graphs[i].size() - 1].push_back(j - m + 1);
				
		}
		
		for(int j = 0; j < m; j++)
			graphs[i][j].push_back(j + m + 1);
			
	}

	set<pair<int,int>> set_pairs;
	for(int i = 1; i <= graphs[n].size(); i++)
	{
		for(int node: graphs[n][i - 1])
		{
			if((set_pairs.find({i, node}) == set_pairs.end()) and 
				(set_pairs.find({node, i}) == set_pairs.end()))
			{
				set_pairs.insert({i, node});
				cout << i - 1 << " " << node - 1 << endl;
			}
		}
	}
	
	
	return 0;
}
