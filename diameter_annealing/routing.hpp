#ifndef _ROUTING_H_
#define _ROUTING_H_

#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include "utils.hpp"

using namespace std;

#define um unordered_map
#define us unordered_set
//default_random_engine generator;
// random_device rd;
// mt19937 generator(rd());



class BaseRouting {
	/**
	 *
	 * @brief Abstract class used on the simulated annealing class/algorithm to control how to find neighborhood_solution
	 */
	public:
		int sample_size;
		/**
		 * This is equivalent to calculate the R[node].size() * charge, is only useful to simplify the use of symmetric paths
		 *
		 * @brief Charge of everynode on the routing or how many paths use this node as internal 
		 */
		vector<int> charge_per_node;

		/**
		 * The data structure is equivalent to adjacency list representation but using arrays/vectors
		 * @brief Graph that represent the connection between every pair of nodes
		 */
		vector<vector<int>> g;

		/**
		 * This is equivalent to the sum of the charge of every internal node on the routing
		 *
		 * @brief Charge of every path
		 */
		vector<int> total_charge_per_path;

		/**
		 * This can be implemented in multiple ways
		 * @brief Find a neighborhood solution
		 */
		virtual void neighborhood_solution() = 0;
		
		/**
		 * It is useful for expensive computations like update the distances using a dijkstra algorithm
		 * @brief Method executed before reduce the temperature on the simulated annealing algorithm
		 */
		virtual void update_routing() = 0;

		/**
		 * This method is called everytime the simulated annealing algorithm decide to not use the new neighborhood solution
		 *
		 * @brief Method used to forget the new solution. 
		 */
		virtual void forget_neighborhood_solution() = 0;

		/**
		 * This method is called everytime the simulated annealing algorithm decide to not use the new neighborhood solution
		 *
		 * @brief Method used to forget the new solution. 
		 */
		virtual vector<vector<vector<int>>> get_standard_format_routing() = 0;


		/**
		 *
		 * @brief Return a pair with the forwarding diameter (1) and any other value that is going to break a tie (2)
		 */
		virtual pair<int, int> evaluate() = 0;

};


class SymmetricRouting: virtual public BaseRouting {
	/**
	 * Implementation of the BaseRounting using symmetric, simple paths and dijkstra to find neighborhood solutions.
	 *
	 * @brief Routing of the graph
	 */
	private:
		/**
		 * @brief Ids of all the path, they are generated as autoincrement indexes.
		 */
		unordered_set<int> paths_ids;

		/**
		 * This is useful for speed up the recalculation of the forwarding diameter and is also useful
		 * on the forget_neighborhood_solution method.
		 * @brief Paths added during the execution of the neighborhood_solution method.
		 */
		vector<unordered_map<int, int>> added_paths;
		/**
		 * This is useful for speed up the recalculation of the forwarding diameter and is also useful
		 * on the forget_neighborhood_solution method.
		 * @brief Paths deleted during the execution of the neighborhood_solution method.
		 */
		vector<unordered_map<int, int>> erased_paths;

		/**
		 * @brief Vector that contains all the paths that are going to be use to modify the actual routing 
		 */
		vector<int> selected_paths;

		/**
		 * @brief Dijkstra distance of every path
		 */
		vector<int> dijkstra_path_distances;

		/**
		 * Useful only for symmetric paths
		 * @brief Represent the charge generated by every use of a node as a internal path
		 */
		int charge;

		/**
		 * @brief Number of node on the graph
		 */
		int num_nodes;

		/**
		 * @brief Probability of use the dijkstra algorithm to find the path instead of the random path method.
		 */
		float regular_path_probability;

		void generate_initial_routing();

		vector<int> build_path(vector<int> &prev, int x, int y);

		pair<vector<int>, vector<int>> dijkstra(int x, int y);

		void modify_routing(vector<unordered_map<int, int>> &new_paths, vector<unordered_map<int, int>> &old_paths);

		void get_bfs_random_paths(int source);

	public:
		/**
		 * Every cell/node of the vector is a set with all the paths (IDs) that use that node as an internal one.
		 * @brief SymmetricRouting of the network
		 */
		vector<unordered_map<int, int>> R;

		
		/**
		 * This attribute is useful to restrict the randomness of the algorithm, an smaller integer (call it N)
		 * means that the algorithm is only going to take into account the top N paths based on the three metrics
		 * described on the neighborhood_solution method, this parameter is combined with the sample_size
		 * 
		 * @brief Determine the number of best possible paths to select randonmly
		 */
		int universe_size;

		/**
		 * Once the numbre of paths is restricted based on the universe_size, they are sorted randomly and then M of them
		 * are selected (this m is the sample_size)
		 * @brief Number of paths that are going to be modified every time the neighborhood_solution method is called
		 */
		int sample_size;

		// uniform_real_distribution<double> uniform_real
	
		SymmetricRouting(){}

		SymmetricRouting(
			vector<vector<int>> &g, 
			int sample_size, 
			int universe_size,
			float regular_path_probability
		);

		SymmetricRouting(const SymmetricRouting &a);

		void neighborhood_solution();

		void update_routing();

		void forget_neighborhood_solution();

		pair<int, int> evaluate();

		vector<vector<vector<int>>> get_standard_format_routing();

		/**
		 * The validation is simple, it takes every pair of nodes call them X and Y and then it takes the path from X to Y
		 * and verify if it possible to reach Y from X using that path, then it verify if the forwarding diameter is the same.
		 * @brief Validate if the routing is valid and has the expected forwarding diameter, useful for debug
		 */
		bool validate_routing(vector<vector<vector<int>>> &routing, int expected_forwarding_diameter);
};


#endif
