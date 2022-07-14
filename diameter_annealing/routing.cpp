#include <random>
#include <math.h>
#include <climits>
#include <set>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <chrono>
#include <map>
#include <iostream>
#include <assert.h>

#include "routing.hpp"

using namespace std;
using namespace std::chrono;

#define um unordered_map
#define us unordered_set
// default_random_engine generator;



template <class bidiiter>
bidiiter BaseRouting::random_choose(bidiiter begin, bidiiter end, size_t num_random)
{
    size_t left = std::distance(begin, end);
    while (num_random--)
    {
        bidiiter r = begin;
        uniform_int_distribution<int> uniform_int_path(0, left - 1);
        int mov = uniform_int_path(this->generator);
        std::advance(r, mov);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

double SymmetricRouting::get_random_probability(){
    uniform_real_distribution<double> uniform_real(0, 1);
    return uniform_real(this->generator);
}

int SymmetricRouting::get_random_int_number(int l, int r){;
    uniform_int_distribution<int> uniform_int_n(l, r);
    return uniform_int_n(this->generator);
}


void SymmetricRouting::update_internal_routing_data()
{
    for (int i = 0; i < this->num_nodes; i++)
    {
        this->charge_per_node[i] = this->R[i].size() * this->charge;
        for (const pair<int, int> path_map : this->R[i])
        {
            this->total_charge_per_path[path_map.first] += this->charge_per_node[i];
            this->paths_ids.insert(path_map.first);
        }
    }
}

void SymmetricRouting::generate_initial_routing()
{
    for (int i = 0; i < this->num_nodes; i++)
        get_bfs_random_paths(i);

    this->update_internal_routing_data();
    this->update_routing();
}
vector<int> SymmetricRouting::build_path(vector<int> &prev, int x, int y)
{
    vector<int> path;
    int act = y;
    while (prev[act] != act)
    {
        if (prev[act] != x)
            path.push_back(prev[act]);
        act = prev[act];
    }
    reverse(path.begin(), path.end());
    return path;
}

pair<vector<int>, vector<int>> SymmetricRouting::dijkstra(int x, int y)
{
    set<tuple<int, float, int>> q;
    q.insert(make_tuple(0, this->get_random_probability(), x));
    vector<int> dist_per_node(this->num_nodes, INT_MAX);
    vector<int> prev(this->num_nodes);
    prev[x] = x;
    dist_per_node[x] = 0;
    while (!q.empty())
    {
        tuple<int, int, int> temp = *q.begin();
        int dist = get<0>(temp), random_value = get<1>(temp), node = get<2>(temp);
        if (dist < 0)
            break;
        dist_per_node[node] = dist;
        q.erase(q.begin());
        if (node == y)
            break;

        for (int &adj : this->g[node])
        {
            if (dist_per_node[adj] > dist + this->charge_per_node[adj] + 2)
            {
                prev[adj] = node;
                q.erase(make_tuple(dist_per_node[adj], random_value, adj));
                dist_per_node[adj] = dist + this->charge_per_node[adj] + 2;
                q.insert(make_tuple(dist_per_node[adj], this->get_random_probability(), adj));
            }
        }
    }
    return {dist_per_node, prev};
}
void SymmetricRouting::modify_routing(vector<unordered_map<int, int>> &new_paths, vector<unordered_map<int, int>> &old_paths)
{
    unordered_map<int, int> new_paths_charge, changed_charge_paths;

    for (int i = 0; i < this->num_nodes; i++)
    {
        int new_charge = this->charge * (this->R[i].size() + new_paths[i].size() - old_paths[i].size());
        int change_value = this->charge * (new_paths[i].size() - old_paths[i].size());
        // calculate the new charge for every node based on the amount of new
        this->charge_per_node[i] = new_charge;

        for (const pair<int, int> &id_path_map : old_paths[i])
            this->R[i].erase(id_path_map.first);

        if (change_value != 0)
            for (const pair<int, int> &path_map : this->R[i])
                changed_charge_paths[path_map.first] += change_value;

        for (const pair<int, int> &id_path_map : new_paths[i])
        {
            this->R[i].insert(id_path_map);
            new_paths_charge[id_path_map.first] += new_charge;
        }
    }
    for (const pair<int, int> &path_charge : new_paths_charge)
        this->total_charge_per_path[path_charge.first] = path_charge.second;

    for (const pair<int, int> &path_charge : changed_charge_paths)
        this->total_charge_per_path[path_charge.first] += path_charge.second;
}
void SymmetricRouting::get_bfs_random_paths(int source)
{
    queue<pair<int, int>> q;
    q.push({source, 0});
    vector<int> distances(this->num_nodes, -1);
    while (!q.empty())
    {
        pair<int, int> temp = q.front();
        int node = temp.first;
        int len = temp.second;
        q.pop();
        if (distances[node] != -1)
            continue;
        distances[node] = len;
        for (int v : this->g[node])
        {
            if (distances[v] != -1)
                continue;
            q.push({v, len + 1});
        }
    }
    for (int i = source + 1; i < this->num_nodes; i++)
    {
        int best_distance = distances[i];
        int act_node = i;
        int path_id = source * this->num_nodes + i;
        for (int dist = best_distance; dist > 1; dist--)
        {
            vector<int> valid;
            for (int &node : this->g[act_node])
                if (distances[node] == dist - 1)
                    valid.push_back(node);

            // modify this
            int pos = this->get_random_int_number(0, valid.size() - 1);
            act_node = valid[pos];
            this->R[act_node].insert({path_id, dist - 2});
        }
    }
}


SymmetricRouting::SymmetricRouting() {}

SymmetricRouting::SymmetricRouting(
    vector<vector<int>> &g,
    int sample_size,
    int universe_size,
    float regular_path_probability)
{
    this->g = g;
    this->num_nodes = g.size();
    this->R = vector<unordered_map<int, int>>(this->num_nodes);
    this->charge_per_node = vector<int>(this->num_nodes);
    this->added_paths = vector<unordered_map<int, int>>(this->num_nodes, unordered_map<int, int>());
    this->erased_paths = vector<unordered_map<int, int>>(this->num_nodes, unordered_map<int, int>());
    this->total_charge_per_path = vector<int>(this->num_nodes * this->num_nodes, 0);
    this->dijkstra_path_distances = vector<int>(this->num_nodes * this->num_nodes, 0);
    this->paths_ids = unordered_set<int>();
    this->charge = 2;
    this->regular_path_probability = regular_path_probability;

    this->generate_initial_routing();

    this->universe_size = min(universe_size, (int)this->paths_ids.size());
    this->sample_size = min(sample_size, this->universe_size);
    this->selected_paths = vector<int>(this->universe_size, -1);
}
SymmetricRouting::SymmetricRouting(const SymmetricRouting &a)
{
    this->num_nodes = a.num_nodes;
    this->R = a.R;
    this->charge_per_node = a.charge_per_node;
    this->g = a.g;
    this->total_charge_per_path = a.total_charge_per_path;
    this->paths_ids = a.paths_ids;
}
void SymmetricRouting::store_routing()
{
    // Store this information on memory to reduce the amount of memory consumption.
    this->best_R = this->R;
}

void SymmetricRouting::set_stored_routing()
{
    this->R = this->best_R;
    this->best_R.clear();
    this->total_charge_per_path = vector<int>(this->num_nodes * this->num_nodes, 0);
    this->charge_per_node = vector<int>(this->num_nodes);
    this->added_paths = vector<unordered_map<int, int>>(this->num_nodes, unordered_map<int, int>());
    this->erased_paths = vector<unordered_map<int, int>>(this->num_nodes, unordered_map<int, int>());
    this->total_charge_per_path = vector<int>(this->num_nodes * this->num_nodes, 0);
    this->paths_ids = unordered_set<int>();
    this->update_internal_routing_data();
    this->update_routing();
}

void SymmetricRouting::neighborhood_solution()
{
    for (int i = 0; i < this->num_nodes; i++)
    {
        this->erased_paths[i].clear();
        this->added_paths[i].clear();
    }

    struct CompareTuple
    {
        bool operator()(const tuple<int, int, int> &x, const tuple<int, int, int> &y) const
        {
            if (get<0>(x) == get<0>(y))
            {
                if (get<1>(x) == get<1>(y))
                {
                    // smaller nodes ID has priority
                    return get<2>(x) < get<2>(y);
                }
                // smaller dijkstra distance has priority
                return get<1>(x) < get<1>(y);
            }
            // biggest path charge has highest priority
            return get<0>(x) > get<0>(y);
        }
    };

    set<tuple<int, int, int>, CompareTuple> ordered_charge_per_path;

    for (int id_path : this->paths_ids)
    {
        /*
            Every path is sorted based on three elements:
            1. The path charge, which is the sum of the charge of the nodes used as internal nodes on the path
            2. The dijkstra distance, explained in the update_routing method
            3. The id of the path, used to break ties
        */
        ordered_charge_per_path.insert(make_tuple(
            this->total_charge_per_path[id_path],
            this->dijkstra_path_distances[id_path],
            id_path
        ));

        // preserve only universe_size elements
        if ((int)ordered_charge_per_path.size() > this->universe_size)
        {
            ordered_charge_per_path.erase(prev(ordered_charge_per_path.end()));
        }
    }
    // Select and sort randomly the paths on the ordered_charge_per_path
    int i = 0;
    for (const tuple<int, int, int> &element : ordered_charge_per_path)
    {
        this->selected_paths[i++] = get<2>(element);
    }
    this->random_choose(this->selected_paths.begin(), this->selected_paths.end(), this->sample_size);

    for (int i = 0; i < this->sample_size; i++)
    {
        int id_path = this->selected_paths[i];
        int x = id_path / this->num_nodes, y = id_path % this->num_nodes;

        // validate always that x is smaller than y, only for debugging
        assert(x < y);

        // remove the previous path first to improve the results of the dijkstra
        for (int j = 0; j < this->num_nodes; j++)
        {
            // verify if the path was using the jth node
            auto r_element = this->R[j].find(id_path);
            if (r_element == this->R[j].end())
                continue;
            // if it was using the node then reduce the charge of that node (only useful for the dijkstra)
            // and add it to the erased_paths set
            this->charge_per_node[j] -= this->charge;
            this->erased_paths[j].insert(*r_element);
        }

        double random_number = this->get_random_probability();
        vector<int> path;

        if (random_number <= this->regular_path_probability)
        {
            // get the dijkstra path from the node X to Ys
            vector<int> prev_path_node = this->dijkstra(x, y).second;
            path = this->build_path(prev_path_node, x, y);
        }
        else
        {
            // create a random path generator using the dijkstra algorithm and random charges per node
            // that are generated shuffeling the charge_per_node vector
            vector<int> copy_charge_per_node = this->charge_per_node;
            shuffle(this->charge_per_node.begin(), this->charge_per_node.end(), generator);
            vector<int> prev_path_node = this->dijkstra(x, y).second;
            path = this->build_path(prev_path_node, x, y);
            this->charge_per_node = copy_charge_per_node;
        }

        for (int i = 0; i < (int)path.size(); i++)
        {
            // for every node on the path increment the charge (only useful for the dijkstra)
            // and add it to the added_paths set
            int v = path[i];
            this->charge_per_node[v] += this->charge;
            this->added_paths[v][id_path] = i;
        }
    }
    // modify the routing based on the added and erased paths
    this->modify_routing(this->added_paths, this->erased_paths);
}

void SymmetricRouting::update_routing()
{
    // Update all the dijkstra distances, only applied every temperture reduction.
    for (int i = 0; i < this->num_nodes; i++)
    {
        vector<int> distances_dijkstra = this->dijkstra(i, -1).first;
        for (int j = 0; j < this->num_nodes; j++)
            this->dijkstra_path_distances[i * this->num_nodes + j] = distances_dijkstra[j];
    }
}
void SymmetricRouting::forget_neighborhood_solution()
{
    this->modify_routing(this->erased_paths, this->added_paths);
}
pair<int, int> SymmetricRouting::evaluate()
{
    return {
        *max_element(this->total_charge_per_path.begin(), this->total_charge_per_path.end()),
        *max_element(this->charge_per_node.begin(), this->charge_per_node.end())
    };
}

vector<vector<vector<int>>> SymmetricRouting::get_standard_format_routing()
{

    vector<vector<vector<int>>> routing = vector<vector<vector<int>>>(
        this->num_nodes,
        vector<vector<int>>(
            this->num_nodes,
            vector<int>()));

    for (int id_path : this->paths_ids)
    {
        int x = id_path / this->num_nodes, y = id_path % this->num_nodes;
        vector<pair<int, int>> elements = vector<pair<int, int>>();

        for (int i = 0; i < this->num_nodes; i++)
        {
            auto r_element = this->R[i].find(id_path);
            if (r_element == this->R[i].end())
                continue;
            elements.push_back({r_element->second, i});
        }
        sort(elements.begin(), elements.end());

        vector<int> path = vector<int>(elements.size());
        for (int i = 0; i < (int)path.size(); i++)
            path[i] = elements[i].second;

        routing[x][y] = path;
        vector<int> inverse_path(path.rbegin(), path.rend());
        routing[y][x] = inverse_path;
    }
    return routing;
}
bool SymmetricRouting::validate_routing(vector<vector<vector<int>>> &routing, int expected_forwarding_diameter)
{

    int n = this->g.size();
    vector<int> charge_per_node = vector<int>(n, 0);
    vector<vector<int>> charge_per_path = vector<vector<int>>(n, vector<int>(n, 0));
    bool incorrect = false;
    vector<unordered_set<int>> graph_find = vector<unordered_set<int>>(this->num_nodes, unordered_set<int>());
    for(int i = 0; i < this->num_nodes; i++)
        for(int &node: this->g[i])
            graph_find[i].insert(node);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            int prev_node = i;
            for (int node : routing[i][j])
            {
                if (graph_find[prev_node].find(node) == graph_find[prev_node].end())
                    incorrect = true;
                
                charge_per_node[node] += 1;
                prev_node = node;
            }
            if (graph_find[prev_node].find(j) == graph_find[prev_node].end())
                incorrect = true;
            
            if (incorrect)
            {
                cout << "The path between the nodes: " << i << " and " << j << " is incorrect" << endl;
                for (int node : routing[i][j])
                    cout << node << " ";
                cout << endl;
                return false;
            }
        }
    }
    int forwarding_diameter = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int node : routing[i][j])
            {
                charge_per_path[i][j] += charge_per_node[node];
            }
            forwarding_diameter = max(forwarding_diameter, charge_per_path[i][j]);
        }
    }

    if (forwarding_diameter != expected_forwarding_diameter)
    {
        cout << "SymmetricRouting forwarding diameter " << forwarding_diameter << ", expected: " << expected_forwarding_diameter << endl;
        return false;
    }
    return true;
}