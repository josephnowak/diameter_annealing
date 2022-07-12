// #include <climits>
// #include <set>
// #include <random>
// #include <math.h>
// #include <numeric>
// #include <unordered_set>
// #include <unordered_map>
// #include <queue>
// #include <algorithm>
// #include <chrono>
// #include <map>
// #include <iostream>
// #include "simulated_annealing.hpp"

// using namespace std;
// using namespace std::chrono;

// #define um unordered_map
// #define us unordered_set

// // default_random_engine generator;
// random_device rd;
// mt19937 generator(rd());
// uniform_real_distribution<double> uniform_real(0, 1);


// template<typename T>
// SimulatedAnnealing::SimulatedAnnealing(
//     T &routing, 
//     double temperature, 
//     double min_temperature, 
//     double beta, 
//     int n_iterations = 100
// )
// {
//     this->temperature = temperature;
//     this->routing = routing;
//     this->best_routing = routing;
//     this->min_temperature = min_temperature;
//     this->beta = beta;
//     this->n_iterations = n_iterations;
// }

// /**
//  *
//  * @brief Simulated annealing algorithm to optimize routings based on any metric
//  */
// template<typename T>
// T SimulatedAnnealing::optimize()
// {
//     double temperature = this->temperature;
//     pair<int, int> act_cost = this->routing.evaluate();
//     pair<int, int> best_cost = act_cost;
//     while (temperature > this->min_temperature)
//     {
//         for (int i = 0; i < this->n_iterations; i++)
//         {
//             this->routing.neighborhood_solution();

//             pair<int, int> new_cost = this->routing.evaluate();

//             float dif_cost = ((act_cost.first - new_cost.first) / 2) / (float)this->routing.sample_size;
//             double x = uniform_real(generator);
//             if (dif_cost > 0)
//                 act_cost = new_cost;
//             else if (dif_cost == 0 and (act_cost.second >= new_cost.second))
//                 act_cost = new_cost;
//             else if (dif_cost < 0 and x < exp(dif_cost / temperature))
//                 act_cost = new_cost;
//             else
//             {
//                 // cout << "act_sol: " << act_cost.first << " new_cost: " << new_cost.first << endl;
//                 this->routing.forget_neighborhood_solution();
//             }

//             if (best_cost.first > new_cost.first or
//                 (new_cost.first == best_cost.first and (best_cost.second > new_cost.second)))
//             {
//                 this->best_routing = this->routing;
//                 best_cost = new_cost;
//             }
//         }
//         this->routing.update_routing();
//         cout << "function value: " << this->routing.evaluate().first << ", temperature: " << temperature << endl;
//         temperature = this->reduce_temperature(temperature);
//     }
//     return this->best_routing;
// }
// double SimulatedAnnealing::reduce_temperature(double temperature)
// {
//     return temperature - this->beta;
// }