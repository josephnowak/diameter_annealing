
#include "simulated_annealing.hpp"
#include "routing.hpp"
#include <climits>
#include <assert.h>
#include <chrono>


SimulatedAnnealing::SimulatedAnnealing(
    BaseRouting &routing,
    double temperature, 
    double min_temperature, 
    double beta,
    int n_iterations,
    bool verbose
)
{
    this->temperature = temperature;
    this->routing = &routing;
    this->min_temperature = min_temperature;
    this->beta = beta;
    this->n_iterations = n_iterations;
    this->verbose = verbose;
}

void SimulatedAnnealing::optimize()
{
    double temperature = this->temperature;
    pair<int, int> act_cost = this->routing->evaluate();
    if(this->verbose)
        cout << "Initial routing created, starting the algorithm" << endl;
    this->initial_routing_diameter = act_cost.first;
    pair<int, int> best_cost = {INT_MAX, INT_MAX};
    while (temperature > this->min_temperature)
    {
        for (int i = 0; i < this->n_iterations; i++)
        {
            this->routing->neighborhood_solution();

            pair<int, int> new_cost = this->routing->evaluate();

            float dif_cost = act_cost.first - new_cost.first;
            double x = this->routing->get_random_probability();

            if (act_cost >= new_cost)
                act_cost = new_cost;
            else if (x < exp(dif_cost / temperature))
                act_cost = new_cost;
            else
                this->routing->forget_neighborhood_solution();
            
            if(best_cost > new_cost)
            {
                this->routing->store_routing();
                best_cost = new_cost;
            }
        }    
        this->routing->update_routing();
        if(this->verbose)
            cout << "function value: " << this->routing->evaluate().first << ", temperature: " << temperature << endl;
        temperature = this->reduce_temperature(temperature);
    }
    this->routing->set_stored_routing();
}

double SimulatedAnnealing::reduce_temperature(double temperature)
{
    return temperature - this->beta;
}
