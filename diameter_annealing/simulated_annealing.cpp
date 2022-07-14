
#include "simulated_annealing.hpp"
#include "routing.hpp"


SimulatedAnnealing::SimulatedAnnealing(
    BaseRouting &routing,
    double temperature, 
    double min_temperature, 
    double beta,
    int n_iterations
)
{
    this->temperature = temperature;
    this->routing = &routing;
    this->min_temperature = min_temperature;
    this->beta = beta;
    this->n_iterations = n_iterations;
}

void SimulatedAnnealing::optimize()
{
    double temperature = this->temperature;
    pair<int, int> act_cost = this->routing->evaluate();
    pair<int, int> best_cost = act_cost;
    while (temperature > this->min_temperature)
    {
        for (int i = 0; i < this->n_iterations; i++)
        {
            this->routing->neighborhood_solution();

            pair<int, int> new_cost = this->routing->evaluate();

            float dif_cost = act_cost.first - new_cost.first;
            double x = this->routing->get_random_probability();

            if (dif_cost > 0)
                act_cost = new_cost;
            else if (dif_cost == 0 and (act_cost.second >= new_cost.second))
                act_cost = new_cost;
            else if (dif_cost < 0 and x < exp(dif_cost / temperature))
                act_cost = new_cost;
            else
            {
                // cout << "act_sol: " << act_cost.first << " new_cost: " << new_cost.first << endl;
                this->routing->forget_neighborhood_solution();
            }

            if (best_cost.first > new_cost.first or
                (new_cost.first == best_cost.first and (best_cost.second > new_cost.second)))
            {
                this->routing->store_routing();
                best_cost = new_cost;
            }
        }
        this->routing->update_routing();
        cout << "function value: " << this->routing->evaluate().first << ", temperature: " << temperature << endl;
        temperature = this->reduce_temperature(temperature);
    }
    this->routing->set_stored_routing();
}

double SimulatedAnnealing::reduce_temperature(double temperature)
{
    return temperature - this->beta;
}
