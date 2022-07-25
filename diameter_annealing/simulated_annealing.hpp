#ifndef _SIMULATED_ANNEALING_H_
#define _SIMULATED_ANNEALING_H_


#include "routing.hpp"


class SimulatedAnnealing
{

public:
	/**
	 * @brief Initial temperature of the system.
	 */
	double temperature;
	/**
	 * @brief Reduction factor explained on reduce_temperature method
	 */
	double beta; 
	/**
	 * @brief Minimum temperature of the system.
	 */
	double min_temperature;
	/**
	 * @brief Number of repetitions to find a neighborhood solution before reduce the temperature
	 */
	int n_iterations;
	/**
	 * @brief Routing
	 */
	BaseRouting * routing;

	/**
	 * @brief Only useful for benchmarking
	 */
	double initial_routing_diameter;

	/**
	 * @brief Indicate if the algorithm should print the progress
	 */
	bool verbose = true;


	SimulatedAnnealing(
        BaseRouting &routing, 
        double temperature, 
        double min_temperature, 
        double beta, 
        int n_iterations = 100,
		bool verbose=true
    );

	/**
	 * 
	 * @brief Simulated annealing algorithm to optimize routings based on any metric
	 */
	void optimize();
    
    double reduce_temperature(double temperature);
};


#endif