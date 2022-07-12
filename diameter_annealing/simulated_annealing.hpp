#ifndef _SIMULATED_ANNEALING_H_
#define _SIMULATED_ANNEALING_H_


using namespace std;


template<typename T>
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
	T routing;

	/**
	 * @brief Best routing found during the execution of the algorithm
	 */
	T best_routing;


	SimulatedAnnealing(T &routing, double temperature, double min_temperature, double beta, int n_iterations = 100);

	/**
	 * 
	 * @brief Simulated annealing algorithm to optimize routings based on any metric
	 */
	T optimize();


	double reduce_temperature(double temperature);
};


#endif
