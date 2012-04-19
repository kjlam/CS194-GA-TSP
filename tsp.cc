/*
 * TSP.cc
 *
 *  Created on: Apr 17, 2012
 *      Author: Kelvin
 */

#include "city.h"
#include "tour.h"

int num_cities;
int populatoin_size;
float greedy_selection_percentage;
int num_closer_way_points;
int group_size;
float mutation_percent;
int termination_step
int** distance_matrix;
tour* population;
	//Tour new_population[]

/*
 * Tour run():  generates initial population, then runs through
 *  a for loop that goes through the genetic algorithm until the
 *  termination_step, return the best population
 */
public void run(){
	num_cities = 10;
	population_size = 10;
	greedy_selection_percentage = 0.8;
	num_closer_way_points = 3;
	group_size = 3;
	mutation_percent = 0.7;
	termination_step = 20;
	distance_matrix = new int[num_cities][num_cities];
	population = new population[population_size];
	generate_distance_matrix();
	generate_initial_population();
}

/*
 * generates a distance matrix in distance_matrix
 */
public generate_distance_matrix(){

}
/*
 * 	void generate_initial_population():  generates population by
 * 	calling generate tour many times,
 */
public void generate_initial_population(){
	for(int i =0; i < population_size; i++){
		population[i] = generate_tour();
	}
}

/*
 * 	Tour generate_tour(): returns a tour based on the distance_matrix
	paper had 2 methods for picking the next waypoint in a tour that
	were influenced by the greedy_selection_percentage
 */
public void generate_tour(){

}

/*
 * 	int greedy_selection(int city): greedy method for selecting next
 * 	waypoint in tour in generateTour
 */
private int greedy_selection(int city){

}

/*
 * int random_selection(int city): randomly pick 2 waypoints among
 *  list of cities (not including the city entered into the method,
 *  then pick the shorter cost of the two
 */
private int random_selection(int city){

}

/*
 * 	Tour[] select_group(int group_size): returns an array of length
 * 	 group_size of the most optimal tours (run sort_population first
 * 	 to figure out best tours)
 */
tour* select_group(int group_size){

}

/*
 * void sort_population(): sorts population so that the GA can
 * organize (quicksort)
 */
void sort_population(){

}


tour* create_children(){
}

/*
 * 	Tour crossover(Tour parent1, Tour Parent2): Crossover of
 * 	 2 parents and then compute the length
 */
tour crossover(tour parent1, tour parent2){

}

/*
 * 	Tour mutate(Tour t ): mutates the tour and returns the
 * 	optimal tour of the mutated or original
 */

tour mutate(tour t ){

}
