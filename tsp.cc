/*
 * TSP.cc
 *
 *  Created on: Apr 17, 2012
 *      Author: Kelvin
 */
#include <stdio.h>
#include <stdlib.h>
#include "city.h"
#include "tour.h"

int num_cities;
int populatoin_size;
float greedy_selection_percentage;
int num_closer_way_points;
int group_size;
float mutation_percent;
int termination_step;
int** distance_matrix;
int ** closest_neighbors;
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
	closest_neightbors = new int[num_cities][num_closer_way_points];
	generate_distance_matrix();
	generate_initial_population();

}

/*
 * partition: used in quick_select
 */
int partition(int* input, int p, int r)
{
    int pivot = input[r];

    while ( p < r )
    {
        while ( input[p] < pivot )
            p++;

        while ( input[r] > pivot )
            r--;

        if ( input[p] == input[r] )
            p++;
        else if ( p < r ) {
            int tmp = input[p];
            input[p] = input[r];
            input[r] = tmp;
        }
    }

    return r;
}

/*
 * quickselect: finds the kth smallest value within the index of p and r of input
 */
int quick_select(int* input, int p, int r, int k)
{
	while(){
		if ( p == r )
			return input[p];
		int j = partition(input, p, r);
		int length = j - p + 1;
		if ( length == k )
			return input[j];
		else if ( k < length )
			r = j-1;
		else{
			k -= length;
			p = j + 1;
		}
	}
}

/*
 * generates a distance matrix in distance_matrix
 */
public generate_distance_matrix(){
	for(int i = 0; i < population_size; i++){

	}

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
 * find_n_closest_neighbors finds the n closest neighbors for all the cities
 * and place them in closest_neighbors 2d array
 */
private void generate_closest_neighbors(){
	for(int i = 0; i < population_size; i++){
		int closest_index = 0;
	    int m = quick_select(distance_matrix[i], 0, population_size - 1, num_closer_way_points);
	    closest_neighbors[i][closest_index] = m;
	    closest_index++;
		for(int k = 0; k < population_size; k++){
			if(distance_matrix[i][k] < m){
				closest_neighbors]i][closest_index] = distance_matrix[i][k];
				closest_index++;
				if(closest_index == num_closer_way_points){
					return;
				}
			}
		}
	}
}



/*
 * 	Tour generate_tour(): returns a tour based on the distance_matrix
	paper had 2 methods for picking the next waypoint in a tour that
	were influenced by the greedy_selection_percentage
 */
private tour generate_tour(){

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
