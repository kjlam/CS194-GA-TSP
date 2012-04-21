/*
 * TSP.cc
 *
 *  Created on: Apr 17, 2012
 *      Author: Kelvin
 */
#include <stdio.h>
#include <stdlib.h>


int num_cities;
int population_size;
int greedy_selection_percentage;
int num_closer_way_points;
int group_size;
float mutation_percent;
int termination_step;
int** distance_matrix;
int ** closest_neighbors;
tour* populationu	//Tour new_population[]

struct tour{
	int fitness;
	int* tour;
	int *tour_lengths;
};

struct city
{
	int cityNum;
	int* closest_neighbors;

};


/*
 * Tour run():  generates initial population, then runs through
 *  a for loop that goes through the genetic algorithm until the
 *  termination_step, return the best population
 */
public void run(){
	num_cities = 10;
	population_size = 10;
	greedy_selection_percentage = 80;
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
*vector input version of partition
*/
int partition(vector<input> input, int p, int r){
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
*vector version of quickselect
*/
int quick_select(vector<int> input, int p, int r, int k)
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
	int *linear_cities = new int[num_cities];
	for(int i = 0; i < num_cities; i++){
		linear_cities[i] = i;
	}
	for(int j =0; j < population_size; j++){
		population[j] = generate_tour(linear_cities);
		
	}
}

/*
 * find_n_closest_neighbors finds the n closest neighbors for all the cities
 * and place them in closest_neighbors 2d array
 */
private void generate_closest_neighbors(){
	for(int i = 0; i < population_size; i++){
		int closest_index = 0;
		//num_closer_way_points + 1 used as quick_select will pick out i as well, since distance is 0 with itself
	    int m = quick_select(distance_matrix[i], 0, population_size - 1, num_closer_way_points + 1);
	    closest_neighbors[i][closest_index] = m;
	    closest_index++;
		for(int k = 0; k < population_size; k++){
			if(distance_matrix[i][k] < m && distance_matrix[i][k] != 0 ){
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
private tour generate_tour(int* linear_cities){
	tour new_tour;
	new_tour.fitness = 0;
	new_tour.tour = new int[num_cities];
	new_tour.tour_lengths = new int[num_cities];
	vector<int> available_cities(linear_cities, linear_cities + sizeof(linear_cities)/sizeof(int));
	new_tour.tour[0] = 0;
	int current_city = 0;
	for(int i = 1; i < num_cities; i ++){
		if(i == 0){
			available_cities.erase(
		}
		int selection = rand() % num_cities;
		//choose greedily 
		if(selection < greedy_selection_percentage && available_cities.size > num_closer_way_points){
			//TODO: check that the ranodm selection producing desired values
			int random_closest= rand() % (num_closer_way_points); 
			int next_city = quickselect(available_cities, 0, availabe_cities.size() - 1, random_closest);
			
		}else // choose next city randomly 
		{
			int next_city = rand() % available_cities.size();
		}
		for(int j = 0; j < available_cities.size; j++){
			if(available_cities[j] == next_city){
				available_cities.erase(available_cities.begin() + j);
				break;
			}
		}
		new_tour.tour[i] = next_city;
		new_tour.tour_lengths[i-1] = distance_matrix[current_city][next_city];
		new_tour.fitness += distance_matrix[current_city][next_city];
		current_city = next_city;
	}
	new_tour.tour_length[i] = distance_matrix[current_city][0];
	new_tour.fitness += distance_matrix[current_city][next_city];
	return tour;
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
 * DEPRECATED, just uses quicksort 
 * 	select_group(int group_size): quicksorts the population 
 * 	 group_size of the most optimal tours (run sort_population first
 * 	 to figure out best tours)
 */
private void select_group(int group_size){

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

/*
* create_new_generation sorts the current population tour array and the children tour array, and then proceeds to replace
the group_size weakest tours in the population array with the children if the fitness of the child is higher
*/
private tour* create_new_generation(tour* population, tour* children, int population_size, int group_size){
	qsort_population(0, population_size -1, population);
	qsort_population(0, group_size -1, children);
	int population_index = population_size - group_size;
	int children_index = 0;
	for(int i = 0; i < group_size && population_index < population_size; i ++){
		if(children[children_index].fitness > population[population_index].fitness){
			population[population_index] = children[children_index];
			
		}
		children_index++;
		population_index++;
	}
}
