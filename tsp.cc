/*
 * pet
 * TSP.cc
 *
 * Created on: Apr 17, 2012
 * Author: Kelvin
 */
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "tsp.h"
#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
#include <time.h>
using std::vector;
using namespace std;

int num_cities = 1;
int population_size = 0;
int greedy_selection_percentage;
int num_closer_way_points;
int group_size;
int mutation_percentage;
int termination_step;
//int** distance_matrix;
float** distance_matrix;
int ** closest_neighbors;
tour* population; //Tour new_population[]




/*
 * partition: used in quick_select
 */
int array_partition(float* input, int p, int r)
{
	float pivot = input[r];
	
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
int array_quick_select(float* input, int p, int r, int k)
{
	while(1){
		if ( p == r )
			return input[p];
		int j = array_partition(input, p, r);
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
 * vector input version of partition
 */
int vector_partition(vector<int> input, int p, int r){
	float pivot = input[r];
	
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
 * vector version of quickselect
 */
int vector_quick_select(vector<int> input, int p, int r, int k)
{
	while(1){
		if ( p == r )
			return input[p];
		int j = vector_partition(input, p, r);
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
void generate_distance_matrix(){
	//TODO
	
}
/*
 * void generate_initial_population(): generates population by
 * calling generate tour many times,
 */
void generate_initial_population(){
	int *linear_cities = new int[num_cities];
	for(int i = 0; i < num_cities; i++){
		linear_cities[i] = i;
	}
	for(int j =0; j < population_size; j++){
		generate_tour(linear_cities, j);
		/*
		 *	cout << "tour number " << j << ": " << endl;
		 *	for (int k = 0; k < num_cities; k++) {
		 *	cout << population[j].tour[k] << " ";
		 }
		 cout << endl;
		 */
	}
}

/*
 * find_n_closest_neighbors finds the n closest neighbors for all the cities
 * and place them in closest_neighbors 2d array NOT WORKING AS INTENDED, QUICK_SELECT 
 * RETURNS VALUE NOT INDEX NEED TO FIND INDEX VALUE
 */
void generate_closest_neighbors(){
	for(int i = 0; i < population_size; i++){
		int closest_index = 0;
		//num_closer_way_points + 1 used as quick_select will pick out i as well, since distance is 0 with itself
		float m = array_quick_select(distance_matrix[i], 0, population_size - 1, num_closer_way_points + 1);
		
		closest_neighbors[i][closest_index] = m;
		closest_index++;
		for(int k = 0; k < population_size; k++){
			if(distance_matrix[i][k] < m && distance_matrix[i][k] != 0 ){
				closest_neighbors[i][closest_index] = distance_matrix[i][k];
				closest_index++;
				if(closest_index == num_closer_way_points){
					return;
				}
			}
		}
	}
}



/*
 * Tour generate_tour(): returns a tour based on the distance_matrix
 * paper had 2 methods for picking the next waypoint in a tour that
 * were influenced by the greedy_selection_percentage
 */
void generate_tour(int* linear_cities, int index){
	tour new_tour;
	new_tour.fitness = 0;
	new_tour.tour = new int[num_cities];
	new_tour.tour_lengths = new float[num_cities];
	//cout << "tour number: " << index << endl;
	//vector of available_cities, cities will get deleted from this resizable array as they get added to the tour
	//vector<int> available_cities(linear_cities, linear_cities + sizeof(linear_cities)/sizeof(int));
	vector<int> available_cities(linear_cities, linear_cities + num_cities);
	
	//first city in tour will always be the first city (doesn't matter where we start as tours will loop through all cities
	
	new_tour.tour[0] = 0;
	int current_city = 0;
	available_cities.erase(available_cities.begin());
	int next_city = 0;
	
	//loop num_cities times to form the tour
	for(int i = 1; i < num_cities; i++){
		int selection = (rand() % 100) + 1;
		//cout << "213 selection percent " << selection << endl;
		//choose greedily
		//if greedy strategy selected and the number of cities available is greater than the number of closer way points
		//perform the greedy strategy
		if((selection < greedy_selection_percentage) and (available_cities.size() > num_closer_way_points)){
			//cout << "217 greedy selection" << endl;
			//TODO: check that the ranodm selection producing desired values
			int random_closest= rand() % num_closer_way_points;
			next_city = vector_quick_select(available_cities, 0, available_cities.size() - 1, random_closest);
		}else // choose next city randomly
		{
		//cout << "224 city chosen randomly" << endl;
		int city_index = rand() % available_cities.size();
		//cout << "227 next city " << city_index << endl;
		next_city = available_cities[city_index];
		}
		//loop through the available_cities vector until u found the city value selected to be added
		//delete the value from the available_cities vector
		for(int j = 0; j < available_cities.size(); j++){
			if(available_cities[j] == next_city){
				available_cities.erase(available_cities.begin() + j);
				break;            
			}
		}
		
		// add the city to the tour, calculate the fitness it adds, and add the tour_length as well
		
		//cout << next_city << endl;
		new_tour.tour[i] = next_city;
		new_tour.tour_lengths[i-1] = distance_matrix[current_city][next_city];
		new_tour.fitness += distance_matrix[current_city][next_city];
		current_city = next_city;
	}
	
	//compute the final fitness and tour_length connecting the final city to the first city
	new_tour.tour_lengths[num_cities-1] = distance_matrix[current_city][0];
	new_tour.fitness += distance_matrix[current_city][next_city];
	cout << endl;
	population[index] = new_tour;
}

/*
 * DEPRECATED, method's funcionality built into generate_tour
 * int greedy_selection(int city): greedy method for selecting next
 * waypoint in tour in generateTour
 */
int greedy_selection(int city){
	return 0;
}

/*
 * DEPRECATED, method's functionality build into generate_tour
 * int random_selection(int city): randomly pick 2 waypoints among
 * list of cities (not including the city entered into the method,
 * then pick the shorter cost of the two
 */
int random_selection(int city){
	return 0;
	
}


/*
 * DEPRECATED, just uses quicksort
 * select_group(int group_size): quicksorts the population
 * group_size of the most optimal tours (run sort_population first
 * to figure out best tours)
 */
void select_group(int group_size){
	
}


//given a tour with initialized tour array, returns fitness of tour
int compute_fitness(tour t){
	int fitness = 0;
	for(int i = 0; i < num_cities -1; i ++){
		fitness = distance_matrix[t.tour[i]][t.tour[i+1]];
	}
	return fitness;
}



/*
 * Tour crossover(Tour parent1, Tour Parent2): Crossover of
 * 2 parents and then compute the length
 */
void crossover(tour parent1, tour parent2, tour* children, int index){
	tour child1;
	tour child2;
	child1.fitness = 0;
	child1.tour = new int[num_cities];
	child1.tour_lengths = new float[num_cities];
	child2.fitness = 0;
	child2.tour = new int[num_cities];
	child2.tour_lengths = new float[num_cities];
	
	for (int k = 0; k < num_cities; k++) {
		child1.tour[k] = -1;
		child2.tour[k] = -1;
	}
	
	
	child1.tour[0] = parent1.tour[0];
	child2.tour[0] = parent2.tour[0];
	int p2city = parent2.tour[0];
	
	
	while(1) {
		bool visited = false;
		for (int j = 0; j < num_cities; j++) {
			if (child1.tour[j] == p2city) {
				visited = true;
			}
		}
		if (visited) {
			break;
		}
		
		//since p2city hasn't yet been visited by child1, find where p2city occurs in parent1 and insert it into child1 at the same index
		for (int j = 0; j < num_cities; j++) {
			if (parent1.tour[j] == p2city) {
				child1.tour[j] = p2city;
				child2.tour[j] = parent2.tour[j];
				p2city = parent2.tour[j];
			}
		}
	}
	
	//fill in the -1 values in the children
	for (int k = 0; k < num_cities; k++) {
		if (child1.tour[k] == -1) {
			child1.tour[k] = parent2.tour[k];
			child2.tour[k] = parent1.tour[k];
		}
	}
	child1.fitness = compute_fitness(child1);
	child2.fitness = compute_fitness(child2);
	children[index] = child1;
	children[index + 1] = child2;
	
}



/*
 * void sort_population(): sorts population so that the GA can
 * organize (quicksort)
 */
void sort_population(){
	qsort_population(0, population_size - 1, population);
}

void qsort_population(int left, int right, tour* population) {
	if (right > left) {
		int pivotIndex = rand() % (right - left + 1);
		tour pivot = population[left + pivotIndex];
		int pivotfitness = pivot.fitness;
		population[left + pivotIndex] = population[right];
		population[right] = pivot;
		
		int i = left - 1;
		int j = right;
		
		do {
			do { i++; } while (population[i].fitness < pivotfitness);
			do { j--; } while (population[j].fitness > pivotfitness && j > left);
			if (i < j) {
				tour ith = population[i];
				population[i] = population[j];
				population[j] = ith;
			}
			
		} while (i < j);
		
		population[right] = population[i];
		population[i] = pivot;
		qsort_population(left, i - 1, population);
		qsort_population(i + 1, right, population);
	}
}



tour* create_children(){
	/*cout << "population before sorting\n";
	for (int i = 0; i < population_size; i++) {
		cout << population[i].fitness << " ";
	}
	cout << endl;
	qsort_population(0, num_cities - 1, population);
	cout << "population after sorting\n";
	for (int i = 0; i < population_size; i++) {
		cout << population[i].fitness << " ";
	}
	cout << endl;
	cout << "\n" << endl;  
	*/
	tour* children = new tour[group_size];
	for (int i = 0; i < group_size; i += 2){
		crossover(population[i], population[group_size-i], children, i);
		int mutate_or_not = rand() % 100 + 1;
		if(mutate_or_not < mutation_percentage){
			mutate(&children[i]);
		}
		children[i].fitness = compute_fitness(children[i]);
	}
	return children;
}


/*
 * Tour mutate(Tour t ): mutates the tour and returns the
 * optimal tour of the mutated or original
 */

void mutate(tour* t ){
	int start;
	int end;
	do {
		start = rand() % num_cities;
		end = rand() % num_cities;
	} while ((start >= end) or (start == 0 and end == (num_cities - 1)));
	
	while (start < end) {
		int temp = t->tour[start];
		t->tour[start] = t->tour[end];
		t->tour[end] = temp;
		start++;
		end--;
	}
}










/*
 * create_new_generation sorts the current population tour array and the children tour array, and then proceeds to replace
 * the group_size weakest tours in the population array with the children if the fitness of the child is higher
 */
tour* create_new_generation(tour* population, tour* children, int population_size, int group_size){
	qsort_population(0, population_size -1, population);
	qsort_population(0, group_size -1, children);
	int population_index = population_size - group_size;
	int children_index = 0;
	/*
	 * this for loop only selects the group_size best between group_size lowest from previous generation and 
	 * the new children
	 */
	for(int i = 0; i < group_size && population_index < population_size; i ++){
		if(children[children_index].fitness > population[population_index].fitness){
			population[population_index] = children[children_index];
			population_index ++;
		}
		children_index++;
	}


}



/*
 * Tour run(): generates initial population, then runs through
 * a for loop that goes through the genetic algorithm until the
 * termination_step, return the best population
 */
void run_genetic_algorithm(){
	/*distance_matrix = new int*[num_cities];
	 *	for(int i = 0; i < num_cities; i ++){
	 *	distance_matrix[i] = new int[num_cities];
	 }
	 */
	
	population = new tour[population_size];
	closest_neighbors = new int*[num_cities];
	for(int i = 0; i < num_cities; i ++){
		closest_neighbors[i] = new int[num_closer_way_points];
		
	}
	
	generate_initial_population();
	for(int j = 0; j < termination_step; j++){
		tour* children = create_children();
		create_new_generation(population, children, population_size, group_size);
	}
	print_best_tour();
}

/*
 * print_best_tour: prints out the best tour along with its fitness, run at the end of
 * run_genetic_algorithm()
 */
void print_best_tour(){
	qsort_population(0, num_cities - 1, population);
	cout << "Best Tour Generated After " << termination_step << "generations \n";
	cout << "Fitness: " << population[0].fitness << endl << "Tour \n";
	for (int i = 0; i < num_cities; i++){
		cout << population[0].tour[i] << endl;
	}
	
}


/*
 * program requires all these parameters
 * program_name tsp_file num_cities population_size greedy_selection_percentage num_closer_way_points group_size
 * mutation_percentage termination_step
 */
int main(int argc, char** argv){
	if(argc != 9){
		cout << "not enough arguments\n" << " program requires 9 arguments "
		<< "program_name tsp_file num_cities population_size greedy_selection_percentage"
		<< "num_closer_way_points group_size mutation_percentage termination_step";
		return 0;
		
		
		//./a.out ./burma14.tsp 14 100 20 5 20 50 10
		
	}
	/*
	 *	num_cities = 10;
	 *	population_size = 10;
	 *	greedy_selection_percentage = 80;
	 *	num_closer_way_points = 3;
	 *	group_size = 3;
	 *	mutation_percentage = 70;
	 *	termination_step = 20;
	 */
	
	char* filename = argv[1];
	const char * num_cities_string = argv[2];
	num_cities = atoi(num_cities_string);
	population_size = atoi(argv[3]);
	greedy_selection_percentage = atoi(argv[4]);
	num_closer_way_points = atoi(argv[5]);
	group_size = atoi(argv[6]);
	mutation_percentage = atoi(argv[7]);
	termination_step = atoi(argv[8]);
	distance_matrix = new float*[num_cities];
	for (int i = 0; i< num_cities; i++) {
		distance_matrix[i] = new float[num_cities];
	}
	//initialize the random seed, ONLY CALL ONCE in program 
	srand(time(NULL));
	
	coordinates cities[num_cities];
	
	char line[256];
	ifstream myfile;
	myfile.open(filename);
	
	const char * lastline = "NODE_COORD_SECTION";
	while (1) {
		myfile.getline(line, 256);
		//cout << "asdf\n";
		if (!strcmp(line, lastline)) {
			break;
		}
	}
	
	for (int i = 0; i < num_cities; i++) {
		myfile.getline(line, 256);
		char * pch;
		pch = strtok(line, " ");
		pch = strtok(NULL, " ");
		cities[i].x = atoi(pch);
		pch = strtok(NULL, " ");
		cities[i].y = atoi(pch);
	}
	
	for (int i = 0; i < num_cities; i++) {
		for (int j = 0; j < num_cities; j++) {
			distance_matrix[i][j] = sqrt(pow(cities[i].x - cities[j].x,2) + pow(cities[i].y - cities[j].y,2));
		}
	}
	
	myfile.close();
	
	cout << "file parsing reached" << endl;
	run_genetic_algorithm();
	for (int i = 0; i < num_cities; i++) {
		delete[] distance_matrix[i];
	}
	delete[] distance_matrix;
}
