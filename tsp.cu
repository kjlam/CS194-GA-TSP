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
#include <sys/time.h>
#include <time.h>
#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/sort.h>

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
float* distance_matrix;
float* d_distance_matrix;
int ** closest_neighbors;
tour* population; //Tour new_population[]

tour *h_children;
tour *h_parent_set1;
tour *h_parent_set2;
tour* h_index_of_parent_set1;
tour* h_index_of_parent_set2;
tour* h_cycles;

tour* d_children;
tour* d_parent_set1;
tour* d_parent_set2;
tour* d_index_of_parent_set1;
tour* d_index_of_parent_set2;
tour* d_cycles;
tour* d_population;

thrust::device_vector<tour> thrust_population; 
thrust::device_vector<tour> thrust_children;


double
timestamp (){
	struct timeval tv;
	gettimeofday (&tv, 0);
	return tv.tv_sec + 1e-6*tv.tv_usec;
}



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
 * DEPRECATED find_n_closest_neighbors finds the n closest neighbors for all the cities
 * and place them in closest_neighbors 2d array NOT WORKING AS INTENDED, QUICK_SELECT 
 * RETURNS VALUE NOT INDEX NEED TO FIND INDEX VALUE
 *
 * void generate_closest_neighbors(){
 *	for(int i = 0; i < population_size; i++){
 *		int closest_index = 0;
 *		//num_closer_way_points + 1 used as quick_select will pick out i as well, since distance is 0 with itself
 *		float m = array_quick_select(distance_matrix[i], 0, population_size - 1, num_closer_way_points + 1);
 *		
 *		closest_neighbors[i][closest_index] = m;
 *		closest_index++;
 *		for(int k = 0; k < population_size; k++){
 *			if(distance_matrix[i][k] < m && distance_matrix[i][k] != 0 ){
 *				closest_neighbors[i][closest_index] = distance_matrix[i][k];
 *				closest_index++;
 *				if(closest_index == num_closer_way_points){
 *					return;
 }
 }
 }
 }
 }
 */


/*
 * Tour generate_tour(): returns a tour based on the distance_matrix
 * paper had 2 methods for picking the next waypoint in a tour that
 * were influenced by the greedy_selection_percentage
 */
void generate_tour(int* linear_cities, int index){
	tour new_tour;
	new_tour.fitness = 0;
	//new_tour.tour = new int[num_cities];
	//new_tour.tour_lengths = new float[num_cities];
	//cout << "tour number: " << index << endl;
	//vector of available_cities, cities will get deleted from this resizable array as they get added to the tour
	//vector<int> available_cities(linear_cities, linear_cities + sizeof(linear_cities)/sizeof(int));
	vector<int> available_cities(linear_cities, linear_cities + num_cities);
	
	//first city in tour will always be the first city (doesn't matter where we start as tours will loop through all cities
	
	new_tour.path[0] = 0;
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
			//cout<< "227 next city " << city_index << endl;
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
		new_tour.path[i] = next_city;
		//new_tour.tour_lengths[i-1] = distance_matrix[current_city*num_cities + next_city];
		new_tour.fitness += distance_matrix[current_city * num_cities + next_city];
		current_city = next_city;
	}
	
	//compute the final fitness and tour_length connecting the final city to the first city
	//new_tour.tour_lengths[num_cities-1] = distance_matrix[current_city * num_cities + 0];
	new_tour.fitness += distance_matrix[current_city * num_cities + next_city];
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


/*//given a tour with initialized tour array, returns fitness of tour
float compute_fitness(tour t){
	float fitness = 0;
	for(int i = 0; i < num_cities -1; i ++){
		//cout << "305 " << t.path[i] << " " << t.path[i+1];
		fitness += distance_matrix[t.path[i] * num_cities + t.path[i+1]];
		if(t.path[i] != -1 and t.path[i+1] != -1){
		//	cout << "299 " <<  t.path[i]  << " " << t.path[i+1] << " " << endl;
		}
	}
	fitness += distance_matrix[t.path[num_cities-1] * num_cities + t.path[0]];
	return fitness;
}
*/

__global__ static void
	compute_fitness(tour* children, float* d_matrix, int num_cities){
		int id = blockIdx.x*blockDim.x + threadIdx.x;
		children[id].fitness = 0;
		for(int i = 0; i < num_cities - 1; i++){
			children[id].fitness += d_matrix[ children[id].path[i] * num_cities + children[id].path[i+1]];
		}
		children[id].fitness += d_matrix[ children[id].path[num_cities - 1] * num_cities + children[id].path[0]];
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



struct tour_pair{
	tour t1; 
	tour t2;
};


/*
 * Tour crossover(Tour parent1, Tour Parent2): Crossover of
 * 2 parents and then compute the length
 */
__global__ static void
	crossover(tour* parent1, tour* parent2, tour* children, tour* indexOfParent1, tour* indexOfParent2, tour* cycles, int num_cities, int group_size){
		int instance = blockIdx.x*blockDim.x + threadIdx.x;
		if(instance >= group_size/2){
			return;
		}

		for (int k = 0; k < num_cities; k++) {
			//children[2*instance].path[k] = -1;
			//children[2*instance + 1].path[k] = -1;
			cycles[instance].path[k] = -1;
		}
		
		int cycle_index = 0;
		
		int p1Index = parent1[instance].path[0];
		cycles[instance].path[cycle_index] = p1Index;
		cycle_index = cycle_index + 1;

		for(int i = 0; i < num_cities; i ++){
			int p1_value_at_i = parent1[instance].path[i];
			int p2_value_at_i = parent2[instance].path[i];
			indexOfParent1[instance].path[p1_value_at_i] = i;
			indexOfParent2[instance].path[p2_value_at_i] = i;
			children[2*instance].path[i] = p1_value_at_i;
			children[2*instance + 1].path[i] = p2_value_at_i;
		}
		
	
		
		children[2*instance].fitness = 0;
		children[2*instance + 1].fitness = 0;
		
		//int p2city = indexOfParent2[instance].path[p1Index];
		int p2city = children[2*instance +1].path[p1Index];
		 p1Index = indexOfParent1[instance].path[p2city];
		
		while(p1Index != cycles[instance].path[0]){
			cycles[instance].path[cycle_index] = p1Index;
			cycle_index++;
			p2city = children[2*instance+1].path[p1Index];
			p1Index = indexOfParent1[instance].path[p2city];
		}
		
		for(int i = 0; i < num_cities; i++){
			if(cycles[instance].path[i] != -1){
				int reverse = cycles[instance].path[i];
				int tmp = children[2*instance].path[reverse];
				children[2*instance].path[reverse] = children[2*instance+1].path[reverse];
				children[2*instance+1].path[reverse] = tmp;
			}
		}
			

		/*
		while(1) {
			bool visited = false;
			for (int j = 0; j < num_cities; j++) {
				if (children[2*instance].path[j] == p2city) {
					visited = true;
				}
			}
			if (visited) {
				break;
			}
			
			//since p2city hasn't yet been visited by child1, find where p2city occurs in parent1 and insert it into child1 at the same index
			for (int j = 0; j < num_cities; j++) {
				if (parent1[instance].path[j] == p2city) {
					children[2*instance].path[j] = p2city;
					children[2*instance + 1].path[j] = parent2[instance].path[j];
					p2city = parent2[instance].path[j];
				}
			}
		}
		
		//fill in the -1 values in the children
		for (int k = 0; k < num_cities; k++) {
			if (children[2*instance].path[k] == -1) {
				children[2*instance].path[k] = parent2[instance].path[k];
				children[2*instance + 1].path[k] = parent1[instance].path[k];
			}
		}
		//can’t compute fitness here as requires distance_matrix
		//child1.fitness = compute_fitness(child1);
		//child2.fitness = compute_fitness(child2);
		//children[index] = child1;
		//children[index + 1] = child2;
		*/
}
	


/*
 * parallalel mutate 
 */

__global__ void parallel_mutate(tour* d_children, int* d_mutate_indices, int group_size, int num_cities){

        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int tid = threadIdx.x;
       
        //tour tours[32];
        //tours[tid] = d_children[i];
       
        int num1 = d_mutate_indices[2*blockIdx.x] % num_cities;
        int num2 = d_mutate_indices[2*blockIdx.x+1] % num_cities;
		
		int start = min(num1, num2);
		int end = max(num1, num2);

        if (i < group_size) {
                
                while (start < end) {
                        int temp = d_children[i].path[start];
                        d_children[i].path[start] = d_children[i].path[end];
                        d_children[i].path[end] = temp;
                        start++;
                        end--;
                }
                
               /* while (start < end) {
                        int temp = tours[tid].path[start];
                        tours[tid].path[start] = tours[tid].path[end];
                        tours[tid].path[end] = temp;
                        start++;
                        end--;
                }
               */
               // d_children[i] = tours[tid];
        }
}

__global__ static void 
	select_group(tour* d_p, tour* d_p1, tour* d_p2){
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		
		d_p1[id]= d_p[id*2];
		d_p2[id] = d_p[id*2 + 1];
		
	}


void create_children(){
	/*cout << "population before sorting\n";
	for (int i = 0; i < population_size; i++) {
		cout << population[i].fitness << " ";
	}
	cout << endl;
	*/
	//qsort population (parallelize?)
	//qsort_population(0, num_cities - 1, population);
	/*thrust::host_vector <tour> h_population(population_size);
	for(int k = 0; k < population_size; k++){
		h_population[k] = population[k];
	}
	thrust::device_vector<tour> d_population = h_population;
	thrust::sort(d_population.begin(), d_population.end());
	thrust::copy(d_population.begin(), d_population.end(), h_population.begin());
	*/
	
	thrust::sort(thrust_population.begin(), thrust_population.end());
	
	/*
	for(int m = 0; m < population_size; m++){
		population[m] = h_population[m];
	}
	*/
	
	int num_threads = group_size/2;
	int num_blocks = (num_threads + 31)/32;
	dim3 block(32, 1);
	dim3 grid(num_blocks, 1);
	
	select_group<<<grid, block>>>(d_population, d_parent_set1, d_parent_set2);
	//cout << "before crossover" << endl;
	//apply crossover on adjacent pairs of elements in the parent set
	//thrust::transform(d_parent_set1.begin(), d_parent_set1.end(), d_parent_set2.begin(), d_children.begin(), crossover_functor(num_cities));
	crossover<<<grid, block>>>(d_parent_set1, d_parent_set2, d_children, d_index_of_parent_set1, d_index_of_parent_set2, d_cycles, num_cities, group_size);
	
	
	num_threads = group_size;
	num_blocks = (num_threads + 31)/32;
	dim3 block1(32, 1);
	dim3 grid1(num_blocks, 1);
	int* h_random_arr = new int[num_blocks*2];
	for(int i = 0; i < num_blocks * 2; i++){
		h_random_arr[i] = rand();
	}
	
	int* d_random_arr = 0;
	
	cudaMalloc((void**)&d_random_arr, sizeof(int) * num_blocks * 2);
	cudaMemcpy(d_random_arr, h_random_arr, sizeof(int) * num_blocks * 2, cudaMemcpyHostToDevice);
	
	parallel_mutate<<<grid1, block1>>> (d_children, d_random_arr, group_size, num_cities);
	
	//cout << "after crossover" << endl;
	//cudaMemcpy(h_children, d_children, sizeof(tour)*group_size, cudaMemcpyDeviceToHost);
	/*
	//cudaFree(d_children);
	cudaFree(d_parent_set1);
	cudaFree(d_parent_set2);
	cudaFree(d_index_of_parent_set1);
	cudaFree(d_index_of_parent_set2);
	cudaFree(d_cycles);
	*/
	
	//cout << "525 after cudaFrees" << endl;
	//serial mutation
	/*
	 for (int i = 0; i < group_size; i ++){
	 		
	 		int mutate_or_not = rand() % 100 + 1;
	 		if(mutate_or_not < mutation_percentage){
	 			mutate(&h_children[i]);
		}
		//h_children[i].fitness = compute_fitness(children[i]);
	 }
	 */

	//cout << "535 after mutaion" << endl;
	//tour* c = new tour[group_size];
	//thrust::host_vector<tour_pair> h_children(group_size);
	//thrust::transform(d_children.begin(), d_children.end(), d_children.begin(), mutate_functor());
	
	//cout << "478 " << h_children[0].path[3] << endl;
	//cout << "478 " << h_children[0].path[40] << endl;
	//thrust::copy(d_children.begin(), d_children.end(), h_children.begin());


	

	//cudaMemcpy(d_children, h_children, sizeof(tour)*group_size, cudaMemcpyHostToDevice);
	
	compute_fitness<<<grid1, block1>>>(d_children, d_distance_matrix, num_cities);
	/*for (int i = 0; i < group_size; i++){
		//cout << "479 " << h_children[i].fitness << endl;
		//cout << "before h_children[i]" << endl;
		h_children[i].fitness = compute_fitness(h_children[i]);
		//cout << "child fitness : " << i << " " << h_children[i].fitness << endl;
	}*/
	
	cudaFree(d_random_arr);
	
	/*cudaMemcpy(h_children, d_children, sizeof(tour)*group_size, cudaMemcpyDeviceToHost);
	cudaFree(d_children);
		//cudaFree(d_children);
	cudaFree(d_parent_set1);
	cudaFree(d_parent_set2);
	cudaFree(d_index_of_parent_set1);
	cudaFree(d_index_of_parent_set2);
	cudaFree(d_cycles);
	
	//cout << "546 after fitness computation" << endl;
	
	return h_children;
	*/
}



/*
 * Tour mutate(Tour t ): mutates the tour and returns the
 * optimal tour of the mutated or original
 */

void mutate(tour* t ){
	int start;
	int end;
	do {
		start = rand() % (num_cities -1) + 1;
		end = rand() % (num_cities - 1) + 1;
	} while ((start >= end) or (start == 0 and end == (num_cities - 1)));
	
	while (start < end) {
		int temp = t->path[start];
		t->path[start] = t->path[end];
		t->path[end] = temp;
		start++;
		end--;
	}
}




/*
 * create_new_generation sorts the current population tour array and the children tour array, and then proceeds to replace
 * the group_size weakest tours in the population array with the children if the fitness of the child is higher
 */
__global__ void
	create_new_generation(tour* population, tour* children, int population_size, int group_size){
	//qsort_population(0, population_size -1, population);
	//qsort_population(0, group_size -1, children);
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int population_index = population_size - group_size;
	//int children_index = 0;
	/*
	 * this for loop only selects the group_size best between group_size lowest from previous generation and 
	 * the new children
	 */
	/*for(int i = 0; i < group_size && population_index < population_size; i ++){

		if(children[children_index].fitness < population[population_index].fitness){
			population[population_index] = children[children_index];
			population_index ++;
		}
		children_index++;
	}*/
	//replaces group_size worst solutions with children regardless of fitness
	population[population_index + id] = children[id];
	
	
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
	
	thrust::host_vector<tour> h_population(population_size);
	for(int i = 0; i < population_size; i++){
		h_population[i] = population[i];
	}
	thrust_population = h_population;
	thrust::sort(thrust_population.begin(), thrust_population.end());
	//move distance_matrix into gpu
	d_distance_matrix = 0;
	
	cudaMalloc((void**)&d_distance_matrix, sizeof(float) * num_cities*num_cities);
	
	int num_threads = group_size/2;
	int num_blocks = (num_threads + 31)/32;
	dim3 block(32, 1);
	dim3 grid(num_blocks, 1);
	h_children = new tour[group_size];
	h_parent_set1 = new tour[group_size/2];
	h_parent_set2 = new tour[group_size/2];
	h_index_of_parent_set1 = new tour[group_size/2];
	h_index_of_parent_set2 = new tour[group_size/2];
	h_cycles = new tour[group_size/2];
	
	
	for(int i = 0 ; i < group_size; i =i+2){
		h_parent_set1[i/2] = population[i];
		h_parent_set2[i/2] = population[i+1];
		//cout << "h_parent_set1 " << i/2 << " " <<  h_parent_set1[i/2].fitness << endl;
		//cout << "h_parent_set2 " << i/2 << " " <<h_parent_set2[i/2].fitness << endl;
	} 
	
	d_children = 0;
	d_parent_set1 = 0;
	d_parent_set2 = 0;
	d_index_of_parent_set1 = 0;
	d_index_of_parent_set2 = 0;
	d_cycles = 0;
	d_population = thrust::raw_pointer_cast(thrust_population.data());
	d_children = thrust::raw_pointer_cast(thrust_children.data());
	
	//cout << "497 before cudaMallocs" << endl;
	cudaMalloc((void**)&d_children, sizeof(tour)*group_size);
	cudaMalloc((void**)&d_parent_set1, sizeof(tour)*group_size/2);
	cudaMalloc((void**)&d_parent_set2, sizeof(tour)*group_size/2);
	cudaMalloc((void**)&d_index_of_parent_set1, sizeof(tour)*group_size/2);
	cudaMalloc((void**)&d_index_of_parent_set2, sizeof(tour)*group_size/2);
	cudaMalloc((void**)&d_cycles, sizeof(tour)*group_size/2);
	//cudaMalloc((void**)&d_population, sizeof(tour) * population_size);
	//cout << "504 before cudaMemcpys " << endl;
	
	cudaMemcpy(d_children, h_children, sizeof(tour)*group_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_parent_set1, h_parent_set1, sizeof(tour)*group_size/2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_parent_set2, h_parent_set2, sizeof(tour)*group_size/2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_index_of_parent_set1, h_index_of_parent_set1, sizeof(tour)*group_size/2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_index_of_parent_set2, h_index_of_parent_set2, sizeof(tour)*group_size/2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_cycles, h_cycles, sizeof(tour)*group_size/2, cudaMemcpyHostToDevice);
	cudaMemcpy(d_distance_matrix, distance_matrix, sizeof(float)*num_cities*num_cities, cudaMemcpyHostToDevice);
	//cudaMemcpy(d_population, population, sizeof(tour) * population_size);
	

	num_threads = group_size;
	num_blocks = (num_threads + 31)/32;
	dim3 block1(32, 1);
	dim3 grid1(num_blocks, 1);

	for(int j = 0; j < termination_step; j++){
		create_children();
		thrust::sort(thrust_children.begin(), thrust_children.end());
		create_new_generation<<<block1, grid1>>>(d_population, d_children, population_size, group_size);

	}
	
	/*for(int m = 0; m < population_size; m++){
		population[m] = h_population[m];
	}
	*/
	thrust::sort(thrust_population.begin(), thrust_population.end());
		
	cudaFree(d_children);
	cudaFree(d_parent_set1);
	cudaFree(d_parent_set2);
	cudaFree(d_index_of_parent_set1);
	cudaFree(d_index_of_parent_set2);
	cudaFree(d_cycles);
	cudaFree(d_distance_matrix);
	

	
	cudaMemcpy(population, d_population, sizeof(tour)*population_size, cudaMemcpyDeviceToHost);
	//cudaMemcpy(best_tour, d_population[0], sizeof(tour), cudaMemcpyDeviceToHost);
	
	cudaFree(d_population);
	
	
	cout << "Fitness: " << population[0].fitness << endl << "Tour \n";
	for (int i = 0; i < num_cities; i++){
		cout << population[0].path[i] << endl;
	}
	//print_best_tour();
}

/*
 * print_best_tour: prints out the best tour along with its fitness, run at the end of
 * run_genetic_algorithm()
 */
void print_best_tour(){
	qsort_population(0, population_size - 1, population);
	cout << "Best Tour Generated After " << termination_step << "generations \n";
	cout << "Fitness: " << population[0].fitness << endl << "Tour \n";
	for (int i = 0; i < num_cities; i++){
		cout << population[0].path[i] << endl;
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
	double start = timestamp();
// your codes here

	
	char* filename = argv[1];
	const char * num_cities_string = argv[2];
	num_cities = atoi(num_cities_string);
	population_size = atoi(argv[3]);
	greedy_selection_percentage = atoi(argv[4]);
	num_closer_way_points = atoi(argv[5]);
	group_size = atoi(argv[6]);
	mutation_percentage = atoi(argv[7]);
	termination_step = atoi(argv[8]);
	distance_matrix = new float[num_cities * num_cities];
	/*for loop no longe rneeded as changed distance_matrix to 1d array
	 * for (int i = 0; i< num_cities; i++) {
	 *		distance_matrix[i] = new float[num_cities];
}
*/
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
			distance_matrix[i * num_cities + j] = sqrt(pow(cities[i].x - cities[j].x,2) + pow(cities[i].y - cities[j].y,2));
		}
	}
	
	myfile.close();
	
	cout << "file parsing reached" << endl;
	run_genetic_algorithm();
	/* no longer needed as distance_matrix is a 1d array
	 *	for (int i = 0; i < num_cities; i++) {
	 *		delete[] distance_matrix[i];
}
*/
	delete[] distance_matrix;
	double end = timestamp();
	cout << "time elapsed" << (start - end) * 1000 << endl;
}
