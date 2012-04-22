#ifndef TSP_H
#define TSP_H

using std::vector;
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

struct coordinates {
float x;
float y;
};

int array_partition(int* input, int p, int r);
int array_quick_select(int* input, int p, int r, int k);
int vector_partition(vector<int> input, int p, int r);
int vector_quick_select(vector<int> input, int p, int r, int k);
void generate_initial_population();
void generate_closest_neighbors();
 void generate_tour(int* linear_cities, int index);
 void crossover(tour parent1, tour parent2, tour* children, int index);
 void sort_population();
 void qsort_population(int left, int right, tour* population) ;
 tour* create_children();
 void mutate(tour* t );
 tour* create_new_generation(tour* population, tour* children, int population_size, int group_size);
 void print_best_tour();
  void run();
  
#endif /* TSP_H_ */




