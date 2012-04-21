 
/*
 *  tsp-parser.cpp
 *  
 *
 *  Created by Preetham Soliman on 4/20/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
using namespace std;


struct city {
    float x;
    float y;
};

int main(int argc, char* argv[]) {
	
	while(argc!= 3){
		cout << "not enough arguments, force close program";
		return 1;
		
	}
	cout << "asdf1\n";

    char* filename = argv[1];
    const char * num_cities_string = argv[2];
    int num_cities = atoi(num_cities_string);
    
    float distance_matrix[num_cities][num_cities];
    city cities[num_cities];
    
    cout << "asdf2" + num_cities;
    
    char line[256];
    ifstream myfile;
	myfile.open(filename);
    
	
	
	cout << "asdf3\n";
	char* comp = "NODE_COORD_SECTION";
    while (1) {
   	 myfile.getline(line, 256);
   	 if (!strcmp(comp, line)) {
		  cout << line << "\n"; 
   		 break;
		
   	 }
   	 	 cout << "t\n " ;
    }
    
    for (int i = 0; i < num_cities; i++) {
   	 myfile.getline(line, 256);
   	 char * pch;

	
	
   	pch = strtok(line, " ");
   	 pch = strtok(NULL, " ");
   	 cities[i].x = atoi(strtok(NULL, " "));
   	 cities[i].y = atoi(strtok(NULL, " "));
   	 
    }
    
    for (int i = 0; i < num_cities; i++) {
   	 for (int j = 0; j < num_cities; j++) {
   		 distance_matrix[i][j] = sqrt(pow(cities[i].x - cities[j].x,2) + pow(cities[i].y - cities[j].y,2));
   	 }
    }

    myfile.close();
    
    
    
    
}
