/*
 * city.h
 *
 *  Created on: Apr 17, 2012
 *      Author: Kelvin
 */

#ifndef CITY_H
#define CITY_H

struct city
{
	int cityNum;
	int* closest_neighbors;

};

int* find_n_closest_neighbors(city* c, int n);


#endif /* CITY_H */
