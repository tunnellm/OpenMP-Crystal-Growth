#include "functions.hpp"
#include <vector>
#include <ostream>
#include <string>
#include <iostream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <climits>
#include <numeric>
#include <random>

using namespace std;

/** @Author: Marc Tunnell
* 	@Date: 9/30/21
* 	CIS677 High-performance computing
*	Project 2: Crystal simulation
*		
*		This file contains all of the necessary
*		functions and utilities to simulate
*		crystal growth via diffusion in an n-
*		dimensional lattice. Represented as abort
*		flattened 1-dimensional vector.
*	*/


/** 
*	Converts coordinates from 1-dimensional to n-dimension.
*
*	@param dimensions : Number of simulated dimensions
*	@param structureSize : Length of cube edge
*	@param flatCoordinate : 1-dimensional coordinate
*/

std::vector<int> inflateCoordinates (const int dimensions, const  int structureSize,
                                    const long long flatCoordinate ) {
	
	long long reduction = 0;
	long long tempFlat = flatCoordinate;
	std::vector<int> result;
	
	for (int i = dimensions - 1; i >= 0; i--){
		long long powNum = (long long) pow((double) structureSize, (double) i);
		reduction = tempFlat / powNum;
		tempFlat = tempFlat % powNum;
		result.push_back(reduction);
	}
	return result;
}

/** 
*	Converts from inflated coordinate to 1-dimension
*
*	@param coords : Vector of coordinates to be converted
*	@param dimensions : Number of simulated dimensions
*	@param structureSize : Length of cube edge
*/

long long flattenCoordinate (const std::vector<int>& coords, const int dimensions, 
                            const int structureSize) {
	long long flatCoord = 0;
	
	/** Must be in order of largest dimension to smallest. */

	for (int i = dimensions - 1; i >= 0; i--) {
		flatCoord += ((long long) coords.at(i) * pow(structureSize, i));
		if (coords.at(i) < 0) {
			return -1;
		}
	}
	return flatCoord;
}

/** 
*	Generates a list of relative movement locations given
*		the simulated dimension.
*
*	@param dimensions : Number of simulated dimensions
* 	@param structureSize : Length of cube edge
*/

std::vector<long long> generateRules(const int dimensions, const int structureSize){
	std::vector <long long> temp;
	std::vector <long long> output;
	
	output.push_back(0);
	temp.push_back(0);
	
	for (int i = 0; i < dimensions; i++) {
		for (auto& it : output){
			long long offset = (long long) pow(structureSize, i);
			temp.push_back(it - offset);
			temp.push_back(it + offset);
		}
		output = temp;
		std::vector <long long> temp; 
	}
	
	return output;
	
}

/** 
*	Checks to see if a new location is out of bounds relative
*		to the current location.
*
*	@param dimension : Number of simulated dimensions
* 	@param structureSize : Length of cube edge
* 	@param old_location : Current location as flattened coordinates
*	@param new_location : Target location 
*/

bool checkBoundaries (const int dimension, const int structureSize, 
                        const long long old_location, const long long new_location){
	
	int difference = 0;
	
	for (auto& it : inflateCoordinates(dimension, structureSize, old_location))
		difference += it;
	for (auto& it : inflateCoordinates(dimension, structureSize, new_location))
		difference -= it;

	return abs(difference) <= dimension;
}

/** 
*	Sets or resets a Floater object to the starting state.
*	
*	@param currentRadius : Length of main starting seed from origin
*	@param dimensionValue : Number of simulated dimensions
*	@param structureSizeValue : Length of cube edge
*/

void Floater::init(const int currentRadius, const int dimensionValue, 
                    const int structureSizeValue){
	structureSize = structureSizeValue;
	dimension = dimensionValue;
	iterations = 0;
	newLocation = -1;
	isSeed = false;
	falseSeed = false;
	location = spawn(currentRadius);
}

/** 
*	Attempts to spawn at a random location outside of the
* 		current radius from the origin.
*
*	@param radius : Length of main starting seed from origin
*/

long long Floater::spawn(const int radius){

	std::vector<int> coords;
	int localOrigin = structureSize / 2;
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distribution(1, (localOrigin - (radius + 1)));	
	std::bernoulli_distribution binomial(0.5);
	
	// Determines a spawn location outside of the current radius.
	// push_back(-1) if no spawn is possible.
	for (int i = 0; i < dimension; i++) {
		int tempVar = 0;
		if (binomial(gen) == 0) 
			if ((tempVar = (distribution(gen))) <= 0) 
				coords.push_back(-1);
			else 
				coords.push_back(tempVar);
		else 
			if ((tempVar = (localOrigin + distribution(gen))) >= structureSize) 
				coords.push_back(-1);
			else
				coords.push_back(tempVar);
	}
	return flattenCoordinate(coords, dimension, structureSize);
}

/** 
*	Returns the current location as a flattened coordinate.
*/

long long Floater::getLocation(){
	return location;
}

/** 
*	Checks whether the current location is a valid extension of the
*		seed radius.
*
*	@param dimension : Number of simulated dimensions
*	@param localOrigin : Half the length of cube edge
* 	@param currentRadius : Length of main starting seed from origin
* 	@param structureSize : Length of cube edge
* 	@param inflatedCoords : Current location as a vector of 
				inflated coordinates.
*/

int validateRadius (const int dimension, const int localOrigin, const int currentRadius, 
                        const int structureSize, const std::vector<int> inflatedCoords) {
	int newRadi = -1;
	int compare = currentRadius;
	
	for (auto& innerIter : inflatedCoords){
		newRadi = abs(localOrigin - innerIter);
		if (newRadi > compare)
			compare = newRadi;
	}
	return compare;
}

/** 
*	Calculates a new move and checks for collision with Floater
*		or seed objects.
*
*	@param inputLattice : 1-dimensional representation of the cube.
*	@param rules : Vector of relative movement rules.
*	@param chance : Percentage chance to stick to adjacent seeds.
*/

void Floater::requestMove(const std::vector<int>& inputLattice, 
                            const std::vector<long long>& rules, const int chance){
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> rulesDistribution(1, rules.size() - 1);
	std::uniform_int_distribution<> chanceDistribution(1, 100);
	
	// Looks for nearby seeds or Floater objects.
	
	newLocation = - 1;
	int current_location = getLocation();
	for (auto it = rules.begin() + 1; it != rules.end(); it++){
		if (current_location + *it >= 0 && 
			current_location + *it < inputLattice.size())
			switch (inputLattice.at(current_location + *it)){	
			case 3:
				if (checkBoundaries(dimension, structureSize, 
					current_location, current_location + *it)){
					falseSeed = true;
					return;
				} else
					break;
			case 2:
				if (checkBoundaries(dimension, structureSize, 
					current_location, current_location + *it)){
					isSeed = true;
					return;
				} else
					break;
			case 1:
				if (chanceDistribution(gen) <= chance) {
					if (checkBoundaries(dimension, structureSize, 
						current_location, current_location + *it)){
						falseSeed = true;
					} else
						break;
				} else
					break;
			default:
				break;
			}
	}
	
	// Checks if out of bounds in order of least expensive
	//	calculation to most.
	
	if ((newLocation = current_location + rules[rulesDistribution(gen)]) < 0)
		newLocation = -1;
	else if (newLocation > inputLattice.size())
		newLocation = -1;
	else if ((!checkBoundaries(dimension, structureSize, current_location, newLocation)))
		newLocation = -1;
	else
		iterations++;
	
}

/** 
*	Spawns a 1-dimensional representation of the n-dimension cube world
*
*	@param dimension : Number of simulated dimensions
*	@param structureSize : Length of cube edge
*/

std::vector<int> spawnLattice (const int dimension, const int structureSize){
	std::vector <int> lat(pow(structureSize, dimension), 0);
	return lat;
}