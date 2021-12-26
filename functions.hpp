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

using namespace std;



/** 
*	Converts coordinates from 1-dimensional to n-dimension.
*
*	@param dimensions : Number of simulated dimensions
*	@param structureSize : Length of cube edge
*	@param flatCoordinate : 1-dimensional coordinate
*/

std::vector<int> inflateCoordinates (const int dimensions, const int structureSize, 
					const long long flatCoordinate );



/** 
*	Converts from inflated coordinate to 1-dimension
*
*	@param coords : Vector of coordinates to be converted
*	@param dimensions : Number of simulated dimensions
*	@param structureSize : Length of cube edge
*/

long long flattenCoordinate (const std::vector<int>& coords, const int dimensions, 
								const int structureSize);



/** 
*	Generates a list of relative movement locations given
*		the simulated dimension.
*
*	@param dimensions : Number of simulated dimensions
* 	@param structureSize : Length of cube edge
*/

std::vector<long long> generateRules(const int dimensions, const int structureSize);



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
			const long long old_location, const long long new_location);



/** 
*	Spawns a 1-dimensional representation of the n-dimension cube world
*
*	@param dimension : Number of simulated dimensions
*	@param structureSize : Length of cube edge
*/

std::vector<int> spawnLattice (const int dimension, const int structureSize);



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
			const int structureSize, const std::vector<int> inflatedCoords);



/** 
*	Contains all of the relevant information and information
*		for traversal and storage of the floating objects.
*
*	@param location : Current location as flattened coordinates
*	@param newLocation : Location of requested move 
						as flattened coordinates
*	@param isSeed : Boolean flag determining whether this object
							is connected to the starter seed.
*	@param falseSeed : Boolean flag determining whether this object
						is stuck from encountering another object.
*	@param iterations : Number of movement requests since initialization
*/

class Floater {
	
	public:
		long long location;
		long long newLocation;
		bool isSeed;
		bool falseSeed;
		int iterations;
		
		/** 
		*	Calculates a new move and checks for collision with Floater
		*		or seed objects.
		*
		*	@param inputLattice : 1-dimensional representation of the cube.
		*	@param rules : Vector of relative movement rules.
		*	@param chance : Percentage chance to stick to adjacent seeds.
		*/
		
		void requestMove(const std::vector<int>& inputLattice, 
				const std::vector<long long>& rules, const int chance);
		
		/** 
		*	Sets or resets a Floater object to 	the starting state.
		*	
		*	@param currentRadius : Length of main starting seed from origin
		*	@param dimensionValue : Number of simulated dimensions
		*	@param structureSizeValue : Length of cube edge
		*/
		
		void init(const int currentRadius, const int dimensionValue, 
					const int structureSizeValue);
	private:
		int structureSize;
		int dimension;
		
		/** 
		*	Returns the current location as a flattened coordinate.
		*/
		
		long long getLocation();

		/** 
		*	Attempts to spawn at a random location outside of the
		* 		current radius from the origin.
		*
		*	@param radius : Length of main starting seed from origin
		*/
		
		long long spawn(const int radius);
		std::vector<bool> calculateRandomDirection();
		
};