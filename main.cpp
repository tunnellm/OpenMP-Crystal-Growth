#include "functions.cpp"
#include <vector>
#include <ostream>
#include <string>
#include <iostream>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <random>

using namespace std;

int g_max;

void sigHandler(int sigNum, siginfo_t *info, void *ucontext);

int main(int argc, char ** argv){
	
	struct sigaction userExit;
	userExit.sa_sigaction = sigHandler;
	sigemptyset(&userExit.sa_mask);
	sigaddset(&userExit.sa_mask, SIGINT);
	userExit.sa_flags = SA_SIGINFO;
	sigaction(SIGINT,&userExit,NULL);
		
	int dimensions = 0;
	int structureSize = 0;
	int initFloaters = 0;
	string options = "";
	string seedOption = "";
	int seedingDifficulty = 0;
	int seedExtension = 0;
	int numThreads;
	
	if (argc == 7) {
		
		/** Skip setup process for ease of automated benchmarking. */
		
		numThreads = atoi (argv[1]);
		dimensions = atoi (argv[2]);
		structureSize = atoi (argv[3]);
		initFloaters = atoi (argv[4]);
		g_max = atoi (argv[5]);
		seedingDifficulty = atoi (argv[6]);
		omp_set_num_threads (numThreads);
		
	} else {
		
		cout << "Number of threads: ";
		cin >> numThreads;
		cout << "Dimensions: ";
		cin >> dimensions;
		cout << "Structure Size: ";
		cin >> structureSize;
		cout << "Maximum Floaters: ";
		cin >> initFloaters;
		cout << "Maximum Iterations: ";
		cin >> g_max;
		cout << "Allow secondary seeds to form (y/n): ";
		cin >> options;
		omp_set_num_threads (numThreads);
			
		bool optionsChoice;
		optionsChoice = (options[0] == 'y' ? true : false);
		
		if (optionsChoice) {
			cout << "Set seeding difficulty." << endl;
			cout << "(int) % chance to create seed on collision: ";
			cin >> seedingDifficulty;
		}
		
		cout << "Set larger starting seed? (y/n): ";
		cin >> seedOption;
		
		if ((seedOption[0] == 'y' ? true : false)) {
			cout << "Set size of new seed: ";
			cin >> seedExtension;
			seedExtension = ((seedExtension > structureSize / 10 
			                   || seedExtension < 0) ? 1 : seedExtension);
		}
		
	}
	
	std::vector<int> lattice = spawnLattice(dimensions, structureSize);
	std::vector<Floater> floatersVector;
	std::vector<long long> seedLocations;
	std::vector<long long> falseSeedLocations;
	std::vector<long long> rules = generateRules(dimensions, structureSize);
	
	long long origin = 0;
	
	
	/** Convert origin to 1D coordinates */
	for (int i = 0; i < dimensions; i++)
		origin += ((long long) (structureSize / 2)) * pow(structureSize, i);
	
	if (origin <= 0) {
		cout << "Maximum dimension reached." << endl;
		cout << "Boost or other library required"
		     << "for larger dimensions than 'long long' can handle" << endl;
		cerr << "Use smaller dimension and/or structure size" << endl;
		exit(-1);
	}
	
	seedLocations.push_back(origin);
	lattice.at(origin) = 2;
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::bernoulli_distribution binomial(0.5);

	long long tempSeed = 0;
	for (int i = 0; i < seedExtension; ++i) {
		for (int j = 0; j < dimensions; j++)
			tempSeed += (long long) ((structureSize / 2) + 
			            (binomial(gen) > 0 ?  i : -i)) * pow(structureSize, j);
		seedLocations.push_back(tempSeed);
		lattice.at(tempSeed) = 2;
		tempSeed = 0;
	}
	
	
	
	for (int i = 0; i < initFloaters; i++){
		Floater floater;
		floatersVector.push_back(floater);
	}
	
	/** Initialize floaters and first iteration of lattice */
	/** CPP guarantees vectors are thread-safe if accessing different elements */
	
	#pragma omp parallel for
	for (int i = 0; i < floatersVector.size(); i++) {
		floatersVector.at(i).init(0, dimensions, structureSize);
		lattice.at(floatersVector.at(i).location) = 1;
	}

	
	bool endLoop = false;
	int currentRadius = 1 + seedExtension;
	int max_threads = omp_get_max_threads();
	int localOrigin = structureSize / 2;
	
	auto startTime = std::chrono::steady_clock::now();
	
	std::vector<int> tempRadi(omp_get_max_threads(), -1);
	
	
	for (int j = 1; j <= g_max; j++) {
		
		cout << j << "/" << g_max << ". Seeds: " << seedLocations.size() 
			 << ", False Seeds:" << falseSeedLocations.size() <<  endl;
		
		if (endLoop == true)
			break;
		
		
		#pragma omp parallel
		{	
			
			#pragma omp for schedule(guided)
			for (int i = 0; i < floatersVector.size(); i++)
				floatersVector.at(i).requestMove(lattice, rules, seedingDifficulty);
			
			int tempRadiusPlaceholder = -1;
		
		        // Directives are on one line in the code.
			// Wrapped here to be readable.
		
			/** The following compiler directive was modified from a SO post:
		        * https://stackoverflow.com/questions/18669296/
		        * c-openmp-parallel-for-loop-alternatives-to-stdvector */
		        
			#pragma omp declare reduction (merge : std::vector<long long>
			: omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
			
			/** Multi-threading */
			
			#pragma omp for reduction(merge: seedLocations)
			reduction(merge: falseSeedLocations)
			private(tempRadiusPlaceholder) schedule(guided)
			for (int i = 0; i < floatersVector.size(); i++)
				if (floatersVector.at(i).isSeed) {
					seedLocations.push_back(floatersVector.at(i).location);
					tempRadiusPlaceholder = validateRadius(dimensions, 
					       localOrigin, currentRadius, structureSize, 
					       inflateCoordinates(dimensions, structureSize,
					       floatersVector.at(i).location));
					
					tempRadi.at(omp_get_thread_num()) = 
							(tempRadiusPlaceholder >
							tempRadi.at(omp_get_thread_num()) 
							? tempRadiusPlaceholder :
							tempRadi.at(omp_get_thread_num()));
					 
					floatersVector.at(i).init(tempRadi.at(omp_get_thread_num()), dimensions, structureSize);
						
				} 
				else if (floatersVector.at(i).falseSeed) {
					falseSeedLocations.push_back(floatersVector.at(i).location);
					floatersVector.at(i).init(currentRadius,
					dimensions, structureSize);
				}
				else if (floatersVector.at(i).iterations >= 100 
				        || floatersVector.at(i).newLocation < 0)
					floatersVector.at(i).init(currentRadius,
					dimensions, structureSize);
				else 
					floatersVector.at(i).location =
					floatersVector.at(i).newLocation;
				
		
		
			#pragma omp single nowait	
			for (auto& it : tempRadi) {
				currentRadius = (it > currentRadius ? it : currentRadius);
				it = -1;
			}
				
		
			/** The following parallelized std::fill was found on SO:
				* https://stackoverflow.com/questions/42044956/
				* parallel-fill-stdvector-with-zero */
			
			
			#pragma omp private (tid, chunksize, begin, end)
			
			auto tid = omp_get_thread_num();
			auto chunksize = lattice.size() / omp_get_num_threads();
			auto begin = lattice.begin() + chunksize * tid;
			auto end = ((tid == omp_get_num_threads() -1) ?
			            lattice.end() : begin + chunksize);
			std::fill(begin, end, 0);
				
			#pragma omp barrier
			
			/** There should realistically be no race condition in this section.
				* there is certainly not one in which a 2 or 3 can be overwritten.
				*/
			
			#pragma omp for nowait
			for (int i = 0; i < floatersVector.size(); i++)
				if (floatersVector.at(i).location != -1)
					lattice.at(floatersVector.at(i).location) = 1;
				else
					endLoop = true; // Does not matter if there is 
					                // a race condition on 
					                // this variable.
			
			#pragma omp for nowait
			for (int i = 0; i < falseSeedLocations.size(); i++) 
				lattice.at(falseSeedLocations.at(i)) = 3;
			
			
			#pragma omp for nowait
			for (int i = 0; i < seedLocations.size(); i++)
				lattice.at(seedLocations.at(i)) = 2;
			
			
		}
		
	}
	
	auto endTime = std::chrono::steady_clock::now();
	
	std::chrono::duration<double> elapsed_seconds = endTime-startTime;
	
	cout << elapsed_seconds.count() << endl;
	
	ofstream file;
	file.open("./output/" + to_string(dimensions) + " " + to_string(structureSize) + " " 
		+ to_string(initFloaters) + " " + to_string(g_max) + " " + to_string(time(NULL)));
	
	for (auto& it : seedLocations) {
		string writeString = "";
		for (auto& innerIter : inflateCoordinates(dimensions, structureSize, it))
			writeString.insert(0, to_string(innerIter) + ",");
		writeString.pop_back();
		file << writeString << endl;
	}
	
	for (auto& it : falseSeedLocations) {
		string writeString = "";
		for (auto& innerIter : inflateCoordinates(dimensions, structureSize, it))
			writeString.insert(0, to_string(innerIter) + ",");
		writeString.pop_back();
		file << writeString << endl;
	}
	
	file.close();
	
	return 0;
}


void sigHandler (int sigNum, siginfo_t *info, void *ucontext) {
    if (g_max > 0)
		cout << endl << "Received an interrupt, saving results." << endl 
			 << "Ctrl^C again to exit immediately." << endl;
	else 
		exit(0);
	g_max = 0;
}