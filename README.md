# OpenMP-Crystal-Growth
OpenMP-Accelerated Simulation of Biological Crystal Growth via Diffusion-Limited Aggregation


Multi-threaded implementation of biological crystal growth through random movement contained within an n-dimensional cubic structure. In essence, we allocate enough elements for the n-dimensional cube, flatten to a 1-dimensional vector, and cast from a 1-dimensional coordinate system into the appropriate coordinate system for the n-dimensional cube. Directional movement is allowed on any adjacent square in a hypersphere around the object, with the relative rules for movement being generated at run-time based on user input.

Compile as: g++ main.cpp functions.cpp -o out.o

Run as:
./out.o
or
./out.o numThreads dimensions structureSize initialFloaters g_max seedingPercentChance
