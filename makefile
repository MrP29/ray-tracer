CC=g++
CFLAGS=-std=c++11

raytrace: raytrace.cpp
	$(CC) $(CFLAGS) -o raytrace raytrace.cpp
