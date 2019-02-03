#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "string.h"
#include "Algorithms.h"
#include <iostream>

int main(int argc, const char* argv[])
{
	// check for memory leaks
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);

	Algorithm* algorithm = nullptr;

// Parse arguments

	if (argc != 5)
	{
		std::cerr << "Usage: <input.csv> <output.csv> <GA|DP|GR> <parameters.txt>" << std::endl;
		return -1;
	}

	const char* inputFilepath = argv[1];
	const char* outputFilepath = argv[2];
	const char* algType = argv[3];
	const char* parametersFilepath = argv[4];

	if (strncmp(algType, "GA", 2) == 0)
	{
		algorithm = new GeneticAlgorithm();
	}
	else if (strncmp(algType, "DP", 2) == 0)
	{
		algorithm = new DynamicProgramming();
	}
	else if (strncmp(algType, "GR", 2) == 0)
	{
		algorithm = new GreedyAlgorithm();
	}
	else
	{
		std::cerr << "Usage: <input.csv> <output.csv> <GA|DP|GR> <parameters.txt>" << std::endl;
		return -1;
	}

// Perform optimization algorithm

	algorithm->run(inputFilepath, parametersFilepath, outputFilepath);

	delete algorithm;

#ifdef _DEBUG 
	std::cin.get();
#endif
}