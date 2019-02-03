#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <cmath>
#include <ctime>
#include "SettingsParser.h"
#include "Item.h"
#include "csv.h"
#include "ctpl_stl.h"
#include <mutex>

// holds a set of items, represented by flags at their indizes (true -> included in knapsack), and their combined benefit
struct BenefitSet
{
	BenefitSet() : benefit(0u), itemFlags(0) {};
	BenefitSet(unsigned int itemsCount) : benefit(0u), itemFlags(itemsCount) {};

	int benefit;
	std::vector<bool> itemFlags;

	/*
	Default assign (and copy) behaviour:
	Each member is copied using their copy constructor or assignment operator (depending on operation).
	This is applied recursively for members and their members.

	-> itemIdxs will get copied, since std::vector copies its members on assignment.
	*/
};

class Algorithm
{
public:

	Algorithm(std::string algorithmName);
	virtual ~Algorithm() {};

	// Runs the algorithm on items, using the parameters in parametersFile, results will be written to outputFile.
	// Fails if files cannot be found or input contains no items.
	bool run(const char* inputFilepath, const char* parametersFilepath, const char* outputFilepath);

protected:

	// returns a random unsigned integer between 0 and inclusive MAX
	unsigned int randomUInt(const unsigned int max);

	// returns a random float between 0 and 1
	double randomProbability();

	// returns a random float between 0 and max
	double randomDouble(const float max);

	// returns a random bool
	bool randomBool();

// multithread random functions since rand() is not threadsafe

	// returns a random unsigned integer between 0 and inclusive MAX, for use with multiple threads
	unsigned int randomUInt(const unsigned int max, int threadId);

	// returns a random float between 0 and 1, for use with multiple threads
	double randomProbability(int threadId);

	// returns a random double between 0 and max
	double randomDouble(const double max, int threadId);

	// returns a random bool, for use with multiple threads
	bool randomBool(int threadId);

	// wait for all threads to finish their tasks
	void waitForThreads();

	// increase nbr of threads currently solving tasks
	void increaseThreadCtr();

	// decrease nbr of threads currently solving tasks
	void decreaseThreadCtr();

	typedef std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> taskBunchOffsetsType;

	std::chrono::steady_clock::time_point mBegin;
	std::mutex threadCtrMutex;
	unsigned int threadCtr;
	SettingsParser mSettingsParser;
	std::ifstream mInputFile;
	std::ofstream mOutputFile;
	std::vector<std::default_random_engine> mRandNbrGens;
	std::default_random_engine* mSeeder;
	ctpl::thread_pool* mThreadPool;
	taskBunchOffsetsType mTaskBunchOffsets; // defines on which parts of a dataset a thread should perform operations
	std::vector<Item*> mItems;
	unsigned int mCapacity;

private:

	// run implementation specific part of the algorithm
	virtual void runImpl() = 0;
	virtual void init() = 0; // call before runImpl!!!

	// report general algorithm information to outputStream
	void reportHeader(std::ostream& outputStream);
	void reportResults(std::ostream& outputStream);
	void reportFooter(std::ostream& outputStream, std::chrono::milliseconds timeDiff);

	// report implementation specific information to outputStream
	virtual void reportHeaderImpl(std::ostream& outputStream) = 0;
	virtual void reportResultsImpl(std::ostream& outputStream) = 0;
	virtual void reportFooterImpl(std::ostream& outputStream) = 0;

	std::string mAlgorithmName;
};

struct Population; //dataObject holding population information
class GeneticAlgorithm : public Algorithm
{
public:

	GeneticAlgorithm();
	virtual ~GeneticAlgorithm();

private:

	virtual void runImpl();
	virtual void init();
	void addTaskBunchOffsets(const unsigned int nbrElements);
	void createInitialPopulation(const unsigned int popSize, const unsigned int chromosomeSize);
	void calcFitnessValues(const unsigned int iteration);
	bool checkOutstandingAccordance();
	void performSelection();
	void performCrossover();
	void performMutation();

	// report implementation specific information to outputStream
	virtual void reportHeaderImpl(std::ostream& outputStream);
	virtual void reportResultsImpl(std::ostream& outputStream);
	virtual void reportFooterImpl(std::ostream& outputStream);
	void reportGenerationHeader(std::ostream& outputStream);
	void reportGeneration(std::ostream& outputStream, const unsigned int iteration);

	Population* mPopulation;
	double mPc;
	double mPm;
	unsigned int mSameValueMargin;
	double mSameFitnessValueAbortPercentage;
	unsigned int mIterations;
	unsigned int mPopSize;
	unsigned int mThreadPoolSize;
};

class DynamicProgramming : public Algorithm
{
public:

	DynamicProgramming();
	virtual ~DynamicProgramming(){};

private:

	virtual void runImpl();
	virtual void init();

	virtual void reportHeaderImpl(std::ostream& outputStream);
	virtual void reportResultsImpl(std::ostream& outputStream);
	virtual void reportFooterImpl(std::ostream& outputStream);

	// an array of weighedItemIdxs to be used for all possible subcapacities d of total capacity c
	typedef std::vector<BenefitSet> weighedBenefitSets;

	// an array of weighedItemIdxsArray subproblems - one array element for each newly added item
	typedef std::vector<weighedBenefitSets> optimalSolutionsType;

	// holds all optimal solution values and item sets as indizes
	optimalSolutionsType mOptimalSolutions;
};

class GreedyAlgorithm : public Algorithm
{
public:

	GreedyAlgorithm();
	virtual ~GreedyAlgorithm(){};

private:

	virtual void runImpl();
	virtual void init();

	// sort efficiencies in n log n time
	void quickSort(std::vector<std::pair<unsigned int, double>>& efficiencies, const unsigned int p, const unsigned int q);
	
	// find partition element for next recursion of quickSort
	unsigned int partition(std::vector<std::pair<unsigned int, double>>& efficiencies, const unsigned int p, const unsigned int q);

	virtual void reportHeaderImpl(std::ostream& outputStream);
	virtual void reportResultsImpl(std::ostream& outputStream);
	virtual void reportFooterImpl(std::ostream& outputStream);

	unsigned int mWeightSum;
	int mValueSum;
	std::string mContent;
};

#endif /* ALGORITHMS_H */