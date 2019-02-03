//#ifndef PRINT
//#define PRINT
//#endif /* PRINT */

//#ifndef PROGRESS
//#define PROGRESS
//#endif /* PROGRESS */

#include <iostream>
#include <cstdlib>
#include <chrono>
#include <thread>
#include "Algorithms.h"
#include "SettingsParser.h"
#include "Item.h"

struct Population
{
public:

	// preallocate memory for population and set all gene bits to false
	Population(const unsigned int popSize, const unsigned int chromosomeSize, const unsigned int threadPoolSize) : 
		chromosomes(popSize, std::vector<bool>(chromosomeSize, false)),
		fitnessValues(popSize),
		sameFitnessValues(),
		sameFitnessValuesParts(threadPoolSize),
		newChromosomes(popSize),
		crossoverSelectedsParts(threadPoolSize),
		crossoverSelecteds(),
		crossoverPairs(),
		bestChromosomeCurrent(chromosomeSize),
		bestFitnessCurrent(0),
		worstFitnessCurrent(0),
		totalFitnessSumCurrent(0),
		bestChromosomeEver(chromosomeSize),
		bestFitnessEver(0),
		bestFitnessIterationEver(0),
		bestDurationEver(),
		avgFitness(0),
		mPopSize(popSize),
		mChromosomeSize(chromosomeSize),
		abortGeneration(0) {};

	unsigned int getPopulationSize() const { return mPopSize; };
	unsigned int getChromosomeSize() const { return mChromosomeSize; };

	std::vector<std::vector<bool>> chromosomes;
	std::vector<int> fitnessValues;
	std::map<int, unsigned int> sameFitnessValues;
	std::vector<std::map<int, unsigned int>> sameFitnessValuesParts;
	std::vector<std::vector<bool>> newChromosomes;
	std::vector<std::vector<unsigned int>> crossoverSelectedsParts;
	std::vector<unsigned int> crossoverSelecteds;
	std::vector<std::pair<unsigned int, unsigned int>> crossoverPairs;

	std::vector<bool> bestChromosomeCurrent;
	int bestFitnessCurrent;
	int worstFitnessCurrent;
	int totalFitnessSumCurrent;
	std::vector<bool> bestChromosomeEver;
	int bestFitnessEver;
	unsigned int bestFitnessIterationEver;
	std::chrono::milliseconds bestDurationEver;
	double avgFitness;
	unsigned int abortGeneration;

private:
	
	const unsigned int mPopSize;
	const unsigned int mChromosomeSize;
};

GeneticAlgorithm::GeneticAlgorithm() : Algorithm("Genetic Algorithm"),
mPopulation(nullptr), mPc(0.), mPm(0.), mSameValueMargin(0u), mSameFitnessValueAbortPercentage(0.),
mIterations(0u), mPopSize(0u), mThreadPoolSize(0u) {};

void GeneticAlgorithm::init()
{
	// set parameters from parametersFile
	mSettingsParser.get("iterations", mIterations);
	mSettingsParser.get("popSize", mPopSize);
	mSettingsParser.get("pc", mPc);
	mSettingsParser.get("pm", mPm);
	mSettingsParser.get("sameFitnessValueAbortPercentage", mSameFitnessValueAbortPercentage);
	mSettingsParser.get("threadPoolSize", mThreadPoolSize);

	// calculate nbr of identic chromosomes in population, at which algorithm shall abort
	mSameValueMargin = static_cast<unsigned int>(mSameFitnessValueAbortPercentage * mPopSize);

	// create threadpool and completing check mechanism
	mThreadPool = new ctpl::thread_pool(mThreadPoolSize);

	// create random number generators
	std::uniform_int_distribution<int> mSeeder_dist(1, RAND_MAX); // range for the seed	
	while (mRandNbrGens.size() < mThreadPoolSize)
	{
		// mSeeder_dist(seeder) delivers different new seeds within range of seeder_dist to initialize PRNGs
		mRandNbrGens.push_back(std::default_random_engine(mSeeder_dist(*mSeeder)));
	}

	createInitialPopulation(mPopSize, mItems.size());
}

void GeneticAlgorithm::runImpl()
{	
	reportGenerationHeader(mOutputFile);
	#ifdef PROGRESS
		reportGenerationHeader(std::cout);
	#endif

	unsigned int generationsCtr = 0;

	for (unsigned int i = 0; i < mIterations; ++i)
	{
		++generationsCtr;

		calcFitnessValues(i);

		reportGeneration(mOutputFile, i);
		#ifdef PROGRESS
				reportGeneration(std::cout, i);
		#endif

		if (checkOutstandingAccordance())
		{
			break;
		}

		performSelection();
		performCrossover();
		performMutation();
	}

	mPopulation->abortGeneration = generationsCtr;

	mOutputFile << std::endl;
	#ifdef PROGRESS
		std::cout << std::endl;
	#endif
}

// add definition, which chromosomes a thread should perform operations on
// if such a definition for the specified nbr of Elements already exists, do nothing
void GeneticAlgorithm::addTaskBunchOffsets(const unsigned int nbrElements)
{
	if (mTaskBunchOffsets.find(nbrElements) == mTaskBunchOffsets.end())
	{
		std::vector<std::pair<unsigned int, unsigned int>> offsets(0);
		unsigned int baseNbr = nbrElements / mThreadPool->size();
		unsigned int remainder = nbrElements % mThreadPool->size();
		unsigned int offset = 0u;

		while (offsets.size() < static_cast<unsigned int>(mThreadPool->size()))
		{
			unsigned int span = baseNbr - 1;
			std::pair<unsigned int, unsigned int> offsetPair(0u, 0u);

			// if popSize is not divisable by threadPoolSize, remainder of tasks also needs to be distributed
			if (remainder > 0)
			{
				++span;
				--remainder;
			}

			offsetPair.first = offset;
			offset += span; // last offset element (inclusive)
			offsetPair.second = offset;
			++offset; // new first offset element 

			offsets.push_back(offsetPair);
		}

		mTaskBunchOffsets.insert(std::pair<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>>(nbrElements, offsets));
	}
}

// creates a random initial population by "coin flip" method
void GeneticAlgorithm::createInitialPopulation(const unsigned int popSize, const unsigned int chromosomeSize)
{
	mPopulation = new Population(popSize, chromosomeSize, mThreadPool->size());

	addTaskBunchOffsets(popSize); // register population Size for multithreaded offsets

	for (unsigned int x = 0; x < mTaskBunchOffsets.at(popSize).size(); ++x)
	{
		increaseThreadCtr();

		// fill chromosomes with random values in different threads
		mThreadPool->push([this, x](int id)
		{
			for (unsigned int i = mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).first; i <= mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).second; ++i)
			{
				for (unsigned int j = 0; j < mPopulation->getChromosomeSize(); j++)
				{
					mPopulation->chromosomes[i][j] = randomBool(id);
				}
			}

			decreaseThreadCtr();
		});
	}

	waitForThreads();

	#ifdef PRINT
	std::cout << std::endl << "Initial population (size = " << mPopulation->getPopulationSize() << "):" << std::endl;
	for (unsigned int i = 0; i < mPopulation->getPopulationSize(); ++i)
	{
		std::cout << i << ": ";
		for (bool bit : mPopulation->chromosomes[i])
		{
			std::cout << bit;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif
}

// calculate fitness of chromosomes, if weight > capacity: fitness = 0
void GeneticAlgorithm::calcFitnessValues(const unsigned int iteration)
{
	mPopulation->totalFitnessSumCurrent = 0;

	for (unsigned int x = 0; x < mTaskBunchOffsets.at(mPopulation->getPopulationSize()).size(); ++x)
	{
		increaseThreadCtr();

		mThreadPool->push([this, x](int id)
		{
			std::map<int, unsigned int> sameFitnessValuesPart;

			for (unsigned int i = mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).first; i <= mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).second; ++i)
			{
				int fitnessSum = 0u;
				unsigned int weightSum = 0u;

				for (unsigned int j = 0; j < mPopulation->getChromosomeSize(); j++)
				{
					// if item has been added, increase fitnessSum and weightSum
					if (mPopulation->chromosomes[i][j])
					{
						weightSum += mItems[j]->getWeight();

						// The more often an additional item's weight exceeds the maximal capacity,
						// the more negative the resulting fitness Value shall be.
						// A too heavy combination of items can never result in a positive fitness value,
						// since it does not yield a valid solution.
						if (weightSum > mCapacity)
						{
							if (fitnessSum >= 0) // item combination just turned too heavy to carry
							{
								fitnessSum = -mItems[j]->getBenefit(); 
							}
							else
							{
								fitnessSum -= mItems[j]->getBenefit();
							}
						} 
						else
						{
							fitnessSum += mItems[j]->getBenefit();
						}
					}
				}

				// bundle nbr of elements with same fitness value (except for fitnessSum of 0)
				if (sameFitnessValuesPart.find(fitnessSum) != sameFitnessValuesPart.end())
				{
					++sameFitnessValuesPart[fitnessSum];
				}
				else
				{
					sameFitnessValuesPart.insert(std::pair<unsigned int, unsigned int>(fitnessSum, 1u));
				}

				mPopulation->fitnessValues[i] = fitnessSum;
			}

			mPopulation->sameFitnessValuesParts[x] = sameFitnessValuesPart;

			decreaseThreadCtr();
		});
	}

	waitForThreads();

	mPopulation->bestFitnessCurrent = mPopulation->fitnessValues[0];
	mPopulation->worstFitnessCurrent = mPopulation->fitnessValues[0];

	// save best chromosome of this generation - even if first generation has only negative benefit chromosomes
	for (unsigned int i = 0; i < mPopulation->getPopulationSize(); ++i)
	{
		mPopulation->totalFitnessSumCurrent += mPopulation->fitnessValues[i];

		if (mPopulation->fitnessValues[i] > mPopulation->bestFitnessCurrent)
		{
			mPopulation->bestFitnessCurrent = mPopulation->fitnessValues[i];
			mPopulation->bestChromosomeCurrent = mPopulation->chromosomes[i];
		}
		if (mPopulation->fitnessValues[i] < mPopulation->worstFitnessCurrent)
		{
			mPopulation->worstFitnessCurrent = mPopulation->fitnessValues[i];
		}
	}

	// save best chromosome of all generations
	if (mPopulation->bestFitnessCurrent > mPopulation->bestFitnessEver || iteration == 0)
	{
		std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
		mPopulation->bestDurationEver = std::chrono::duration_cast<std::chrono::milliseconds>(now - mBegin);
		mPopulation->bestFitnessEver = mPopulation->bestFitnessCurrent;
		mPopulation->bestFitnessIterationEver = iteration;
		mPopulation->bestChromosomeEver = mPopulation->bestChromosomeCurrent;
	}

	// calc avg fitness of generation
	mPopulation->avgFitness = static_cast<double>(mPopulation->totalFitnessSumCurrent) / static_cast<double>(mPopulation->getPopulationSize());

	#ifdef PRINT
	std::cout << std::endl << "Fitness values: " << std::endl;
	for (unsigned int i = 0; i < mPopulation->getPopulationSize(); ++i)
	{
		std::cout << i << ": " << mPopulation->fitnessValues[i] << std::endl;
	}
	std::cout << std::endl;
	#endif
}

// checks if 90% of population have reached the same fitness value (very likely that maximum is found)
bool GeneticAlgorithm::checkOutstandingAccordance()
{
	for (std::map<int, unsigned int> sameFitnessValuesPart : mPopulation->sameFitnessValuesParts)
	{
		typedef std::map<int, unsigned int>::iterator it_type;
		for (it_type iterator = sameFitnessValuesPart.begin(); iterator != sameFitnessValuesPart.end(); iterator++) 
		{
			if (mPopulation->sameFitnessValues.find(iterator->first) != mPopulation->sameFitnessValues.end())
			{
				mPopulation->sameFitnessValues[iterator->first] += iterator->second;

				// do not abort if all fitness values are 0 (if weightSums are too high this will happen)
				if (mPopulation->sameFitnessValues[iterator->first] >= mSameValueMargin && iterator->first > 0u)
				{
					mPopulation->sameFitnessValues.clear();
					return true;
				}
			}
			else
			{
				mPopulation->sameFitnessValues.insert(std::pair<unsigned int, unsigned int>(iterator->first, iterator->second));
			}
		}
	}

	mPopulation->sameFitnessValues.clear();
	return false;
}

// select those chromosomes which will be part of a new generation with higher chances for fitter chromosomes
void GeneticAlgorithm::performSelection()
{
	std::vector<unsigned int> rouletteWheel(0);
	unsigned int projectedWheelOffset = 0u;

	// allow negativ fitnesses small chances of reproduction too
	for (int fitnessValue : mPopulation->fitnessValues)
	{
		projectedWheelOffset += (fitnessValue - mPopulation->worstFitnessCurrent); // results in 0 additional offset for worst Fitness value
		rouletteWheel.push_back(projectedWheelOffset);
	}

	unsigned int span = rouletteWheel[rouletteWheel.size() - 1];

	for (unsigned int x = 0; x < mTaskBunchOffsets.at(mPopulation->getPopulationSize()).size(); ++x)
	{
		increaseThreadCtr();

		mThreadPool->push([this, rouletteWheel, span, x](int id)
		{
			for (unsigned int i = mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).first;
				i <= mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).second; ++i)
			{
				unsigned int randUInt = randomUInt(span, id);

				for (unsigned int j = 0; j < rouletteWheel.size(); ++j)
				{
					if (randUInt < rouletteWheel[j])
					{
						mPopulation->newChromosomes[i] = mPopulation->chromosomes[j];
						break;
					}
				}
			}

			decreaseThreadCtr();
		});
	}

	waitForThreads();

	mPopulation->chromosomes = mPopulation->newChromosomes;

	#ifdef PRINT
	std::cout << std::endl << "New population:" << std::endl;
	for (unsigned int i = 0; i < mPopulation->getPopulationSize(); ++i)
	{
		std::cout << i << ": ";
		for (bool bit : mPopulation->chromosomes[i])
		{
			std::cout << bit;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif
}

// performs crossover between pairs of chromosomes selected with a probability of pc
void GeneticAlgorithm::performCrossover()
{
	mPopulation->crossoverPairs.clear();

	// select chromosomes for crossover
	for (unsigned int x = 0; x < mTaskBunchOffsets.at(mPopulation->getPopulationSize()).size(); ++x)
	{
		increaseThreadCtr();

		mThreadPool->push([this, x](int id)
		{
			std::vector<unsigned int> crossoverSelectedsPart(0);

			for (unsigned int i = mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).first; i <= mTaskBunchOffsets.at(mPopulation->getPopulationSize()).at(x).second; ++i)
			{
				if (randomProbability(id) < mPc)
				{
					crossoverSelectedsPart.push_back(i);
				}
			}

			//std::cout << "crossoverSelectedsPart.size(): " << crossoverSelectedsPart.size() << std::endl;

			mPopulation->crossoverSelectedsParts[x] = crossoverSelectedsPart;

			decreaseThreadCtr();
		});
	}

	waitForThreads();

	// add all selected chromosomes into one single vector
	for (std::vector<unsigned int> crossoverSelectedPart : mPopulation->crossoverSelectedsParts)
	{
		// std::cout << "crossoverSelectedPart.size(): " << crossoverSelectedPart.size() << std::endl;
		mPopulation->crossoverSelecteds.insert(mPopulation->crossoverSelecteds.end(), crossoverSelectedPart.begin(), crossoverSelectedPart.end());
		// std::cout << "cSelecteds.size(): " << mPopulation->crossoverSelecteds.size() << std::endl;
	}

	unsigned int reservedChromosomeIdx = mPopulation->crossoverSelecteds[randomUInt(mPopulation->crossoverSelecteds.size() - 1)];

	// remove chromosomes from pool of those selected for crossover by forming random crossover pairs
	while (mPopulation->crossoverSelecteds.size() > 1) 
	{
		unsigned int idx1 = randomUInt(mPopulation->crossoverSelecteds.size() - 1);
		unsigned int chromosomeIdx1 = mPopulation->crossoverSelecteds[idx1];
		mPopulation->crossoverSelecteds.erase(mPopulation->crossoverSelecteds.begin() + idx1);
		unsigned int idx2 = randomUInt(mPopulation->crossoverSelecteds.size() - 1);
		unsigned int chromosomeIdx2 = mPopulation->crossoverSelecteds[idx2];
		mPopulation->crossoverSelecteds.erase(mPopulation->crossoverSelecteds.begin() + idx2);

		mPopulation->crossoverPairs.push_back(std::pair<unsigned int, unsigned int>(chromosomeIdx1, chromosomeIdx2));
	}

	// if the number of chromosomes for crossover is uneven, add one additional chromosome from selected chromosomes
	if (mPopulation->crossoverSelecteds.size() == 1)
	{
		unsigned int chromosomeIdx1 = mPopulation->crossoverSelecteds[0];
		mPopulation->crossoverPairs.push_back(std::pair<unsigned int, unsigned int>(chromosomeIdx1, reservedChromosomeIdx));
	}

	addTaskBunchOffsets(mPopulation->crossoverPairs.size()); // register population Size for multithreaded offsets

	#ifdef PRINT
	mPopulation->newChromosomes = mPopulation->chromosomes;
	#endif

	// perform actual crossover
	for (unsigned int x = 0; x < mTaskBunchOffsets.at(mPopulation->crossoverPairs.size()).size(); ++x)
	{
		mThreadPool->push([this, x](int id)
		{
			for (unsigned int i = mTaskBunchOffsets.at(mPopulation->crossoverPairs.size()).at(x).first; i <= mTaskBunchOffsets.at(mPopulation->crossoverPairs.size()).at(x).second; ++i)
			{
				unsigned int randPosition = randomUInt(mPopulation->chromosomes[0].size() - 2, id);

				std::vector<bool> chromosome1 = mPopulation->chromosomes[mPopulation->crossoverPairs[i].first];
				std::vector<bool> chromosome2 = mPopulation->chromosomes[mPopulation->crossoverPairs[i].second];

				for (unsigned int j = randPosition; j < chromosome1.size(); ++j)
				{
					mPopulation->chromosomes[mPopulation->crossoverPairs[i].first][j] = chromosome2[j];
					mPopulation->chromosomes[mPopulation->crossoverPairs[i].second][j] = chromosome1[j];
				}
			}
		});
	}

	waitForThreads();

	#ifdef PRINT
	std::cout << std::endl << "Crossedover population:" << std::endl;
	for (unsigned int i = 0; i < mPopulation->getPopulationSize(); ++i)
	{
		std::cout << i << ": ";
		for (bool bit : mPopulation->newChromosomes[i])
		{
			std::cout << bit;
		}
		std::cout << " -> ";
		for (bool bit : mPopulation->chromosomes[i])
		{
			std::cout << bit;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif
}

// change bit values for all bits of all chromosomes with probability of pm
void GeneticAlgorithm::performMutation()
{
	#ifdef PRINT
	mPopulation->newChromosomes = mPopulation->chromosomes;
	#endif

	for (unsigned int x = 0; x < mTaskBunchOffsets.at(mPopulation->chromosomes.size()).size(); ++x)
	{
		increaseThreadCtr();

		mThreadPool->push([this, x](int id)
		{
			for (unsigned int i = mTaskBunchOffsets.at(mPopulation->chromosomes.size()).at(x).first; i <= mTaskBunchOffsets.at(mPopulation->chromosomes.size()).at(x).second; ++i)
			{
				for (unsigned int j = 0; j < mPopulation->chromosomes[0].size(); ++j)
				{
					if (randomProbability(id) < mPm)
					{
						mPopulation->chromosomes[i][j] = !mPopulation->chromosomes[i][j];
					}
				}
			}

			decreaseThreadCtr();
		});
	}

	waitForThreads();

	#ifdef PRINT
	std::cout << std::endl << "mutation results: " << std::endl;
	for (unsigned int i = 0; i < mPopulation->chromosomes.size(); ++i)
	{
		bool wasMutated = false;
		for (unsigned int j = 0; j < mPopulation->chromosomes[0].size(); ++j)
		{
			if (mPopulation->chromosomes[i][j] != mPopulation->newChromosomes[i][j])
			{
				wasMutated = true;
			}
		}

		if (wasMutated)
		{
			std::cout << i << ": ";
			for (bool bit : mPopulation->newChromosomes[i])
			{
				std::cout << bit;
			}
			std::cout << " -> ";
			for (bool bit : mPopulation->chromosomes[i])
			{
				std::cout << bit;
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	#endif
}

// writes header for outputFile
void GeneticAlgorithm::reportGenerationHeader(std::ostream& outputStream)
{

	if (&outputStream == &std::cout)
	{
		outputStream << "iteration" << ": "
			<< "aValue" << ", "
			<< "bValue" << ", "
			<< "bContent"
			<< std::endl;
	}
	else
	{
		outputStream << "iteration" << ","
			<< "aValue" << ","
			<< "cValue" << ","
			<< "cWeight" << ","
			<< "cContent" << ","
			<< "bValue" << ","
			<< "bWeight" << ","
			<< "bContent"
			<< std::endl;
	}
}

// write fittest chromosome to outputFile
void GeneticAlgorithm::reportGeneration(std::ostream& outputStream, const unsigned int iteration)
{
	std::string knapsackContentCurrent = "";
	std::string knapsackContentEver = "";
	unsigned int weightSumCurrent = 0u;
	unsigned int weightSumEver = 0u;

	for (unsigned int i = 0; i < mPopulation->getChromosomeSize(); ++i)
	{
		if (mPopulation->bestChromosomeCurrent[i])
		{
			weightSumCurrent += mItems[i]->getWeight();
			knapsackContentCurrent += mItems[i]->getName();
			knapsackContentCurrent += ";";
		}
		if (mPopulation->bestChromosomeEver[i])
		{
			weightSumEver += mItems[i]->getWeight();
			knapsackContentEver += mItems[i]->getName();
			knapsackContentEver += ";";
		}
	}

	if (&outputStream == &std::cout)
	{
		outputStream << (iteration + 1) << ": "
			<< mPopulation->avgFitness << ", "
			<< mPopulation->bestFitnessEver << ", "
			<< knapsackContentEver.c_str()
			<< std::endl;
	}
	else
	{
		outputStream << (iteration + 1) << ","
			<< mPopulation->avgFitness << ","
			<< mPopulation->bestFitnessCurrent << ","
			<< weightSumCurrent << ","
			<< knapsackContentCurrent.c_str() << ","
			<< mPopulation->bestFitnessEver << ","
			<< weightSumEver << ","
			<< knapsackContentEver.c_str()
			<< std::endl;
	}
}

void GeneticAlgorithm::reportHeaderImpl(std::ostream& outputStream)
{
	outputStream
		<< "capacity = " << mCapacity
		<< "; maxIterations = " << mIterations
		<< "; popSize = " << mPopulation->getPopulationSize()
		<< "; pc = " << mPc
		<< "; pm = " << mPm
		<< "; sameFitnessValueAbortPercentage = " << mSameFitnessValueAbortPercentage
		<< "; threadPoolSize = " << mThreadPoolSize
		<< std::endl;
};

void GeneticAlgorithm::reportResultsImpl(std::ostream& outputStream)
{
	unsigned int bestWeight = 0;
	std::string bestContent = "";
	std::string delimiter = ",";

	if (&outputStream == &std::cout)
	{
		delimiter = ", ";
	}

	for (unsigned int i = 0; i < mPopulation->bestChromosomeEver.size(); ++i)
	{
		if (mPopulation->bestChromosomeEver[i])
		{
			bestWeight += mItems[i]->getWeight();
			bestContent += mItems[i]->getName() + ";";
		}
	}

	outputStream << mPopulation->bestFitnessEver << delimiter << bestWeight << delimiter << bestContent << std::endl;
}

void GeneticAlgorithm::reportFooterImpl(std::ostream& outputStream)
{
	if (&outputStream == &std::cout)
	{
		outputStream << "bestExecTime: " << mPopulation->bestDurationEver.count() << "ms" << std::endl;
		outputStream << "bestFitnessIteration: " << (mPopulation->bestFitnessIterationEver + 1) << std::endl;
		outputStream << "abortIteration: " << mPopulation->abortGeneration << std::endl;
	}
	else
	{
		outputStream << "bestExecTime:," << mPopulation->bestDurationEver.count() << ",ms" << std::endl;
		outputStream << "bestFitnessIteration:," << (mPopulation->bestFitnessIterationEver + 1) << std::endl;
		outputStream << "abortIteration:," << mPopulation->abortGeneration << std::endl;
	}

	outputStream << "Result is an approximation." << std::endl;
}

GeneticAlgorithm::~GeneticAlgorithm()
{
	delete mThreadPool;
	delete mPopulation;
};