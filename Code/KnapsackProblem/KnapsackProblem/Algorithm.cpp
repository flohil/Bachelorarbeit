#include "Algorithms.h"

#ifndef CONSOLE_REPORT
#define CONSOLE_REPORT
#endif /* CONSOLE_REPORT */

Algorithm::Algorithm(std::string algorithmName) : mBegin(), threadCtrMutex(), threadCtr(0u), 
mSettingsParser(), mInputFile(), mOutputFile(), mRandNbrGens(), mSeeder(nullptr), mThreadPool(nullptr),
mTaskBunchOffsets(), mItems(), mCapacity(0), mAlgorithmName(algorithmName)
{
	// seed random number generator with current time
	srand(static_cast <unsigned> (time(0)));
};

bool Algorithm::run(const char* inputFilepath, const char* parametersFilepath, const char* outputFilepath)
{
	// try to open inputFile for parsing
	FILE *inputFile;
	if (inputFile = fopen(inputFilepath, "r")) {

		// parse inputFile and store items
		io::CSVReader<3> in(inputFilepath, inputFile);
		in.read_header(io::ignore_extra_column, "name", "weight", "benefit");

		if (!in.has_column("name"))
		{
			std::cerr << "Column 'name' is missing in " << inputFilepath << std::endl;
		}
		else if (!in.has_column("weight"))
		{
			std::cerr << "Column 'weight' is missing in " << inputFilepath << std::endl;
		}
		else if (!in.has_column("benefit"))
		{
			std::cerr << "Column 'benefit' is missing in " << inputFilepath << std::endl;
		}

		std::string name;
		unsigned int weight;
		int value;

		while (in.read_row(name, weight, value))
		{
			mItems.push_back(new Item(name, weight, value));
		}

		fclose(inputFile);
	}
	else
	{
		std::cerr << "Could not read from " << inputFilepath << std::endl;
		return false;
	}

	// check if there are items at all
	if (mItems.size() < 1)
	{
		std::cerr << "There must be at least one item" << std::endl;
		return false;
	}

	// try to load parameters
	if (!mSettingsParser.loadFromFile(parametersFilepath))
	{
		std::cerr << "Could not read parameters from " << parametersFilepath << std::endl;
		for (Item* item : mItems)
		{
			delete item;
		}
		return false;
	}

	// try to open outputFile for writing
	mOutputFile.open(outputFilepath);
	if (!mOutputFile.is_open())
	{
		std::cerr << "Could not write to " << outputFilepath << std::endl;
		for (Item* item : mItems)
		{
			delete item;
		}
		return false;
	}

	std::random_device rd; // Seed with a real random value, if available
	mSeeder = new std::default_random_engine(rd()); // use one Pseudo Random Number Generator to seed other PRNGs for each thread
	std::uniform_int_distribution<int> mSeeder_dist(1, RAND_MAX); // range for the seed

	// make sure there is at least one random number generator available
	mRandNbrGens.push_back(std::default_random_engine(mSeeder_dist(*mSeeder)));

	// set parameters from parametersFile
	mSettingsParser.get("capacity", mCapacity);

	// read implementation specific parameters and initialize member variables
	init();

	// report header
	reportHeader(mOutputFile);

	#ifdef CONSOLE_REPORT
		reportHeader(std::cout);
	#endif

	// perform specific algorithm and measure ellapsed time
	mBegin = std::chrono::steady_clock::now();
	runImpl();
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	// report results and footer
	reportResults(mOutputFile);
	reportFooter(mOutputFile, std::chrono::duration_cast<std::chrono::milliseconds>(end - mBegin));

	#ifdef CONSOLE_REPORT
		reportResults(std::cout);
		reportFooter(std::cout, std::chrono::duration_cast<std::chrono::milliseconds>(end - mBegin));
	#endif
	
	// cleanup
	mOutputFile.close();
	delete mSeeder;
	for (Item* item : mItems)
	{
		delete item;
	}

	return true;
}

unsigned int Algorithm::randomUInt(const unsigned int max)
{
	return randomUInt(max, 0);
};

double Algorithm::randomProbability()
{
	return randomProbability(0);
};

double Algorithm::randomDouble(const float max)
{
	return randomDouble(max, 0);
}

bool Algorithm::randomBool()
{
	return randomBool(0);
};

unsigned int Algorithm::randomUInt(const unsigned int max, int threadId)
{
	std::uniform_int_distribution<unsigned int> uniform_dist(0, max);
	return uniform_dist(mRandNbrGens[threadId]); // uses different PRNGs for each thread
};

double Algorithm::randomProbability(int threadId)
{
	return static_cast<float>(randomUInt(RAND_MAX, threadId)) / static_cast<float>(RAND_MAX);
};

double Algorithm::randomDouble(const double max, int threadId)
{
	return static_cast<double>(randomUInt(RAND_MAX, threadId)) / static_cast<double>(RAND_MAX)* max;
}

bool Algorithm::randomBool(int threadId)
{
	if (randomProbability(threadId) < 0.5f)
	{
		return true;
	}
	else
	{
		return false;
	}
};

void Algorithm::reportHeader(std::ostream& outputStream)
{
	outputStream << "Executing " << mAlgorithmName << " with parameters: " << std::endl;
	reportHeaderImpl(outputStream);
	outputStream << std::endl;
};

void Algorithm::reportResults(std::ostream& outputStream)
{
	if (&outputStream == &std::cout)
	{
		outputStream << "rValue, rWeight, rContent" << std::endl;
	}
	else
	{
		outputStream << "rValue,rWeight,rContent" << std::endl;
	}
	reportResultsImpl(outputStream);
	outputStream << std::endl;
};

void Algorithm::reportFooter(std::ostream& outputStream, std::chrono::milliseconds timeDiff)
{
	if (&outputStream == &std::cout)
	{
		outputStream << "execTime: " << timeDiff.count() << "ms" << std::endl;
	}
	else
	{
		outputStream << "execTime:," << timeDiff.count() << ",ms" << std::endl;
	}
	reportFooterImpl(outputStream);
};

// wait for all threads to finish current tasks
void Algorithm::waitForThreads()
{
	bool finished = false;
	while (!finished)
	{
		threadCtrMutex.lock();
		if (threadCtr == 0)
		{
			finished = true;
		}
		threadCtrMutex.unlock();
	}
}

void Algorithm::increaseThreadCtr()
{
	threadCtrMutex.lock();
	std::thread::id this_id = std::this_thread::get_id();
	++threadCtr;
	threadCtrMutex.unlock();
}

void Algorithm::decreaseThreadCtr()
{
	threadCtrMutex.lock();
	std::thread::id this_id = std::this_thread::get_id();
	--threadCtr;
	threadCtrMutex.unlock();
}