#include "Algorithms.h"

DynamicProgramming::DynamicProgramming() : Algorithm("Dynamic Programming"), mOptimalSolutions()
{

};

void DynamicProgramming::init()
{

}

// Solves smaller subproblems (less items, less capacity than original) of knapsack problem.
// By iteratively adding one new item and increasing the possible capacity until the original parameters are reached,
// the optimal solution for the knapsack problem can be found.
void DynamicProgramming::runImpl()
{
	unsigned int i = 0u; // iterator for items' indizes
	unsigned int d = 0u; // iterator for capacities, can be seen as remaining capacity in a knapsack

	// reserve memory
	mOptimalSolutions.resize(mItems.size() + 1);

	for (unsigned int i = 0; i < mOptimalSolutions.size(); ++i)
	{
		mOptimalSolutions[i].resize(mCapacity + 1);

		for (unsigned int j = 0; j < mOptimalSolutions[i].size(); ++j)
		{
			mOptimalSolutions[i][j].itemFlags.resize(mItems.size());
		}
	}

	// Build optimalSolutions in bottom up manner (starting with smallest optimal subproblem)
	for (i = 0; i <= mItems.size(); ++i) // try adding a new item to the previous optimal subset
	{
		/*
			At first, this part of the algorithm might seem a little unintuitive:

			There is no check if the weight of a new item, added to the sum of weights of all already packed items
			would exceeds the knapsack's total capacity c.	However, this check happens implicitly!

			Whenever an item is added to the current subset of the knapsack, 
			the capacity index at which the item is added is decreased by the items weight.

			So, if current capacity d = 6, and the added item has a weight w = 3, the item will be added at 
			capacities[d-w] = capacities[6-3] = capacities[3] instead of capacities[d] = capacities[6].

			Following this concept, checking whether an item's weight exceeds the current capacity d
			suffices to ensure that the weight of all packed items never surpasses the knapsacks total capacity c,
			because the weight of the other packed items is already considered through decreasing the index where 
			the items are added by the item's weight.
		*/

		for (d = 0; d <= mCapacity; ++d) // compute possible change of the optimal solutions again for all capacities
		{
			if (i == 0 || d == 0) // no items in knapsack or no capacity -> no items can be packed
			{
				mOptimalSolutions[i][d].benefit = 0; // ...hence set optimal solution value to 0
			}
			else if (d >= mItems[i - 1]->getWeight()) // new item's weight does not exceed remaining capacity d
			{
				unsigned int newTotalBenefit = 
					mItems[i - 1]->getBenefit() + mOptimalSolutions[i - 1][d - mItems[i - 1]->getWeight()].benefit;
				unsigned int previousTotalBenefit = mOptimalSolutions[i - 1][d].benefit;

				/* 
					If packing a new item into the knapsack would result in a higher optimal solution value,
					pack the item and increase the solution value of the previous optimal subset by the item's value.
					Otherwise, do not pack the new item and leave the previous optimal solution value unchanged.
				*/
				if (newTotalBenefit > previousTotalBenefit)
				{
					mOptimalSolutions[i][d].benefit = newTotalBenefit; // total benefit increases
					mOptimalSolutions[i][d].itemFlags = mOptimalSolutions[i - 1][d - mItems[i - 1]->getWeight()].itemFlags; // in addition to previous items
					mOptimalSolutions[i][d].itemFlags[i - 1] = true; // ... new item is added to optimal subset of items
				}
				else
				{
					mOptimalSolutions[i][d] = mOptimalSolutions[i - 1][d]; // item not added -> benefit stays unchanged
				}
			}
			else // new item is too heavy to fit into knapsack with remaining capacity d
			{
				mOptimalSolutions[i][d] = mOptimalSolutions[i - 1][d]; // ... hence knapsack stays unchanged
			}
		}
	}
}

void DynamicProgramming::reportHeaderImpl(std::ostream& outputStream)
{
	outputStream << "capacity = " << mCapacity << std::endl;
};

void DynamicProgramming::reportResultsImpl(std::ostream& outputStream)
{
	int optimalSolution = mOptimalSolutions[mItems.size()][mCapacity].benefit;
	unsigned int optimalWeight = 0;
	std::string optimalContent = "";
	std::string delimiter = ",";

	if (&outputStream == &std::cout)
	{
		delimiter = ", ";
	}

	for (unsigned int i = 0; i < mOptimalSolutions[mItems.size()][mCapacity].itemFlags.size(); ++i)
	{
		if (mOptimalSolutions[mItems.size()][mCapacity].itemFlags[i])
		{
			optimalWeight += mItems[i]->getWeight();
			optimalContent += mItems[i]->getName() + ";";
		}
	}
	
	outputStream << optimalSolution << delimiter << optimalWeight << delimiter << optimalContent << std::endl;
}

void DynamicProgramming::reportFooterImpl(std::ostream& outputStream)
{
	outputStream << "Result is an optimal solution." << std::endl;
}