#include "Algorithms.h"

GreedyAlgorithm::GreedyAlgorithm() : Algorithm("Greedy Algorithm"), mWeightSum(0u), mValueSum(0u), mContent()
{

}

void GreedyAlgorithm::init()
{

}

void GreedyAlgorithm::runImpl()
{
	std::vector<std::pair<unsigned int, double>> efficiencies(mItems.size(), std::pair<unsigned int, double>(0, -1));
	std::vector<bool> usedItems(mItems.size(), false);

	// calculate items' efficiency
	for (unsigned int i = 0; i < mItems.size(); ++i)
	{
		double efficiency = static_cast<double>(mItems[i]->getBenefit()) / mItems[i]->getWeight();
		efficiencies[i].first = i;
		efficiencies[i].second = efficiency;
	}

	// sort efficiencies from smallest to largest
	quickSort(efficiencies, 0, efficiencies.size());

	// fill backpack with greedy strategy

	for (unsigned int i = 0; i < efficiencies.size(); ++i)
	{
		unsigned int newWeightSum = mWeightSum + mItems[efficiencies[i].first]->getWeight();

		if (newWeightSum <= mCapacity)
		{
			mWeightSum = newWeightSum;
			mValueSum += mItems[efficiencies[i].first]->getBenefit();
			usedItems[efficiencies[i].first] = true;
		}
	}

	for (unsigned int i = 0; i < usedItems.size(); ++i)
	{
		if (usedItems[i])
		{
			mContent += mItems[i]->getName() + ";";
		}
	}
}

void GreedyAlgorithm::quickSort(std::vector<std::pair<unsigned int, double>>& efficiencies, const unsigned int p, const unsigned int q)
{
	unsigned int r = 0u; // pivot element

	if (p < q)
	{
		r = partition(efficiencies, p, q);
		quickSort(efficiencies, p, r);
		quickSort(efficiencies, r + 1, q);
	}
}

/*Pivot element is chosen as leftmost element initially.
Imagining all array elements to be cards, flip all card except for pivot element down.
Compare Pivot element with all other cards:
- if other card is <= pivot element, swap it with first open card and then close it
- if other card is > pivat element, leave it open
*/
unsigned int GreedyAlgorithm::partition(std::vector<std::pair<unsigned int, double>>& efficiencies, const unsigned int p, const unsigned int q)
{
	double efficiency = efficiencies[p].second;
	unsigned int i = p;

	for (unsigned int j = p + 1; j < q; j++)
	{
		double otherEfficiency = efficiencies[j].second;

		if (otherEfficiency <= efficiency)
		{
			i = i + 1;
			std::swap(efficiencies[i], efficiencies[j]);
		}
	}

	std::swap(efficiencies[i], efficiencies[p]);
	return i;
}

void GreedyAlgorithm::reportHeaderImpl(std::ostream& outputStream)
{
	outputStream << "capacity = " << mCapacity << std::endl;
};

void GreedyAlgorithm::reportResultsImpl(std::ostream& outputStream)
{
	std::string delimiter = ",";

	if (&outputStream == &std::cout)
	{
		delimiter = ", ";
	}

	outputStream << mValueSum << delimiter << mWeightSum << delimiter << mContent << std::endl;
}

void GreedyAlgorithm::reportFooterImpl(std::ostream& outputStream)
{
	outputStream << "Result is a greedy solution." << std::endl;
}

