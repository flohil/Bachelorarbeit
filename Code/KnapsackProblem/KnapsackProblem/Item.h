#ifndef ITEM_H
#define ITEM_H

class Item
{
public:

	Item(std::string name, const unsigned int weight, const unsigned int benefit) : mName(name), mWeight(weight), mBenefit(benefit) {};

	inline std::string getName() { return mName; };
	inline const unsigned int getWeight() { return mWeight; };
	inline const int getBenefit() { return mBenefit; };

private:

	std::string mName;
	const unsigned int mWeight;
	const int mBenefit;
};

#endif /* ITEM_H */