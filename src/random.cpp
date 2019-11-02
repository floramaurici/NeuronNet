#include "random.h"

RandomNumbers::RandomNumbers(unsigned long int s) : seed(s) {
	if (seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    rng = std::mt19937(seed);
	}
void RandomNumbers::uniform_double(std::vector<double>& tab, double lower, double upper) {
	for (auto i = tab.begin(); i != tab.end(); i++) *i = uniform_double(lower, upper);
	}

double RandomNumbers::uniform_double(double lower, double upper) {
	std::uniform_real_distribution<> unif(lower, upper);
	return unif (rng);
	}
	
void RandomNumbers::normal(std::vector<double>& tab, double mean, double sd) {
	for (auto i = tab.begin(); i != tab.end(); i++) *i = normal(mean, sd);
	}

double RandomNumbers::normal(double mean, double sd) {
	std::normal_distribution<> norm(mean, sd);
	return norm(rng);
	}

void RandomNumbers::poisson(std::vector<int>& tab, double mean) {
	for (auto i = tab.begin(); i != tab.end(); i++) *i = poisson(mean);
	}

int RandomNumbers::poisson(double mean) {
	 std::poisson_distribution<int> poiss(mean);
	 return poiss(rng);
	}
