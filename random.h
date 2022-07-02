using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>

std::uint32_t calc_seed(int ndata=1);

std::vector<int> calc_random_int(size_t ndata, unsigned int seedv=1, int minv=0, int maxv=100);

std::vector<int> calc_random_int2(size_t ndata, int seedv=1, int minv=0, int maxv=100);

std::vector<double> calc_random_double(size_t ndata, unsigned int seedv=1, double minv=0.0, double maxv=1.0);

