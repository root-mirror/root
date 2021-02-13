#include "runvector.h"

std::vector<int> global_vector;
double global_effect = 0.;

static const int N = 100000000; // 10^8

namespace {
    struct Initializer {
        Initializer() {
            global_vector.reserve(N);
            for (int i=0; i < N; ++i) global_vector.push_back(i);
        }
    } _init_global_vector;
}
