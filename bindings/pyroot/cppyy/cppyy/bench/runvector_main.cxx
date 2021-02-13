#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <unistd.h>

#include "runvector.h"


static const int N = 100000000; // 10^8, i.e. per call is in 10 nanoseconds

int cpp_loop_offset() {
    int i = 0;
    for ( ; i < N; ++i)
        ;
    return i;
}


//- group: stl-vector --------------------------------------------------------
void cpp_runvector() {
    for (auto i: global_vector)
        global_effect += i;
}


//----------------------------------------------------------------------------
void run_bench(void (*cpp_bench)(), const char* label) {
    clock_t t1 = clock();
    cpp_loop_offset();
    clock_t t2 = clock();
    cpp_bench();
    clock_t t3 = clock();

    std::cout << label << ": " << std::setprecision(8)
              << 10.*((t3-t2) - (t2-t1))/((double)CLOCKS_PER_SEC) << " nanoseconds" << std::endl;
}

int main() {
    run_bench(cpp_runvector, "cpp_runvector");
    return 0;
}
