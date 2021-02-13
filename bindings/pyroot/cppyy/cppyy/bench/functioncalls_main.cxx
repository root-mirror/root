#include <iostream>
#include <iomanip>
#include <time.h>
#include <unistd.h>

#include "functioncalls.h"


static const int N = 1000000000; // 10^9, i.e. per call is in nanoseconds

int cpp_loop_offset() {
    int i = 0;
    for ( ; i < N; ++i)
        ;
    return i;
}

//- group: empty-free --------------------------------------------------------
void cpp_free_empty_call() {
    for (int i=0; i < N; ++i)
        empty_call();
}

//- group: empty-inst --------------------------------------------------------
void cpp_inst_empty_call() {
    EmptyCall e;
    for (int i=0; i < N; ++i)
        e.empty_call();
}


//- group: builtin-args-inst -------------------------------------------------
void cpp_inst_pass_int() {
    TakeAValue t;
    for (int i=0; i < N; ++i)
        t.pass_int(i);
}


//- group: do-work -----------------------------------------------------------
void cpp_inst_do_work() {
    DoWork d;
    for (int i=0; i < N; ++i)
        d.do_work(i);
}


//- group: do-overload -------------------------------------------------------
void cpp_inst_overload() {
    OverloadedCall o;
    for (int i=0; i < N; ++i)
        o.add_it((float)i);
}


//----------------------------------------------------------------------------
void run_bench(void (*cpp_bench)(), const char* label) {
    clock_t t1 = clock();
    cpp_loop_offset();
    clock_t t2 = clock();
    cpp_bench();
    clock_t t3 = clock();

    std::cout << label << ": " << std::setprecision(8)
              << ((t3-t2) - (t2-t1))/((double)CLOCKS_PER_SEC) << " nanoseconds" << std::endl;
}

int main() {
    run_bench(cpp_free_empty_call, "cpp_free_empty_call");
    run_bench(cpp_inst_empty_call, "cpp_inst_empty_call");
    run_bench(cpp_inst_pass_int,   "cpp_inst_pass_int");
    run_bench(cpp_inst_do_work,    "cpp_do_work");
    run_bench(cpp_inst_overload,   "cpp_inst_overload");
    return 0;
}
