#include "functioncalls.h"

#include <exception>
#include <math.h>


//- group: empty -------------------------------------------------------------
void empty_call() {
/* empty, to measure pure call overhead */
}

void EmptyCall::empty_call() {
/* empty, to measure pure call overhead */
}


//- group: builtin-args-free -------------------------------------------------
void take_an_int(int /* unused */) {
/* empty, to measure pure call overhead */
}

void take_a_double(double /* unused */) {
/* empty, to measure pure call overhead */
}

void take_a_struct(Value /* unused */) {
/* empty, to measure pure call overhead */
}

//- group: builtin-args-inst -------------------------------------------------
void TakeAValue::take_an_int(int /* unused */) {
/* empty, to measure pure call overhead */
}

void TakeAValue::take_a_double(double /* unused */) {
/* empty, to measure pure call overhead */
}

void TakeAValue::take_a_struct(Value /* unused */) {
/* empty, to measure pure call overhead */
}

//- group: builtin-args-pass -------------------------------------------------
int TakeAValue::pass_int(int a) {
    return a + 42;
}


//- group: do-work -----------------------------------------------------------
double do_work(double arg) {
    return atan(arg);
}

double DoWork::do_work(double arg) {
    return atan(arg);
}


//- group: overload-inst -----------------------------------------------------
double OverloadedCall::add_it(int a, int b)        { std::terminate(); }
double OverloadedCall::add_it(short a)             { std::terminate(); }
double OverloadedCall::add_it(long a)              { std::terminate(); }
double OverloadedCall::add_it(int a, int b, int c) { std::terminate(); }
double OverloadedCall::add_it(double a)            { return 3.1415 + a; }
double OverloadedCall::add_it(float a)             { return 3.1415 + a; }
double OverloadedCall::add_it(int a)               { std::terminate(); }
