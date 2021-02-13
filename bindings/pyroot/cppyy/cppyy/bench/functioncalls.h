#ifndef CPPYY_FUNCTIONCALLS_H
#define CPPYY_FUNCTIONCALLS_H

//- group: empty-free --------------------------------------------------------
void empty_call();

//- group: empty-inst --------------------------------------------------------
class EmptyCall {
public:
    void empty_call();
};


//- group: builtin-args-free -------------------------------------------------
struct Value { int m_int; };

void take_an_int(int);
void take_a_double(double);
void take_a_struct(Value);

//- group: builtin-args-free -------------------------------------------------
class TakeAValue {
public:
    void take_an_int(int);
    void take_a_double(double);
    void take_a_struct(Value);

    int pass_int(int);
};


//- group: do-work-free ------------------------------------------------------
double do_work(double);

//- group: do-work-inst ------------------------------------------------------
class DoWork {
public:
    double do_work(double);
};


//- group: overload-inst -----------------------------------------------------
class OverloadedCall {
public:
    double add_it(int a, int b);
    double add_it(short a);
    double add_it(long a);
    double add_it(int a, int b, int c);
    double add_it(double a);
    double add_it(float a);
    double add_it(int a);
};

#endif // !CPPYY_FUNCTIONCALLS_H
