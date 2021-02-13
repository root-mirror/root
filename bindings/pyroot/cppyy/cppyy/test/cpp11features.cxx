#if __cplusplus >= 201103L

#include "cpp11features.h"


// for std::shared_ptr<> testing
int TestSharedPtr::s_counter = 0;

std::shared_ptr<TestSharedPtr> create_shared_ptr_instance() {
    return std::shared_ptr<TestSharedPtr>(new TestSharedPtr);
}

int TestSharedPtr::get_value() {
    return 17;
}

int DerivedTestSharedPtr::get_value() {
    return m_int + 76;
}

int pass_shared_ptr(std::shared_ptr<TestSharedPtr> p) {
    return p->get_value();
}


// for move ctors etc.
int TestMoving1::s_move_counter = 0;
int TestMoving1::s_instance_counter = 0;
int TestMoving2::s_move_counter = 0;
int TestMoving2::s_instance_counter = 0;

void implicit_converion_move(TestMoving2&&) {
    /* empty */
}


// for std::function testing
std::function<int(const FNTestStruct& t)> FNCreateTestStructFunc() { return [](const FNTestStruct& t) { return t.t; }; }
std::function<int(const FunctionNS::FNTestStruct& t)> FunctionNS::FNCreateTestStructFunc() { return [](const FNTestStruct& t) { return t.t; }; }

#endif // c++11 and later
