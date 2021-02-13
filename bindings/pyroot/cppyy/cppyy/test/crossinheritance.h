#ifndef CPPYY_TEST_CROSSINHERITANCE_H
#define CPPYY_TEST_CROSSINHERITANCE_H

#include <string>


//===========================================================================
namespace CrossInheritance {

class Base1 {                // for overridden method checking
public:
    Base1() : m_int(42) {}
    Base1(int i) : m_int(i) {}
    virtual ~Base1();

    virtual int get_value() { return m_int; }
    static int call_get_value(Base1* b);

    virtual int sum_value(int i) { return m_int + i; }
    static int call_sum_value(Base1* b, int);

    virtual int sum_all(int i) { return m_int + i; }
    virtual int sum_all(int i, int j) { return m_int + i + j; }
    static int call_sum_all(Base1* b, int);
    static int call_sum_all(Base1* b, int, int);

    virtual int pass_value1(int a) { return a; }
    virtual int pass_value2(int& a) { return a; }
    virtual int pass_value3(const int& a) { return a; }
    virtual int pass_value4(const Base1& b) { return b.m_int; }
    virtual int pass_value5(Base1& b) { return b.m_int; }
    static int sum_pass_value(Base1* b);

public:
    int m_int;
};

class IBase2 {
public:
    IBase2() {}
    virtual ~IBase2() {}
    virtual int get_value() = 0;
    static int call_get_value(IBase2* b);
};

class IBase3 : IBase2 {
public:
    IBase3(int);
    int m_int;
};

class CBase2 : public IBase2 {
public:
    int get_value();
};

class IBase4 {
public:
    IBase4() {}
    virtual ~IBase4() {}
    virtual int get_value() const = 0;      // <- const, as opposed to IBase2
    static int call_get_value(IBase4* b);
};

class CBase4 : public IBase4 {
public:
    int get_value() const;
};

template<typename T>
class TBase1 {
public:
    virtual ~TBase1() {}
    virtual int get_value() {
        return 42;
    }
};

class TDerived1 : public TBase1<int> {
public:
    int get_value();
};

using TBase1_I = TBase1<int>;

class CountableBase {
public:
    CountableBase();
    CountableBase(const CountableBase&);
    CountableBase& operator=(const CountableBase&);
    virtual ~CountableBase();

    virtual int call();

    static int s_count;
};

} // namespace CrossInheritance

#endif // !CPPYY_TEST_CROSSINHERITANCE_H
