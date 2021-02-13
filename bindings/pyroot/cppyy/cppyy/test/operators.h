#ifndef CPPYY_TEST_OPERATORS_H
#define CPPYY_TEST_OPERATORS_H

class number  {
public:
    number() { m_int = 0; }
    number(int i) { m_int = i; }

    number operator+(const number& n) const { return number(m_int + n.m_int); }
    number operator+(int n) const { return number(m_int + n); }
    number operator-(const number& n) const { return number(m_int - n.m_int); }
    number operator-(int n) const { return number(m_int - n); }
    number operator*(const number& n) const { return number(m_int * n.m_int); }
    number operator*(int n) const { return number(m_int * n); }
    number operator/(const number& n) const { return number(m_int / n.m_int); }
    number operator/(int n) const { return number(m_int / n); }
    number operator%(const number& n) const { return number(m_int % n.m_int); }
    number operator%(int n) const { return number(m_int % n); }

    number& operator+=(const number& n) { m_int += n.m_int; return *this; }
    number& operator-=(const number& n) { m_int -= n.m_int; return *this; }
    number& operator*=(const number& n) { m_int *= n.m_int; return *this; }
    number& operator/=(const number& n) { m_int /= n.m_int; return *this; }
    number& operator%=(const number& n) { m_int %= n.m_int; return *this; }

    number operator-() { return number( -m_int ); }

    bool operator<(const number& n) const { return m_int < n.m_int; }
    bool operator>(const number& n) const { return m_int > n.m_int; }
    bool operator<=(const number& n) const { return m_int <= n.m_int; }
    bool operator>=(const number& n) const { return m_int >= n.m_int; }
    bool operator!=(const number& n) const { return m_int != n.m_int; }
    bool operator==(const number& n) const { return m_int == n.m_int; }

    operator bool() { return m_int != 0; }

    number operator&(const number& n) const { return number(m_int & n.m_int); }
    number operator|(const number& n) const { return number(m_int | n.m_int); }
    number operator^(const number& n) const { return number(m_int ^ n.m_int); }

    number& operator&=(const number& n) { m_int &= n.m_int; return *this; }
    number& operator|=(const number& n) { m_int |= n.m_int; return *this; }
    number& operator^=(const number& n) { m_int ^= n.m_int; return *this; }

    number operator<<(int i) const { return number(m_int << i); }
    number operator>>(int i) const { return number(m_int >> i); }

private:
    int m_int;
};


//----------------------------------------------------------------------------
struct operator_char_star {       // for testing user-defined implicit casts
    operator_char_star() : m_str((char*)"operator_char_star") {}
    operator char*() { return m_str; }
    char* m_str;
};

struct operator_const_char_star {
    operator_const_char_star() : m_str("operator_const_char_star" ) {}
    operator const char*() { return m_str; }
    const char* m_str;
};

struct operator_int {
    operator int() { return m_int; }
    int m_int;
};

struct operator_long {
    operator long() { return m_long; }
    long m_long;
};

struct operator_double {
    operator double() { return m_double; }
    double m_double;
};

struct operator_short {
    operator short() { return m_short; }
    unsigned short m_short;
};

struct operator_unsigned_int {
    operator unsigned int() { return m_uint; }
    unsigned int m_uint;
};

struct operator_unsigned_long {
    operator unsigned long() { return m_ulong; }
    unsigned long m_ulong;
};

struct operator_float {
    operator float() { return m_float; }
    float m_float;
};


//----------------------------------------------------------------------------
class v_opeq_base {
public:
    v_opeq_base(int val);
    virtual ~v_opeq_base();

    virtual bool operator==(const v_opeq_base& other);

protected:
    int m_val;
};

class v_opeq_derived : public v_opeq_base {
public:
    v_opeq_derived(int val);
    virtual ~v_opeq_derived();

    virtual bool operator==(const v_opeq_derived& other);
};


//----------------------------------------------------------------------------
class YAMatrix1 {        // YetAnotherMatrix class for indexing tests
public:
    YAMatrix1() : m_val(42) {}

    int& operator() (int i, int j);
    const int& operator() (int i, int j) const;

    int m_val;
};

class YAMatrix2 {
public:
    YAMatrix2() : m_val(42) {}

    int& operator[] (int i);
    const int& operator[] (int i) const;

    int m_val;
};

class YAMatrix3 {
public:
    YAMatrix3() : m_val(42) {}

    int& operator() (int i, int j);
    const int& operator() (int i, int j) const;

    int& operator[] (int i);
    const int& operator[] (int i) const;

    int m_val;
};

class YAMatrix4 {
public:
    YAMatrix4() : m_val(42) {}

// opposite order of method declarations from YAMatrix3
    int& operator[] (int i);
    const int& operator[] (int i) const;

    int& operator() (int i, int j);
    const int& operator() (int i, int j) const;

    int m_val;
};

class YAMatrix5 {
public:
    YAMatrix5() : m_val(42) {}

    int& operator[] (int i);
    int operator[] (int i) const;

    int& operator() (int i, int j);
    int operator() (int i, int j) const;

    int m_val;
};

class YAMatrix6 {
public:
    YAMatrix6() : m_val(42) {}

    int& operator[] (int i);
    int operator[] (int i) const;

    int& operator() (int i, int j);

    int m_val;
};

class YAMatrix7 {
public:
    YAMatrix7() : m_val(42) {}

    int& operator[] (int i);
    int& operator() (int i, int j);

    int m_val;
};


//----------------------------------------------------------------------------
// templated operator
class TOIClass {
public:
    template<typename T>
    bool operator<(T t) { return true; }
};


//----------------------------------------------------------------------------
// __radd__/__rmul__, (non-)associative
class AssocADD {
public:
    AssocADD(double d) : fx(d) {}
    bool operator==(double d) { return d == fx; }

    friend AssocADD operator+(int n, AssocADD& a);
    friend AssocADD operator+(AssocADD& a, int n);

private:
    double fx;
};

AssocADD operator+(int n, AssocADD& a);
AssocADD operator+(AssocADD& a, int n);

class NonAssocRADD {
public:
    NonAssocRADD(double d) : fx(d) {}
    bool operator==(double d) { return d == fx; }

    friend NonAssocRADD operator+(int n, NonAssocRADD& a);

private:
    double fx;
};

NonAssocRADD operator+(int n, NonAssocRADD& a);

class AssocMUL {
public:
    AssocMUL(double d) : fx(d) {}
    bool operator==(double d) { return d == fx; }

    friend AssocMUL operator*(int n, AssocMUL& m);
    friend AssocMUL operator*(AssocMUL& m, int n);

private:
    double fx;
};

AssocMUL operator*(int n, AssocMUL& m);
AssocMUL operator*(AssocMUL& m, int n);

class NonAssocRMUL {
public:
    NonAssocRMUL(double d) : fx(d) {}
    bool operator==(double d) { return d == fx; }

    friend NonAssocRMUL operator*(int n, NonAssocRMUL& m);

private:
    double fx;
};

NonAssocRMUL operator*(int n, NonAssocRMUL& m);


//----------------------------------------------------------------------------
// multi-lookup
namespace MultiLookup {

struct Vector2 {
    Vector2(double d1, double d2) : x(d1), y(d2) {}
    double x, y;
};

double operator*(const Vector2&, const Vector2&);
Vector2 operator*(const Vector2&, double);

double operator/(const Vector2&, const Vector2&);
Vector2 operator/(const Vector2&, double);

double operator+(const Vector2&, const Vector2&);
Vector2 operator+(const Vector2&, double);

double operator-(const Vector2&, const Vector2&);
Vector2 operator-(const Vector2&, double);

} // namespace MultiLookup


//----------------------------------------------------------------------------
// unary functions
struct SomeGlobalNumber {
    SomeGlobalNumber(int n) : i(n) {}
    int i;
};

SomeGlobalNumber operator-(const SomeGlobalNumber&);
SomeGlobalNumber operator+(const SomeGlobalNumber&);
SomeGlobalNumber operator~(const SomeGlobalNumber&);

namespace Unary {

struct SomeNumber {
    SomeNumber(int n) : i(n) {}
    int i;
};

SomeNumber operator-(const SomeNumber&);
SomeNumber operator+(const SomeNumber&);
SomeNumber operator~(const SomeNumber&);

}


//----------------------------------------------------------------------------
// comma operator (mapped to __comma__)
struct CommaOperator {
    CommaOperator(int i) : fInt(i) {}
    CommaOperator& operator,(int i) { fInt += i; return *this; }
    int fInt;
};

#endif // !CPPYY_TEST_OPERATORS_H
