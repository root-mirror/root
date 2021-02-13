#include "operators.h"

// for testing the case of virtual operator==
v_opeq_base::v_opeq_base(int val) : m_val(val) {}
v_opeq_base::~v_opeq_base() {}

bool v_opeq_base::operator==(const v_opeq_base& other) {
   return m_val == other.m_val;
}

v_opeq_derived::v_opeq_derived(int val) : v_opeq_base(val) {}
v_opeq_derived::~v_opeq_derived() {}

bool v_opeq_derived::operator==(const v_opeq_derived& other) {
   return m_val != other.m_val;
}


// for indexing tests
int& YAMatrix1::operator() (int, int) {
    return m_val;
}

const int& YAMatrix1::operator() (int, int) const {
    return m_val;
}

//-
int& YAMatrix2::operator[] (int) {
    return m_val;
}

const int& YAMatrix2::operator[] (int) const {
    return m_val;
}

//-
int& YAMatrix3::operator() (int, int) {
    return m_val;
}

const int& YAMatrix3::operator() (int, int) const {
    return m_val;
}

int& YAMatrix3::operator[] (int) {
    return m_val;
}

const int& YAMatrix3::operator[] (int) const {
    return m_val;
}

//-
int& YAMatrix4::operator[] (int) {
    return m_val;
}

const int& YAMatrix4::operator[] (int) const {
    return m_val;
}

int& YAMatrix4::operator() (int, int) {
    return m_val;
}

const int& YAMatrix4::operator() (int, int) const {
    return m_val;
}

//-
int& YAMatrix5::operator[] (int) {
    return m_val;
}

int YAMatrix5::operator[] (int) const {
    return m_val;
}

int& YAMatrix5::operator() (int, int) {
    return m_val;
}

int YAMatrix5::operator() (int, int) const {
    return m_val;
}

//-
int& YAMatrix6::operator[] (int) {
    return m_val;
}

int YAMatrix6::operator[] (int) const {
    return m_val;
}

int& YAMatrix6::operator() (int, int) {
    return m_val;
}

//-
int& YAMatrix7::operator[] (int) {
    return m_val;
}

int& YAMatrix7::operator() (int, int) {
    return m_val;
}


//- for __radd__/__rmul__, non-associative
AssocADD operator+(int n, AssocADD& a) {
    return AssocADD(n+a.fx);
}

AssocADD operator+(AssocADD& a, int n) {
    return AssocADD(a.fx+n);
}

NonAssocRADD operator+(int n, NonAssocRADD& a) {
    return NonAssocRADD(n+a.fx);
}

AssocMUL operator*(int n, AssocMUL& m) {
    return AssocMUL(n*m.fx);
}

AssocMUL operator*(AssocMUL& m, int n) {
    return AssocMUL(m.fx*n);
}

NonAssocRMUL operator*(int n, NonAssocRMUL& m) {
    return NonAssocRMUL(n*m.fx);
}


//- for multi-lookup
double MultiLookup::operator*(const Vector2& v1, const Vector2& v2) {
    return v1.x*v2.x + v1.y*v2.y;
}

MultiLookup::Vector2 MultiLookup::operator*(const Vector2& v, double a) {
    return Vector2{v.x*a, v.y*a};
}

double MultiLookup::operator/(const Vector2& v1, const Vector2& v2) {
    return v1.x/v2.x + v1.y/v2.y;
}

MultiLookup::Vector2 MultiLookup::operator/(const Vector2& v, double a) {
    return Vector2{v.x/a, v.y/a};
}

double MultiLookup::operator+(const Vector2& v1, const Vector2& v2) {
    return v1.x+v2.x + v1.y+v2.y;
}

MultiLookup::Vector2 MultiLookup::operator+(const Vector2& v, double a) {
    return Vector2{v.x+a, v.y+a};
}

double MultiLookup::operator-(const Vector2& v1, const Vector2& v2) {
    return v1.x-v2.x + v1.y-v2.y;
}

MultiLookup::Vector2 MultiLookup::operator-(const Vector2& v, double a) {
    return Vector2{v.x-a, v.y-a};
}


//- for unary functions
SomeGlobalNumber operator-(const SomeGlobalNumber& n) {
    return SomeGlobalNumber{-n.i};
}

SomeGlobalNumber operator+(const SomeGlobalNumber& n) {
    return SomeGlobalNumber{+n.i};
}

SomeGlobalNumber operator~(const SomeGlobalNumber& n) {
    return SomeGlobalNumber{~n.i};
}

Unary::SomeNumber Unary::operator-(const SomeNumber& n) {
    return SomeNumber{-n.i};
}

Unary::SomeNumber Unary::operator+(const SomeNumber& n) {
    return SomeNumber{+n.i};
}

Unary::SomeNumber Unary::operator~(const SomeNumber& n) {
    return SomeNumber{~n.i};
}
