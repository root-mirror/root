#include <memory>
#include <string>
#include <vector>

// Python
struct _object;
typedef _object PyObject;


namespace pyzables {

//===========================================================================
class SomeDummy1 {};
class SomeDummy2 {};


//===========================================================================
class NakedBuffers {
public:
    NakedBuffers(int size, double valx, double valy);
    NakedBuffers(const NakedBuffers&) = delete;
    NakedBuffers& operator=(const NakedBuffers&) = delete;
    ~NakedBuffers();

public:
    int GetN();
    double* GetX();
    double* GetY();

private:
    double* m_Xbuf;
    double* m_Ybuf;
    int m_size;
};

template<typename C>
class NakedBuffers2 {
public:
    NakedBuffers2(int size, double valx, double valy) : m_Xbuf(size), m_Ybuf(size) {
        for (int i=0; i<size; ++i) {
            m_Xbuf[i] = valx*i;
            m_Ybuf[i] = valy*i;
        }
    }

public:
    int GetN() { return m_Xbuf.size(); }
    double* GetX() { return m_Xbuf.data(); }
    double* GetY() { return m_Ybuf.data(); }

private:
    C m_Xbuf;
    C m_Ybuf;
    int m_size;
};

class Vector : public std::vector<double> {
public:
    Vector(int size) : std::vector<double>(size) {}
};


//===========================================================================
class MyBase {
public:
    virtual ~MyBase();
};
class MyDerived : public MyBase {
public:
    virtual ~MyDerived();
};

MyBase* GimeDerived();


//===========================================================================
class Countable {
public:
    Countable() { ++sInstances; }
    Countable(const Countable&) { ++sInstances; }
    Countable& operator=(const Countable&) { return *this; }
    ~Countable() { --sInstances; }

public:
    virtual const char* say_hi() { return "Hi!"; }

public:
    unsigned int m_check = 0xcdcdcdcd;

public:
    static int sInstances;
};

typedef std::shared_ptr<Countable> SharedCountable_t; 
extern SharedCountable_t mine;

void renew_mine();

SharedCountable_t gime_mine();
SharedCountable_t* gime_mine_ptr();
SharedCountable_t& gime_mine_ref();

unsigned int pass_mine_sp(SharedCountable_t p);
unsigned int pass_mine_sp_ref(SharedCountable_t& p);
unsigned int pass_mine_sp_ptr(SharedCountable_t* p);

unsigned int pass_mine_rp(Countable);
unsigned int pass_mine_rp_ref(const Countable&);
unsigned int pass_mine_rp_ptr(const Countable*);

Countable* gime_naked_countable();


//===========================================================================
class unknown_iterator;
class IndexableBase {
public:
    unknown_iterator* begin() { return nullptr; }
    unknown_iterator* end() { return (unknown_iterator*)1; }
    int operator[](int) { return 42; }
    int size() { return 1; }
};

class IndexableDerived : public IndexableBase {};


//===========================================================================
class WithCallback1 {
public:
    WithCallback1(int i);

public:
    int get_int();
    void set_int(int i);

private:
    int m_int;

public:
    static void __cppyy_explicit_pythonize__(PyObject* klass, const std::string&);
    static std::string klass_name;
};

class WithCallback2 {
public:
    WithCallback2(int i);

public:
    int get_int();
    void set_int(int i);

protected:
    int m_int;

public:
    static void __cppyy_pythonize__(PyObject* klass, const std::string&);
    static std::string klass_name;
};

class WithCallback3 : public WithCallback2 {
public:
    using WithCallback2::WithCallback2;

public:
    int get_int();
    void set_int(int i);
};

} // namespace pyzables
