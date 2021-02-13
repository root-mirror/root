File features.h
===============

.. code-block:: c++
    :linenos:

    #include <cmath>
    #include <iostream>
    #include <vector>


    //-----
    unsigned int gUint = 0;

    //-----
    class Abstract {
    public:
        virtual ~Abstract() {}
        virtual void abstract_method() = 0;
        virtual void concrete_method() = 0;
    };

    void Abstract::concrete_method() {
        std::cout << "called Abstract::concrete_method" << std::endl;
    }

    //-----
    class Concrete : Abstract {
    public:
        Concrete(int n=42) : m_int(n), m_const_int(17) {}
        ~Concrete() {}

        virtual void abstract_method() {
            std::cout << "called Concrete::abstract_method" << std::endl;
        }

        virtual void concrete_method() {
            std::cout << "called Concrete::concrete_method" << std::endl;
        }

        void array_method(int* ad, int size) {
            for (int i=0; i < size; ++i)
                std::cout << ad[i] << ' ';
            std::cout << '\n';
        }

        void array_method(double* ad, int size) {
            for (int i=0; i < size; ++i)
                std::cout << ad[i] << ' ';
            std::cout << '\n';
        }

        void uint_ref_assign(unsigned int& target, unsigned int value) {
            target = value;
        }

        Abstract* show_autocast() {
            return this;
        }

        operator const char*() {
            return "Hello operator const char*!";
        }

    public:
        double m_data[4];
        int m_int;
        const int m_const_int;

        static int s_int;
    };

    typedef Concrete Concrete_t;

    int Concrete::s_int = 321;

    void call_abstract_method(Abstract* a) {
        a->abstract_method();
    }

    //-----
    int global_function(int) {
        return 42;
    }

    double global_function(double) {
        return std::exp(1);
    }

    int call_int_int(int (*f)(int, int), int i1, int i2) {
        return f(i1, i2);
    }

    template<class A, class B, class C = A>
    C multiply(A a, B b) {
        return C{a*b};
    }

    //-----
    namespace Namespace {

        class Concrete {
        public:
            class NestedClass {
            public:
                std::vector<int> m_v;
            };

        };

        int global_function(int i) {
            return 2*::global_function(i);
        }

        double global_function(double d) {
            return 2*::global_function(d);
        }

    } // namespace Namespace

    //-----
    enum EFruit {kApple=78, kBanana=29, kCitrus=34};
    enum class NamedClassEnum { E1 = 42 };

    //-----
    void throw_an_error(int i);

    class SomeError : public std::exception {
    public:
        explicit SomeError(const std::string& msg) : fMsg(msg) {}
        const char* what() const throw() override { return fMsg.c_str(); }

    private:
        std::string fMsg;
    };

    class SomeOtherError : public SomeError {
    public:
        explicit SomeOtherError(const std::string& msg) : SomeError(msg) {}
        SomeOtherError(const SomeOtherError& s) : SomeError(s) {}
    };
