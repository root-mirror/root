#include <vector>


namespace CNS {

//===========================================================================
double sumit(const std::vector<double>&);
double sumit(const std::vector<double>&, const std::vector<double>&);

class Counter {
public:
    Counter();
    Counter(const Counter&);
    Counter& operator=(const Counter&);
    ~Counter();

    static int s_count;
};

double myhowmany(const std::vector<Counter>&);
double myhowmany(const std::vector<Counter>&, const std::vector<Counter>&);

int sumints(const std::vector<int>&);
int sumints(const std::vector<int>&, const std::vector<int>&);

double notallowed(std::vector<double>&);

} // namespace CNS
