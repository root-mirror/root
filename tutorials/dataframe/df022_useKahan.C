/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Implement a custom action that evaluates a Kahan sum.
///
/// This tutorial shows how to implement a Kahan summation custom action.
///
/// \macro_code
/// \macro_output
///
/// \date July 2018
/// \authors Enrico Guiraud, Danilo Piparo (CERN), Massimo Tumolo (Politecnico di Torino)

template <typename T>
class KahanSum final : public ROOT::Detail::RDF::RActionImpl<class KahanSum<T>>  {
public:
   /// This type is a requirement for every helper.
   using Result_t = T;

private:
   std::vector<T> fPartialSums;
   std::vector<T> fCompensations;
   int fNSlots;

   std::shared_ptr<T> fResultSum;

   void KahanAlgorithm(const T &x, T &sum, T &compensation){
      T y = x - compensation;
      T t = sum + y;
      compensation = (t - sum) - y;
      sum = t;
   }

public:
   KahanSum(KahanSum &&) = default;
   KahanSum(const KahanSum &) = delete;

   KahanSum(const std::shared_ptr<T> &r) : fResultSum(r)
   {
      static_assert(std::is_floating_point<T>::value, "Kahan sum makes sense only on floating point numbers");

      fNSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
      fPartialSums.resize(fNSlots, 0.);
      fCompensations.resize(fNSlots, 0.);
   }

   std::shared_ptr<Result_t> GetResultPtr() const { return fResultSum; }

   void Initialize() {}
   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int slot, T x)
   {
      KahanAlgorithm(x, fPartialSums[slot], fCompensations[slot]);
   }

   template <typename V=T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<V>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      for (auto &&v : vs) {
         Exec(slot, v);
      }
   }

   void Finalize()
   {
      T sum(0) ;
      T compensation(0);
      for (int i = 0; i < fNSlots; ++i) {
         KahanAlgorithm(fPartialSums[i], sum, compensation);
      }
      *fResultSum = sum;
   }

   std::string GetActionName(){
      return "THnHelper";
   }

};

void df022_useKahan()
{
   // We enable implicit parallelism
   ROOT::EnableImplicitMT(2);

   ROOT::RDataFrame d(20);
   auto dd = d.Define("x", "(rdfentry_ %2 == 0) ? 0.00000001 : 100000000.");

   auto ptr = std::make_shared<double>();
   KahanSum<double> helper(ptr);

   auto kahanResult = dd.Book<double>(std::move(helper), {"x"});
   auto plainResult = dd.Sum<double>({"x"});

   std::cout << std::setprecision(24) << "Kahan: " << *kahanResult << " Classical: " << *plainResult << std::endl;
   // Outputs: Kahan: 1000000000.00000011920929 Classical: 1000000000
}
