/// \file
/// \ingroup tutorial_vecops
/// \notebook -nodraw
/// In this tutorial we learn how the RVec class can be used to
/// express easily mathematical operations involving arrays and scalars.
///
/// \macro_code
/// \macro_output
///
/// \date May 2018
/// \author Danilo Piparo

using namespace ROOT::VecOps;

void vo002_VectorCalculations()
{

   // Operations on RVec instances are made to be *fast* (vectorisation is exploited)
   // and easy to use.
   RVec<float> v1{1., 2., 3.};
   RVec<float> v2{4., 5., 6.};

   // Arithmetic operations are to be intended on pairs of elements with identical index
   auto v_sum = v1 + v2;
   auto v_mul = v1 * v2;

   // Easy to inspect:
   std::cout << "v1 = " << v1 << "\n"
             << "v2 = " << v2 << "\n"
             << "v1 + v2 = " << v_sum << "\n"
             << "v1 * v2 = " << v_mul << std::endl;

   // It's also possible to mix scalars and RVecs
   auto v_diff_s_0 = v1 - 2;
   auto v_diff_s_1 = 2 - v1;
   auto v_div_s_0 = v1 / 2.;
   auto v_div_s_1 = 2. / v1;

   std::cout << v1 << " - 2 = " << v_diff_s_0 << "\n"
             << "2 - " << v1 << " = " << v_diff_s_1 << "\n"
             << v1 << " / 2 = " << v_div_s_0 << "\n"
             << "2 / " << v1 << " = " << v_div_s_1 << std::endl;

   // Dot product and the extraction of quantities such as Mean, Min and Max
   // are also easy to express (see here for the full list:
   // https://root.cern/doc/master/namespaceROOT_1_1VecOps.html)
   auto v1_mean = Mean(v1);
   auto v1_dot_v2 = Dot(v1, v2);

   std::cout << "Mean of " << v1 << " is " << v1_mean << "\n"
             << "Dot product of " << v1 << " and " << v2 << " is " << v1_dot_v2 << std::endl;

   // Most used mathematical functions are supported
   auto v_exp = exp(v1);
   auto v_log = log(v1);
   auto v_sin = sin(v1);

   std::cout << "exp(" << v1 << ") = " << v_exp << "\n"
             << "log(" << v1 << ") = " << v_log << "\n"
             << "sin(" << v1 << ") = " << v_sin << std::endl;

   // Even an optimised version of the functions is available
   // provided that VDT is not disabled during the configuration
#ifdef R__HAS_VDT
   auto v_fast_exp = fast_exp(v1);
   auto v_fast_log = fast_log(v1);
   auto v_fast_sin = fast_sin(v1);

   std::cout << "fast_exp(" << v1 << ") = " << v_fast_exp << "\n"
             << "fast_log(" << v1 << ") = " << v_fast_log << "\n"
             << "fast_sin(" << v1 << ") = " << v_fast_sin << std::endl;

   // It may happen that a custom operation needs to be applied to the RVec.
   // In this case, the Map utitlity can be used:
   auto v_transf = Map(v1, [](double x) { return x * 2 / 3; });

   std::cout << "Applying [](double x){return x * 2 / 3;} to " << v1 << " leads to " << v_transf << "\n";
#endif
}
