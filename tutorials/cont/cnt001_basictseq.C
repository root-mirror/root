/// \file
/// \ingroup tutorial_cont
/// \notebook -nodraw
/// Example showing possible usages of the TSeq class.
///
/// \macro_code
/// \macro_output
///
/// \author Danilo Piparo

using namespace ROOT;

void cnt001_basictseq()
{
   cout << "Loop on sequence of integers from 0 to 10" << endl;
   for (auto i : TSeqI(10)) {
      cout << "Element " << i << endl;
   }
   //
   cout << "Loop on sequence of integers from -5 to 29 in steps of 6" << endl;
   for (auto i : TSeqI(-5, 29, 6)) {
      cout << "Element " << i << endl;
   }
   //
   cout << "Loop backwards on sequence of integers from 50 to 30 in steps of 3" << endl;
   for (auto i : TSeqI(50, 30, -3)) {
      cout << "Element " << i << endl;
   }
   //
   cout << "stl algorithm, for_each" << endl;
   TSeqUL ulSeq(2,30,3);
   std::for_each(std::begin(ulSeq),std::end(ulSeq),[](ULong_t i){cout << "For each: " << i <<endl;});

   cout << "Random access: 3rd element is " << ulSeq[2] << endl;
   //
   cout << "Loop using MakeSeq" << endl;
   for (auto i : MakeSeq(1000000000000UL, 1000000000003UL)) {
      cout << "Element " << i << endl;
   }
}

