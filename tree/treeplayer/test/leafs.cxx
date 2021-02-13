#include "ROOT/RMakeUnique.hxx"
#include "TInterpreter.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "gtest/gtest.h"

#include "data.h"

#include "RErrorIgnoreRAII.hxx"

TEST(TTreeReaderLeafs, LeafListCaseA) {
   // From "Case A" of the TTree class doc:
   const char* str = "This is a null-terminated string literal";
   signed char SChar = 2;
   unsigned char UChar = 3;
   signed short SShort = 4;
   unsigned short UShort = 5;
   signed int SInt = 6;
   unsigned int UInt = 7;
   float Float = 8.;
   double Double = 9.;
   long long SLL = 10;
   unsigned long long ULL = 11;
   long SL = 12;
   unsigned long UL = 13;
   bool Bool = true;

   auto tree = std::make_unique<TTree>("T", "In-memory test tree");
   tree->Branch("C", &str, "C/C");
   tree->Branch("B", &SChar, "B/B");
   tree->Branch("b", &UChar, "b/b");
   tree->Branch("S", &SShort, "S/S");
   tree->Branch("s", &UShort, "s/s");
   tree->Branch("I", &SInt, "I/I");
   tree->Branch("i", &UInt, "i/i");
   tree->Branch("F", &Float, "F/F");
   tree->Branch("D", &Double, "D/D");
   tree->Branch("L", &SLL, "L/L");
   tree->Branch("l", &ULL, "l/l");
   tree->Branch("G", &SL, "G/G");
   tree->Branch("g", &UL, "g/g");
   tree->Branch("O", &Bool, "O/O");

   tree->Fill();
   tree->Fill();
   tree->Fill();

   TTreeReader TR(tree.get());
   //TTreeReaderValue<const char*> trStr(TR, "C");
   TTreeReaderValue<signed char> trSChar(TR, "B");
   TTreeReaderValue<unsigned char> trUChar(TR, "b");
   TTreeReaderValue<signed short> trSShort(TR, "S");
   TTreeReaderValue<unsigned short> trUShort(TR, "s");
   TTreeReaderValue<signed int> trSInt(TR, "I");
   TTreeReaderValue<unsigned int> trUInt(TR, "i");
   TTreeReaderValue<float> trFloat(TR, "F");
   TTreeReaderValue<double> trDouble(TR, "D");
   TTreeReaderValue<signed long long> trSLL(TR, "L");
   TTreeReaderValue<unsigned long long> trULL(TR, "l");
   TTreeReaderValue<signed long> trSL(TR, "G");
   TTreeReaderValue<unsigned long> trUL(TR, "g");
   TTreeReaderValue<bool> trBool(TR, "O");

   TR.SetEntry(1);

   //EXPECT_STREQ(str, *trStr);
   EXPECT_EQ(SChar, *trSChar);
   EXPECT_EQ(UChar, *trUChar);
   EXPECT_EQ(SShort, *trSShort);
   EXPECT_EQ(UShort, *trUShort);
   EXPECT_EQ(SInt, *trSInt);
   EXPECT_EQ(UInt, *trUInt);
   EXPECT_FLOAT_EQ(Float, *trFloat);
   EXPECT_DOUBLE_EQ(Double, *trDouble);
   EXPECT_EQ(SLL, *trSLL);
   EXPECT_EQ(ULL, *trULL);
   EXPECT_EQ(SL, *trSL);
   EXPECT_EQ(UL, *trUL);
   EXPECT_EQ(Bool, *trBool);
}



std::unique_ptr<TTree> CreateTree() {
   TInterpreter::EErrorCode error = TInterpreter::kNoError;
   gInterpreter->ProcessLine("#include \"data.h\"", &error);
   if (error != TInterpreter::kNoError)
      return {};

   Data data;
   auto tree = std::make_unique<TTree>("T", "test tree");
   tree->Branch("Data", &data);
   data.fArray = new double[4]{12., 13., 14., 15.};
   data.fSize = 4;
   data.fUArray = new float[2]{42., 43.};
   data.fUSize = 2;
   data.fVec = { 17., 18., 19., 20., 21., 22.};
   data.fDouble32 = 17.;
   data.fFloat16 = 44.;
   tree->Fill();

   data.fVec.clear();
   data.fVec.resize(3210, 1001.f); // ROOT-8747
   tree->Fill();

   data.fVec.clear();
   data.fVec.resize(2, 42.f); // ROOT-8747
   tree->Fill();

   tree->ResetBranchAddresses();
   return tree;
}

TEST(TTreeReaderLeafs, LeafList) {
   auto tree = CreateTree();
   ASSERT_NE(nullptr, tree.get());

   TTreeReader tr(tree.get());
   TTreeReaderArray<double> arr(tr, "fArray");
   TTreeReaderArray<float> arrU(tr, "fUArray");
   TTreeReaderArray<double> vec(tr, "fVec");
   TTreeReaderValue<Double32_t> d32(tr, "fDouble32");
   TTreeReaderValue<Float16_t> f16(tr, "fFloat16");

   tr.Next();
   EXPECT_EQ(4u, arr.GetSize());
   EXPECT_EQ(2u, arrU.GetSize());
   EXPECT_EQ(6u, vec.GetSize());
   //FAILS EXPECT_FLOAT_EQ(13., arr[1]);
   //FAILS EXPECT_DOUBLE_EQ(43., arrU[1]);
   EXPECT_DOUBLE_EQ(19., vec[2]);
   EXPECT_DOUBLE_EQ(17., vec[0]);
   // T->Scan("fUArray") claims fUArray only has one instance per row.

   EXPECT_FLOAT_EQ(17, *d32);
   EXPECT_FLOAT_EQ(44, *f16);

   tr.Next();
   EXPECT_FLOAT_EQ(1001.f, vec[1]); // ROOT-8747
   EXPECT_EQ(3210u, vec.GetSize());

   tr.Next();
   EXPECT_FLOAT_EQ(42.f, vec[1]); // ROOT-8747
   EXPECT_EQ(2u, vec.GetSize());
   tree.release();
}

TEST(TTreeReaderLeafs, TArrayD) {
   // https://root-forum.cern.ch/t/tarrayd-in-ttreereadervalue/24495
   auto tree = std::make_unique<TTree>("TTreeReaderLeafsTArrayD", "In-memory test tree");
   TArrayD arrD(7);
   for (int i = 0; i < arrD.GetSize(); ++i)
      arrD.SetAt(i + 2., i);
   tree->Branch("arrD", &arrD);

   tree->Fill();
   tree->Fill();
   tree->Fill();

   tree->ResetBranchAddresses();

   TTreeReader tr(tree.get());
   TTreeReaderValue<TArrayD> arr(tr, "arrD");

   tr.SetEntry(1);

   EXPECT_EQ(7, arr->GetSize());
   EXPECT_DOUBLE_EQ(4., (*arr)[2]);
   EXPECT_DOUBLE_EQ(2., (*arr)[0]);
}

TEST(TTreeReaderLeafs, ArrayWithReaderValue)
{
   // reading a float[] with a TTreeReaderValue should cause an error

   auto tree = std::make_unique<TTree>("arraywithreadervaluetree", "test tree");
   std::vector<double> arr = {42., 84.};
   tree->Branch("arr", arr.data(), "arr[2]/D");
   tree->Fill();
   tree->ResetBranchAddresses();

   TTreeReader tr(tree.get());
   TTreeReaderValue<double> valueOfArr(tr, "arr");
   {
      RErrorIgnoreRAII errorIgnRAII;
      tr.Next();
      *valueOfArr;
   }
   EXPECT_FALSE(valueOfArr.IsValid());
}


TEST(TTreeReaderLeafs, NamesWithDots)
{
   gInterpreter->ProcessLine(".L data.h+");

   TTree tree("t", "t");
   V v;
   v.a = 64;
   tree.Branch("v", &v, "a/I");
   W w;
   w.v.a = 132;
   tree.Branch("w", &w);
   tree.Fill();

   TTreeReader tr(&tree);
   TTreeReaderValue<int> rv(tr, "v.a");
   tr.Next();
   EXPECT_EQ(*rv, 64) << "The wrong leaf has been read!";
}

// This is ROOT-9743
struct Event {
   float bla = 3.14f;
   int truth_type = 1;
};

TEST(TTreeReaderLeafs, MultipleReaders) {
   TTree t("t","t");
   Event event;
   t.Branch ("event", &event, "bla/F:truth_type/I");
   t.Fill();

   TTreeReader r(&t);
   TTreeReaderValue<int> v1(r, "event.truth_type");
   TTreeReaderValue<int> v2(r, "event.truth_type");
   TTreeReaderValue<int> v3(r, "event.truth_type");

   r.Next();
   EXPECT_EQ(*v1, 1) << "Wrong value read for rv1!";
   EXPECT_EQ(*v2, 1) << "Wrong value read for rv2!";
   EXPECT_EQ(*v3, 1) << "Wrong value read for rv3!";
}
