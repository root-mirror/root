/****** Run RDataFrame tests both with and without IMT enabled *******/
#include <gtest/gtest.h>
#include <ROOT/RDataFrame.hxx>

using namespace ROOT;
using namespace ROOT::RDF;
using namespace ROOT::VecOps;

static const std::string DisplayPrintDefaultRows(
   "b1 | b2  | b3        | \n0  | 1   | 2.0000000 | \n   | ... |           | \n   | 3   |           | \n0  | 1   | "
   "2.0000000 | \n   | ... |           | \n   | 3   |           | \n0  | 1   | 2.0000000 | \n   | ... |           | \n "
   "  | 3   |           | \n0  | 1   | 2.0000000 | \n   | ... |           | \n   | 3   |           | \n0  | 1   | "
   "2.0000000 | \n   | ... |           | \n   | 3   |           | \n");

static const std::string DisplayAsStringDefaultRows(
   "b1 | b2  | b3        | \n0  | 1   | 2.0000000 | \n   | 2   |           | \n   | 3   |           | \n0  | 1   | "
   "2.0000000 | \n   | 2   |           | \n   | 3   |           | \n0  | 1   | 2.0000000 | \n   | 2   |           | \n "
   "  | 3   |           | \n0  | 1   | 2.0000000 | \n   | 2   |           | \n   | 3   |           | \n0  | 1   | "
   "2.0000000 | \n   | 2   |           | \n   | 3   |           | \n   |     |           | \n");

TEST(RDFDisplayTests, DisplayNoJitDefaultRows)
{
   RDataFrame rd1(10);
   auto dd = rd1.Define("b1", []() { return 0; })
                .Define("b2",
                        []() {
                           return std::vector<int>({1, 2, 3});
                        })
                .Define("b3", []() { return 2.; })
                .Display<int, std::vector<int>, double>({"b1", "b2", "b3"});

   // Testing the std output printing
   std::cout << std::flush;
   // Redirect cout.
   std::streambuf *oldCoutStreamBuf = std::cout.rdbuf();
   std::ostringstream strCout;
   std::cout.rdbuf(strCout.rdbuf());
   dd->Print();
   // Restore old cout.
   std::cout.rdbuf(oldCoutStreamBuf);

   EXPECT_EQ(strCout.str(), DisplayPrintDefaultRows);

   // Testing the string returned
   EXPECT_EQ(dd->AsString(), DisplayAsStringDefaultRows);
}

TEST(RDFDisplayTests, DisplayJitDefaultRows)
{
   RDataFrame rd1(10);
   auto dd = rd1.Define("b1", []() { return 0; })
                .Define("b2",
                        []() {
                           return std::vector<int>({1, 2, 3});
                        })
                .Define("b3", []() { return 2.; })
                .Display({"b1", "b2", "b3"});

   // Testing the std output printing
   std::cout << std::flush;
   // Redirect cout.
   std::streambuf *oldCoutStreamBuf = std::cout.rdbuf();
   std::ostringstream strCout;
   std::cout.rdbuf(strCout.rdbuf());
   dd->Print();
   // Restore old cout.
   std::cout.rdbuf(oldCoutStreamBuf);

   EXPECT_EQ(strCout.str(), DisplayPrintDefaultRows);

   // Testing the string returned
   EXPECT_EQ(dd->AsString(), DisplayAsStringDefaultRows);
}

TEST(RDFDisplayTests, DisplayRegexDefaultRows)
{
   RDataFrame rd1(10);
   auto dd = rd1.Define("b1", []() { return 0; })
                .Define("b2",
                        []() {
                           return std::vector<int>({1, 2, 3});
                        })
                .Define("b3", []() { return 2.; })
                .Display("");

   // Testing the std output printing
   std::cout << std::flush;
   // Redirect cout.
   std::streambuf *oldCoutStreamBuf = std::cout.rdbuf();
   std::ostringstream strCout;
   std::cout.rdbuf(strCout.rdbuf());
   dd->Print();
   // Restore old cout.
   std::cout.rdbuf(oldCoutStreamBuf);

   EXPECT_EQ(strCout.str(), DisplayPrintDefaultRows);

   // Testing the string returned
   EXPECT_EQ(dd->AsString(), DisplayAsStringDefaultRows);
}

static const std::string
   DisplayPrintTwoRows("b1 | b2  | b3        | \n0  | 1   | 2.0000000 | \n   | ... |           | \n   | 3   |          "
                       " | \n0  | 1   | 2.0000000 | \n   | ... |           | \n   | 3   |           | \n");

static const std::string DisplayAsStringTwoRows(
   "b1 | b2  | b3        | \n0  | 1   | 2.0000000 | \n   | 2   |           | \n   | 3   |           | \n0  | 1   | "
   "2.0000000 | \n   | 2   |           | \n   | 3   |           | \n   |     |           | \n");

TEST(RDFDisplayTests, DisplayJitTwoRows)
{
   RDataFrame rd1(10);
   auto dd = rd1.Define("b1", []() { return 0; })
                .Define("b2",
                        []() {
                           return std::vector<int>({1, 2, 3});
                        })
                .Define("b3", []() { return 2.; })
                .Display({"b1", "b2", "b3"}, 2);

   // Testing the std output printing
   std::cout << std::flush;
   // Redirect cout.
   std::streambuf *oldCoutStreamBuf = std::cout.rdbuf();
   std::ostringstream strCout;
   std::cout.rdbuf(strCout.rdbuf());
   dd->Print();
   // Restore old cout.
   std::cout.rdbuf(oldCoutStreamBuf);

   EXPECT_EQ(strCout.str(), DisplayPrintTwoRows);

   // Testing the string returned
   EXPECT_EQ(dd->AsString(), DisplayAsStringTwoRows);
}

static const std::string DisplayAsStringOneColumn("b1 | \n0  | \n0  | \n0  | \n0  | \n0  | \n   | \n");
static const std::string DisplayAsStringTwoColumns(
   "b1 | b2  | \n0  | 1   | \n   | 2   | \n   | 3   | \n0  | 1   | \n   | 2   | \n   | 3   | \n0  | 1   | \n   | 2   | "
   "\n   | 3   | \n0  | 1   | \n   | 2   | \n   | 3   | \n0  | 1   | \n   | 2   | \n   | 3   | \n   |     | \n");

TEST(RDFDisplayTests, DisplayAmbiguity)
{
   // This test verifies that the correct method is called and there is no ambiguity between the JIT call to Display
   // using a column list as a parameter and the JIT call to Display using the Regexp.
   RDataFrame rd1(10);
   auto dd = rd1.Define("b1", []() { return 0; }).Define("b2", []() { return std::vector<int>({1, 2, 3}); });

   auto display_1 = dd.Display({"b1"});
   auto display_2 = dd.Display({"b1", "b2"});

   EXPECT_EQ(display_1->AsString(), DisplayAsStringOneColumn);
   EXPECT_EQ(display_2->AsString(), DisplayAsStringTwoColumns);
}

static const std::string DisplayAsStringString("b1    | \n\"foo\" | \n\"foo\" | \n      | \n");

TEST(RDFDisplayTests, DisplayPrintString)
{
   RDataFrame rd1(2);
   auto dd = rd1.Define("b1", []() { return std::string("foo"); })
                .Display({"b1"});

   // Testing the std output printing
   std::cout << std::flush;
   // Redirect cout.
   std::streambuf *oldCoutStreamBuf = std::cout.rdbuf();
   std::ostringstream strCout;
   std::cout.rdbuf(strCout.rdbuf());
   dd->Print();
   // Restore old cout.
   std::cout.rdbuf(oldCoutStreamBuf);

   // Testing the string returned
   EXPECT_EQ(dd->AsString(), DisplayAsStringString);
}
