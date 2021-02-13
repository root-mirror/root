/// \file TTabComTests.cxx
///
/// \brief The file contain unit tests which test the TTabCom class.
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date Jul, 2020
///
/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOTUnitTestSupport.h"
#include <ROOT/RConfig.hxx>

#include "TTabCom.h"

#include <string>

TEST(TTabComTests, Sanity)
{
   ASSERT_FALSE(gTabCom);
}

static std::string SortCompletions(std::string result, std::set<std::string> ignore)
{
   std::replace(result.begin(), result.end(), '\n', ' ');
   std::istringstream iss(result);
   std::vector<std::string> completions{std::istream_iterator<std::string>{iss},
         std::istream_iterator<std::string>{}};
   std::sort(completions.begin(), completions.end());
   result = "";
   for (size_t i = 0, e = completions.size(); i < e; ++i) {
      if (ignore.count(completions[i]))
         continue;
      result += completions[i];
      if (i != e-1)
         result += ' ';
   }

   return result;
}

static std::string GetCompletions(const std::string& pattern, std::set<std::string> ignore = {})
{
   static auto ttc = new TTabCom;
   const size_t lineBufSize = 2*1024;  // must be equal to/larger than BUF_SIZE in TTabCom.cxx
   std::unique_ptr<char[]> completed(new char[lineBufSize]);
   strncpy(completed.get(), pattern.c_str(), lineBufSize);
   completed[lineBufSize-1] = '\0';
   int pLoc = strlen(completed.get());
   std::ostringstream oss;
   ttc->Hook(completed.get(), &pLoc, oss);
   return SortCompletions(oss.str(), ignore);
}

TEST(TTabComTests, CompleteTH1)
{
   // FIXME: The first call is unsuccessful due to a bug in the TTabCom::Hook
   // on some systems.
   GetCompletions("TH1");
   std::string expected = "TH1 TH1C TH1D"
#if defined(R__USE_CXXMODULES) && defined(R__HAS_DATAFRAME)
      // FIXME: See ROOT-10989
      " TH1DModel"
#endif
      " TH1Editor TH1F TH1I TH1K TH1S";

   ASSERT_STREQ(expected.c_str(), GetCompletions("TH1").c_str());
}

TEST(TTabComTests, CompleteTProfile)
{
   std::string expected = "TProfile"
#if defined(R__USE_CXXMODULES) && defined(R__HAS_DATAFRAME)
      // FIXME: See ROOT-10989
      " TProfile1DModel"
#endif
      " TProfile2D"
#if defined(R__USE_CXXMODULES) && defined(R__HAS_DATAFRAME)
      // FIXME: See ROOT-10989
      " TProfile2DModel"
#endif
      " TProfile2Poly TProfile2PolyBin TProfile3D";

   ASSERT_STREQ(expected.c_str(), GetCompletions("TProfile").c_str());
}

TEST(TTabComTests, CompleteTObj)
{
   std::string expected = "TObjArray TObjArrayIter TObjLink TObjOptLink TObjString"
      " TObject TObjectRefSpy TObjectSpy TObjectTable";
   // FIXME: See ROOT-10989
   ASSERT_STREQ(expected.c_str(), GetCompletions("TObj",
                                                 /*ignore=*/{"TObjectDisplayItem", "TObjectDrawable", "TObjectHolder",
                                                             "TObjectItem", "TObjectElement", "TObject::EDeprecatedStatusBits", "TObject::EStatusBits"}).c_str());
}
