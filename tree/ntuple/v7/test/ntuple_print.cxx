#include "ntuple_test.hxx"

TEST(RNtuplePrint, FullString)
{
   FileRaii fileGuard("test_ntuple_print_fullstring.root");
   {
      auto model = RNTupleModel::Create();
      auto fieldPx = model->MakeField<float>("px");
      auto fieldPy = model->MakeField<float>("py");
      auto fieldPz = model->MakeField<float>("pz");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), "Staff", fileGuard.GetPath());
      *fieldPx = 1.0;
      *fieldPy = 1.0;
      *fieldPz = 1.0;
      ntuple->Fill();
   }
   auto ntuple2 = RNTupleReader::Open("Staff", fileGuard.GetPath());
   std::ostringstream os;
   ntuple2->PrintInfo(ROOT::Experimental::ENTupleInfo::kSummary, os);
   std::string fString{std::string("")
       + "************************************ NTUPLE ************************************\n"
       + "* N-Tuple : Staff                                                              *\n"
       + "* Entries : 1                                                                  *\n"
       + "********************************************************************************\n"
       + "* Field 1   : px (float)                                                       *\n"
       + "* Field 2   : py (float)                                                       *\n"
       + "* Field 3   : pz (float)                                                       *\n"
       + "********************************************************************************\n"};
   EXPECT_EQ(fString, os.str());
}

TEST(RNtuplePrint, Int)
{
   std::stringstream os;
   RPrintSchemaVisitor visitor(os);
   RField<int> testField("intTest");
   testField.AcceptVisitor(visitor);
   std::string expected{std::string("")
       + "* Field 1   : intTest (std::int32_t)                                           *\n"};
   EXPECT_EQ(expected, os.str());
}

TEST(RNtuplePrint, Float)
{
   std::stringstream os;
   RPrintSchemaVisitor visitor(os, 'a');
   RField<float> testField("floatTest");
   testField.AcceptVisitor(visitor);
   std::string expected{std::string("")
       + "a Field 1   : floatTest (float)                                                a\n"};
   EXPECT_EQ(expected, os.str());
}

TEST(RNtuplePrint, Vector)
{
   std::stringstream os;
   RPrepareVisitor prepVisitor;
   RField<std::vector<float>> testField("floatVecTest");
   testField.AcceptVisitor(prepVisitor);
   RPrintSchemaVisitor visitor(os, '$');
   visitor.SetDeepestLevel(prepVisitor.GetDeepestLevel());
   visitor.SetNumFields(prepVisitor.GetNumFields());
   testField.AcceptVisitor(visitor);
   std::string expected{std::string("")
       + "$ Field 1       : floatVecTest (std::vector<float>)                            $\n"
       + "$   Field 1.1   : float (float)                                                $\n"};
   EXPECT_EQ(expected, os.str());
}

TEST(RNtuplePrint, VectorNested)
{
   std::stringstream os;
   RPrepareVisitor prepVisitor;
   RField<std::vector<std::vector<float>>> testField("floatVecVecTest");
   testField.AcceptVisitor(prepVisitor);
   RPrintSchemaVisitor visitor(os, 'x');
   visitor.SetDeepestLevel(prepVisitor.GetDeepestLevel());
   visitor.SetNumFields(prepVisitor.GetNumFields());
   testField.AcceptVisitor(visitor);
   std::string expected{std::string("")
       + "x Field 1           : floatVecVecTest (std::vector<std::vector<float>>)        x\n"
       + "x   Field 1.1       : std::vector<float> (std::vector<float>)                  x\n"
       + "x     Field 1.1.1   : float (float)                                            x\n"};
   EXPECT_EQ(expected, os.str());
}

TEST(RNtuplePrint, NarrowManyEntriesVecVecTraverse)
{
   std::stringstream os;
   RPrepareVisitor prepVisitor;
   RField<std::vector<std::vector<float>>> testField("floatVecVecTest");
   testField.AcceptVisitor(prepVisitor);
   RPrintSchemaVisitor visitor(os, ' ', 25);
   visitor.SetDeepestLevel(prepVisitor.GetDeepestLevel());
   visitor.SetNumFields(prepVisitor.GetNumFields());
   testField.AcceptVisitor(visitor);
   std::string expected{std::string("")
       + "  Field 1    : floatV... \n"
       + "    Field... : std::v... \n"
       + "      Fie... : float ... \n"};
   EXPECT_EQ(expected, os.str());
}

/* Currently the width can't be set by PrintInfo(). This test will be enabled when this feature is added.
TEST(RNTuplePrint, TooShort)
{
FileRaii fileGuard("test.root");
{
   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float>("pt");
   auto ntuple = RNTupleWriter::Recreate(std::move(model), "Staff", "test.root");
}
auto ntuple2 = RNTupleReader::Open("Staff", "test.root");
std::ostringstream os;
ntuple2->PrintInfo(ROOT::Experimental::ENTupleInfo::kSummary, os, '+', 29);
std::string fString{"The width is too small! Should be at least 30.\n"};
EXPECT_EQ(fString, os.str());
}
*/
