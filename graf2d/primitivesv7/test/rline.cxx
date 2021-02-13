// @(#)root/graf2d:$Id$
// Author: Sergey Linev <s.linev@gsi.de>, 2020-02-18


#include "gtest/gtest.h"

#include <ROOT/RCanvas.hxx>
#include <ROOT/RColor.hxx>
#include <ROOT/RLine.hxx>
#include <ROOT/RPadPos.hxx>

#include "TMemFile.h"

#include <string>

using namespace ROOT::Experimental;

using namespace std::string_literals;

TEST(RCanvas, SimpleIO)
{
   TMemFile* f = new TMemFile("testio","recreate");

   {
      // Create a canvas to be displayed.
      auto canvas = RCanvas::Create("Canvas Title");

      // draw several RLine objects

      for (int n=0;n<10;++n) {
         RColor col( 10 + n*10, 50 + n*5, 90 + n);

         auto line = canvas->Draw<RLine>();

         line->SetP1({0.0_normal, 0.0_normal});
         line->SetP2({1.0_normal, 1.0_normal});

         line->AttrLine().SetColor(col);

         line->SetId("line"s + std::to_string(n));
      }

      f->WriteObject(canvas.get(), "canvas");
   }


   std::shared_ptr<RCanvas> canvas_read(f->Get<RCanvas>("canvas"));

   delete f;

   ASSERT_NE(canvas_read, nullptr);

   for (int n=0;n<10;++n) {

      auto line = std::dynamic_pointer_cast<RLine>(canvas_read->FindPrimitive("line"s + std::to_string(n)));

      ASSERT_NE(line, nullptr);

      RColor col( 10 + n*10, 50 + n*5, 90 + n);

      EXPECT_EQ(line->AttrLine().GetColor(), col);
   }

}
