/// \file
/// \ingroup tutorial_eve
/// Demonstrates usage of REveBoxSet class.
///
/// \image html eve_boxset.png
/// \macro_code
///
/// \author Matevz Tadel

#include "TRandom.h"
#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveBoxSet.hxx>

using namespace ROOT::Experimental;


REveBoxSet* boxset(Int_t num=100)
{
   auto eveMng = REveManager::Create();

   TRandom r(0);

   auto pal = new REveRGBAPalette(0, 130);

   auto q = new REveBoxSet("BoxSet");
   q->SetPalette(pal);
   q->Reset(REveBoxSet::kBT_FreeBox, kFALSE, 64);

#define RND_BOX(x) (Float_t)r.Uniform(-(x), (x))

   Float_t verts[24];
   for (Int_t i=0; i<num; ++i) {
      Float_t x = RND_BOX(10);
      Float_t y = RND_BOX(10);
      Float_t z = RND_BOX(10);
      Float_t a = r.Uniform(0.2, 0.5);
      Float_t d = 0.05;
      Float_t verts[24] = {
                           x - a + RND_BOX(d), y - a + RND_BOX(d), z - a + RND_BOX(d),
                           x - a + RND_BOX(d), y + a + RND_BOX(d), z - a + RND_BOX(d),
                           x + a + RND_BOX(d), y + a + RND_BOX(d), z - a + RND_BOX(d),
                           x + a + RND_BOX(d), y - a + RND_BOX(d), z - a + RND_BOX(d),
                           x - a + RND_BOX(d), y - a + RND_BOX(d), z + a + RND_BOX(d),
                           x - a + RND_BOX(d), y + a + RND_BOX(d), z + a + RND_BOX(d),
                           x + a + RND_BOX(d), y + a + RND_BOX(d), z + a + RND_BOX(d),
                           x + a + RND_BOX(d), y - a + RND_BOX(d), z + a + RND_BOX(d) };
      q->AddBox(verts);
      q->DigitValue(r.Uniform(0, 130));
   }
   q->RefitPlex();

#undef RND_BOX

   // Uncomment these two lines to get internal highlight / selection.
   q->SetPickable(1);
   q->SetAlwaysSecSelect(1);

   eveMng->GetEventScene()->AddElement(q);

   eveMng->Show();
   return q;
}

REveBoxSet* boxset_axisaligned(Float_t x=0, Float_t y=0, Float_t z=0,
                   Int_t num=100, Bool_t registerSet=kTRUE)
{
   auto eveMng = REveManager::Create();

   TRandom r(0);

   auto pal = new REveRGBAPalette(0, 130);

   auto frm = new REveFrameBox();
   frm->SetAABoxCenterHalfSize(0, 0, 0, 12, 12, 12);
   frm->SetFrameColor(kCyan);
   frm->SetBackColorRGBA(120,120,120,20);
   frm->SetDrawBack(kTRUE);

   auto q = new REveBoxSet("BoxSet");
   q->SetPalette(pal);
   q->SetFrame(frm);
   q->Reset(REveBoxSet::kBT_AABox, kFALSE, 64);
   for (Int_t i=0; i<num; ++i) {
      q->AddBox(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
                r.Uniform(0.2, 1),  r.Uniform(0.2, 1),  r.Uniform(0.2, 1));
      q->DigitValue(r.Uniform(0, 130));
   }
   q->RefitPlex();

   REveTrans& t = q->RefMainTrans();
   t.SetPos(x, y, z);

   // Uncomment these two lines to get internal highlight / selection.
   q->SetPickable(1);
   q->SetAlwaysSecSelect(1);

   if (registerSet)
   {
      eveMng->GetEventScene()->AddElement(q);
      eveMng->Show();
   }

   return q;
}

REveBoxSet* boxset_colisval(Float_t x=0, Float_t y=0, Float_t z=0,
                            Int_t num=100, Bool_t registerSet=kTRUE)
{
   auto eveMng = REveManager::Create();

   TRandom r(0);

   auto q = new REveBoxSet("BoxSet");
   q->Reset(REveBoxSet::kBT_AABox, kTRUE, 64);
   for (Int_t i=0; i<num; ++i) {
      q->AddBox(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
                r.Uniform(0.2, 1),  r.Uniform(0.2, 1),  r.Uniform(0.2, 1));
      q->DigitColor(r.Uniform(20, 255), r.Uniform(20, 255),
                    r.Uniform(20, 255), r.Uniform(20, 255));
   }
   q->RefitPlex();

   REveTrans& t = q->RefMainTrans();
   t.SetPos(x, y, z);

   if (registerSet)
   {
      eveMng->GetEventScene()->AddElement(q);
      eveMng->Show();
   }

   return q;
}

REveBoxSet* boxset_single_color(Float_t x=0, Float_t y=0, Float_t z=0,
                                Int_t num=100, Bool_t registerSet=kTRUE)
{
   auto eveMng = REveManager::Create();

   TRandom r(0);

   auto q = new REveBoxSet("BoxSet");
   q->UseSingleColor();
   q->SetMainColorPtr(new Color_t);
   q->SetMainColor(kRed);
   q->SetMainTransparency(50);
   q->Reset(REveBoxSet::kBT_AABox, kFALSE, 64);

   for (Int_t i=0; i<num; ++i) {
      q->AddBox(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
                r.Uniform(0.2, 1),  r.Uniform(0.2, 1),  r.Uniform(0.2, 1));
   }
   q->RefitPlex();

   REveTrans& t = q->RefMainTrans();
   t.SetPos(x, y, z);

   if (registerSet) {
      eveMng->GetEventScene()->AddElement(q);
      eveMng->Show();
   }

   return q;
}

/*
REveBoxSet* boxset_hex(Float_t x=0, Float_t y=0, Float_t z=0,
                       Int_t num=100, Bool_t registerSet=kTRUE)
{
   auto eveMng = REveManager::Create();

   TRandom r(0);

   auto q = new REveBoxSet("BoxSet");
   q->Reset(REveBoxSet::kBT_Hex, kTRUE, 64);

   for (Int_t i=0; i<num; ++i) {
      q->AddHex(REveVector(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10)),
                r.Uniform(0.2, 1), r.Uniform(0, 60), r.Uniform(0.2, 5));
      q->DigitColor(r.Uniform(20, 255), r.Uniform(20, 255),
                    r.Uniform(20, 255), r.Uniform(20, 255));
   }
   q->RefitPlex();

   q->SetPickable(true);
   q->SetAlwaysSecSelect(true);

   REveTrans& t = q->RefMainTrans();
   t.SetPos(x, y, z);

   if (registerSet)
   {
      eveMng->GetEventScene()->AddElement(q);
      eveMng->Show();
   }

   return q;
}
*/
