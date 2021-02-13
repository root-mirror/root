/// \file
/// \ingroup tutorial_eve
/// Demonstrates usage of TEveTrackPRopagator with different magnetic
/// field configurations.
/// Needs to be run in compiled mode.
/// root
/// ~~~{.cpp}
///   .L track.C+
///   track(3, kTRUE)
/// ~~~
/// void track(Int_t mode = 5, Bool_t isRungeKutta = kTRUE)
/// Modes are
///  0. B = 0, no difference between signed and charge particles;
///  1. constant B field (along z, but could have arbitrary direction);
///  2. variable B field, sign change at  R = 200 cm;
///  3. magnetic field with a zero-field region;
///  4. CMS magnetic field - simple track;
///  5. CMS magnetic field - track with different path-marks.
///  6. Conceptual ILC detector, problematic track
///
/// \image html eve_track.png
/// \macro_code
///
/// \author Alja Mrak-Tadel


#include "TEveTrackPropagator.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TEveViewer.h"
#include "TSystem.h"
#include "TGLViewer.h"
#include "TMath.h"

#include "TEveViewer.h"
#include "TEvePointSet.h"

#include <iostream>

TEveTrackPropagator* g_prop = 0;

class GappedField : public TEveMagField
{
public:
   GappedField():TEveMagField(){}
   ~GappedField(){};

   virtual Double_t    GetMaxFieldMagD() { return 4; }

   virtual TEveVectorD GetFieldD(Double_t /*x*/, Double_t /*y*/, Double_t z) const
   {
      if (TMath::Abs(z) < 300) return TEveVectorD(0, 0, -4);
      if (TMath::Abs(z) < 600) return TEveVectorD(0, 0, 0);
      return TEveVectorD(0, 0, 4);
   }
};

//==============================================================================

class CmsMagField: public TEveMagField
{
   bool m_magnetIsOn;
   bool m_reverse;
   bool m_simpleModel;

public:
   CmsMagField():
      m_magnetIsOn(true),
      m_reverse(false),
      m_simpleModel(true){}

   virtual ~CmsMagField(){}

   void setMagnetState( bool state )
   {
      if (state != m_magnetIsOn)
      {
         if ( state )
            std::cout << "Magnet state is changed to ON" << std::endl;
         else
            std::cout << "Magnet state is changed to OFF" << std::endl;
      }
      m_magnetIsOn = state;
   }

   bool isMagnetOn() const               { return m_magnetIsOn;}
   void setReverseState(bool state)      { m_reverse = state; }
   bool isReverse() const                { return m_reverse;}
   void setSimpleModel(bool simpleModel) { m_simpleModel = simpleModel; }
   bool isSimpleModel() const            { return m_simpleModel;}

   virtual Double_t    GetMaxFieldMagD() const { return m_magnetIsOn ? 3.8 : 0.0; }

   virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const
   {
      double R = sqrt(x*x+y*y);
      double field = m_reverse ? -GetMaxFieldMagD() : GetMaxFieldMagD();
      //barrel
      if ( TMath::Abs(z) < 724 )
      {
         //inside solenoid
         if ( R < 300) return TEveVectorD(0,0,field);
         // outside solenoid
         if ( m_simpleModel ||
              ( R>461.0 && R<490.5 ) ||
              ( R>534.5 && R<597.5 ) ||
              ( R>637.0 && R<700.0 ) )
            return TEveVectorD(0,0,-field/3.8*1.2);

      } else {
         // endcaps
         if (m_simpleModel)
         {
            if ( R < 50 ) return TEveVectorD(0,0,field);
            if ( z > 0 )
               return TEveVectorD(x/R*field/3.8*2.0, y/R*field/3.8*2.0, 0);
            else
               return TEveVectorD(-x/R*field/3.8*2.0, -y/R*field/3.8*2.0, 0);
         }
         // proper model
         if ( ( TMath::Abs(z)>724 && TMath::Abs(z)<786  ) ||
              ( TMath::Abs(z)>850 && TMath::Abs(z)<910  ) ||
              ( TMath::Abs(z)>975 && TMath::Abs(z)<1003 ) )
         {
            if ( z > 0 )
               return TEveVectorD(x/R*field/3.8*2.0, y/R*field/3.8*2.0, 0);
            else
               return TEveVectorD(-x/R*field/3.8*2.0, -y/R*field/3.8*2.0, 0);
         }
      }
      return TEveVectorD(0,0,0);
   }
};


//==============================================================================
//==============================================================================

//______________________________________________________________________________
TEveTrack* make_track(TEveTrackPropagator* prop, Int_t sign)
{
  // Make track with given propagator.
  // Add to math-marks to test fit.

  auto rc = new TEveRecTrackD();
  rc->fV.Set(0.028558, -0.000918, 3.691919);
  rc->fP.Set(0.767095, -2.400006, -0.313103);
  rc->fSign = sign;

  auto track = new TEveTrack(rc, prop);
  track->SetName(Form("Charge %d", sign));
  // daughter 0
  auto pm1 = new TEvePathMarkD(TEvePathMarkD::kDaughter);
  pm1->fV.Set(1.479084, -4.370661, 3.119761);
  track->AddPathMark(*pm1);
  // daughter 1
  auto pm2 = new TEvePathMarkD(TEvePathMarkD::kDaughter);
  pm2->fV.Set(57.72345, -89.77011, -9.783746);
  track->AddPathMark(*pm2);

  return track;
}


void track(Int_t mode = 1, Bool_t isRungeKutta = kTRUE)
{

   gSystem->IgnoreSignal(kSigSegmentationViolation, true);
   TEveManager::Create();

   auto list = new TEveTrackList();
   auto prop = g_prop = list->GetPropagator();
   prop->SetFitDaughters(kFALSE);
   prop->SetMaxZ(1000);

   if (isRungeKutta)
   {
      prop->SetStepper(TEveTrackPropagator::kRungeKutta);
      list->SetName("RK Propagator");
   }
   else
   {
      list->SetName("Heix Propagator");
   }

   TEveTrack *track = 0;
   switch (mode)
   {
      case 0:
      {
         // B = 0 no difference between signed and charge particles
         prop->SetMagField(0.);
         list->SetElementName(Form("%s, zeroB", list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }

      case 1:
      {
         // constant B field, const angle
         prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., -3.8));
         list->SetElementName(Form("%s, constB", list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }

      case 2:
      {
         // variable B field, sign change at  R = 200 cm
         prop->SetMagFieldObj(new TEveMagFieldDuo(200, -4.4, 2));
         list->SetElementName(Form("%s, duoB", list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }

      case 3:
      {
         // gapped field
         prop->SetMagFieldObj(new GappedField());
         list->SetElementName(Form("%s, gappedB", list->GetElementName()));


         auto rc = new TEveRecTrackD();
         rc->fV.Set(0.028558, -0.000918, 3.691919);
         rc->fP.Set(0.767095, -0.400006, 2.313103);
         rc->fSign = 1;
         track = new TEveTrack(rc, prop);

         auto marker = new TEvePointSet(2);
         marker->SetElementName("B field break points");
         marker->SetPoint(0, 0., 0., 300.f);
         marker->SetPoint(1, 0., 0., 600.f);
         marker->SetMarkerColor(3);
         gEve->AddElement(marker);
         break;
      }

      case 4:
      {
         // Magnetic field of CMS I.
         auto mf = new CmsMagField;
         mf->setReverseState(true);

         prop->SetMagFieldObj(mf);
         prop->SetMaxR(1000);
         prop->SetMaxZ(1000);
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         list->SetElementName(Form("%s, CMS field", list->GetElementName()));


         auto rc = new TEveRecTrackD();
         rc->fV.Set(0.027667, 0.007919, 0.895964);
         rc->fP.Set(3.903134, 2.252232, -3.731366);
         rc->fSign = -1;
         track = new TEveTrack(rc, prop);

         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(3.576755e+00, 2.080579e+00, -2.507230e+00)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(8.440379e+01, 6.548286e+01, -8.788129e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.841321e+02, 3.915693e+02, -3.843072e+02)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.946167e+02, 4.793932e+02, -4.615060e+02)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDecay,
                  TEveVectorD(2.249656e+02, 5.835767e+02, -5.565275e+02)));

         track->SetRnrPoints(kTRUE);
         track->SetMarkerStyle(4);

         break;
      }

      case 5:
      {
         // Magnetic field of CMS I.
         auto mf = new CmsMagField;
         mf->setReverseState(true);
         mf->setSimpleModel(false);

         prop->SetMagFieldObj(mf);
         prop->SetMaxR(1000);
         prop->SetMaxZ(1000);
         prop->SetRnrReferences(kTRUE);
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         list->SetElementName(Form("%s, CMS field", list->GetElementName()));

         auto rc = new TEveRecTrackD();
         rc->fV.Set(-16.426592, 16.403185, -19.782692);
         rc->fP.Set(3.631100, 3.643450, 0.682254);
         rc->fSign = -1;
         track = new TEveTrack(rc, prop);

         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                  TEveVectorD(-1.642659e+01, 1.640318e+01, -1.978269e+01),
                  TEveVectorD(3.631100, 3.643450, 0.682254)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                  TEveVectorD(-1.859987e+00, 3.172243e+01, -1.697866e+01),
                  TEveVectorD(3.456056, 3.809894, 0.682254)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                  TEveVectorD(4.847579e+01, 9.871711e+01, -5.835719e+00),
                  TEveVectorD(2.711614, 4.409945, 0.687656)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.342045e+02, 4.203950e+02, 3.846268e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.483827e+02, 5.124750e+02, 5.064311e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.674676e+02, 6.167731e+02, 6.517403e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDecay,
                  TEveVectorD(1.884976e+02, 7.202000e+02, 7.919290e+01)));

         track->SetRnrPoints(kTRUE);
         track->SetMarkerStyle(4);

         break;
      }

      case 6:
      {
         // Problematic track from Druid
         prop->SetMagFieldObj(new TEveMagFieldDuo(350, -3.5, 2.0));
         prop->SetMaxR(1000);
         prop->SetMaxZ(1000);
         prop->SetRnrReferences(kTRUE);
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         list->SetElementName(Form("%s, Some ILC Detector field",
                                   list->GetElementName()));

         auto rc = new TEveRecTrackD();
         rc->fV.Set(57.1068, 31.2401, -7.07629);
         rc->fP.Set(4.82895, 2.35083, -0.611757);
         rc->fSign = 1;
         track = new TEveTrack(rc, prop);

         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.692235e+02, 7.047929e+01, -2.064785e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(5.806180e+02, 6.990633e+01, -6.450000e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDecay,
                  TEveVectorD(6.527213e+02, 1.473249e+02, -8.348498e+01)));

         track->SetRnrPoints(kTRUE);
         track->SetMarkerStyle(4);

         break;
      }
   };

   if (isRungeKutta)
      list->SetLineColor(kMagenta);
   else
      list->SetLineColor(kCyan);

   track->SetLineColor(list->GetLineColor());

   gEve->AddElement(list);
   list->AddElement(track);

   track->MakeTrack();

   TEveViewer *ev = gEve->GetDefaultViewer();
   TGLViewer  *gv = ev->GetGLViewer();
   gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

   gEve->Redraw3D(kTRUE);
   gSystem->ProcessEvents();

   gv->CurrentCamera().RotateRad(-0.5, 1.4);
   gv->RequestDraw();
}
