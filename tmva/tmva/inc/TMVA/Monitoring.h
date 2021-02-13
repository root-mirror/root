#pragma once

#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

// FIXME: We should not rely on std::cout but on the ROOT printing facilities or
// MsgLogger!
#include <iostream> // for std::cout
#include <map>
#include <string>

namespace TMVA
{
   class  Monitoring
   {

   public:
      /* Monitoring (int argc, char* /\*argv[]*\/) */
      /* { */
      /* }     */

      Monitoring ()
         : fCanvas (NULL)
         {
         }    

      ~Monitoring () 
         { 
            delete fCanvas; 
            //            delete fApplication;
         }

      void Start ()
      {
         /*             std::cout << "start monitoring" << std::endl; */
         /*             std::cout << "  new tapp " <<  std::endl; */
         /*             fApplication = new TApplication ("TMVA Monitoring", 0, 0); */
         /*             std::cout << "  set return from run" << std::endl; */
         /* //            fApplication->SetReturnFromRun (true); */

         std::cout << "  new tcanvas" << std::endl;
         fCanvas = new TCanvas ("TMVA Monitoring", "Monitoring", 1000, 500);
         std::cout << "  draw" << std::endl;
         fCanvas->Draw ();
         std::cout << "  update" << std::endl;
         GetCanvas ()->Update();
         std::cout << "  process events" << std::endl;
         gSystem->ProcessEvents(); //canvas can be edited during the loop
         std::cout << "  run app" << std::endl;
         //            fApplication->Run ();
         std::cout << "  run app executed" << std::endl;

         gStyle->SetOptStat (0);
      }


      void ProcessEvents ()
      {
         GetCanvas ()->Modified();
         GetCanvas ()->Update();
         gSystem->ProcessEvents(); //canvas can be edited during the loop
      }

      TCanvas* GetCanvas () { return fCanvas; }

      void pads (int numPads);
      void create (std::string histoName, int bins, double min, double max);
      void create (std::string histoName, int bins, double min, double max, int bins2, double min2, double max2);
      void addPoint (std::string histoName, double x);
      void addPoint (std::string histoName, double x, double y);
      void plot (std::string histoName, std::string options = "L", int pad = 0, EColor color = kBlue);
      void clear (std::string histoName);
      bool exists (std::string histoName);
      bool exists (TH1F* dummy, std::string histoName);
      bool exists (TH2F* dummy, std::string histoName);

   protected:

      TH1F* getHistogram (const TH1F* dummy, std::string histoName, int bins = 0, double min = 0, double max = 0);
      TH2F* getHistogram (const TH2F* dummy, std::string histoName, int bins = 0, double min = 0, double max = 0, int bins2 = 0, double min2 = 0, double max2 = 0);


   private:
      TCanvas* fCanvas;

      //        TApplication* fApplication;


      std::map<std::string, TH1F*> m_histos1D;
      std::map<std::string, TH2F*> m_histos2D;
   };



   inline bool Monitoring::exists (TH1F* /*dummy*/, std::string histoName)
      {
         auto it = m_histos1D.find (histoName);
         if (it != m_histos1D.end ())
            return true;
         return false;
      }

   inline bool Monitoring::exists (TH2F* /*dummy*/, std::string histoName)
      {
         auto it2 = m_histos2D.find (histoName);
         if (it2 != m_histos2D.end ())
            return true;
         return false;
      }


   inline bool Monitoring::exists (std::string histoName)
      {
         TH1F* dummy1D (NULL);
         TH2F* dummy2D (NULL);
         return exists (dummy1D, histoName) || exists (dummy2D, histoName);
      }

   inline void Monitoring::pads (int numPads)
   {
      TCanvas* canvas = GetCanvas ();
      canvas->Clear ();
      std::cout << "divide canvas " << canvas << " into " << numPads << "numPads" << std::endl;
      GetCanvas ()->DivideSquare (numPads);
   }


   inline void Monitoring::create (std::string histoName, int bins, double min, double max)
   {
      TH1F* dummy (NULL);
      getHistogram (dummy, histoName, bins, min, max);
   }

   inline void Monitoring::create (std::string histoName, int bins, double min, double max, int bins2, double min2, double max2)
   {
      TH2F* dummy (NULL);
      getHistogram (dummy, histoName, bins, min, max, bins2, min2, max2);
   }



   inline TH1F* Monitoring::getHistogram (const TH1F* /*dummy*/, std::string histoName, int bins, double min, double max)
   {
      auto it = m_histos1D.find (histoName);
      if (it != m_histos1D.end ())
         return it->second;
      std::cout << "new 1D histogram " << histoName << std::endl;
      TH1F* histogram = m_histos1D.insert (std::make_pair (histoName, new TH1F (histoName.c_str (), histoName.c_str (), bins, min, max))).first->second;
      //    int numPads = m_histos1D.size () + m_histos2D.size ();
      return histogram;
   }

   inline TH2F* Monitoring::getHistogram (const TH2F* /*dummy*/, std::string histoName, int bins, double min, double max, int bins2, double min2, double max2)
   {
      // 2D histogram
      auto it = m_histos2D.find (histoName);
      if (it != m_histos2D.end ())
         return it->second;
      std::cout << "new 2D histogram " << histoName << std::endl;
      TH2F* histogram = m_histos2D.insert (std::make_pair (histoName, new TH2F (histoName.c_str (), histoName.c_str (), bins, min, max, bins2, min2, max2))).first->second;
      //    int numPads = m_histos1D.size () + m_histos2D.size ();
      return histogram;
   }

   inline void Monitoring::addPoint (std::string histoName, double x)
   {
      TH1F* dummy (NULL);
      TH1F* hist = getHistogram (dummy, histoName, 100, 0, 1);
      hist->Fill (x);
   }

   inline void Monitoring::addPoint (std::string histoName, double x, double y)
   {
      TH2F* dummy (NULL);
      TH2F* hist = getHistogram (dummy, histoName, 100, 0, 1, 100, 0, 1);
      hist->Fill (x, y);
   }

   inline void Monitoring::clear (std::string histoName)
   {
      //    std::cout << "clear histo " << histoName << std::endl;
      if (!exists (histoName))
         return;

      //    std::cout << "clear histo which exists " << histoName << std::endl;
      TH1F* hist1D (NULL);
      TH2F* hist2D (NULL);
      if (exists (hist1D, histoName))
         {
            hist1D = getHistogram (hist1D, histoName, 100, 0,1);
            hist1D->Reset ();
            return;
         }

      if (exists (hist2D, histoName))
         {
            hist2D = getHistogram (hist2D, histoName, 100, 0,1,100,0,1);
            hist2D->Reset ();
         }
   }


   inline void Monitoring::plot (std::string histoName, std::string options, int pad, EColor color)
   {
      TCanvas* canvas = GetCanvas ();
      canvas->cd (pad);
      auto it1D = m_histos1D.find (histoName);
      if (it1D != m_histos1D.end ())
         {
            TH1F* dummy (NULL);
            TH1F* histogram = getHistogram (dummy, histoName);
            //        histogram->SetBit (TH1::kCanRebin);
            histogram->SetLineColor (color);
            histogram->SetMarkerColor (color);
            //        std::cout << "draw " << histoName << " 1D on canvas " << canvas << " on pad " << pad << " with  options " << options << std::endl;
            histogram->Draw (options.c_str ());
            canvas->Modified ();
            canvas->Update ();
            return;
         }
      auto it2D = m_histos2D.find (histoName);
      if (it2D != m_histos2D.end ())
         {
            TH2F* dummy (NULL);
            TH2F* histogram = getHistogram (dummy, histoName);
            //        histogram->SetBit (TH1::kCanRebin);
            histogram->SetLineColor (color);
            histogram->SetMarkerColor (color);
            //        std::cout << "draw " << histoName << " 2D on canvas " << canvas << " on pad " << pad << " with  options " << options << std::endl;
            histogram->Draw (options.c_str ());
            canvas->Modified ();
            canvas->Update ();
         }
   }



} // namespace TMVA
