/// \file
/// \ingroup tutorial_dataframe
/// Header file with functions needed to execute the Python version
/// of the NanoAOD Higgs tutorial. The header is declared to the
/// ROOT C++ interpreter prior to the start of the analysis via the
/// `ROOT.gInterpreter.Declare()` function.
///
/// \date June 2019
/// \author Stefan Wunsch (KIT, CERN), Vincenzo Eduardo Padulano (UniMiB, CERN)

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
const auto z_mass = 91.2;

// Reconstruct two Z candidates from four leptons of the same kind
RVec<RVec<size_t>> reco_zz_to_4l(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_i charge)
{
   RVec<RVec<size_t>> idx(2);
   idx[0].reserve(2); idx[1].reserve(2);
   // Find first lepton pair with invariant mass closest to Z mass
   auto idx_cmb = Combinations(pt, 2);
   auto best_mass = -1;
   size_t best_i1 = 0; size_t best_i2 = 0;
   for (size_t i = 0; i < idx_cmb[0].size(); i++) {
      const auto i1 = idx_cmb[0][i];
      const auto i2 = idx_cmb[1][i];
      if (charge[i1] != charge[i2]) {
         TLorentzVector p1, p2;
         p1.SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], mass[i1]);
         p2.SetPtEtaPhiM(pt[i2], eta[i2], phi[i2], mass[i2]);
         const auto this_mass = (p1 + p2).M();
         if (std::abs(z_mass - this_mass) < std::abs(z_mass - best_mass)) {
            best_mass = this_mass;
            best_i1 = i1;
            best_i2 = i2;
         }
      }
   }
   idx[0].emplace_back(best_i1);
   idx[0].emplace_back(best_i2);
   // Reconstruct second Z from remaining lepton pair
   for (size_t i = 0; i < 4; i++) {
      if (i != best_i1 && i != best_i2) {
         idx[1].emplace_back(i);
      }
   }
   // Return indices of the pairs building two Z bosons
   return idx;
}

// Compute mass of two Z candidates from two electrons and two muons and sort ascending in distance to Z mass
RVec<float> compute_z_masses_2el2mu(rvec_f el_pt, rvec_f el_eta, rvec_f el_phi, rvec_f el_mass, rvec_f mu_pt,
                                  rvec_f mu_eta, rvec_f mu_phi, rvec_f mu_mass)
{
   TLorentzVector p1, p2, p3, p4;
   p1.SetPtEtaPhiM(mu_pt[0], mu_eta[0], mu_phi[0], mu_mass[0]);
   p2.SetPtEtaPhiM(mu_pt[1], mu_eta[1], mu_phi[1], mu_mass[1]);
   p3.SetPtEtaPhiM(el_pt[0], el_eta[0], el_phi[0], el_mass[0]);
   p4.SetPtEtaPhiM(el_pt[1], el_eta[1], el_phi[1], el_mass[1]);
   auto mu_z = (p1 + p2).M();
   auto el_z = (p3 + p4).M();
   RVec<float> z_masses(2);
   if (std::abs(mu_z - z_mass) < std::abs(el_z - z_mass)) {
      z_masses[0] = mu_z;
      z_masses[1] = el_z;
   } else {
      z_masses[0] = el_z;
      z_masses[1] = mu_z;
   }
   return z_masses;
}

// Compute Higgs mass from two electrons and two muons
float compute_higgs_mass_2el2mu(rvec_f el_pt, rvec_f el_eta, rvec_f el_phi, rvec_f el_mass, rvec_f mu_pt, rvec_f mu_eta,
                                rvec_f mu_phi, rvec_f mu_mass)
{
   TLorentzVector p1, p2, p3, p4;
   p1.SetPtEtaPhiM(mu_pt[0], mu_eta[0], mu_phi[0], mu_mass[0]);
   p2.SetPtEtaPhiM(mu_pt[1], mu_eta[1], mu_phi[1], mu_mass[1]);
   p3.SetPtEtaPhiM(el_pt[0], el_eta[0], el_phi[0], el_mass[0]);
   p4.SetPtEtaPhiM(el_pt[1], el_eta[1], el_phi[1], el_mass[1]);
   return (p1 + p2 + p3 + p4).M();
}

// Compute Z masses from four leptons of the same kind and sort ascending in distance to Z mass
RVec<float> compute_z_masses_4l(const RVec<RVec<size_t>> &idx, rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass)
{
   RVec<float> z_masses(2);
   for (size_t i = 0; i < 2; i++) {
      TLorentzVector p1, p2;
      const auto i1 = idx[i][0]; const auto i2 = idx[i][1];
      p1.SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], mass[i1]);
      p2.SetPtEtaPhiM(pt[i2], eta[i2], phi[i2], mass[i2]);
      z_masses[i] = (p1 + p2).M();
   }
   if (std::abs(z_masses[0] - z_mass) < std::abs(z_masses[1] - z_mass)) {
      return z_masses;
   } else {
      return Reverse(z_masses);
   }
}

// Compute mass of Higgs from four leptons of the same kind
float compute_higgs_mass_4l(const RVec<RVec<size_t>> &idx, rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass)
{
   TLorentzVector p1, p2, p3, p4;
   const auto i1 = idx[0][0]; const auto i2 = idx[0][1];
   const auto i3 = idx[1][0]; const auto i4 = idx[1][1];
   p1.SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], mass[i1]);
   p2.SetPtEtaPhiM(pt[i2], eta[i2], phi[i2], mass[i2]);
   p3.SetPtEtaPhiM(pt[i3], eta[i3], phi[i3], mass[i3]);
   p4.SetPtEtaPhiM(pt[i4], eta[i4], phi[i4], mass[i4]);
   return (p1 + p2 + p3 + p4).M();
}

bool filter_z_dr(const RVec<RVec<size_t>> &idx, rvec_f eta, rvec_f phi)
{
   for (size_t i = 0; i < 2; i++) {
      const auto i1 = idx[i][0];
      const auto i2 = idx[i][1];
      const auto dr = sqrt(pow(eta[i1] - eta[i2], 2) + pow(phi[i1] - phi[i2], 2));
      if (dr < 0.02) {
         return false;
      }
   }
   return true;
};

bool pt_cuts(rvec_f mu_pt, rvec_f el_pt)
{
   auto mu_pt_sorted = Reverse(Sort(mu_pt));
   if (mu_pt_sorted[0] > 20 && mu_pt_sorted[1] > 10) {
      return true;
   }
   auto el_pt_sorted = Reverse(Sort(el_pt));
   if (el_pt_sorted[0] > 20 && el_pt_sorted[1] > 10) {
      return true;
   }
   return false;
}

bool dr_cuts(rvec_f mu_eta, rvec_f mu_phi, rvec_f el_eta, rvec_f el_phi)
{
   auto mu_dr = sqrt(pow(mu_eta[0] - mu_eta[1], 2) + pow(mu_phi[0] - mu_phi[1], 2));
   auto el_dr = sqrt(pow(el_eta[0] - el_eta[1], 2) + pow(el_phi[0] - el_phi[1], 2));
   if (mu_dr < 0.02 || el_dr < 0.02) {
      return false;
   }
   return true;
}

float add_weight_higgs_sig_4mu()
{
   const auto luminosity = 11580.0;            // Integrated luminosity of the data samples
   const auto xsec_SMHiggsToZZTo4L = 0.0065;   // H->4l: Standard Model cross-section
   const auto nevt_SMHiggsToZZTo4L = 299973.0; // H->4l: Number of simulated events

   return luminosity * xsec_SMHiggsToZZTo4L / nevt_SMHiggsToZZTo4L;
}

float add_weight_higgs_bkg_4mu()
{
   const auto luminosity = 11580.0;            // Integrated luminosity of the data samples
   const auto scale_ZZTo4l = 1.386;     // ZZ->4mu: Scale factor for ZZ to four leptons
   const auto xsec_ZZTo4mu = 0.077;     // ZZ->4mu: Standard Model cross-section
   const auto nevt_ZZTo4mu = 1499064.0; // ZZ->4mu: Number of simulated events

   return luminosity * xsec_ZZTo4mu * scale_ZZTo4l / nevt_ZZTo4mu;
}

float add_weight_higgs_sig_4el()
{
   const auto luminosity = 11580.0;            // Integrated luminosity of the data samples
   const auto xsec_SMHiggsToZZTo4L = 0.0065;   // H->4l: Standard Model cross-section
   const auto nevt_SMHiggsToZZTo4L = 299973.0; // H->4l: Number of simulated events

   return luminosity * xsec_SMHiggsToZZTo4L / nevt_SMHiggsToZZTo4L;
}

float add_weight_higgs_bkg_4el()
{
   const auto luminosity = 11580.0;            // Integrated luminosity of the data samples
   const auto xsec_ZZTo4el = 0.077;     // ZZ->4mu: Standard Model cross-section
   const auto scale_ZZTo4l = 1.386;     // ZZ->4mu: Scale factor for ZZ to four leptons
   const auto nevt_ZZTo4el = 1499093.0;    // ZZ->4el: Number of simulated events

   return luminosity * xsec_ZZTo4el * scale_ZZTo4l / nevt_ZZTo4el;
}

float add_weight_higgs_sig_2el2mu()
{
   const auto luminosity = 11580.0;            // Integrated luminosity of the data samples
   const auto xsec_SMHiggsToZZTo4L = 0.0065;   // H->4l: Standard Model cross-section
   const auto nevt_SMHiggsToZZTo4L = 299973.0; // H->4l: Number of simulated events

   return luminosity * xsec_SMHiggsToZZTo4L / nevt_SMHiggsToZZTo4L;
}

float add_weight_higgs_bkg_2el2mu()
{
   const auto luminosity = 11580.0;            // Integrated luminosity of the data samples
   const auto xsec_ZZTo2el2mu = 0.18;      // ZZ->2el2mu: Standard Model cross-section
   const auto scale_ZZTo4l = 1.386;     // ZZ->4mu: Scale factor for ZZ to four leptons
   const auto nevt_ZZTo2el2mu = 1497445.0; // ZZ->2el2mu: Number of simulated events

   return luminosity * xsec_ZZTo2el2mu * scale_ZZTo4l / nevt_ZZTo2el2mu;
}

