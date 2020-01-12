#ifndef ntpl005_convert_hxx
#define ntpl005_convert_hxx

#include "h1event.hxx"

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// Import classes from experimental namespace for the time being
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;

const std::string fileLocation = "http://root.cern.ch/files/h1/";

void convert(const std::string fileName)
{
   const std::string fullPath = fileLocation + fileName;

   // Opens the .root file containing a TTree. By using a unique_ptr the statement f->Close() is not needed.
   std::unique_ptr<TFile> f(TFile::Open(fullPath.c_str()));
   assert(f.get() && ! f->IsZombie());

   // h42 is the name of the TTree
   auto tree = f->Get<TTree>("h42");
   TTreeReader reader("h42", f.get());

   TTreeReaderValue<std::int32_t>   nrun(reader, "nrun"); // 0
   TTreeReaderValue<std::int32_t>   nevent(reader, "nevent"); // 1
   TTreeReaderValue<std::int32_t>   nentry(reader, "nentry"); // 2
   TTreeReaderArray<unsigned char>  trelem(reader, "trelem"); // 3
   TTreeReaderArray<unsigned char>  subtr(reader, "subtr"); // 4
   TTreeReaderArray<unsigned char>  rawtr(reader, "rawtr"); // 5
   TTreeReaderArray<unsigned char>  L4subtr(reader, "L4subtr"); // 6
   TTreeReaderArray<unsigned char>  L5class(reader, "L5class"); // 7
   TTreeReaderValue<float>          E33(reader, "E33"); // 8
   TTreeReaderValue<float>          de33(reader, "de33"); // 9
   TTreeReaderValue<float>          x33(reader, "x33"); // 10
   TTreeReaderValue<float>          dx33(reader, "dx33"); // 11
   TTreeReaderValue<float>          y33(reader, "y33"); // 12
   TTreeReaderValue<float>          dy33(reader, "dy33"); // 13
   TTreeReaderValue<float>          E44(reader, "E44"); // 14
   TTreeReaderValue<float>          de44(reader, "de44"); // 15
   TTreeReaderValue<float>          x44(reader, "x44"); // 16
   TTreeReaderValue<float>          dx44(reader, "dx44"); // 17
   TTreeReaderValue<float>          y44(reader, "y44"); // 18
   TTreeReaderValue<float>          dy44(reader, "dy44"); // 19
   TTreeReaderValue<float>          Ept(reader, "Ept"); // 20
   TTreeReaderValue<float>          dept(reader, "dept"); // 21
   TTreeReaderValue<float>          xpt(reader, "xpt"); // 22
   TTreeReaderValue<float>          dxpt(reader, "dxpt"); // 23
   TTreeReaderValue<float>          ypt(reader, "ypt"); // 24
   TTreeReaderValue<float>          dypt(reader, "dypt"); // 25
   TTreeReaderArray<float>          pelec(reader, "pelec"); // 26, size = 4
   TTreeReaderValue<std::int32_t>   flagelec(reader, "flagelec"); // 27
   TTreeReaderValue<float>          xeelec(reader, "xeelec"); // 28
   TTreeReaderValue<float>          yeelec(reader, "yeelec"); // 29
   TTreeReaderValue<float>          Q2eelec(reader, "Q2eelec"); // 30
   TTreeReaderValue<std::int32_t>   nelec(reader, "nelec"); // 31
   
   TTreeReaderArray<float>          Eelec(reader, "Eelec"); // 32, size = nelec
   TTreeReaderArray<float>          thetelec(reader, "thetelec"); // 33, size = nelec
   TTreeReaderArray<float>          phielec(reader, "phielec"); // 34, size = nelec
   TTreeReaderArray<float>          xelec(reader, "xelec"); // 35, size = nelec
   TTreeReaderArray<float>          Q2elec(reader, "Q2elec"); // 36, size = nelec
   TTreeReaderArray<float>          xsigma(reader, "xsigma"); // 37, size = nelec
   TTreeReaderArray<float>          Q2sigma(reader, "Q2sigma"); // 38, size = nelec
   
   TTreeReaderArray<float>          sumc(reader, "sumc"); // 39, size = 4
   TTreeReaderValue<float>          sumetc(reader, "sumetc"); // 40
   TTreeReaderValue<float>          yjbc(reader, "yjbc"); // 41
   TTreeReaderValue<float>          Q2jbc(reader, "Q2jbc"); // 42
   TTreeReaderArray<float>          sumct(reader, "sumct"); // 43, size = 4
   TTreeReaderValue<float>          sumetct(reader, "sumetct"); // 44
   TTreeReaderValue<float>          yjbct(reader, "yjbct"); // 45
   TTreeReaderValue<float>          Q2jbct(reader, "Q2jbct"); // 46
   TTreeReaderValue<float>          Ebeamel(reader, "Ebeamel"); // 47
   TTreeReaderValue<float>          Ebeampr(reader, "Ebeampr"); // 48
   
   TTreeReaderArray<float>          pvtx_d(reader, "xelec"); // 49, size = 3
   TTreeReaderArray<float>          cpvtx_d(reader, "cpvtx_d"); // 50, size = 6
   TTreeReaderArray<float>          pvtx_t(reader, "pvtx_t"); // 51, size = 3
   TTreeReaderArray<float>          cpvtx_t(reader, "cpvtx_t"); // 52, size = 6
   
   TTreeReaderValue<std::int32_t>   ntrkxy_t(reader, "ntrkxy_t"); // 53
   TTreeReaderValue<float>          prbxy_t(reader, "prbxy_t"); // 54
   TTreeReaderValue<std::int32_t>   ntrkz_t(reader, "ntrkz_t"); // 55
   TTreeReaderValue<float>          prbz_t(reader, "prbz_t"); // 56
   TTreeReaderValue<std::int32_t>   nds(reader, "nds"); // 57
   TTreeReaderValue<std::int32_t>   rankds(reader, "rankds"); // 58
   TTreeReaderValue<std::int32_t>   qds(reader, "qds"); // 59
   TTreeReaderArray<float>          pds_d(reader, "pds_d"); // 60, size = 4
   TTreeReaderValue<float>          ptds_d(reader, "ptds_d"); // 61
   TTreeReaderValue<float>          etads_d(reader, "etads_d"); // 62
   TTreeReaderValue<float>          dm_d(reader, "dm_d"); // 63
   TTreeReaderValue<float>          ddm_d(reader, "ddm_d"); // 64
   TTreeReaderArray<float>          pds_t(reader, "pds_t"); // 65, size = 4
   TTreeReaderValue<float>          dm_t(reader, "dm_t"); // 66
   TTreeReaderValue<float>          ddm_t(reader, "ddm_t"); // 67
   
   TTreeReaderValue<std::int32_t>   ik(reader, "ik"); // 68
   TTreeReaderValue<std::int32_t>   ipi(reader, "ipi"); // 69
   TTreeReaderValue<std::int32_t>   ipis(reader, "ipis"); // 70
   TTreeReaderArray<float>          pd0_d(reader, "pd0_d"); // 71, size = 4
   TTreeReaderValue<float>          ptd0_d(reader, "ptd0_d"); // 72
   TTreeReaderValue<float>          etad0_d(reader, "etad0_d"); // 73
   TTreeReaderValue<float>          md0_d(reader, "md0_d"); // 74
   TTreeReaderValue<float>          dmd0_d(reader, "dmd0_d"); // 75
   TTreeReaderArray<float>          pd0_t(reader, "pd0_t"); // 76
   TTreeReaderValue<float>          md0_t(reader, "md0_t"); // 77
   TTreeReaderValue<float>          dmd0_t(reader, "dmd0_t"); // 78
   TTreeReaderArray<float>          pk_r(reader, "pk_r"); // 79, size = 4
   TTreeReaderArray<float>          ppi_r(reader, "ppi_r"); // 80, size = 4
   TTreeReaderArray<float>          pd0_r(reader, "pd0_r"); // 81, size = 4
   
   TTreeReaderValue<float>          md0_r(reader, "md0_r"); // 82
   TTreeReaderArray<float>          Vtxd0_r(reader, "Vtxd0_r"); // 83, size = 3
   TTreeReaderArray<float>          cvtxd0_r(reader, "cvtxd0_r"); // 84, size = 6
   TTreeReaderValue<float>          dxy_r(reader, "dxy_r"); // 85
   TTreeReaderValue<float>          dz_r(reader, "dz_r"); // 86
   TTreeReaderValue<float>          psi_r(reader, "psi_r"); // 87
   TTreeReaderValue<float>          rd0_d(reader, "rd0_d"); // 88
   TTreeReaderValue<float>          drd0_d(reader, "drd0_d"); // 89
   TTreeReaderValue<float>          rpd0_d(reader, "rpd0_d"); // 90
   TTreeReaderValue<float>          drpd0_d(reader, "drpd0_d"); // 91
   TTreeReaderValue<float>          rd0_t(reader, "rd0_t"); // 92
   TTreeReaderValue<float>          drd0_t(reader, "drd0_t"); // 93
   TTreeReaderValue<float>          rpd0_t(reader, "rpd0_t"); // 94
   TTreeReaderValue<float>          drpd0_t(reader, "drpd0_t"); // 95
   TTreeReaderValue<float>          rd0_dt(reader, "rd0_dt"); // 96
   TTreeReaderValue<float>          drd0_dt(reader, "drd0_dt"); // 97
   TTreeReaderValue<float>          prbr_dt(reader, "prbr_dt"); // 98
   TTreeReaderValue<float>          prbz_dt(reader, "prbz_dt"); // 99
   TTreeReaderValue<float>          rd0_tt(reader, "rd0_tt"); // 100
   TTreeReaderValue<float>          drd0_tt(reader, "drd0_tt"); // 101
   TTreeReaderValue<float>          prbr_tt(reader, "prbr_tt"); // 102
   TTreeReaderValue<float>          prbz_tt(reader, "prbz_tt"); // 103
   TTreeReaderValue<std::int32_t>   ijetd0(reader, "ijetd0"); // 104
   TTreeReaderValue<float>          ptr3d0_j(reader, "ptr3d0_j"); // 105
   TTreeReaderValue<float>          ptr2d0_j(reader, "ptr2d0_j"); // 106
   TTreeReaderValue<float>          ptr3d0_3(reader, "ptr3d0_3"); // 107
   TTreeReaderValue<float>          ptr2d0_3(reader, "ptr2d0_3"); // 108
   TTreeReaderValue<float>          ptr2d0_2(reader, "ptr2d0_2"); // 109
   TTreeReaderValue<float>          Mimpds_r(reader, "Mimpds_r"); // 110
   TTreeReaderValue<float>          Mimpbk_r(reader, "Mimpbk_r"); // 111
   
   TTreeReaderValue<std::int32_t>   ntracks(reader, "ntracks"); // 112
   TTreeReaderArray<float>          pt(reader, "pt"); // 113
   TTreeReaderArray<float>          kappa(reader, "kappa"); // 114
   TTreeReaderArray<float>          phi(reader, "phi"); // 115
   TTreeReaderArray<float>          theta(reader, "theta"); // 116
   TTreeReaderArray<float>          dca(reader, "dca"); // 117
   TTreeReaderArray<float>          z0(reader, "z0"); // 118
   //TTreeReaderArray<float>        covar(reader, "covar"); // 119
   TTreeReaderArray<std::int32_t>   nhitrp(reader, "nhitrp"); // 120
   TTreeReaderArray<float>          prbrp(reader, "prbrp"); // 121
   TTreeReaderArray<std::int32_t>   nhitz(reader, "nhitz"); // 122
   TTreeReaderArray<float>          prbz(reader, "prbz"); // 123
   TTreeReaderArray<float>          rstart(reader, "rstart"); // 124
   TTreeReaderArray<float>          rend(reader, "rend"); // 125
   TTreeReaderArray<float>          lhk(reader, "lhk"); // 126
   TTreeReaderArray<float>          lhpi(reader, "lhpi"); // 127
   TTreeReaderArray<float>          nlhk(reader, "nlhk"); // 128
   TTreeReaderArray<float>          nlhpi(reader, "nlhpi"); // 129
   TTreeReaderArray<float>          dca_d(reader, "dca_d"); // 130
   TTreeReaderArray<float>          ddca_d(reader, "ddca_d"); // 131
   TTreeReaderArray<float>          dca_t(reader, "dca_t"); // 132
   TTreeReaderArray<float>          ddca_t(reader, "ddca_t"); // 133
   TTreeReaderArray<std::int32_t>   muqual(reader, "muqual"); // 134
   
   TTreeReaderValue<std::int32_t>   imu(reader, "imu"); // 135
   TTreeReaderValue<std::int32_t>   imufe(reader, "imufe"); // 136
   
   TTreeReaderValue<std::int32_t>   njets(reader, "njets"); // 137
   TTreeReaderArray<float>          E_j(reader, "E_j"); // 138
   TTreeReaderArray<float>          pt_j(reader, "pt_j"); // 139
   TTreeReaderArray<float>          theta_j(reader, "theta_j"); // 140
   TTreeReaderArray<float>          eta_j(reader, "eta_j"); // 141
   TTreeReaderArray<float>          phi_j(reader, "phi_j"); // 142
   TTreeReaderArray<float>          m_j(reader, "m_j"); // 143
   
   TTreeReaderValue<float>          thrust(reader, "thrust"); // 144
   TTreeReaderArray<float>          pthrust(reader, "pthrust"); // 145, size = 4
   TTreeReaderValue<float>          thrust2(reader, "thrust2"); // 146
   TTreeReaderArray<float>          pthrust2(reader, "pthrust2"); // 147, size = 4
   TTreeReaderValue<float>          spher(reader, "spher"); // 148
   TTreeReaderValue<float>          aplan(reader, "aplan"); // 149
   TTreeReaderValue<float>          plan(reader, "plan"); // 150
   TTreeReaderArray<float>          nnout(reader, "nnout"); // 151
   
   

   // Loads the dictionary for the H1event struct. Without the dictionary MakeField<H1Event> would fail.
   //gSystem->Load("./libH1event.so");
   {
      // Create a ntuple model with a single field
      auto model = RNTupleModel::Create();
      auto ev = model->MakeField<H1Event>("event");

      // h42 refers to the name of the ntuple.
      auto ntuple = RNTupleWriter::Recreate(std::move(model), "h42", fileName);
      int count = 0;

      // Fills the ntuple with entries from the TTree.
      while(reader.Next()) {
         std::array<bool, 192> trelemNTuple;
         for (int i = 0; i < 192; ++i) {
            trelemNTuple.at(i) = trelem[i];
         }
         
         std::array<bool, 128> subtrNTuple;
         for (int i = 0; i < 128; ++i) {
            subtrNTuple.at(i) = subtr[i];
         }
         
         std::array<bool, 128> rawtrNTuple;
         for (int i = 0; i < 128; ++i) {
            rawtrNTuple.at(i) = rawtr[i];
         }
         
         std::array<bool, 128> L4subtrNTuple;
         for (int i = 0; i < 128; ++i) {
            L4subtrNTuple.at(i) = L4subtr[i];
         }
         
         std::array<bool, 32> L5classNTuple;
         for (int i = 0; i < 32; ++i) {
            L5classNTuple.at(i) = L5class[i];
         }
         
         std::array<float, 4> pelecNTuple;
         for (int i = 0; i < 4; ++i) {
            pelecNTuple.at(i) = pelec[i];
         }
         
         std::vector<H1Event::Electron> nelecNTuple(*nelec);
         for (int i = 0; i < *nelec; ++i) {
            nelecNTuple.at(i) = H1Event::Electron{Eelec[i], thetelec[i], phielec[i], xelec[i], Q2elec[i], xsigma[i], Q2sigma[i]};
         }
         
         std::array<float, 4> sumcNTuple;
         for (int i = 0; i < 4; ++i) {
            sumcNTuple.at(i) = sumc[i];
         }
         
         std::array<float, 4> sumctNTuple;
         for (int i = 0; i < 4; ++i) {
            sumctNTuple.at(i) = sumct[i];
         }
         
         std::array<float, 3> pvtx_dNTuple;
         for (int i = 0; i < 3; ++i) {
            pvtx_dNTuple.at(i) = pvtx_d[i];
         }
         
         std::array<float, 6> cpvtx_dNTuple;
         for (int i = 0; i < 6; ++i) {
            cpvtx_dNTuple.at(i) = cpvtx_d[i];
         }
         
         std::array<float, 3> pvtx_tNTuple;
         for (int i = 0; i < 3; ++i) {
            pvtx_tNTuple.at(i) = pvtx_t[i];
         }
         
         std::array<float, 6> cpvtx_tNTuple;
         for (int i = 0; i < 6; ++i) {
            cpvtx_tNTuple.at(i) = cpvtx_t[i];
         }
         
         std::array<float, 4> pds_dNTuple;
         for (int i = 0; i < 4; ++i) {
            pds_dNTuple.at(i) = pds_d[i];
         }
         
         std::array<float, 4> pds_tNTuple;
         for (int i = 0; i < 4; ++i) {
            pds_tNTuple.at(i) = pds_t[i];
         }
         
         std::array<float, 4> pd0_dNTuple;
         for (int i = 0; i < 4; ++i) {
            pd0_dNTuple.at(i) = pd0_d[i];
         }
         
         std::array<float, 4> pd0_tNTuple;
         for (int i = 0; i < 4; ++i) {
            pd0_tNTuple.at(i) = pd0_t[i];
         }
         
         std::array<float, 4> pk_rNTuple;
         for (int i = 0; i < 4; ++i) {
            pk_rNTuple.at(i) = pk_r[i];
         }
         
         std::array<float, 4> ppi_rNTuple;
         for (int i = 0; i < 4; ++i) {
            ppi_rNTuple.at(i) = ppi_r[i];
         }
         
         std::array<float, 4> pd0_rNTuple;
         for (int i = 0; i < 4; ++i) {
            pd0_rNTuple.at(i) = pd0_r[i];
         }
         
         std::array<float, 3> Vtxd0_rNTuple;
         for (int i = 0; i < 3; ++i) {
            Vtxd0_rNTuple.at(i) = Vtxd0_r[i];
         }
         
         std::array<float, 6> cvtxd0_rNTuple;
         for (int i = 0; i < 6; ++i) {
            cvtxd0_rNTuple.at(i) = cvtxd0_r[i];
         }
         
         
         static float covar[200][15];
         tree->SetBranchAddress("covar", covar);
         tree->GetEntry(count++);
         std::vector<std::array<float, 15>> covarVec;
         for(int i = *ntracks; i > 0; --i) {
            std::array<float, 15> ar;
            for(int j = 0; j < 15; ++j) {
               ar.at(j) = covar[i][j];
            }
            covarVec.emplace_back(ar);
         }
         
         std::vector<H1Event::Track> ntrackNTuple(*ntracks);
         for (int i = 0; i < *ntracks; ++i) {
            ntrackNTuple.at(i) = H1Event::Track{ pt[i], kappa[i], phi[i], theta[i], dca[i], z0[i], covarVec.at(i), nhitrp[i], prbrp[i], nhitz[i], prbz[i], rstart[i], rend[i], lhk[i], lhpi[i], nlhk[i], nlhpi[i], dca_d[i], ddca_d[i], dca_t[i], ddca_t[i], muqual[i]};
         }

         std::vector<H1Event::Jet> njetNTuple(*njets);
         for (int i = 0; i < *njets; ++i) {
            njetNTuple.at(i) = H1Event::Jet{E_j[i], pt_j[i], theta_j[i], eta_j[i], phi_j[i], m_j[i]};
         }

         std::array<float, 4> pthrustNTuple;
         for (int i = 0; i < 4; ++i) {
            pthrustNTuple.at(i) = pthrust[i];
         }

         std::array<float, 4> pthrust2NTuple;
         for (int i = 0; i < 4; ++i) {
            pthrust2NTuple.at(i) = pthrust2[i];
         }

         H1Event eventEntry{/*0-9*/ *nrun, *nevent, *nentry, std::move(trelemNTuple), std::move(subtrNTuple), std::move(rawtrNTuple), std::move(L4subtrNTuple), std::move(L5classNTuple), *E33, *de33, /*10-19*/ *x33, *dx33, *y33, *dy33, *E44, *de44, *x44, *dx44, *y44, *dy44, /*20-29*/ *Ept, *dept, *xpt, *dxpt, *ypt, *dypt, std::move(pelecNTuple), *flagelec, *xeelec, *yeelec, /*30-39*/ *Q2eelec, /* *nelec,*/ std::move(nelecNTuple), sumcNTuple, /*40-49*/ *sumetc, *yjbc, *Q2jbc, std::move(sumctNTuple), *sumetct, *yjbct, *Q2jbct, *yjbct, *Q2jbct, std::move(pvtx_dNTuple), /*50-59*/ std::move(cpvtx_dNTuple), std::move(pvtx_tNTuple), std::move(cpvtx_tNTuple), *ntrkxy_t, *prbxy_t, *ntrkz_t, *prbz_t, *nds, *rankds, *qds, /*60-69*/ std::move(pds_dNTuple), *ptds_d, *etads_d, *dm_d, *ddm_d, std::move(pds_tNTuple), *dm_t, *ddm_t, *ik, *ipi, /*70-79*/ *ipis, std::move(pd0_dNTuple), *ptd0_d, *etad0_d, *md0_d, *dmd0_d, std::move(pd0_tNTuple), *md0_t, *dmd0_t, std::move(pk_rNTuple), /*80-89*/ std::move(ppi_rNTuple), std::move(pd0_rNTuple), *md0_r, std::move(Vtxd0_rNTuple), std::move(cvtxd0_rNTuple), *dxy_r, *dz_r, *psi_r, *rd0_d, *drd0_d, /*90-99*/ *rpd0_d, *drpd0_d, *rd0_t, *drd0_t, *rpd0_t, *drpd0_t, *rd0_dt, *drd0_dt, *prbr_dt, *prbz_dt, /*100-109*/ *rd0_tt, *drd0_tt, *prbr_tt, *prbz_tt, *ijetd0, *ptr3d0_j, *ptr2d0_j, *ptr3d0_3, *ptr2d0_3, *ptr2d0_2, /*110-134*/ *Mimpds_r, *Mimpbk_r, /* *ntracks,*/ std::move(ntrackNTuple), /*135-143*/ *imu, *imufe, /* *njets,*/ std::move(njetNTuple), /*144-151*/ *thrust, std::move(pthrustNTuple), *thrust2, std::move(pthrust2NTuple), *spher, *aplan, *plan, {nnout[0]}};
         *ev = eventEntry;
         ntuple->Fill();
      }
   } // Upon destruction of the RNTupleWriter object, the data is flushed to a ntuple file.
}

#endif
