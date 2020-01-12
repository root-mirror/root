#ifndef ntpl005_h1event_hxx
#define ntpl005_h1event_hxx

#include <array>
#include <cstdint>
#include <vector>


// comment on the right side tells branch ordering.
struct H1Event {
   std::int32_t nrun; // 0
   std::int32_t nevent; // 1
   std::int32_t nentry; // 2
   std::array<bool, 192> trelem; // 3
   std::array<bool, 128> subtr; // 4
   std::array<bool, 128> rawtr; // 5
   std::array<bool, 128> L4subtr; // 6
   std::array<bool, 32> L5class; // 7
   float E33; // 8
   float de33; // 9
   float x33; // 10
   float dx33; // 11
   float y33; // 12
   float dy33; // 13
   float E44; // 14
   float de44; // 15
   float x44; // 16
   float dx44; // 17
   float y44; // 18
   float dy44; // 19
   float Ept; // 20
   float dept; // 21
   float xpt; // 22
   float dxpt; // 23
   float ypt; // 24
   float dypt; // 25
   std::array<float, 4> pelec; // 26
   std::int32_t flagelec; // 27
   float xeelec; // 28
   float yeelec; // 29
   float Q2eelec; // 30
   //std::int32_t nelec; // 31
   
   struct Electron {
      float Eelec; // 32
      float thetelec; // 33
      float phielec; // 34
      float xelec; // 35
      float Q2elec; // 36
      float xsigma; // 37
      float Q2sigma; // 38
   };
   
   std::vector<Electron> electrons; // elecInfo.size() = nelec
   std::array<float, 4> sumc; // 39
   float sumetc; // 40
   float yjbc; // 41
   float Q2jbc; // 42
   std::array<float, 4> sumct; // 43
   float sumetct; // 44
   float yjbct; // 45
   float Q2jbct; // 46
   float Ebeamel; // 47
   float Ebeampr; // 48
   std::array<float, 3> pvtx_d; // 49
   std::array<float, 6> cpvtx_d; // 50
   std::array<float, 3> pvtx_t; // 51
   std::array<float, 6> cpvtx_t; // 52
   std::int32_t ntrkxy_t; // 53
   float prbxy_t; // 54
   std::int32_t ntrkz_t; // 55
   float prbz_t; // 56
   std::int32_t nds; // 57
   std::int32_t rankds; // 58
   std::int32_t qds; // 59
   std::array<float, 4> pds_d; // 60
   float ptds_d; // 61
   float etads_d; // 62
   float dm_d; // 63
   float ddm_d; // 64
   std::array<float, 4> pds_t; // 65
   float dm_t; // 66
   float ddm_t; // 67
   std::int32_t ik; // 68
   std::int32_t ipi; // 69
   std::int32_t ipis; // 70
   std::array<float, 4> pd0_d; // 71
   float ptd0_d; // 72
   float etad0_d; // 73
   float md0_d; // 74
   float dmd0_d; // 75
   std::array<float, 4> pd0_t; // 76
   float md0_t; // 77
   float dmd0_t; // 78
   std::array<float, 4> pk_r; // 79
   std::array<float, 4> ppi_r; // 80
   std::array<float, 4> pd0_r; // 81
   float md0_r; // 82
   std::array<float, 3> Vtxd0_r; // 83
   std::array<float, 6> cvtxd0_r; // 84
   float dxy_r; // 85
   float dz_r; // 86
   float psi_r; // 87
   float rd0_d; // 88
   float drd0_d; // 89
   float rpd0_d; // 90
   float drpd0_d; // 91
   float rd0_t; // 92
   float drd0_t; // 93
   float rpd0_t; // 94
   float drpd0_t; // 95
   float rd0_dt; // 96
   float drd0_dt; // 97
   float prbr_dt; // 98
   float prbz_dt; // 99
   float rd0_tt; // 100
   float drd0_tt; // 101
   float prbr_tt; // 102
   float prbz_tt; // 103
   std::int32_t ijetd0; // 104
   float ptr3d0_j; // 105
   float ptr2d0_j; // 106
   float ptr3d0_3; // 107
   float ptr2d0_3; // 108
   float ptr2d0_2; // 109
   float Mimpds_r; // 110
   float Mimpbk_r; // 111
   //std::int32_t ntracks; // 112 -> almost always has value 3
   
   struct Track {
      float pt; // 113
      float kappa; // 114
      float phi; // 115
      float theta; // 116
      float dca; // 117
      float z0; // 118
      std::array<float, 15> covar; // 119
      std::int32_t nhitrp; //120
      float prbrp; // 121
      std::int32_t nhitz; // 122
      float prbz; // 123
      float rstart; // 124
      float rend; // 125
      float lhk; // 126
      float lhpi; // 127
      float nlhk; // 128
      float nlhpi; // 129
      float dca_d; // 130
      float ddca_d; // 131
      float dca_t; // 132
      float ddca_t; // 133
      std::int32_t muqual; // 134
   };
   
   std::vector<Track> tracks; // trackInfo.size() = ntracks
   std::int32_t imu; // 135
   std::int32_t imufe; // 136
   //std::int32_t njets; // 137
   
   struct Jet {
      float E_j; // 138
      float pt_j; // 139
      float theta_j; // 140
      float eta_j; // 141
      float phi_j; // 142
      float m_j; // 143
   };
   
   std::vector<Jet> jets; // jetInfo.size() = njets
   float thrust; // 144
   std::array<float, 4> pthrust; // 145
   float thrust2; // 146
   std::array<float, 4> pthrust2; // 147
   float spher; // 148
   float aplan; // 149
   float plan; // 150
   std::array<float, 1> nnout; // 151
};

#endif
