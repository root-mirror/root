#include <vector>
#include <map>

typedef std::multimap<TObject*, TF1*>::iterator fPrevFitIter;
typedef std::vector<TF1*>::iterator             fSystemFuncIter;

enum EFitPanel {
   kFP_FLIST = 0, kFP_GAUS = 1,  kFP_GAUSN, kFP_EXPO,  kFP_LAND,  kFP_LANDN,
   kFP_POL0,  kFP_POL1,  kFP_POL2,  kFP_POL3,  kFP_POL4,  kFP_POL5,
   kFP_POL6,  kFP_POL7,  kFP_POL8,  kFP_POL9,
   kFP_CHEB0,  kFP_CHEB1,  kFP_CHEB2,  kFP_CHEB3,  kFP_CHEB4,  kFP_CHEB5,
   kFP_CHEB6,  kFP_CHEB7,  kFP_CHEB8,  kFP_CHEB9,
   kFP_XYGAUS, kFP_BIGAUS, kFP_XYEXP, kFP_XYLAN, kFP_XYLANN,
// Above here -> All editable formulaes!
   kFP_USER,
   kFP_NONE,  kFP_ADD,   kFP_NORMADD,   kFP_CONV,  kFP_FILE,  kFP_PARS,  kFP_RBUST, kFP_EMPW1,
   kFP_INTEG, kFP_IMERR, kFP_USERG, kFP_ADDLS, kFP_ALLW1, kFP_IFITR, kFP_NOCHI,
   kFP_MLIST, kFP_MCHIS, kFP_MBINL, kFP_MUBIN, kFP_MUSER, kFP_MLINF, kFP_MUSR,
   kFP_DSAME, kFP_DNONE, kFP_DADVB, kFP_DNOST, kFP_PDEF,  kFP_PVER,  kFP_PQET,
   kFP_XMIN,  kFP_XMAX,  kFP_YMIN,  kFP_YMAX,  kFP_ZMIN,  kFP_ZMAX,

   kFP_LMIN,  kFP_LMIN2, kFP_LFUM,  kFP_LGSL,  kFP_LGAS,  kFP_MIGRAD,kFP_SIMPLX,
   kFP_FUMILI,kFP_COMBINATION,      kFP_MINMETHOD,
   kFP_GSLFR, kFP_GSLPR, kFP_BFGS,  kFP_BFGS2, kFP_GSLLM, kFP_GSLSA,
   kFP_SCAN,  kFP_TMVAGA,kFP_GALIB,

   kFP_MERR,  kFP_MTOL,  kFP_MITR,

   kFP_UPDATE, kFP_FIT,   kFP_RESET, kFP_CLOSE,

   // New GUI elements from here!
   kFP_TLIST, kFP_PRED1D, kFP_PRED2D, kFP_PRED3D, kFP_UFUNC, kFP_PREVFIT, kFP_ROOFIT,
   kFP_DATAS,

   kFP_NOSEL = 8000,
   kFP_ALTFUNC = 10000

};

enum EParStruct {
   PAR_VAL = 0,
   PAR_MIN = 1,
   PAR_MAX = 2
};
