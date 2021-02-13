// @(#)root/pythia6:$Id$
// Author: Rene Brun   19/10/99

/** \class TPythia6
    \ingroup pythia6

TPythia is an interface class to F77 version of Pythia 6.2          

To use this class you must install a version of pythia6.            
See the installation instructions at                              
  http://root.cern.ch/root/Install.html                         
                                                                    
CERNLIB event generators, written by T.Sjostrand.                   
For the details about these generators look at Pythia/Jetset manual:

\verbatim
******************************************************************************
**                                                                          **
**                                                                          **
**  PPP  Y   Y TTTTT H   H III   A        JJJJ EEEE TTTTT  SSS  EEEE TTTTT  **
**  P  P  Y Y    T   H   H  I   A A          J E      T   S     E      T    **
**  PPP    Y     T   HHHHH  I  AAAAA         J EEE    T    SSS  EEE    T    **
**  P      Y     T   H   H  I  A   A      J  J E      T       S E      T    **
**  P      Y     T   H   H III A   A       JJ  EEEE   T    SSS  EEEE   T    **
**                                                                          **
**                                                                          **
**              *......*                  Welcome to the Lund Monte Carlo!  **
**         *:::!!:::::::::::*                                               **
**      *::::::!!::::::::::::::*            This is PYTHIA version 5.720    **
**    *::::::::!!::::::::::::::::*        Last date of change: 29 Nov 1995  **
**   *:::::::::!!:::::::::::::::::*                                         **
**   *:::::::::!!:::::::::::::::::*         This is JETSET version 7.408    **
**    *::::::::!!::::::::::::::::*!       Last date of change: 23 Aug 1995  **
**      *::::::!!::::::::::::::* !!                                         **
**      !! *:::!!:::::::::::*    !!                 Main author:            **
**      !!     !* -><- *         !!              Torbjorn Sjostrand         **
**      !!     !!                !!        Dept. of theoretical physics 2   **
**      !!     !!                !!              University of Lund         **
**      !!                       !!                Solvegatan 14A           **
**      !!        ep             !!             S-223 62 Lund, Sweden       **
**      !!                       !!          phone: +46 - 46 - 222 48 16    **
**      !!                 pp    !!          E-mail: torbjorn@thep.lu.se    **
**      !!   e+e-                !!                                         **
**      !!                       !!         Copyright Torbjorn Sjostrand    **
**      !!                                     and CERN, Geneva 1993        **
**                                                                          **
**                                                                          **
** The latest program versions and documentation is found on WWW address    **
** http://thep.lu.se/tf2/staff/torbjorn/Welcome.html                        **
**                                                                          **
** When you cite these programs, priority should always be given to the     **
** latest published description. Currently this is                          **
** T. Sjostrand, Computer Physics Commun. 82 (1994) 74.                     **
** The most recent long description (unpublished) is                        **
** T. Sjostrand, LU TP 95-20 and CERN-TH.7112/93 (revised August 1995).     **
** Also remember that the programs, to a large extent, represent original   **
** physics research. Other publications of special relevance to your        **
** studies may therefore deserve separate mention.                          **
**                                                                          **
**                                                                          **
******************************************************************************
\endverbatim
*/

#include "TPythia6.h"

#include "TClonesArray.h"
#include "TMCParticle.h"
#include "TParticle.h"
#include "snprintf.h"

TPythia6*  TPythia6::fgInstance = nullptr;


#ifndef WIN32
# define pydiff pydiff_
# define pyevnt pyevnt_
# define pyinit pyinit_
# define pychge pychge_
# define pycomp pycomp_
# define pyedit pyedit_
# define pyexec pyexec_
# define pyhepc pyhepc_
# define pygive pygive_
# define pylist pylist_
# define pymass pymass_
# define pyname pyname_
# define pyr    pyr_
# define pyrget pyrget_
# define pyrset pyrset_
# define pystat pystat_
# define pytest pytest_
# define pytune pytune_
# define pyupda pyupda_
# define py1ent py1ent_
# ifdef PYTHIA6_DOUBLE_UNDERSCORE
#  define tpythia6_open_fortran_file tpythia6_open_fortran_file__
#  define tpythia6_close_fortran_file tpythia6_close_fortran_file__
#  define pythia6_common_address pythia6_common_address__
# elif PYTHIA6_SINGLE_UNDERSCORE
#  define tpythia6_open_fortran_file tpythia6_open_fortran_file_
#  define tpythia6_close_fortran_file tpythia6_close_fortran_file_
#  define pythia6_common_address pythia6_common_address
# else
#  define pythia6_common_address pythia6_common_address
#  define tpythia6_open_fortran_file tpythia6_open_fortran_file_
#  define tpythia6_close_fortran_file tpythia6_close_fortran_file_
# endif
# define type_of_call
#else
# define pydiff PYDIFF
# define pyevnt PYEVNT
# define pyinit PYINIT
# define pychge PYCHGE
# define pycomp PYCOMP
# define pyedit PYEDIT
# define pyexec PYEXEC
# define pygive PYGIVE
# define pyhepc PYHEPC
# define pylist PYLIST
# define pymass PYMASS
# define pyname PYNAME
# define pyr    PYR
# define pyrget PYRGET
# define pyrset PYRSET
# define pystat PYSTAT
# define pytest PYTEST
# define pytune PYTUNE
# define pyupda PYUPDA
# define py1ent PY1ENT
# define tpythia6_open_fortran_file TPYTHIA6_OPEN_FORTRAN_FILE
# define tpythia6_close_fortran_file TPYTHIA6_CLOSE_FORTRAN_FILE
# define type_of_call _stdcall
#endif


extern "C" void type_of_call pyevnt();
extern "C" void type_of_call pystat(int *key);
extern "C" void type_of_call pylist(int *key);
extern "C" void type_of_call pyedit(int *medit);
extern "C" void type_of_call pydiff();
extern "C" void type_of_call pyexec();
extern "C" void type_of_call pygive(const char *param, Long_t lparam);
extern "C" void type_of_call pyhepc(int *mconv);
extern "C" void type_of_call pylist(int *flag);
extern "C" int  type_of_call pychge(int *kf);
extern "C" int  type_of_call pycomp(int *kf);
extern "C" double  type_of_call pymass(int *kf);
extern "C" void type_of_call pyname(int *kf, char *name, Long_t l_name);
extern "C" int  type_of_call pyr(int *dummy);
extern "C" int  type_of_call pyrget(int *lun, int *move);
extern "C" int  type_of_call pyrset(int *lun, int *move);
extern "C" int  type_of_call pytest(int *flag);
extern "C" int  type_of_call pytune(int *itune);
extern "C" int  type_of_call pyupda(int *mupda, int *lun);
extern "C" void type_of_call py1ent(Int_t&, Int_t&, Double_t&, Double_t&, Double_t&);

#ifndef WIN32
extern "C" void type_of_call pyinit(char *frame, char *beam, char *target,
                                     double *win, Long_t l_frame, Long_t l_beam,
                                     Long_t l_target);
#else
extern "C" void type_of_call pyinit(char *frame,  Long_t l_frame,
                                    char *beam,   Long_t l_beam,
                                    char *target, Long_t l_target,
                                    double *win
                                    );
#endif

extern "C" {
   void*  pythia6_common_address(const char*);
   void   type_of_call tpythia6_open_fortran_file(int* lun, char* name, int);
   void   type_of_call tpythia6_close_fortran_file(int* lun);
}

ClassImp(TPythia6);

/** \class TPythia6::TPythia6Cleaner
    \ingroup pythia6

Utility class to manage the TPythia6 instance
*/

TPythia6::TPythia6Cleaner::TPythia6Cleaner() {
}

////////////////////////////////////////////////////////////////////////////////
///delete the TPythia6 insntance

TPythia6::TPythia6Cleaner::~TPythia6Cleaner() {
   if (TPythia6::fgInstance) {
      delete TPythia6::fgInstance;
      TPythia6::fgInstance = 0;
   }
}

//------------------------------------------------------------------------------
//  constructor is not supposed to be called from the outside - only
//  Initialize() method
////////////////////////////////////////////////////////////////////////////////
/// TPythia6 constructor: creates a TClonesArray in which it will store all
/// particles. Note that there may be only one functional TPythia6 object
/// at a time, so it's not use to create more than one instance of it.

TPythia6::TPythia6() : TGenerator("TPythia6","TPythia6") {
   // Protect against multiple objects.   All access should be via the
   // Instance member function.
   if (fgInstance)
      Fatal("TPythia6", "There's already an instance of TPythia6");

   delete fParticles; // was allocated as TObjArray in TGenerator

   fParticles = new TClonesArray("TMCParticle",50);

   // initialize common-blocks
   // the functions/subroutines referenced by TPythia6 can be found
   // at ftp://root.cern.ch/root/pythia6.tar.gz

   fPyjets = (Pyjets_t*) pythia6_common_address("PYJETS");
   fPydat1 = (Pydat1_t*) pythia6_common_address("PYDAT1");
   fPydat2 = (Pydat2_t*) pythia6_common_address("PYDAT2");
   fPydat3 = (Pydat3_t*) pythia6_common_address("PYDAT3");
   fPydat4 = (Pydat4_t*) pythia6_common_address("PYDAT4");
   fPydatr = (Pydatr_t*) pythia6_common_address("PYDATR");
   fPysubs = (Pysubs_t*) pythia6_common_address("PYSUBS");
   fPypars = (Pypars_t*) pythia6_common_address("PYPARS");
   fPyint1 = (Pyint1_t*) pythia6_common_address("PYINT1");
   fPyint2 = (Pyint2_t*) pythia6_common_address("PYINT2");
   fPyint3 = (Pyint3_t*) pythia6_common_address("PYINT3");
   fPyint4 = (Pyint4_t*) pythia6_common_address("PYINT4");
   fPyint5 = (Pyint5_t*) pythia6_common_address("PYINT5");
   fPyint6 = (Pyint6_t*) pythia6_common_address("PYINT6");
   fPyint7 = (Pyint7_t*) pythia6_common_address("PYINT7");
   fPyint8 = (Pyint8_t*) pythia6_common_address("PYINT8");
   fPyint9 = (Pyint9_t*) pythia6_common_address("PYINT9");
   fPymssm = (Pymssm_t*) pythia6_common_address("PYMSSM");
   fPyssmt = (Pyssmt_t*) pythia6_common_address("PYSSMT");
   fPyints = (Pyints_t*) pythia6_common_address("PYINTS");
   fPybins = (Pybins_t*) pythia6_common_address("PYBINS");
}

////////////////////////////////////////////////////////////////////////////////

TPythia6::TPythia6(const TPythia6& p6) :
  TGenerator(p6),
  fPyjets(p6.fPyjets),
  fPydat1(p6.fPydat1),
  fPydat2(p6.fPydat2),
  fPydat3(p6.fPydat3),
  fPydat4(p6.fPydat4),
  fPydatr(p6.fPydatr),
  fPysubs(p6.fPysubs),
  fPypars(p6.fPypars),
  fPyint1(p6.fPyint1),
  fPyint2(p6.fPyint2),
  fPyint3(p6.fPyint3),
  fPyint4(p6.fPyint4),
  fPyint5(p6.fPyint5),
  fPyint6(p6.fPyint6),
  fPyint7(p6.fPyint7),
  fPyint8(p6.fPyint8),
  fPyint9(p6.fPyint9),
  fPymssm(p6.fPymssm),
  fPyssmt(p6.fPyssmt),
  fPyints(p6.fPyints),
  fPybins(p6.fPybins)
{ }

////////////////////////////////////////////////////////////////////////////////
/// Destroys the object, deletes and disposes all TMCParticles currently on list.

TPythia6::~TPythia6()
{
   if (fParticles) {
      fParticles->Delete();
      delete fParticles;
      fParticles = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// model of automatic memory cleanup suggested by Jim Kowalkovski:
/// destructor for local static variable `cleaner' is  always called in the end
/// of the job thus deleting the only TPythia6 instance

TPythia6* TPythia6::Instance() {
   static TPythia6::TPythia6Cleaner cleaner;
   return fgInstance ? fgInstance : (fgInstance=new TPythia6()) ;
}





////////////////////////////////////////////////////////////////////////////////
///  generate event and copy the information from /HEPEVT/ to fPrimaries

void TPythia6::GenerateEvent() {
   pyevnt();
   ImportParticles();
}

////////////////////////////////////////////////////////////////////////////////
///interface with fortran i/o

void TPythia6::OpenFortranFile(int lun, char* name) {
   tpythia6_open_fortran_file(&lun, name, strlen(name));
}

////////////////////////////////////////////////////////////////////////////////
///interface with fortran i/o

void TPythia6::CloseFortranFile(int lun) {
   tpythia6_close_fortran_file(&lun);
}


////////////////////////////////////////////////////////////////////////////////
/// Fills TObjArray fParticles list with particles from common LUJETS
/// Old contents of a list are cleared. This function should be called after
/// any change in common LUJETS, however GetParticles() method  calls it
/// automatically - user don't need to care about it. In case you make a call
/// to LuExec() you must call this method yourself to transfer new data from
/// common LUJETS to the fParticles list.

TObjArray *TPythia6::ImportParticles(Option_t *)
{
   fParticles->Clear();
   Int_t numpart   = fPyjets->N;
   TClonesArray &a = *((TClonesArray*)fParticles);
   for (Int_t i = 0; i<numpart; i++) {
      new(a[i]) TMCParticle(fPyjets->K[0][i] ,
                            fPyjets->K[1][i] ,
                            fPyjets->K[2][i] ,
                            fPyjets->K[3][i] ,
                            fPyjets->K[4][i] ,
                            fPyjets->P[0][i] ,
                            fPyjets->P[1][i] ,
                            fPyjets->P[2][i] ,
                            fPyjets->P[3][i] ,
                            fPyjets->P[4][i] ,
                            fPyjets->V[0][i] ,
                            fPyjets->V[1][i] ,
                            fPyjets->V[2][i] ,
                            fPyjets->V[3][i] ,
                            fPyjets->V[4][i]);
   }
   return fParticles;
}

////////////////////////////////////////////////////////////////////////////////
///  Default primary creation method. It reads the /HEPEVT/ common block which
///  has been filled by the GenerateEvent method. If the event generator does
///  not use the HEPEVT common block, This routine has to be overloaded by
///  the subclasses.
///  The function loops on the generated particles and store them in
///  the TClonesArray pointed by the argument particles.
///  The default action is to store only the stable particles (ISTHEP = 1)
///  This can be demanded explicitly by setting the option = "Final"
///  If the option = "All", all the particles are stored.
///

Int_t TPythia6::ImportParticles(TClonesArray *particles, Option_t *option)
{
   if (particles == 0) return 0;
   TClonesArray &clonesParticles = *particles;
   clonesParticles.Clear();
   Int_t numpart = fPyjets->N;
   Int_t nparts=0;
   if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (Int_t i = 0; i<numpart; i++) {

      if (fPyjets->K[0][i] == 1) {
         //
         //  Use the common block values for the TParticle constructor
         //
         new(clonesParticles[nparts]) TParticle(
                            fPyjets->K[1][i] ,
                            fPyjets->K[0][i] ,
                            fPyjets->K[2][i] ,
                            -1,
                            fPyjets->K[3][i] ,
                            fPyjets->K[4][i] ,
                            fPyjets->P[0][i] ,
                            fPyjets->P[1][i] ,
                            fPyjets->P[2][i] ,
                            fPyjets->P[3][i] ,
                            fPyjets->V[0][i] ,
                            fPyjets->V[1][i] ,
                            fPyjets->V[2][i] ,
                            fPyjets->V[3][i]);

         //     if(gDebug) printf("%d %d %d! ",i,fPyjets->K[1][i],numpart);
         nparts++;
         }
      }
   } else if (!strcmp(option,"All")) {
      for (Int_t i = 0; i<numpart; i++) {
         new(clonesParticles[i]) TParticle(
                            fPyjets->K[1][i] ,
                            fPyjets->K[0][i] ,
                            fPyjets->K[2][i] ,
                            -1,
                            fPyjets->K[3][i] ,
                            fPyjets->K[4][i] ,
                            fPyjets->P[0][i] ,
                            fPyjets->P[1][i] ,
                            fPyjets->P[2][i] ,
                            fPyjets->P[3][i] ,
                            fPyjets->V[0][i] ,
                            fPyjets->V[1][i] ,
                            fPyjets->V[2][i] ,
                            fPyjets->V[3][i]);
      }
      nparts=numpart;
   }

   return nparts;
}

////////////////////////////////////////////////////////////////////////////////
/// Calls PyInit with the same parameters after performing some checking,
/// sets correct title. This method should preferably be called instead of PyInit.
/// PURPOSE: to initialize the generation procedure.
/// ARGUMENTS: See documentation for details.
/// -  frame:  - specifies the frame of the experiment:
///                "CMS","FIXT","USER","FOUR","FIVE","NONE"
/// -  beam,
/// -  target: - beam and target particles (with additionaly charges, tildes or "bar":
///              e,nu_e,mu,nu_mu,tau,nu_tau,gamma,pi,n,p,Lambda,Sigma,Xi,Omega,
///              pomeron,reggeon
/// -  win:    - related to energy system:
///              - for frame=="CMS" - total energy of system
///              - for frame=="FIXT" - momentum of beam particle
///              - for frame=="USER" - dummy - see documentation.

void TPythia6::Initialize(const char *frame, const char *beam, const char *target, float win)
{
   char  cframe[4];
   strlcpy(cframe,frame,4);
   char  cbeam[10];
   strlcpy(cbeam,beam,10);
   char  ctarget[10];
   strlcpy(ctarget,target,10);

   // For frames "3MOM", "4MOM" and "5MOM" see p. 181-182 of the version 6 manual,
   // http://home.thep.lu.se/~torbjorn/pythia/lutp0613man2.pdf
   // their usage may depend on the version of Pythia6 used
   if ( (!strncmp(frame, "CMS"  ,3)) &&
        (!strncmp(frame, "FIXT" ,4)) &&
        (!strncmp(frame, "USER" ,4)) &&
        (!strncmp(frame, "FOUR" ,4)) &&
        (!strncmp(frame, "FIVE" ,4)) &&
        (!strncmp(frame, "3MOM" ,4)) &&
        (!strncmp(frame, "4MOM" ,4)) &&
        (!strncmp(frame, "5MOM" ,4)) &&
        (!strncmp(frame, "NONE" ,4)) ) {
      printf("WARNING! In TPythia6:Initialize():\n");
      printf(" specified frame=%s is neither of CMS,FIXT,USER,FOUR,FIVE,NONE,3MOM,4MOM,5MOM\n",frame);
      printf(" resetting to \"CMS\" .");
      snprintf(cframe,4,"CMS");
   }

   if ( (!strncmp(beam, "e"       ,1)) &&
        (!strncmp(beam, "nu_e"    ,4)) &&
        (!strncmp(beam, "mu"      ,2)) &&
        (!strncmp(beam, "nu_mu"   ,5)) &&
        (!strncmp(beam, "tau"     ,3)) &&
        (!strncmp(beam, "nu_tau"  ,6)) &&
        (!strncmp(beam, "gamma"   ,5)) &&
        (!strncmp(beam, "pi"      ,2)) &&
        (!strncmp(beam, "n"       ,1)) &&
        (!strncmp(beam, "p"       ,1)) &&
        (!strncmp(beam, "Lambda"  ,6)) &&
        (!strncmp(beam, "Sigma"   ,5)) &&
        (!strncmp(beam, "Xi"      ,2)) &&
        (!strncmp(beam, "Omega"   ,5)) &&
        (!strncmp(beam, "pomeron" ,7)) &&
        (!strncmp(beam, "reggeon" ,7)) ) {
      printf("WARNING! In TPythia6:Initialize():\n");
      printf(" specified beam=%s is unrecognized .\n",beam);
      printf(" resetting to \"p+\" .");
      snprintf(cbeam,8,"p+");
   }

   if ( (!strncmp(target, "e"       ,1)) &&
        (!strncmp(target, "nu_e"    ,4)) &&
        (!strncmp(target, "mu"      ,2)) &&
        (!strncmp(target, "nu_mu"   ,5)) &&
        (!strncmp(target, "tau"     ,3)) &&
        (!strncmp(target, "nu_tau"  ,6)) &&
        (!strncmp(target, "gamma"   ,5)) &&
        (!strncmp(target, "pi"      ,2)) &&
        (!strncmp(target, "n"       ,1)) &&
        (!strncmp(target, "p"       ,1)) &&
        (!strncmp(target, "Lambda"  ,6)) &&
        (!strncmp(target, "Sigma"   ,5)) &&
        (!strncmp(target, "Xi"      ,2)) &&
        (!strncmp(target, "Omega"   ,5)) &&
        (!strncmp(target, "pomeron" ,7)) &&
        (!strncmp(target, "reggeon" ,7)) ){
      printf("WARNING! In TPythia6:Initialize():\n");
      printf(" specified target=%s is unrecognized.\n",target);
      printf(" resetting to \"p+\" .");
      snprintf(ctarget,8,"p+");
   }

   Pyinit(cframe, cbeam ,ctarget, win);

   char atitle[64];
   snprintf(atitle, sizeof(atitle)," %s-%s at %g GeV", cbeam, ctarget, win);
   SetTitle(atitle);
}


void TPythia6::Pyinit(char* frame, char* beam, char* target, double win)
{
   // Calls Pyinit with the same parameters after performing some checking,
   // sets correct title. This method should preferably be called instead of PyInit.
   // PURPOSE: to initialize the generation procedure.
   // ARGUMENTS: See documentation for details.
   //    frame:  - specifies the frame of the experiment:
   //                "CMS","FIXT","USER","FOUR","FIVE","NONE"
   //    beam,
   //    target: - beam and target particles (with additionaly charges,
   //              tildes or "bar":
   //              e,nu_e,mu,nu_mu,tau,nu_tau,gamma,pi,n,p,Lambda,Sigma,Xi,Omega,
   //              pomeron,reggeon
   //    win:    - related to energy system:
   //              for frame=="CMS" - total energy of system
   //              for frame=="FIXT" - momentum of beam particle
   //              for frame=="USER" - dummy - see documentation.

   Double_t lwin = win;
   Long_t  s1    = strlen(frame);
   Long_t  s2    = strlen(beam);
   Long_t  s3    = strlen(target);
#ifndef WIN32
   pyinit(frame,beam,target,&lwin,s1,s2,s3);
#else
   pyinit(frame, s1, beam , s2, target, s3, &lwin);
#endif
}


int TPythia6::Pycomp(int kf) {
   //interface with fortran routine pycomp
   return pycomp(&kf);
}

void TPythia6::Pyedit(int medit) {
   //interface with fortran routine pyedit
   pyedit(&medit);
   ImportParticles();
}

void TPythia6::Pydiff() {
   //interface with fortran routine pydiff
   pydiff();
}

void TPythia6::Pyevnt() {
   //interface with fortran routine pyevnt
   pyevnt();
}

void TPythia6::Pyexec() {
   //interface with fortran routine pyexec
   pyexec();
}

void TPythia6::Pygive(const char *param) {
   //interface with fortran routine pygive
   Long_t lparam = strlen(param);
   pygive(param,lparam);
}

void TPythia6::Pyhepc(int mconv) {
   //interface with fortran routine pyhepc
   pyhepc(&mconv);
}

void TPythia6::Pylist(int flag) {
   //interface with fortran routine pylist
   pylist(&flag);
}

void TPythia6::Pyname(int kf, char* name) {
   //Note that the array name must be dimensioned in the calling program
   //to at least name[16]

   pyname(&kf,name,15);
   // cut trailing blanks to get C string
   name[15] = 0;
   //for (int i=15; (i>=0) && (name[i] == ' '); i--) {
   //  name[i] = 0;
   // }
}

double TPythia6::Pyr(int idummy) {
   //interface with fortran routine pyr
   return pyr(&idummy);
}

void TPythia6::Pyrget(int lun, int move) {
   //interface with fortran routine pyrget
   pyrget(&lun,&move);
}

void TPythia6::Pyrset(int lun, int move) {
   //interface with fortran routine pyrset
   pyrset(&lun,&move);
}

void TPythia6::Pystat(int flag) {
   //interface with fortran routine pystat
   pystat(&flag);
}

void TPythia6::Pytest(int flag) {
   //interface with fortran routine pytest
   pytest(&flag);
}

void TPythia6::Pytune(int itune) {
   //interface with fortran routine pytune
   pytune(&itune);
}

void TPythia6::Pyupda(int mupda, int lun) {
   //interface with fortran routine pyupda
   pyupda(&mupda,&lun);
}

double TPythia6::Pymass(int kf) {
   //interface with fortran routine pymass
   return pymass(&kf);
}

int TPythia6::Pychge(int kf) {
   //interface with fortran routine pychge
   return pychge(&kf);
}

////////////////////////////////////////////////////////////////////////////////
/// Add one entry to the event record, i.e. either a parton or a
/// particle.
///
/// - IP:   normally line number for the parton/particle. There are two
///         exceptions:
///         - If IP = 0: line number 1 is used and PYEXEC is called.
///         - If IP < 0: line -IP is used, with status code K(-IP,2)=2
///                   rather than 1; thus a parton system may be built
///                   up by filling all but the last parton of the
///                   system with IP < 0.
///  - KF:   parton/particle flavour code (PDG code)
///  - PE:   parton/particle energy. If PE is smaller than the mass,
///          the parton/particle is taken to be at rest.
///  - THETA:
///  - PHI:  polar and azimuthal angle for the momentum vector of the
///          parton/particle.

void TPythia6::Py1ent(Int_t ip, Int_t kf, Double_t pe, Double_t theta, Double_t phi)
{
   py1ent(ip, kf, pe, theta, phi);
}


////////////////////////////////////////////////////////////////////////////////
/// Exemplary setup of Pythia parameters:
/// Switches on processes 102,123,124 (Higgs generation) and switches off
/// interactions, fragmentation, ISR, FSR...

void TPythia6::SetupTest()
{
   SetMSEL(0);            // full user controll;

   SetMSUB(102,1);        // g + g -> H0
   SetMSUB(123,1);        // f + f' -> f + f' + H0
   SetMSUB(124,1);        // f + f' -> f" + f"' + H0


   SetPMAS(6,1,175.0);   // mass of TOP
   SetPMAS(25,1,300);    // mass of Higgs


   SetCKIN(1,290.0);     // range of allowed mass
   SetCKIN(2,310.0);

   SetMSTP(61,  0);      // switch off ISR
   SetMSTP(71,  0);      // switch off FSR
   SetMSTP(81,  0);      // switch off multiple interactions
   SetMSTP(111, 0);      // switch off fragmentation/decay
}
