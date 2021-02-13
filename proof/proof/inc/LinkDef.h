/* @(#)root/proof:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedtypedefs;
#pragma link C++ nestedclasses;

#pragma link C++ global gProof;
#pragma link C++ global gProofServ;
#pragma link C++ global gProofDebugMask;
#pragma link C++ global gProofDebugLevel;

#pragma link C++ class TDSet-;
#pragma link C++ class TDSetElement-;
#pragma link C++ class TLockPath+;
#pragma link C++ class TPackMgr+;
#pragma link C++ class TProofChain+;
#pragma link C++ class TProofMgr;
#pragma link C++ class TProofMgrLite;
#pragma link C++ class TProofDesc;
#pragma link C++ class TProof;
#pragma link C++ class TProofCondor;
#pragma link C++ class TProofLite;
#pragma link C++ class TProofSuperMaster;
#pragma link C++ class TSlaveInfo+;
#pragma link C++ class TProofServ;
#pragma link C++ class TProofServLite;
#pragma link C++ class TProofDebug;
#pragma link C++ class TProofLog;
#pragma link C++ class TProofLogElem;
#pragma link C++ class TSlave;
#pragma link C++ class TSlaveLite;
#pragma link C++ class TVirtualPacketizer+;
#pragma link C++ class TVirtualProofPlayer+;
#pragma link C++ class TProofQueryResult+;
#pragma link C++ class TQueryResultManager+;
#pragma link C++ class TDSetProxy+;
#pragma link C++ class TCondor+;
#pragma link C++ class TCondorSlave+;
#pragma link C++ class TProofNodeInfo;
#pragma link C++ class TProofResources;
#pragma link C++ class TProofResourcesStatic;
#pragma link C++ class TProofProgressStatus;
#pragma link C++ class TMergerInfo;
#pragma link C++ class TProofProgressInfo;
#pragma link C++ class TProofOutputList+;
#pragma link C++ class TProofOutputFile+;

#pragma link C++ class TDataSetManager;
#pragma link C++ class TDataSetManagerFile;

// Dictionary for TDataSetManagerAliEn only if requested
#ifdef ALIENDSMGR
#pragma link C++ class TDataSetManagerAliEn;
#pragma link C++ class TAliEnFind;
#endif

#pragma link C++ class TSelVerifyDataSet+;

// For backward compatibility with old client / masters
#pragma link C++ class std::pair<TDSetElement*, TString>;
#pragma link C++ class std::list<std::pair<TDSetElement*,TString> >;

#endif
