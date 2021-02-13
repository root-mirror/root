// @(#)root/alien:$Id: db925cce3e3dc34e650f92e7e84655fb84ef04e9 $
// Author: Andreas Peters   5/5/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAlien                                                               //
//                                                                      //
// Class defining interface to TAlien GRID services.                    //
//                                                                      //
// To force to connect to a running API Service, use                    //
//   - TGrid::Connect("alien://<apihosturl>/","user");                  //
//                                                                      //
// If you want to use TGrid::Connect("alien://");                       //
// you have to set the environment variables                            //
// export alien_API_HOST=<host>                                         //
// export alien_API_PORT=<port>                                         //
// export alien_API_USER=<user>                                         //
//                                                                      //
//                                                                      //
// some examples:                                                       //
// -------------------------------------------------------------------- //
// connect in an initialized API environemnt (see above)                //
// > TGrid::Connect("alien://");                                        //
// -> this is not thread safe                                           //
// connect using an already established shell token:                    //
// > TGrid::Connect("alien://",0,0,"t");                                //
// -> this is thread safe                                               //
// -------------------------------------------------------------------- //
// change the working directory                                         //
// > gGrid->Cd("/alice"); //=> returns 1 for success, 0 for error       //
//                                                                      //
//                                                                      //
// -------------------------------------------------------------------- //
// get the working directory                                            //
// > printf("Working Directory is %s\n",gGrid->Pwd());                  //
//                                                                      //
//                                                                      //
// -------------------------------------------------------------------- //
// get a filename list of the working directory and print it            //
// > TGridResult* result = gGrid->Ls("/alice");                         //
// > Int_t i=0;                                                         //
// > while (result->GetFileName(i))\                                    //
// printf("File %s\n",result->GetFileName(i++));                        //
// -------------------------------------------------------------------- //
//                                                                      //
// get all file permissions                                             //
// > TGridResult* result = gGrid->Ls("/alice","-la");                   //
// > while (result->GetFileName(i))\                                    //
// printf("The permissions are %s\n",result->GetKey(i++,"permissions")  //
//                                                                      //
// => the defined keys for GetKey can be seen using result->Print();    //
// -------------------------------------------------------------------- //
//////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include "TAlien.h"
#include "TUrl.h"
#include "Riostream.h"
#include "TString.h"
#include "TObjString.h"
#include "TMap.h"
#include "TSystem.h"
#include "TAlienJDL.h"
#include "TAlienResult.h"
#include "TAlienJob.h"
#include "TAlienCollection.h"
#include "TAlienMasterJobStatus.h"
#include "TAlienJobStatusList.h"
#include "TAlienJobStatus.h"
#include "gapiUI.h"
#include "gapi_job_operations.h"

using namespace std;



ClassImp(TAlien);

////////////////////////////////////////////////////////////////////////////////

TAlien::TAlien(const char *gridurl, const char *uid, const char * passwd,
               const char *options)
{
   gSystem->Unsetenv("GCLIENT_EXTRA_ARG");

   // Connect to the AliEn grid.
   TUrl *gurl = new TUrl(gridurl);

   fGridUrl = gridurl;
   fGrid    = "alien";
   if (!strlen(gurl->GetHost())) {
      if (gSystem->Getenv("alien_API_HOST"))
         fHost = gSystem->Getenv("alien_API_HOST");
      else
         fHost = "";
   } else {
      fHost = gurl->GetHost();
   }

   if (gurl->GetPort()<=0) {
      if (gSystem->Getenv("alien_API_PORT"))
         fPort = atoi(gSystem->Getenv("alien_API_PORT"));
      else
         fPort = 0;
   } else {
      fPort = gurl->GetPort();
   }

   if (!uid[0]) {
      if (gSystem->Getenv("alien_API_USER")) {
         fUser = gSystem->Getenv("alien_API_USER");
      } else {
         if (gSystem->Getenv("LOGNAME")) {
            // we try the LOGNAME env
            fUser = gSystem->Getenv("LOGNAME");
         } else {
            // we set the USER env
            fUser = gSystem->Getenv("USER");
         }
      }
   } else {
      fUser = uid;
   }

   fOptions = options;

   if (gDebug > 1)
      Info("TAlien", "%s => %s port: %d user: %s",gridurl,fHost.Data(),fPort,fUser.Data());

   Bool_t fstoken;
   fstoken = kFALSE;

   if (options && (options[0] == 't')) {
       fGc = GapiUI::MakeGapiUI(kTRUE);
       fstoken = kTRUE;
   } else {
       fGc = GapiUI::MakeGapiUI(kFALSE);
   }

   if (!fGc) {
      Error("TAlien", "could not connect to a alien service at:");
      Error("TAlien", "host: %s port: %d user: %s", fHost.Data(), fPort, fUser.Data());
      MakeZombie();
   } else {
      if (passwd) {
         if (!passwd[0]) {
            passwd = 0;
         }
      }
      if (fstoken) {
         if (!Command("motd")) {
            Error("TAlien", "we have no valid connection ... I try to connect ...");
            fGc->Connect(fHost, fPort, fUser, passwd);
         }
      } else {
         fGc->Connect(fHost, fPort, fUser, passwd);
      }

      char cwd[1024];
      const char* cwdp=0;
      cwdp=Pwd();
      fHome=TString("/");
      if (cwdp) {
         sprintf(cwd,"%s",cwdp);
         if (Cd("")) {
            cwdp=Pwd();
            if (cwdp) {
               fHome = TString(cwdp);
            }
         }
      }

      if (!fGc->Connected()) {
         Error("TAlien", "could not authenticate at:");
         Error("TAlien", "host: %s port: %d user: %s",fHost.Data(),fPort,fUser.Data());
         MakeZombie();
      } else {
         gGrid = this;
         Command("motd");
         Stdout();
      }
      // export this UI to all POSIX functions like gapi_dir_xxx gapi_job_xxx
      fGc->SetGapiUI(fGc);
   }
   delete gurl;
}

////////////////////////////////////////////////////////////////////////////////
/// do we need to delete fGc ? (rdm)

TAlien::~TAlien()
{
   if (gDebug > 1)
      Info("~TAlien", "destructor called");
}

////////////////////////////////////////////////////////////////////////////////
/// Start an interactive AliEn shell.

void TAlien::Shell()
{
   fGc->Shell();
}

////////////////////////////////////////////////////////////////////////////////
/// Escape \" by \\\".

TString TAlien::Escape(const char *input)
{
   if (!input)
      return TString();

   TString output(input);
   output.ReplaceAll("\"", "\\\"");

   return output;
}

////////////////////////////////////////////////////////////////////////////////
/// Submit a command to AliEn. Returns 0 in case of error.

TGridJob *TAlien::Submit(const char *jdl)
{
   if (!jdl)
      return 0;

   TString command("submit =< ");
   //command += Escape(jdl);
   command += jdl;

   cout << command << endl;

   TGridResult* result = Command(command,kFALSE,kOUTPUT);
   TAlienResult* alienResult = dynamic_cast<TAlienResult*>(result);
   TList* list = dynamic_cast<TList*>(alienResult);
   if (!list) {
      if (result)
         delete result;
      return 0;
   }

   alienResult->DumpResult();

   TString jobID = "0";

   TIterator* iter = list->MakeIterator();
   TObject* object = 0;
   while ((object = iter->Next()) != 0) {
      TMap* map = dynamic_cast<TMap*>(object);

      TObject* jobIDObject = map->GetValue("jobId");
      TObjString* jobIDStr = dynamic_cast<TObjString*>(jobIDObject);
      if (jobIDStr) {
         jobID = jobIDStr->GetString();
      }
   }
   delete iter;
   delete result;

   if (jobID == "0") {
      Error("Submit", "error submitting job");
      return 0;
   }

   Info("Submit", "your job was submitted with the ID = %s", jobID.Data());

   return dynamic_cast<TGridJob*>(new TAlienJob(jobID));
}

////////////////////////////////////////////////////////////////////////////////
/// Create a Alien JDL.

TGridJDL *TAlien::GetJDLGenerator()
{
   return new TAlienJDL();
}

////////////////////////////////////////////////////////////////////////////////
/// Get job status list.

TGridJobStatusList *TAlien::Ps(const char* /*options*/, Bool_t verbose)
{
   GAPI_JOBARRAY* gjobarray = gapi_queryjobs("-", gGrid->GetUser(), "-", "-", "-",
                                             "-", "-", "-", "-");
   if (!gjobarray)
      return 0;

   if (gjobarray->size() == 0) {
      delete gjobarray;
      return 0;
   }

   TGridJobStatusList* joblist = new TAlienJobStatusList();

   for (UInt_t i = 0; i < gjobarray->size(); i++ ) {
      TMap* statusmap = new TMap();
      GAPI_JOB gjob = gjobarray->at(i);
      std::map<std::string, std::string>::const_iterator iter = gjob.gapi_jobmap.begin();
      for (; iter != gjob.gapi_jobmap.end(); ++iter) {
         //      printf("Adding %s = %s\n", iter->first.c_str(),iter->second.c_str());
         statusmap->Add(new TObjString(iter->first.c_str()), new TObjString(iter->second.c_str()));
      }
      TGridJobStatus *status = (TGridJobStatus*)new TAlienJobStatus(statusmap);
      joblist->Add((TGridJobStatus*)status);
      if (verbose) {
      }
   }
   delete gjobarray;
   return joblist;
}

////////////////////////////////////////////////////////////////////////////////
/// Kill a specific job.

Bool_t TAlien::KillById(TString jobid)
{
   if (gapi_kill(atoi(jobid.Data())))
      return kFALSE;
   else
      return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Resubmit a specific job.

Bool_t TAlien::ResubmitById(TString jobid)
{
  if (gapi_resubmit(atoi(jobid.Data())))
      return kFALSE;
   else
      return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute AliEn command. Returns 0 in case or error.

TGridResult *TAlien::Command(const char *command, bool interactive, UInt_t stream)
{
   if (fGc) {
      if (fGc->Command(command)) {
         // command successful
         TAlienResult* gresult = new TAlienResult();

         for (Int_t column = 0 ; column < (fGc->GetStreamColumns(stream)); column++) {
            TMap *gmap = new TMap();
            gmap->SetOwner(kTRUE);
            for (Int_t row=0; row < fGc->GetStreamRows(stream,column); row++) {
               gmap->Add((TObject*)(new TObjString(fGc->GetStreamFieldKey(stream,column,row))),
                         (TObject*)(new TObjString(fGc->GetStreamFieldValue(stream,column,row))));
            }
            gresult->Add(gmap);
         }

         if (interactive) {
            // write also stdout/stderr to the screen
            fGc->DebugDumpStreams();
         }
         return gresult;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

TGridResult *TAlien::LocateSites()
{
   return Command("locatesites");
}

////////////////////////////////////////////////////////////////////////////////

void TAlien::Stdout()
{
   if (fGc) {
      fGc->PrintCommandStdout();
   }
}

////////////////////////////////////////////////////////////////////////////////

void TAlien::Stderr()
{
   if (fGc) {
      fGc->PrintCommandStderr();
   }
}

////////////////////////////////////////////////////////////////////////////////

TGridResult *TAlien::Query(const char *path, const char *pattern,
                           const char *conditions, const char *options)
{
   TString cmdline = TString("find -z ") + TString(options) + TString(" ") + TString(path) + TString(" ")  + TString(pattern) + TString(" ") + TString(conditions);
   return Command(cmdline);
}

////////////////////////////////////////////////////////////////////////////////

TGridResult *TAlien::OpenDataset(const char *lfn, const char *options)
{
   TString cmdline = TString("getdataset") + TString(" ") + TString(options) + TString(" ") + TString(lfn);
   return Command(cmdline, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////

TMap *TAlien::GetColumn(UInt_t stream, UInt_t column)
{
   TMap *gmap = new TMap();
   for (Int_t row = 0; row < fGc->GetStreamRows(stream,column); row++) {
      gmap->Add((TObject*)(new TObjString(fGc->GetStreamFieldKey(stream,column,row))),
                (TObject*)(new TObjString(fGc->GetStreamFieldValue(stream,column,row))));
   }
   return gmap;
}

////////////////////////////////////////////////////////////////////////////////

const char *TAlien::GetStreamFieldValue(UInt_t stream, UInt_t column, UInt_t row)
{
   return fGc->GetStreamFieldValue(stream,column,row);
}

////////////////////////////////////////////////////////////////////////////////

const char *TAlien::GetStreamFieldKey(UInt_t stream, UInt_t column, UInt_t row)
{
   return fGc->GetStreamFieldKey(stream,column,row);
}

////////////////////////////////////////////////////////////////////////////////

UInt_t TAlien::GetNColumns(UInt_t stream)
{
   return fGc->GetStreamColumns(stream);
}

//--- Catalogue Interface

////////////////////////////////////////////////////////////////////////////////

TGridResult *TAlien::Ls(const char* ldn, Option_t* options, Bool_t verbose)
{
   TString cmdline = TString("ls") + TString(" ") + TString(options) + TString(" ") + TString(ldn);

   return Command(cmdline, verbose);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TAlien::Cd(const char* ldn, Bool_t verbose)
{
   TString cmdline = TString("cd") + TString(" ") + TString(ldn);

   Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   const char* result = (GetStreamFieldValue(kOUTPUT,0,0));
   if (result) {
      if (strlen(result) > 0) {
         if (atoi(result) == 1) {
            return kTRUE;
         }
      }
   }

   Error("Cd","Cannot change to directory %s\n",ldn);
   if (!verbose) Stdout();
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

const char* TAlien::Pwd(Bool_t verbose)
{
   TString cmdline = TString("pwd");

   TGridResult* result = Command(cmdline, kFALSE, kENVIR);

   if (verbose) {
      Stdout();
      Stderr();
   }

   if (result) {
      TMap* resultmap = ((TMap*)result->At(0));
      if (resultmap) {
         TObjString* pwd = (TObjString*)resultmap->GetValue("pwd");
         if (pwd) {
            fPwd = pwd->GetName();
            delete resultmap;
            return fPwd;
         } else {
            delete resultmap;
            return 0;
         }
      }
   }

   Error("pwd","Cannot get current working directory\n");
   if (!verbose) Stdout();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TAlien::Mkdir(const char* ldn, Option_t* options, Bool_t verbose)
{
   TString cmdline = TString("mkdir");
   if (strlen(options)) {
      cmdline +=  (TString(" ") + TString(options));
   }

   cmdline += (TString(" -s ") + TString(ldn));

   Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   const char* result = (GetStreamFieldValue(kOUTPUT,0,0));
   if (result) {
      if (strlen(result) > 0) {
         if (atoi(result) > 0) {
            return atoi(result);
         }
      }
   }

   Error("Mkdir","Cannot create directory %s\n",ldn);
   if (!verbose) Stdout();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TAlien::Rmdir(const char* ldn, Option_t* options, Bool_t verbose)
{
   TString cmdline = TString("rmdir");
   if (strlen(options)) {
      cmdline +=  (TString(" ") + TString(options));
   } else {
      cmdline += (TString(" -s ") + TString(ldn));
   }

   Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   const char* result = (GetStreamFieldValue(kOUTPUT,0,0));
   if (result) {
      if (strlen(result) > 0) {
         if (atoi(result) == 1) {
            return kTRUE;
         }
      }
   }

   Error("Rmdir","Cannot remove directory %s\n",ldn);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TAlien::Register(const char* lfn, const char* turl, Long_t size, const char* se, const char* guid, Bool_t verbose)
{
   TString cmdline = TString("register") + TString(" ") + TString(lfn) + TString(" ") + TString(turl);
   if (se) {
      cmdline += TString(" ");
      cmdline += size;
      cmdline += TString(" ");
      cmdline += TString(se);
      if (guid) {
         cmdline += TString(" ");
         cmdline += TString(guid);
      }
   }

   Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   const char* result = (GetStreamFieldValue(kOUTPUT,0,0));
   if (result) {
      if (strlen(result) > 0) {
         if (atoi(result) == 1) {
            return kTRUE;
         }
      }
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TAlien::Rm(const char* lfn, Option_t* options, Bool_t verbose)
{
   TString cmdline = TString("rm") + TString(" -s ") + TString(options) + TString(" ") + TString(lfn);

   Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   const char* result = (GetStreamFieldValue(kOUTPUT,0,0));
   if (result) {
      if (strlen(result) > 0) {
         if (atoi(result) == 1) {
            return kTRUE;
         }
      }
   }

   Error("Rm","Cannot remove %s\n",lfn);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Factory function for a TAlienCollection based on an XML file.

TGridCollection *TAlien::OpenCollection(const char *collectionfile,
                                        UInt_t maxentries)
{
   TString path(collectionfile);
   if (path.BeginsWith("alien://", TString::kIgnoreCase)) {
      TAlien* alien = dynamic_cast<TAlien*> (gGrid);
      if (!alien) {
        Error("OpenCollection",
              "Trying to read a collection, but gGrid not initialized with AliEn");
        return 0;
      }
      TString lfn = path(strlen("alien://"), path.Length());
      if (alien->Type(lfn) == kCollection) {
        // it is a collection
        TGridResult* gridResult = alien->GetCollection(lfn);
        if (!gridResult) {
          Error("OpenCollection",
                "Could not retrieve collection %s from the catalog", collectionfile);
          return 0;
        }

        return TAlienCollection::OpenAlienCollection(gridResult);
      }
   }

   return TAlienCollection::Open(collectionfile, maxentries);
}

////////////////////////////////////////////////////////////////////////////////
/// Factory function fo a TAlienCollection based on a gGrid Query.

TGridCollection *TAlien::OpenCollectionQuery(TGridResult *queryresult,
                                             Bool_t nogrouping)
{
   return (TGridCollection*)TAlienCollection::OpenQuery(queryresult, nogrouping);
}

////////////////////////////////////////////////////////////////////////////////
/// returns the type of the given lfn

TAlien::CatalogType TAlien::Type(const char* lfn, Option_t* option, Bool_t verbose)
{
   TString cmdline = TString("type -z") + TString(" ") + TString(option) + TString(" ") + TString(lfn);

   TGridResult* gridResult = Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   if (!gridResult) {
      Error("Type", "Did not receive TGridResult from query %s", cmdline.Data());
      return kFailed;
   }

   const char* typeStr = gridResult->GetKey(0, "type");
   if (!typeStr || !typeStr[0]) {
      Error("Type", "Could not get type of %s", lfn);
      delete gridResult;
      return kFailed;
   }

   TAlien::CatalogType type = kFailed;

   if (strcmp(typeStr, "file") == 0) {
      type =  kFile;
   } else if (strcmp(typeStr, "directory") == 0) {
      type = kDirectory;
   } else if (strcmp(typeStr, "collection") == 0) {
      type = kCollection;
   } else
     Error("Type", "Unknown type %s", typeStr);

   delete gridResult;

   return type;
}

////////////////////////////////////////////////////////////////////////////////
///

TGridResult* TAlien::GetCollection(const char* lfn, Option_t* option, Bool_t verbose)
{
   TString cmdline = TString("listFilesFromCollection -z -v") + TString(" ") + TString(option) + TString(" ") + TString(lfn);

   TGridResult* gridResult = Command(cmdline, kFALSE);

   if (verbose) {
      Stdout();
      Stderr();
   }

   return gridResult;
}

////////////////////////////////////////////////////////////////////////////////
/// List packages in the specified directory.

TGridResult* TAlien::ListPackages(const char* alienpackagedir)
{
   if (!alienpackagedir) {
      alienpackagedir = "/alice/packages";
   }

   TGridResult* gr = (TGridResult*) new TAlienResult();

   TGridResult* result = Ls(alienpackagedir);
   if (result) {
      Int_t i=0;
      while ( result->GetFileName(i)) {
         TString pname=result->GetFileName(i);
         TGridResult* version = Ls(Form("%s/%s",alienpackagedir,pname.Data()));
         if (version) {
            Int_t j=0;
            while (version->GetFileName(j)) {
               TString pversion=version->GetFileName(j);
               if ((pversion.Contains("post_"))) {
               } else {
                  TGridResult* platform = Ls(Form("%s/%s/%s",alienpackagedir,pname.Data(),pversion.Data()));
                  if (platform) {
                     Int_t k=0;
                     TString allplatform="";
                     while (platform->GetFileName(k)) {
                        TString pplatform=platform->GetFileName(k);
                        allplatform+=pplatform;
                        allplatform+=" ";
                        TMap* grmap = new TMap();
                        grmap->Add((TObject*) new TObjString("name"),
                           (TObject*) new TObjString(pplatform.Data()));
                        grmap->Add((TObject*) new TObjString("path"),
                           new TObjString( Form ( "%s/%s/%s/%s" , alienpackagedir,pname.Data(),pversion.Data(),pplatform.Data())));
                        gr->Add(grmap);
                        k++;
                     }
                     Info("ListPackages","Package: %-16s Version: %-20s Platform:  [ %s ]",pname.Data(),pversion.Data(),allplatform.Data());
                     delete platform;
                  }
               }
               j++;
            }
            delete version;
         }
         i++;
      }
      delete result;
   }
   return gr;
}
