// @(#)root/netxng:$Id$
/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TNetXNGFile                                                                //
//                                                                            //
// Authors: Justin Salmon, Lukasz Janyst                                      //
//          CERN, 2013                                                        //
//                                                                            //
// Enables access to XRootD files using the new client.                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TArchiveFile.h"
#include "TNetXNGFile.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TVirtualPerfStats.h"
#include "TVirtualMonitoring.h"
#include <XrdCl/XrdClURL.hh>
#include <XrdCl/XrdClFile.hh>
#include <XrdCl/XrdClXRootDResponses.hh>
#include <XrdCl/XrdClDefaultEnv.hh>
#include <XrdVersion.hh>
#include <iostream>

//------------------------------------------------------------------------------
// Open handler for async open requests
////////////////////////////////////////////////////////////////////////////////

class TAsyncOpenHandler: public XrdCl::ResponseHandler
{
   public:
      //------------------------------------------------------------------------
      // Constructor
      //////////////////////////////////////////////////////////////////////////

      TAsyncOpenHandler(TNetXNGFile *file)
      {
         fFile = file;
         fFile->SetAsyncOpenStatus(TFile::kAOSInProgress);
      }

      //------------------------------------------------------------------------
      // Called when a response to open arrives
      //////////////////////////////////////////////////////////////////////////

      virtual void HandleResponse(XrdCl::XRootDStatus *status,
                                  XrdCl::AnyObject    *response)
      {
         if (status->IsOK())
         {
            fFile->SetAsyncOpenStatus(TFile::kAOSSuccess);
         }
         else
         {
            fFile->SetAsyncOpenStatus(TFile::kAOSFailure);
         }

         delete response;
         delete status;
         delete this;
      }

   private:
      TNetXNGFile *fFile;
};

//------------------------------------------------------------------------------
// Async readv handler
////////////////////////////////////////////////////////////////////////////////

class TAsyncReadvHandler: public XrdCl::ResponseHandler
{
   public:
      //------------------------------------------------------------------------
      // Constructor
      //////////////////////////////////////////////////////////////////////////

      TAsyncReadvHandler(std::vector<XrdCl::XRootDStatus*> *statuses,
                         Int_t                              statusIndex,
                         TSemaphore                        *semaphore):
         fStatuses(statuses), fStatusIndex(statusIndex), fSemaphore(semaphore) {}


      //------------------------------------------------------------------------
      // Handle readv response
      //////////////////////////////////////////////////////////////////////////

      virtual void HandleResponse(XrdCl::XRootDStatus *status,
                                  XrdCl::AnyObject    *response)
      {
         fStatuses->at(fStatusIndex) = status;
         fSemaphore->Post();
         delete response;
         delete this;
      }

   private:
      std::vector<XrdCl::XRootDStatus*> *fStatuses;    // Pointer to status vector
      Int_t                              fStatusIndex; // Index into status vector
      TSemaphore                        *fSemaphore;   // Synchronize the responses
};


ClassImp(TNetXNGFile);

////////////////////////////////////////////////////////////////////////////////
/// Constructor
///
/// param url:          URL of the entry-point server to be contacted
/// param mode:         initial file access mode
/// param title:        title of the file (shown by ROOT browser)
/// param compress:     compression level and algorithm
/// param netopt:       TCP window size in bytes (unused)
/// param parallelopen: open asynchronously

TNetXNGFile::TNetXNGFile(const char *url,
                         Option_t   *mode,
                         const char *title,
                         Int_t       compress,
                         Int_t       netopt,
                         Bool_t      parallelopen) :
	TNetXNGFile(url,0,mode,title,compress,netopt,parallelopen){}

TNetXNGFile::TNetXNGFile(const char *url,
		         const char *lurl,
                         Option_t   *mode,
                         const char *title,
                         Int_t       compress,
                         Int_t       /*netopt*/,
                         Bool_t      parallelopen) :
   TFile((lurl ? lurl : url), "NET", title, compress)
{
   using namespace XrdCl;

   // Set the log level
   TString val = gSystem->Getenv("XRD_LOGLEVEL");
   if (val.IsNull()) val = gEnv->GetValue("NetXNG.Debug",     "");
   if (!val.IsNull()) XrdCl::DefaultEnv::SetLogLevel(val.Data());

   // Remove any anchor from the url. It may have been used by the base TFile
   // constructor to setup a TArchiveFile but we should not pass it to the xroot
   // client as a part of the filename
   {
     TUrl urlnoanchor(url);
     urlnoanchor.SetAnchor("");
     fUrl = new URL(std::string(urlnoanchor.GetUrl()));
   }

   fFile        = new File();
   fInitCondVar = new XrdSysCondVar();
   fUrl->SetProtocol(std::string("root"));
   fQueryReadVParams = 1;
   fReadvIorMax = 2097136;
   fReadvIovMax = 1024;

   if (ParseOpenMode(mode, fOption, fMode, kTRUE)<0) {
      Error("Open", "could not parse open mode %s", mode);
      MakeZombie();
      return;
   }

   // Map ROOT and xrootd environment
   SetEnv();

   // Init the monitoring system
   if (gMonitoringWriter) {
      if (!fOpenPhases) {
         fOpenPhases = new TList;
         fOpenPhases->SetOwner();
      }
      gMonitoringWriter->SendFileOpenProgress(this, fOpenPhases, "xrdopen",
                                              kFALSE);
   }

   XRootDStatus status;
   if (parallelopen) {
      // Open the file asynchronously
      TAsyncOpenHandler *handler = new TAsyncOpenHandler(this);
      status = fFile->Open(fUrl->GetURL(), fMode, Access::None, handler);
      if (!status.IsOK()) {
         Error("Open", "%s", status.ToStr().c_str());
         MakeZombie();
      }
      return;
   }

   // Open the file synchronously
   status = fFile->Open(fUrl->GetURL(), fMode);
   if (!status.IsOK()) {
#if XrdVNUMBER >= 40000
      if( status.code == errRedirect )
         fNewUrl = status.GetErrorMessage().c_str();
      else
         Error("Open", "%s", status.ToStr().c_str());
#else
      Error("Open", "%s", status.ToStr().c_str());
#endif
      MakeZombie();
      return;
   }

   if( (fMode & OpenFlags::New) || (fMode & OpenFlags::Delete) ||
       (fMode & OpenFlags::Update) )
      fWritable = true;

   // Initialize the file
   bool create = false;
   if( (fMode & OpenFlags::New) || (fMode & OpenFlags::Delete) )
      create = true;
   TFile::Init(create);

   // Get the vector read limits
   GetVectorReadLimits();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TNetXNGFile::~TNetXNGFile()
{
   if (IsOpen())
      Close();
   delete fFile;
   delete fUrl;
   delete fInitCondVar;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the file. Makes sure that the file is really open before
/// calling TFile::Init. It may block.

void TNetXNGFile::Init(Bool_t create)
{
   using namespace XrdCl;

   if (fInitDone) {
      if (gDebug > 1) Info("Init", "TFile::Init already called once");
      return;
   }

   // If the async open didn't return yet, wait for it
   if (!IsOpen() && fAsyncOpenStatus == kAOSInProgress) {
      fInitCondVar->Wait();
   }

   // Notify the monitoring system
   if (gMonitoringWriter)
      gMonitoringWriter->SendFileOpenProgress(this, fOpenPhases, "rootinit",
                                              kFALSE);

   // Initialize the file
   TFile::Init(create);

   // Notify the monitoring system
   if (gMonitoringWriter)
      gMonitoringWriter->SendFileOpenProgress(this, fOpenPhases, "endopen",
                                              kTRUE);

   // Get the vector read limits
   GetVectorReadLimits();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the file size. Returns -1 in the case that the file could not be
/// stat'ed.

Long64_t TNetXNGFile::GetSize() const
{
   if (fArchive && fArchive->GetMember()) {
     return fArchive->GetMember()->GetDecompressedSize();
   }

   using namespace XrdCl;

   // Check the file isn't a zombie or closed
   if (!IsUseable())
      return -1;

   bool forceStat = true;
   if( fMode == XrdCl::OpenFlags::Read )
      forceStat = false;

   StatInfo *info = 0;
   if( !fFile->Stat( forceStat, info ).IsOK() )
    return -1;
   Long64_t size = info->GetSize();
   delete info;
   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the file is open

Bool_t TNetXNGFile::IsOpen() const
{
   return fFile->IsOpen();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the status of an asynchronous file open

void TNetXNGFile::SetAsyncOpenStatus(EAsyncOpenStatus status)
{
   fAsyncOpenStatus = status;
   // Unblock Init() if it is waiting
   fInitCondVar->Signal();
}

////////////////////////////////////////////////////////////////////////////////
/// Close the file
///
/// param option: if == "R", all TProcessIDs referenced by this file are
///               deleted (is this valid in xrootd context?)

void TNetXNGFile::Close(const Option_t */*option*/)
{
   TFile::Close();

   XrdCl::XRootDStatus status = fFile->Close();
   if (!status.IsOK()) {
      Error("Close", "%s", status.ToStr().c_str());
      MakeZombie();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reopen the file with the new access mode
///
/// param mode: the new access mode
/// returns:    0 in case the mode was successfully modified, 1 in case
///             the mode did not change (was already as requested or wrong
///             input arguments) and -1 in case of failure, in which case
///             the file cannot be used anymore

Int_t TNetXNGFile::ReOpen(Option_t *modestr)
{
   using namespace XrdCl;
   TString newOpt;
   OpenFlags::Flags mode;

   Int_t parseres = ParseOpenMode(modestr, newOpt, mode, kFALSE);

   // Only Read and Update are valid modes
   if (parseres<0 || (mode != OpenFlags::Read && mode != OpenFlags::Update)) {
      Error("ReOpen", "mode must be either READ or UPDATE, not %s", modestr);
      return 1;
   }

   // The mode is not really changing
   if (mode == fMode || (mode == OpenFlags::Update
                         && fMode == OpenFlags::New)) {
      return 1;
   }

   XRootDStatus st = fFile->Close();
   if (!st.IsOK()) {
      Error("ReOpen", "%s", st.ToStr().c_str());
      return 1;
   }
   fOption = newOpt;
   fMode = mode;

   st = fFile->Open(fUrl->GetURL(), fMode);
   if (!st.IsOK()) {
      Error("ReOpen", "%s", st.ToStr().c_str());
      return 1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read a data chunk of the given size
///
/// param buffer: a pointer to a buffer big enough to hold the data
/// param length: number of bytes to be read
/// returns:      kTRUE in case of failure

Bool_t TNetXNGFile::ReadBuffer(char *buffer, Int_t length)
{
   return ReadBuffer(buffer, GetRelOffset(), length);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a data chunk of the given size, starting from the given offset
///
/// param buffer:   a pointer to a buffer big enough to hold the data
/// param position: offset from the beginning of the file
/// param length:   number of bytes to be read
/// returns:        kTRUE in case of failure

Bool_t TNetXNGFile::ReadBuffer(char *buffer, Long64_t position, Int_t length)
{
   using namespace XrdCl;
   if (gDebug > 0)
      Info("ReadBuffer", "offset: %lld length: %d", position, length);

   // Check the file isn't a zombie or closed
   if (!IsUseable())
      return kTRUE;

   // Try to read from cache
   SetOffset(position);
   Int_t status;
   if ((status = ReadBufferViaCache(buffer, length))) {
      if (status == 2)
         return kTRUE;
      return kFALSE;
   }

   Double_t start = 0;
   if (gPerfStats) start = TTimeStamp();

   // Read the data
   uint32_t bytesRead = 0;
   XRootDStatus st = fFile->Read(fOffset, length, buffer, bytesRead);
   if (gDebug > 0)
      Info("ReadBuffer", "%s bytes read: %u", st.ToStr().c_str(), bytesRead);

   if (!st.IsOK()) {
      Error("ReadBuffer", "%s", st.ToStr().c_str());
      return kTRUE;
   }

   if ((Int_t)bytesRead != length) {
      Error("ReadBuffer", "error reading all requested bytes, got %u of %d",
            bytesRead, length);
      return kTRUE;
   }

   // Bump the globals
   fOffset     += bytesRead;
   fBytesRead  += bytesRead;
   fgBytesRead += bytesRead;
   fReadCalls  ++;
   fgReadCalls ++;

   if (gPerfStats)
      gPerfStats->FileReadEvent(this, (Int_t)bytesRead, start);

   if (gMonitoringWriter)
      gMonitoringWriter->SendFileReadProgress(this);

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Read scattered data chunks in one operation
///
/// param buffer:   a pointer to a buffer big enough to hold all of the
///                 requested data
/// param position: position[i] is the seek position of chunk i of len
///                 length[i]
/// param length:   length[i] is the length of the chunk at offset
///                 position[i]
/// param nbuffs:   number of chunks
/// returns:        kTRUE in case of failure

Bool_t TNetXNGFile::ReadBuffers(char *buffer, Long64_t *position, Int_t *length,
      Int_t nbuffs)
{
   using namespace XrdCl;

   // Check the file isn't a zombie or closed
   if (!IsUseable())
      return kTRUE;

   std::vector<ChunkList>      chunkLists;
   ChunkList                   chunks;
   std::vector<XRootDStatus*> *statuses;
   TSemaphore                 *semaphore;
   Int_t                       totalBytes = 0;
   Long64_t                    offset     = 0;
   char                       *cursor     = buffer;

   Double_t start = 0;
   if (gPerfStats) start = TTimeStamp();

   if (fArchiveOffset)
      for (Int_t i = 0; i < nbuffs; i++)
         position[i] += fArchiveOffset;

   // Build a list of chunks. Put the buffers in the ChunkInfo's
   for (Int_t i = 0; i < nbuffs; ++i) {
      totalBytes += length[i];

      // If the length is bigger than max readv size, split into smaller chunks
      if (length[i] > fReadvIorMax) {
         Int_t nsplit = length[i] / fReadvIorMax;
         Int_t rem    = length[i] % fReadvIorMax;
         Int_t j;

         // Add as many max-size chunks as are divisible
         for (j = 0; j < nsplit; ++j) {
            offset = position[i] + (j * fReadvIorMax);
            chunks.push_back(ChunkInfo(offset, fReadvIorMax, cursor));
            cursor += fReadvIorMax;
         }

         // Add the remainder
         offset = position[i] + (j * fReadvIorMax);
         chunks.push_back(ChunkInfo(offset, rem, cursor));
         cursor += rem;
      } else {
         chunks.push_back(ChunkInfo(position[i], length[i], cursor));
         cursor += length[i];
      }

      // If there are more than or equal to max chunks, make another chunk list
      if ((Int_t) chunks.size() == fReadvIovMax) {
         chunkLists.push_back(chunks);
         chunks = ChunkList();
      } else if ((Int_t) chunks.size() > fReadvIovMax) {
         chunkLists.push_back(ChunkList(chunks.begin(),
                                        chunks.begin() + fReadvIovMax));
         chunks = ChunkList(chunks.begin() + fReadvIovMax, chunks.end());
      }
   }

   // Push back the last chunk list
   if( !chunks.empty() )
      chunkLists.push_back(chunks);

   TAsyncReadvHandler *handler;
   XRootDStatus        status;
   semaphore = new TSemaphore(0);
   statuses  = new std::vector<XRootDStatus*>(chunkLists.size());

   // Read asynchronously but wait for all responses
   std::vector<ChunkList>::iterator it;
   for (it = chunkLists.begin(); it != chunkLists.end(); ++it)
   {
      handler = new TAsyncReadvHandler(statuses, it - chunkLists.begin(),
                                       semaphore);
      status = fFile->VectorRead(*it, 0, handler);

      if (!status.IsOK()) {
         Error("ReadBuffers", "%s", status.ToStr().c_str());
         return kTRUE;
      }
   }

   // Wait for all responses
   for (it = chunkLists.begin(); it != chunkLists.end(); ++it) {
      semaphore->Wait();
   }

   // Check for errors
   for (it = chunkLists.begin(); it != chunkLists.end(); ++it) {
      XRootDStatus *st = statuses->at(it - chunkLists.begin());

      if (!st->IsOK()) {
         Error("ReadBuffers", "%s", st->ToStr().c_str());
         for( ; it != chunkLists.end(); ++it )
         {
            st = statuses->at( it - chunkLists.begin() );
            delete st;
         }
         delete statuses;
         delete semaphore;

         return kTRUE;
      }
      delete st;
   }

   // Bump the globals
   fBytesRead  += totalBytes;
   fgBytesRead += totalBytes;
   fReadCalls  ++;
   fgReadCalls ++;

   if (gPerfStats) {
      fOffset = position[0];
      gPerfStats->FileReadEvent(this, totalBytes, start);
   }

   if (gMonitoringWriter)
      gMonitoringWriter->SendFileReadProgress(this);

   delete statuses;
   delete semaphore;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a data chunk
///
/// param buffer: the data to be written
/// param length: the size of the buffer
/// returns:      kTRUE in case of failure

Bool_t TNetXNGFile::WriteBuffer(const char *buffer, Int_t length)
{
   using namespace XrdCl;

   // Check the file isn't a zombie or closed
   if (!IsUseable())
      return kTRUE;

   if (!fWritable) {
      if (gDebug > 1)
         Info("WriteBuffer", "file not writable");
      return kTRUE;
   }

   // Check the write cache
   Int_t status;
   if ((status = WriteBufferViaCache(buffer, length))) {
      if (status == 2)
         return kTRUE;
      return kFALSE;
   }

   // Write the data
   XRootDStatus st = fFile->Write(fOffset, length, buffer);
   if (!st.IsOK()) {
      Error("WriteBuffer", "%s", st.ToStr().c_str());
      return kTRUE;
   }

   // Bump the globals
   fOffset      += length;
   fBytesWrite  += length;
   fgBytesWrite += length;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

void TNetXNGFile::Flush()
{
   if (!IsUseable())
      return;

   if (!fWritable) {
      if (gDebug > 1)
         Info("Flush", "file not writable - do nothing");
      return;
   }

   FlushWriteCache();

   //
   // Flush via the remote xrootd
   XrdCl::XRootDStatus status = fFile->Sync();
   if( !status.IsOK() )
      Error("Flush", "%s", status.ToStr().c_str());

   if (gDebug > 1)
      Info("Flush", "XrdClient::Sync succeeded.");
}

////////////////////////////////////////////////////////////////////////////////
/// Set the position within the file
///
/// param offset:   the new offset relative to position
/// param position: the relative position, either kBeg, kCur or kEnd

void TNetXNGFile::Seek(Long64_t offset, ERelativeTo position)
{
   SetOffset(offset, position);
}

////////////////////////////////////////////////////////////////////////////////
/// Parse a file open mode given as a string into a canonically formatted
/// output mode string and an integer code that the xroot client can use
///
/// param    in:         the file open mode as a string (in)
///          modestr:    open mode string after parsing (out)
///          mode:       correctly parsed option mode code (out)
///          assumeRead: if the open mode is not recognised assume read (in)
/// returns:             0 in case the mode was successfully parsed,
///                     -1 in case of failure

Int_t TNetXNGFile::ParseOpenMode(Option_t *in, TString &modestr,
                                 XrdCl::OpenFlags::Flags &mode,
                                 Bool_t assumeRead)
{
   using namespace XrdCl;
   modestr = ToUpper(TString(in));

   if (modestr == "NEW" || modestr == "CREATE")  mode = OpenFlags::New;
   else if (modestr == "RECREATE")               mode = OpenFlags::Delete;
   else if (modestr == "UPDATE")                 mode = OpenFlags::Update;
   else if (modestr == "READ")                   mode = OpenFlags::Read;
   else {
      if (!assumeRead) {
         return -1;
      }
      modestr = "READ";
      mode = OpenFlags::Read;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the file is open and isn't a zombie

Bool_t TNetXNGFile::IsUseable() const
{
   if (IsZombie()) {
      Error("TNetXNGFile", "Object is in 'zombie' state");
      return kFALSE;
   }

   if (!IsOpen()) {
      Error("TNetXNGFile", "The remote file is not open");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the server-specific readv config params. Returns kFALSE in case of
/// error, kTRUE otherwise.

Bool_t TNetXNGFile::GetVectorReadLimits()
{
   using namespace XrdCl;

   // Check the file isn't a zombie or closed
   if (!IsUseable())
      return kFALSE;

   if (!fQueryReadVParams)
      return kTRUE;

#if XrdVNUMBER >= 40000
   std::string lasturl;
   fFile->GetProperty("LastURL",lasturl);
   URL lrl(lasturl);
   //local redirect will split vector reads into multiple local reads anyway,
   // so we are fine with the default values
   if(lrl.GetProtocol().compare("file") == 0 &&
      lrl.GetHostId().compare("localhost") == 0){
       if (gDebug >= 1)
          Info("GetVectorReadLimits","Local redirect, using default values");
       return kTRUE;
   }

   std::string dataServerStr;
   if( !fFile->GetProperty( "DataServer", dataServerStr ) )
      return kFALSE;
   URL dataServer(dataServerStr);
#else
   URL dataServer(fFile->GetDataServer());
#endif
   FileSystem fs(dataServer);
   Buffer  arg;
   Buffer *response;
   arg.FromString(std::string("readv_ior_max readv_iov_max"));

   XRootDStatus status = fs.Query(QueryCode::Config, arg, response);
   if (!status.IsOK())
      return kFALSE;

   Ssiz_t from = 0;
   TString token;

   std::vector<TString> resps;
   while (TString(response->ToString()).Tokenize(token, from, "\n"))
      resps.push_back(token);

   if (resps.size() != 2)
      return kFALSE;

   if (resps[0].IsDigit())
      fReadvIorMax = resps[0].Atoi();

   if (resps[1].IsDigit())
      fReadvIovMax = resps[1].Atoi();

   delete response;

   // this is to workaround a dCache bug reported here:
   // https://sft.its.cern.ch/jira/browse/ROOT-6639
   if( fReadvIovMax == 0x7FFFFFFF )
   {
     fReadvIovMax = 1024;
     fReadvIorMax = 2097136;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Map ROOT and xrootd environment variables

void TNetXNGFile::SetEnv()
{
   XrdCl::Env *env  = XrdCl::DefaultEnv::GetEnv();
   const char *cenv = 0;
   TString     val;

   val = gEnv->GetValue("NetXNG.ConnectionWindow",     "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_CONNECTIONWINDOW"))
                            || strlen(cenv) <= 0))
      env->PutInt("ConnectionWindow", val.Atoi());

   val = gEnv->GetValue("NetXNG.ConnectionRetry",      "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_CONNECTIONRETRY"))
                            || strlen(cenv) <= 0))
      env->PutInt("RequestTimeout", val.Atoi());

   val = gEnv->GetValue("NetXNG.RequestTimeout",       "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_REQUESTTIMEOUT"))
                            || strlen(cenv) <= 0))
      env->PutInt("RequestTimeout", val.Atoi());

   val = gEnv->GetValue("NetXNG.SubStreamsPerChannel", "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_SUBSTREAMSPERCHANNEL"))
                            || strlen(cenv) <= 0))
      env->PutInt("SubStreamsPerChannel", val.Atoi());

   val = gEnv->GetValue("NetXNG.TimeoutResolution",    "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_TIMEOUTRESOLUTION"))
                            || strlen(cenv) <= 0))
      env->PutInt("TimeoutResolution", val.Atoi());

   val = gEnv->GetValue("NetXNG.StreamErrorWindow",    "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_STREAMERRORWINDOW"))
                            || strlen(cenv) <= 0))
      env->PutInt("StreamErrorWindow", val.Atoi());

   val = gEnv->GetValue("NetXNG.RunForkHandler",       "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_RUNFORKHANDLER"))
                            || strlen(cenv) <= 0))
      env->PutInt("RunForkHandler", val.Atoi());

   val = gEnv->GetValue("NetXNG.RedirectLimit",        "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_REDIRECTLIMIT"))
                            || strlen(cenv) <= 0))
      env->PutInt("RedirectLimit", val.Atoi());

   val = gEnv->GetValue("NetXNG.WorkerThreads",        "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_WORKERTHREADS"))
                            || strlen(cenv) <= 0))
      env->PutInt("WorkerThreads", val.Atoi());

   val = gEnv->GetValue("NetXNG.CPChunkSize",          "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_CPCHUNKSIZE"))
                            || strlen(cenv) <= 0))
      env->PutInt("CPChunkSize", val.Atoi());

   val = gEnv->GetValue("NetXNG.CPParallelChunks",     "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_CPPARALLELCHUNKS"))
                            || strlen(cenv) <= 0))
      env->PutInt("CPParallelChunks", val.Atoi());

   val = gEnv->GetValue("NetXNG.PollerPreference",     "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_POLLERPREFERENCE"))
                            || strlen(cenv) <= 0))
      env->PutString("PollerPreference", val.Data());

   val = gEnv->GetValue("NetXNG.ClientMonitor",        "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_CLIENTMONITOR"))
                            || strlen(cenv) <= 0))
      env->PutString("ClientMonitor", val.Data());

   val = gEnv->GetValue("NetXNG.ClientMonitorParam",   "");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XRD_CLIENTMONITORPARAM"))
                            || strlen(cenv) <= 0))
      env->PutString("ClientMonitorParam", val.Data());

   fQueryReadVParams = gEnv->GetValue("NetXNG.QueryReadVParams", 1);
   env->PutInt( "MultiProtocol", gEnv->GetValue("TFile.CrossProtocolRedirects", 1));

   // Old style netrc file
   TString netrc;
   netrc.Form("%s/.rootnetrc", gSystem->HomeDirectory());
   gSystem->Setenv("XrdSecNETRC", netrc.Data());

   // For authentication
   val = gEnv->GetValue("XSec.Pwd.ALogFile",     "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecPWDALOGFILE",     val.Data());

   val = gEnv->GetValue("XSec.Pwd.ServerPuk",    "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecPWDSRVPUK",       val.Data());

   val = gEnv->GetValue("XSec.GSI.CAdir",        "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSICADIR",        val.Data());

   val = gEnv->GetValue("XSec.GSI.CRLdir",       "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSICRLDIR",       val.Data());

   val = gEnv->GetValue("XSec.GSI.CRLextension", "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSICRLEXT",       val.Data());

   val = gEnv->GetValue("XSec.GSI.UserCert",     "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSIUSERCERT",     val.Data());

   val = gEnv->GetValue("XSec.GSI.UserKey",      "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSIUSERKEY",      val.Data());

   val = gEnv->GetValue("XSec.GSI.UserProxy",    "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSIUSERPROXY",    val.Data());

   val = gEnv->GetValue("XSec.GSI.ProxyValid",   "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSIPROXYVALID",   val.Data());

   val = gEnv->GetValue("XSec.GSI.ProxyKeyBits", "");
   if (val.Length() > 0)
      gSystem->Setenv("XrdSecGSIPROXYKEYBITS", val.Data());

   val = gEnv->GetValue("XSec.GSI.ProxyForward", "0");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XrdSecGSIPROXYDEPLEN"))
                            || strlen(cenv) <= 0))
      gSystem->Setenv("XrdSecGSIPROXYDEPLEN",  val.Data());

   val = gEnv->GetValue("XSec.GSI.CheckCRL",     "1");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XrdSecGSICRLCHECK"))
                            || strlen(cenv) <= 0))
      gSystem->Setenv("XrdSecGSICRLCHECK",     val.Data());

   val = gEnv->GetValue("XSec.GSI.DelegProxy",   "0");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XrdSecGSIDELEGPROXY"))
                            || strlen(cenv) <= 0))
      gSystem->Setenv("XrdSecGSIDELEGPROXY",   val.Data());

   val = gEnv->GetValue("XSec.GSI.SignProxy",    "1");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XrdSecGSISIGNPROXY"))
                            || strlen(cenv) <= 0))
      gSystem->Setenv("XrdSecGSISIGNPROXY",    val.Data());

   val = gEnv->GetValue("XSec.Pwd.AutoLogin",    "1");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XrdSecPWDAUTOLOG"))
                            || strlen(cenv) <= 0))
      gSystem->Setenv("XrdSecPWDAUTOLOG",      val.Data());

   val = gEnv->GetValue("XSec.Pwd.VerifySrv",    "1");
   if (val.Length() > 0 && (!(cenv = gSystem->Getenv("XrdSecPWDVERIFYSRV"))
                            || strlen(cenv) <= 0))
      gSystem->Setenv("XrdSecPWDVERIFYSRV",    val.Data());
}
