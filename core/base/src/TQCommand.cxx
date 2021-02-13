// @(#)root/base:$Id$
// Author: Valeriy Onuchin 04/27/2004

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TQCommand
\ingroup Base

The Command design pattern is based on the idea, that all editing
in an application is done by creating instances of command objects.
Command objects apply changes to the edited object and then are
stored  on a command stack. Furthermore, each command knows how to
undo its changes to bring the edited object back to its previous
state. As long as the application only uses command objects to
change the state of the edited object, it is possible to undo
a sequence of commands by traversing the command stack downwards and
calling the "undo" method of each command in turn. It is also
possible to redo a sequence of commands by traversing the command
stack upwards and calling the "redo" method of each command.


Examples:

1. Create a new command
~~~ {.cpp}
   TQCommand *com = new TQCommand("TH1", hpx, "SetFillColor(Color_t)"
                                "SetFillColor(Color_t)");
~~~
   - 1st parameter - the name of class
   - 2nd parameter - object
   - 3rd parameter - the name of do/redo method
   - 4th parameter - the name of undo method
   Since redo,undo methods are the same, undo name can be omitted, e.g.
~~~ {.cpp}
   TQCommand *com = new TQCommand("TH1", hpx, "SetFillColor(Color_t)");
~~~
   For objects derived from TObject class name can be omitted, e.g.
~~~ {.cpp}
   TQCommand *com = new TQCommand(hpx, "SetFillColor(Color_t)");
~~~

2. Setting undo, redo parameters.
~~~ {.cpp}
   Color_t old_color = hpx->GetFillColor();
   Color_t new_color = 4;  // blue color

   com->SetRedoArgs(1, new_color);
   com->SetUndoArgs(1, old_color);
~~~
   1st argument - the number of undo, redo parameters
   the other arguments - undo, redo values
   Since the number of undo,redo parameters is the same one can use
~~~ {.cpp}
   com->SetArgs(1, new_color, old_color);
~~~

3. Undo, redo method execution
~~~ {.cpp}
   com->Redo(); // execute redo method
   com->Undo(); // execute undo method
~~~

4. Merged commands
   It possible to group several commands together so an end user
   can undo and redo them with one command.
~~~ {.cpp}
   TQCommand *update = new TQCommand(gPad, "Modified()");
   com->Add(update);
~~~

5. Macro commands
   "Merging" allows to create macro commands, e.g.
~~~ {.cpp}
   TQCommand *macro = new TQCommand("my macro");
   macro->Add(com1);
   macro->Add(com2);
   ...
~~~
   During Redo operation commands composing macro command are executed
   sequentially in direct  order (first in first out). During Undo,
   they are executed in reverse order (last in first out).

6. Undo manager.
   TQUndoManager is recorder of undo and redo operations. This is
   command history list which can be traversed backwards and upwards
   performing undo and redo operations.
   To register command TQUndoManager::Add(TObject*) method is used.
~~~ {.cpp}
   TQUndoManager *history = new TQUndoManager();
   history->Add(com);
~~~
   TQUndoManager::Add automatically invokes execution of command's Redo method.

Use TQUndoManager::Undo to undo commands in  history list.
Redo is Undo for undo action. Use TQUndoManager::Redo method for that
*/

#include "TQCommand.h"
#include "TQConnection.h"
#include "TDataType.h"
#include "stdarg.h"
#include "TROOT.h"
#include "ThreadLocalStorage.h"
#include "TVirtualRWMutex.h"

ClassImp(TQCommand);
ClassImp(TQUndoManager);

static TQCommand *gActiveCommand = 0;

////////////////////////////////////////////////////////////////////////////////
/// Common protected method used in several constructors

void TQCommand::Init(const char *clname, void *obj, const char *redo, const char *undo)
{
   TString credo( CompressName(redo) );
   TString cundo( CompressName(undo) );

   fNRargs = fNUargs = -1;
   fNewDelete = kFALSE;
   fObject = obj;

   fRedo = redo ? new TQConnection(clname, obj, credo) : 0;
   fUndo = undo ? new TQConnection(clname, obj, cundo) : fRedo;

   fRedoArgs = 0;
   fUndoArgs = 0;
   fStatus = 0;
   fState = 0;

   if (!obj && !redo && !undo) { // macros
      fName = clname;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
///
///   Input parameters:
///     1. clname - class name.
///     2. obj - an object
///     3. redo - method or function to do/redo operation
///     4. undo - method or function to undo operation
///
/// Comments:
///    - if either clname or obj is NULL that means that redo/undo is function
///    - to specify default arguments for redo/undo method/function
///       '=' must precede to argument value.
///
///  Example:
/// ~~~ {.cpp}
///   TQCommand("TPad", gPad, "SetEditable(=kTRUE)", "SetEditable(=kFALSE)");
/// ~~~
///   undo method can be same as redo one. In that case undo parameter
///   can be omitted.
///
///  Example:
/// ~~~ {.cpp}
///   TQCommand("TPad", gPad, "SetFillStyle(Style_t)");
/// ~~~

TQCommand::TQCommand(const char *clname, void *obj, const char *redo,
                     const char *undo) : TList(), TQObject()
{
   Init(clname, obj, redo, undo);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
///
///   Input parameters:
///     1. obj - an object
///     2. redo - method or function to do/redo operation
///     3. undo - method or function to undo operation
///
/// Comments:
///   to specify default arguments for redo/undo method/function
///   '=' must precede to argument value.
///
///  Example:
/// ~~~ {.cpp}
///    TQCommand(gPad, "SetEditable(=kTRUE)", "SetEditable(=kFALSE)");
/// ~~~
///
///   undo method can be same as redo one. In that case "undo"
///   can parameter be omitted.
///
///  Example:
/// ~~~ {.cpp}
///    TQCommand(gPad, "SetFillStyle(Style_t)");
/// ~~~

TQCommand::TQCommand(TObject *obj, const char *redo, const char *undo) :
           TList(), TQObject()
{
   if (obj) Init(obj->ClassName(), obj, redo, undo);
   else Init(0, 0, redo, undo);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TQCommand::TQCommand(const TQCommand &com) : TList(), TQObject()
{
   fRedo = new TQConnection(*(com.fRedo));
   fUndo = new TQConnection(*(com.fUndo));

   fRedoArgs = 0;
   fUndoArgs = 0;
   fNRargs = com.fNRargs;
   fNUargs = com.fNUargs;

   if (fNRargs > 0) {
      fRedoArgs = new Long_t[fNRargs];
      for (int i = 0; i< fNRargs; i++) {
         fRedoArgs[i] = com.fRedoArgs[i];
      }
   }
   if (fNUargs > 0) {
      fUndoArgs = new Long_t[fNUargs];
      for (int i = 0; i < fNUargs; i++) {
         fUndoArgs[i] = com.fUndoArgs[i];
      }
   }
   fStatus = com.fStatus;
   fNewDelete = com.fNewDelete;
   fName = com.fName;
   fTitle = com.fTitle;
   fObject = com.fObject;
   fState = com.fState;

   // copy merged commands
   TIter next(&com);
   TQCommand *obj;
   while ((obj = (TQCommand*)next())) {
      TList::Add(new TQCommand(*obj));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// dtor.

TQCommand::~TQCommand()
{
   if (fRedo != fUndo) delete fUndo;

   delete fRedo;
   delete [] fRedoArgs;
   delete [] fUndoArgs;

   Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Return a command which is doing redo/undo action.
///
/// This static method allows to set undo parameters dynamically, i.e.
/// during execution of Redo function.
///
/// Example:
///   For redo actions like TGTextEdit::DelChar() it is not possible to
///   know ahead what character will be deleted.
///   To set arguments for undo action ( e.g. TGTextEdit::InsChar(char)),
///   one needs to call TQCommand::SetUndoArgs(1, character) from
///   inside of TGTextEdit::DelChar() method, i.e.
/// ~~~ {.cpp}
///    TQCommand::GetCommand()->SetUndoArgs(1, somechar);
/// ~~~

TQCommand *TQCommand::GetCommand()
{
   return gActiveCommand;
}

////////////////////////////////////////////////////////////////////////////////
/// If "opt" is not zero delete every merged command which option string is
/// equal to "opt". If "opt" is zero - delete all merged commands.

void TQCommand::Delete(Option_t *opt)
{
   if (!opt) {
      TList::Delete();
      return;
   }

   auto lnk = fFirst;
   decltype(lnk) sav;

   while (lnk) {
      sav = lnk->NextSP();
      TString ostr = lnk->GetOption();
      if (ostr.Contains(opt)) {   // remove command
         TObject *obj = lnk->GetObject();
         lnk->SetObject(nullptr);
         delete obj;
         Remove(lnk);
      }
      lnk = sav;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Two commands can be merged if they can be composed into
/// a single command (Macro command).
///
/// To allow merging commands user might override this function.

Bool_t TQCommand::CanMerge(TQCommand *) const
{
   return (!fRedo && !fUndo);
}

////////////////////////////////////////////////////////////////////////////////
/// Add command to the list of merged commands.
/// This make it possible to group complex actions together so an end user
/// can undo and redo them with one command. Execution of TQUndoManager::Undo(),
/// TQUndoManager::Redo() methods only invokes the top level command as a whole.
///
/// Merge method is analogous to logical join operation.
///
/// Note:  Merge method invokes redo action.

void TQCommand::Merge(TQCommand *c)
{
   Add(c, "merge");
}

////////////////////////////////////////////////////////////////////////////////
/// Merge a collection of TQCommand.

Long64_t TQCommand::Merge(TCollection *collection,TFileMergeInfo*)
{
   TIter next(collection);
   while (TObject* o = next()) {
      TQCommand *command = dynamic_cast<TQCommand*> (o);
      if (!command) {
         Error("Merge",
               "Cannot merge - an object which doesn't inherit from TQCommand found in the list");
         return -1;
      }
      Merge(command);
   }
   return GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Add command to the list of merged commands.
///
/// Option string can contain substrings:
///  - "compress" - try to compress input command
///  - "radd" - execute redo action of input command
///  - "uadd" - execute undo action of input command

void TQCommand::Add(TObject *obj, Option_t *opt)
{
   if (!obj->InheritsFrom(TQCommand::Class())) return;

   TQCommand *o = (TQCommand *)obj;
   TQCommand *c = (TQCommand *)Last();
   TString ostr = opt;

   if (c) {
      if (c->CanCompress(o) || (c->IsEqual(o) && ostr.Contains("compress"))) {
         c->Compress(o);
         return;
      }
   }
   TList::AddLast(o, opt);
   if (o->CanRedo() && ostr.Contains("radd")) o->Redo();
   if (o->CanUndo() && ostr.Contains("uadd")) o->Undo();
}

////////////////////////////////////////////////////////////////////////////////
/// By default, commands can be compressed if they are:
///
///  - equal
///  - setter commands
///
/// More complicated commands might want to override this function.

Bool_t TQCommand::CanCompress(TQCommand *c) const
{
   return (IsEqual(c) && IsSetter());
}

////////////////////////////////////////////////////////////////////////////////
/// Compress command. Compression is analogous to arithmetic "addition operation".
///
/// Note:
///   - The compressed command will be deleted.
///   - Execution Compress method invokes Redo action with new redo arguments
///     inherited from compressed command.
///
/// More complicated commands might want to override this function.

void TQCommand::Compress(TQCommand *c)
{
   for (int i = 0; i < fNRargs; i++) {
      fRedoArgs[i] = c->fRedoArgs[i];
   }
   Redo();
   fStatus--;   //do not change the state of command
   delete c;
}

////////////////////////////////////////////////////////////////////////////////
/// Equal comparison. The commands are equal if they are
/// applied to the same object and have the same Redo/Undo actions
///
/// More complicated commands might want to override this function.

Bool_t TQCommand::IsEqual(const TObject* obj) const
{
   if (!obj->InheritsFrom(TQCommand::Class())) return kFALSE;
   TQCommand *c = (TQCommand *)obj;
   if (!fRedo || !fUndo || (c->GetObject() != fObject)) return kFALSE;

   TString cname = fRedo->GetClassName();
   TString rname = fRedo->GetName();

   return ((cname == c->GetRedo()->GetClassName()) &&
           (rname == c->GetRedo()->GetName()));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE is command if Redo is the same as Undo function
/// and is the setter action.
///
/// By default, all functions with names like "SetXXX" or "setXXX"
/// considered as setters. Another type of setters are Move, Resize operations
///
/// More complicated commands might want to override this function.

Bool_t TQCommand::IsSetter() const
{
   TString redo = GetRedoName();
   TString undo = GetUndoName();

   if (!redo || !undo || (redo != undo)) return kFALSE;

   return (redo.BeginsWith("Set") ||
           redo.BeginsWith("set") ||
           redo.BeginsWith("Move") ||
           redo.BeginsWith("move") ||
           redo.BeginsWith("Resize") ||
           redo.BeginsWith("resize"));
}

////////////////////////////////////////////////////////////////////////////////
/// Set do/redo and undo parameters. The format is
///    SetArgs(number_of_params, redo_params, undo_params)
///
/// Example:
/// ~~~ {.cpp}
///     move_command->SetArgs(2, 100, 100, 200, 200);
/// ~~~
///      2 params, (100,100) - do/redo position, (200,200) - undo position

void TQCommand::SetArgs(Int_t narg, ...)
{
   if (narg < 0) {
      return;
   } else if (!narg) {  // no arguments
      fNRargs = fNUargs = narg;
      return;
   }

   va_list ap;
   va_start(ap, narg);

   if (fNRargs != narg ) {
      delete fRedoArgs;
   }
   fRedoArgs = new Long_t[narg];

   if (fNUargs != narg ) {
      delete fUndoArgs;
   }
   fUndoArgs = new Long_t[narg];

   fNRargs = fNUargs = narg;

   Int_t i;
   for (i = 0; i < fNRargs; i++) {
      fRedoArgs[i] = va_arg(ap, Long_t);
   }
   for (i = 0; i < fNUargs; i++) {
      fUndoArgs[i] = va_arg(ap, Long_t);
   }
   va_end(ap);
}

////////////////////////////////////////////////////////////////////////////////
/// Set redo parameters. The format is
///    SetRedoArgs(number_of_params, params)
///
/// Example:
/// ~~~ {.cpp}
///     move_command->SetRedoArgs(2, 100, 100);
/// ~~~

void TQCommand::SetRedoArgs(Int_t narg, ...)
{
   if (narg < 0) {
      return;
   } else if (!narg) {  // no arguments
      fNRargs = 0;
      return;
   }

   va_list ap;
   va_start(ap, narg);

   if (fNRargs != narg ) {
      delete fRedoArgs;
   }
   fRedoArgs = new Long_t[narg];

   fNRargs = narg;

   for (int i = 0; i < fNRargs; i++) {
      fRedoArgs[i] = va_arg(ap, Long_t);
   }
   va_end(ap);
}

////////////////////////////////////////////////////////////////////////////////
/// Set undo parameters. The format is
///    SetUndoArgs(number_of_params, params)
///
/// Example:
/// ~~~ {.cpp}
///     move_command->SetUndoArgs(2, 200, 200);
/// ~~~

void TQCommand::SetUndoArgs(Int_t narg, ...)
{
   if (narg < 0) {
      return;
   } else if (!narg) {  // no arguments
      fNUargs = narg;
      return;
   }

   va_list ap;
   va_start(ap, narg);

   if (fNUargs != narg ) {
      delete fUndoArgs;
   }
   fUndoArgs = new Long_t[narg];

   fNUargs = narg;

   for (int i = 0; i < fNUargs; i++) {
      fUndoArgs[i] = va_arg(ap, Long_t);
   }
   va_end(ap);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if Redo action is possible, kFALSE if it's not.
/// By default, only single sequential redo action is possible.

Bool_t TQCommand::CanRedo() const
{
   return (fStatus <= 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if Undo action is possible, kFALSE if it's not.
/// By default, only single trial undo action is possible.

Bool_t TQCommand::CanUndo() const
{
   return (fStatus > 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute command and then merge commands

void TQCommand::Redo(Option_t *)
{
   Bool_t done = kFALSE;
   fState = 1;

   gActiveCommand = this;

   if (fNRargs > 0) {
      if (fRedo) {
         fRedo->ExecuteMethod(fRedoArgs, fNRargs);
         done = kTRUE;
      }
   } else if (!fNRargs) {
      if (fRedo) {
         fRedo->ExecuteMethod();
         done = kTRUE;
      }
   }

   // execute merged commands
   auto lnk = fFirst;
   while (lnk) {
      TQCommand *c = (TQCommand *)lnk->GetObject();
      c->Redo();
      done = kTRUE;
      lnk = lnk->NextSP();
   }

   if (done) Emit("Redo()");
   fStatus++;
   fState = 0;
   gActiveCommand = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Un-execute all merged commands and the command.
/// Merged commands are executed in reverse order.

void TQCommand::Undo(Option_t *)
{
   Bool_t done = kFALSE;
   fState = -1;

   gActiveCommand = this;

   // unexecute merged commands
   auto lnk = fLast;
   while (lnk) {
      TQCommand *c = (TQCommand *)lnk->GetObject();
      TString opt = lnk->GetOption();
      auto sav = lnk->PrevSP();
      c->Undo();
      done = kTRUE;
      if (opt.Contains("remove")) {   // remove  command
         delete lnk->GetObject();
         Remove(lnk);
      }
      lnk = sav;
   }
   if (fNUargs > 0) {
      if (fUndo) {
         fUndo->ExecuteMethod(fUndoArgs, fNUargs);
         done = kTRUE;
      }
   } else if (!fNUargs) {
      if (fUndo) {
         fUndo->ExecuteMethod();
         done = kTRUE;
      }
   }

   if (done) Emit("Undo()");
   fStatus--;
   fState = 0;
   gActiveCommand = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the command name. Default name is "ClassName::RedoName(args)"
/// If list of merged commands is not empty the name is
/// "ClassName::RedoName(args):cname1:cname2 ..."

const char *TQCommand::GetName() const
{
   const Int_t maxname = 100;

   if (!fName.IsNull())
      return fName.Data();

   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

   // In case another thread already did the work while
   // we were waiting.
   if (!fName.IsNull())
      return fName.Data();

   TString name;

   if (fRedo) {
      if (fRedo->GetClassName()) {
         name = fRedo->GetClassName();
      }
      name += "::";
      name += fRedo->GetName();
   }
   TQCommand *c;
   TObjLink *lnk = fFirst.get();

   while (lnk && (fName.Length() < maxname)) {
      c = (TQCommand *)lnk->GetObject();
      name += ":";
      name += c->GetName();
      lnk = lnk->Next();
   }

   TQCommand *m = const_cast<TQCommand*>(this);
   m->fName = name;

   return fName;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns command description.
/// By default, "ClassName::RedoName(args)_ClassName::UndoName(args)"

const char *TQCommand::GetTitle() const
{
   if (!fTitle.IsNull())
      return fTitle.Data();

   if (fUndo) {
      TTHREAD_TLS_DECL_ARG(TString, title, GetName());

      title += "_";
      title += fUndo->GetClassName();
      title += "::";
      if (fUndo->GetName())
         title += fUndo->GetName();

      return title.Data();
   } else {
      return GetName();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the name of redo command

const char *TQCommand::GetRedoName() const
{
   return (fRedo ? fRedo->GetName() : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the name of undo command

const char *TQCommand::GetUndoName() const
{
   return (fUndo ? fUndo->GetName() : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to array of redo arguments

Long_t *TQCommand::GetRedoArgs() const
{
   return fRedoArgs;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to array of undo arguments

Long_t *TQCommand::GetUndoArgs() const
{
   return fUndoArgs;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a number of redo arguments

Int_t TQCommand::GetNRargs() const
{
   return fNRargs;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a number of undo arguments

Int_t TQCommand::GetNUargs() const
{
   return fNUargs;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an object for which undo redo actions are applied

void *TQCommand::GetObject() const
{
   return fObject;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a number of sequential undo or redo operations

Int_t TQCommand::GetStatus() const
{
   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if neither redo nor undo action specified

Bool_t TQCommand::IsMacro() const
{
   return (!fRedo && !fUndo);
}

////////////////////////////////////////////////////////////////////////////////
/// Undo action is in progress

Bool_t TQCommand::IsUndoing() const
{
   return (fState < 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Redo action is in progress

Bool_t TQCommand::IsRedoing() const
{
   return (fState > 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if command execution is in progress

Bool_t TQCommand::IsExecuting() const
{
   return fState;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets name of the command

void TQCommand::SetName(const char *name)
{
   fName = name;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets description of the command

void TQCommand::SetTitle(const char *title)
{
   fTitle = title;
}

////////////////////////////////////////////////////////////////////////////////
/// ls this command and merged commands

void TQCommand::ls(Option_t *) const
{
   TString name = GetName();
   printf("%d %s\n", fStatus, name.Data());

   TObjLink *lnk = fFirst.get();
   while (lnk) {
      printf("\t");
      lnk->GetObject()->ls();
      lnk = lnk->Next();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print collection header.

void TQCommand::PrintCollectionHeader(Option_t* /*option*/) const
{
   TROOT::IndentLevel();
   printf("%d %s\n", fStatus, GetName());
}

/** \class TQUndoManager
Recorder of operations for undo and redo
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQUndoManager::TQUndoManager() : TQCommand(0, 0, 0, 0)
{
   fCursor = 0;
   fLimit = kMaxUInt;   // maximum value for UInt_t
   fLogging = kFALSE;
   fLogBook = 0;
   fCurrent = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TQUndoManager::~TQUndoManager()
{
   Delete();

   if (fLogBook) {
      delete fLogBook;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Lists all commands in stack

void TQUndoManager::ls(Option_t *option) const
{
   if (!IsEmpty()) {
      TObjLink *lnk = fFirst.get();
      while (lnk) {
         if (lnk == fCursor) {
            printf("->");
         } else {
            printf("  ");
         }
         TQCommand *com = (TQCommand*)lnk->GetObject();
         com->ls(option);
         lnk = lnk->Next();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print collection entry.

void TQUndoManager::PrintCollectionEntry(TObject* entry, Option_t* option,
                                         Int_t /*recurse*/) const
{
   TQCommand *com = (TQCommand*) entry;
   TROOT::IndentLevel();
   if (fCursor && fCursor->GetObject() == entry) {
      printf("->");
   } else {
      printf("  ");
   }
   com->ls(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Start logging. Delete all previous log records
/// Note: logging is not implemented yet

void  TQUndoManager::SetLogging(Bool_t on)
{
   fLogging = on;

   if (fLogging) {
      if (fLogBook) {
         fLogBook->Delete();
      } else {
         fLogBook = new TList();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add command to the stack of commands.
/// Command's redo action will be executed.
///
/// option string can contain the following substrings:
///  - "merge" - input command will be merged
///  - "compress" - input command will be compressed

void TQUndoManager::Add(TObject *obj, Option_t *opt)
{
   if (!obj->InheritsFrom(TQCommand::Class())) return;

   TQCommand *o = (TQCommand *)obj;
   TQCommand *c;
   Bool_t onredo = fCursor && fCursor->Next();
   TString ostr = onredo ? "1radd" : "0radd"; // execute redo on add
   if (opt) ostr += opt;

   if (fState) { // undo/redo in progress
      c = fCurrent;
      if (c) {
         fCurrent = o;
         c->Add(o, "remove");   // add nested command
      }
      return;
   }

   // delete all commands after cursor position
   if (fCursor && fCursor->Next()) {
      TObjLink *lnk = fCursor->Next();
      TObjLink *sav;
      while (lnk) {
         sav = lnk->Next();
         delete lnk->GetObject();
         Remove(lnk);
         lnk = sav;
      }
   }

   c = GetCursor();
   if (c) {
      if (c->CanCompress(o) || c->CanMerge(o) ||
          ostr.Contains("merge") || ostr.Contains("compress")) {
         fState = 1;
         c->Add(o, ostr.Data());
         fState = 0;
         return;
      }
   }

   TList::AddLast(obj, ostr.Data());
   fCursor = fLast.get();
   Redo(ostr.Data());

   if ((fSize > 0) && ((UInt_t)fSize > fLimit)) {
      Remove(fFirst);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// emit signal

void TQUndoManager::CurrentChanged(TQCommand *c)
{
   Emit("CurrentChanged(TQCommand*)", (long)c);
}

////////////////////////////////////////////////////////////////////////////////
/// Performs undo action. Move cursor position backward in history stack

void TQUndoManager::Undo(Option_t *option)
{
   Bool_t done = kFALSE;
   if (!CanUndo()) return;

   TQCommand *sav = fCurrent;
   TQCommand *c = (TQCommand*)fCursor->GetObject();

   if (c->CanUndo()) {
      fState = -1;
      fCurrent = c;
      fCurrent->Undo(option);
      fState = 0;
      done = kTRUE;
      fCursor = fCursor->Prev() ? fCursor->Prev() : fFirst.get();
   } else {
      fCursor = fCursor->Prev();
      fCurrent = (TQCommand*)fCursor->GetObject();
      fState = -1;
      fCurrent->Undo(option);
      fState = 0;
      done = kTRUE;
   }
   if (done && fLogging && fLogBook) {
      fLogBook->Add(new TQCommand(*fCurrent));
   }
   if (sav != fCurrent) CurrentChanged(fCurrent);
}

////////////////////////////////////////////////////////////////////////////////
/// Performs redo action. Move cursor position forward in history stack

void TQUndoManager::Redo(Option_t *option)
{
   Bool_t done = kFALSE;
   if (!CanRedo()) return;

   TQCommand *sav = fCurrent;
   TQCommand *c = (TQCommand*)fCursor->GetObject();

   if (c->CanRedo()) {
      fState = 1;
      fCurrent = c;
      fCurrent->Redo(option);
      fState = 0;
      done = kTRUE;
      fCursor = fCursor->Next() ? fCursor->Next() : fLast.get();
   } else {
      fCursor = fCursor->Next();
      fCurrent = (TQCommand*)fCursor->GetObject();
      fState = 1;
      fCurrent->Redo(option);
      fState = 0;
      done = kTRUE;
   }
   if (done && fLogging && fLogBook) {
      fLogBook->Add(new TQCommand(*fCurrent));
   }
   if (sav != fCurrent) CurrentChanged(fCurrent);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if redo action is possible

Bool_t TQUndoManager::CanRedo() const
{
   if (!fCursor) return kFALSE;

   TQCommand *c = (TQCommand*)fCursor->GetObject();
   if (c->CanRedo()) return kTRUE;

   c = fCursor->Next() ? (TQCommand*)fCursor->Next()->GetObject() : 0;
   return (c && c->CanRedo());
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if undo action is possible

Bool_t TQUndoManager::CanUndo() const
{
   if (!fCursor) return kFALSE;

   TQCommand *c = (TQCommand*)fCursor->GetObject();
   if (c->CanUndo()) return kTRUE;

   c = fCursor->Prev() ? (TQCommand*)fCursor->Prev()->GetObject() : 0;
   return (c && c->CanUndo());
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if logging is ON

Bool_t TQUndoManager::IsLogging() const
{
   return fLogging;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the last executed command

TQCommand *TQUndoManager::GetCurrent() const
{
   return fCurrent;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a command correspondent to the current cursor position in stack

TQCommand *TQUndoManager::GetCursor() const
{
   return (TQCommand*)(fCursor ? fCursor->GetObject() : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a maximum number of commands which could be located in stack

void TQUndoManager::SetLimit(UInt_t limit)
{
   fLimit = limit;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a maximum number of commands which  could be located in stack

UInt_t TQUndoManager::GetLimit() const
{
   return fLimit;
}
