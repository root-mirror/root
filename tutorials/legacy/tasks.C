/// \file
/// \ingroup tutorial_legacy
/// Example of TTasks.
/// Create a hierarchy of objects derived from TTask in library Mytasks
/// Show the tasks in a browser.
/// To execute a Task, use the context context menu and select
/// the item "ExecuteTask"
/// see also other functions in the TTask context menu, such as
///  - setting a breakpoint in one or more tasks
///  - enabling/disabling one task, etc
///
/// As of 2021, tbb is a fine replacement for ROOT's task system. It
/// offers better interplay with non-ROOT multi-threading use cases.
///
/// \macro_code
///
/// \author Rene Brun

#ifndef __RUN_TASKS__

void tasks()
{
   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("tasks.C","");
   dir.ReplaceAll("/./","/");
   gROOT->LoadMacro(dir +"MyTasks.cxx+");

   gROOT->ProcessLine("#define __RUN_TASKS__ 1");
   gROOT->ProcessLine(TString("#include \"") + dir + "tasks.C\"");
   gROOT->ProcessLine("runtasks()");
   gROOT->ProcessLine("#undef __RUN_TASKS__");
}

#else

void runtasks()
//void tasks()
{
   TTask *run      = new MyRun("run","Process one run");
   TTask *event    = new MyEvent("event","Process one event");
   TTask *geomInit = new MyGeomInit("geomInit","Geometry Initialisation");
   TTask *matInit  = new MyMaterialInit("matInit","Materials Initialisation");
   TTask *tracker  = new MyTracker("tracker","Tracker manager");
   TTask *tpc      = new MyRecTPC("tpc","TPC Reconstruction");
   TTask *its      = new MyRecITS("its","ITS Reconstruction");
   TTask *muon     = new MyRecMUON("muon","MUON Reconstruction");
   TTask *phos     = new MyRecPHOS("phos","PHOS Reconstruction");
   TTask *rich     = new MyRecRICH("rich","RICH Reconstruction");
   TTask *trd      = new MyRecTRD("trd","TRD Reconstruction");
   TTask *global   = new MyRecGlobal("global","Global Reconstruction");

   run->Add(geomInit);
   run->Add(matInit);
   run->Add(event);
   event->Add(tracker);
   event->Add(global);
   tracker->Add(tpc);
   tracker->Add(its);
   tracker->Add(muon);
   tracker->Add(phos);
   tracker->Add(rich);
   tracker->Add(trd);

   gROOT->GetListOfTasks()->Add(run);
   gROOT->GetListOfBrowsables()->Add(run);
   new TBrowser;
}

#endif
