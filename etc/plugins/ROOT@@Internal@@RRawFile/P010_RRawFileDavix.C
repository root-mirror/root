void P010_RRawFileDavix()
{
   TString configfeatures = gROOT->GetConfigFeatures();

   // only if ROOT was compiled with davix enabled do we configure a handler
   if (configfeatures.Contains("davix") &&
       !gEnv->GetValue("Davix.UseOldClient", 0)) {

      gPluginMgr->AddHandler(
         "ROOT::Internal::RRawFile",
         "^http[s]?:",
         "ROOT::Internal::RRawFileDavix",
         "RDAVIX",
         "RRawFileDavix(std::string_view, ROOT::Internal::RRawFile::ROptions)");
   }
}
