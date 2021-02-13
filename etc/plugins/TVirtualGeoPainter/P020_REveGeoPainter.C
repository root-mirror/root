void P020_REveGeoPainter()
{
   gPluginMgr->AddHandler("TVirtualGeoPainter", "web", "ROOT::Experimental::REveGeoPainter",
      "ROOTEve", "REveGeoPainter(TGeoManager*)");
}
