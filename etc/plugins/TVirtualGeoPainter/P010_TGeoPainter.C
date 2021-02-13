void P010_TGeoPainter()
{
   gPluginMgr->AddHandler("TVirtualGeoPainter", "root", "TGeoPainter",
      "GeomPainter", "TGeoPainter(TGeoManager*)");
}
