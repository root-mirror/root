/// \file
/// \ingroup tutorial_eve
/// Plays back event-recording of a root session running geom_cms.C tutorial.
/// [ Recorded using "new TGRecorder" command. ]
///
/// Script:
/// - type: .x geom_cms.C
/// - demonstrate rotation (left-mouse), zoom (right-mouse left-right)
/// - show GL window Help Window
/// - show wireframe (w), smooth (r, default) and outline (t) render modes
/// - show flip of background color dark-light-dark (e pressed twice)
/// - disable clipping plane in GL-viewer panel
/// - open "Scene" list-tree and further "Geometry scene"
/// - disable drawing of muon system and then calorimeters
/// - select tracker geometry top-node and increase drawing depth
/// - re-enable clipping plane and zoom into pixel detector.
///
/// \macro_code
///
/// \author Matevz Tadel

void geom_cms_playback()
{
   auto r = new TRecorder("http://root.cern/files/geom_cms_recording.root");
}
