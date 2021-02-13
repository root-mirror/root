\defgroup webgui6 ROOT 6 web display
\ingroup webgui
\brief To display %ROOT 6 canvases in the web browser

This group contains TWebCanvas class which provides web-based TCanvasImp
and allows display of **%ROOT 6 TCanvas** in the web browser.

This is fully reimplements TVirtualX and TVirtualPadPainter classes,
supporting majority of existing ROOT classes. Implementation does not
provide some interactive features - like custom mouse events handling.
Object changes performed in the browser (histogram color change)
are not reflected in the C++ objects -
WebGui provides READ-ONLY display capability

