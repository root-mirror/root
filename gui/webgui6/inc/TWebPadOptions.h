// Author:  Sergey Linev, GSI  29/06/2017

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TWebPadOptions
#define ROOT_TWebPadOptions

#include <string>
#include <vector>

/** \class TWebObjectOptions
\ingroup webgui6

Class used to transport drawing options from the client

*/

class TWebObjectOptions {
public:
   std::string snapid;       ///< id of the object
   std::string opt;          ///< drawing options
   std::string fcust;        ///< custom string
   std::vector<double> fopt; ///< custom float array
};

/// Class used to transport ranges from JSROOT canvas
class TWebPadOptions {
public:
   std::string snapid;                        ///< id of pad
   bool active{false};                        ///< if pad selected as active
   int logx{0}, logy{0}, logz{0};             ///< pad log properties
   int gridx{0}, gridy{0};                    ///< pad grid properties
   int tickx{0}, ticky{0};                    ///< pad ticks properties
   float mleft{0}, mright{0}, mtop{0}, mbottom{0}; ///< frame margins
   bool ranges{false};                        ///< if true, pad has ranges
   double px1{0}, py1{0}, px2{0}, py2{0};     ///< pad range
   double ux1{0}, uy1{0}, ux2{0}, uy2{0};     ///< pad axis range
   unsigned bits{0};                          ///< canvas status bits like tool editor
   double zx1{0}, zx2{0}, zy1{0}, zy2{0}, zz1{0}, zz2{0}; ///< zooming ranges
   std::vector<TWebObjectOptions> primitives; ///< drawing options for primitives
};

/// Class used to transport pad click events
class TWebPadClick {
public:
   std::string padid;                         ///< id of pad
   std::string objid;                         ///< id of clicked object, "null" when not defined
   int x{-1};                                 ///< x coordinate of click event
   int y{-1};                                 ///< y coordinate of click event
   bool dbl{false};                           ///< when double-click was performed
};

#endif
