/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsable_RElement
#define ROOT7_Browsable_RElement

#include <ROOT/Browsable/RHolder.hxx>

#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {
namespace Browsable {

using RElementPath_t = std::vector<std::string>;

class RLevelIter;

/** \class RElement
\ingroup rbrowser
\brief Basic element of browsable hierarchy. Provides access to data, creates iterator if any
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RElement {
public:

   enum EContentKind {
      kNone,      ///< not recognized
      kText,      ///< "text" - plain text for code editor
      kImage,     ///< "image64" - base64 for supported image formats (png/gif/gpeg)
      kPng,       ///< "png" - plain png binary code, returned inside std::string
      kJpeg,      ///< "jpg" or "jpeg" - plain jpg binary code, returned inside std::string
      kJson,      ///< "json" representation of object, can be used in code editor
      kFileName   ///< "filename" - file name if applicable
   };

   static EContentKind GetContentKind(const std::string &kind);

   /** Possible actions on double-click */
   enum EActionKind {
      kActNone,    ///< do nothing
      kActBrowse,  ///< just browse (expand) item
      kActEdit,    ///< can provide data for text editor
      kActImage,   ///< can be shown in image viewer, can provide image
      kActDraw6,   ///< can be drawn inside ROOT6 canvas
      kActDraw7,   ///< can be drawn inside ROOT7 canvas
      kActGeom     ///< can be shown in geometry viewer
   };

   virtual ~RElement() = default;

   /** Name of browsable, must be provided in derived classes */
   virtual std::string GetName() const = 0;

   /** Checks if element name match to provided value */
   virtual bool MatchName(const std::string &name) const { return name == GetName(); }

   /** Title of browsable (optional) */
   virtual std::string GetTitle() const { return ""; }

   /** Create iterator for childs elements if any */
   virtual std::unique_ptr<RLevelIter> GetChildsIter();

   virtual int GetNumChilds();

   /** Returns element content, depends from kind. Can be "text" or "image64" or "json" */
   virtual std::string GetContent(const std::string & = "text");

   /** Access object */
   virtual std::unique_ptr<RHolder> GetObject() { return nullptr; }

   /** Get default action */
   virtual EActionKind GetDefaultAction() const { return kActNone; }

   /** Check if want to perform action */
   virtual bool IsCapable(EActionKind action) const { return action == GetDefaultAction(); }

   static std::shared_ptr<RElement> GetSubElement(std::shared_ptr<RElement> &elem, const RElementPath_t &path);

   static RElementPath_t ParsePath(const std::string &str);

   static int ComparePaths(const RElementPath_t &path1, const RElementPath_t &path2);

   static std::string GetPathAsString(const RElementPath_t &path);

   static int ExtractItemIndex(std::string &name);
};

} // namespace Browsable
} // namespace Experimental
} // namespace ROOT

#endif
