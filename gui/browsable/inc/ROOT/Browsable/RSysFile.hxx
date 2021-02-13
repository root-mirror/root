/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsable_RSysFile
#define ROOT7_Browsable_RSysFile

#include <ROOT/Browsable/RElement.hxx>
#include <ROOT/Browsable/RGroup.hxx>

#include "TSystem.h"
#include <string>

namespace ROOT {
namespace Experimental {
namespace Browsable {

class RSysDirLevelIter;

/** \class RSysFile
\ingroup rbrowser
*/

class RSysFile : public RElement {

   friend class RSysDirLevelIter;

   FileStat_t fStat;       ///<! file stat object
   std::string fDirName;   ///<! fully-qualified directory name
   std::string fFileName;  ///<! file name in current dir

   std::string GetFullName() const;

public:
   RSysFile(const std::string &filename);

   RSysFile(const FileStat_t& stat, const std::string &dirname, const std::string &filename);

   virtual ~RSysFile() = default;

   /** Name of RElement - file name in this case */
   std::string GetName() const override;

   /** Checks if element name match to provided value */
   bool MatchName(const std::string &name) const override;

   /** Title of RElement - full file name  */
   std::string GetTitle() const override { return GetFullName(); }

   std::unique_ptr<RLevelIter> GetChildsIter() override;

   std::string GetContent(const std::string &kind) override;

   EActionKind GetDefaultAction() const override;

   static std::string GetFileIcon(const std::string &fname);

   static RElementPath_t ProvideTopEntries(std::shared_ptr<RGroup> &comp, const std::string &workdir = "");

};

} // namespace Browsable
} // namespace Experimental
} // namespace ROOT


#endif
