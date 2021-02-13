// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2018-10-17
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RWebDisplayHandle
#define ROOT7_RWebDisplayHandle

#include <ROOT/RWebDisplayArgs.hxx>

#include <string>
#include <map>
#include <memory>

namespace ROOT {
namespace Experimental {

class RWebDisplayHandle {

   std::string fUrl; ///!< URL used to launch display

   std::string fContent; ///!< page content

protected:
   class Creator {
   public:
      virtual std::unique_ptr<RWebDisplayHandle> Display(const RWebDisplayArgs &args) = 0;
      virtual bool IsActive() const { return true; }
      virtual ~Creator() = default;
   };

   class BrowserCreator : public Creator {
   protected:
      std::string fProg;  ///< browser executable
      std::string fExec;  ///< standard execute line
      std::string fBatchExec; ///< batch execute line

      void TestProg(const std::string &nexttry, bool check_std_paths = false);

      virtual void ProcessGeometry(std::string &, const RWebDisplayArgs &) {}
      virtual std::string MakeProfile(std::string &, bool) { return ""; }

   public:

      BrowserCreator(bool custom = true, const std::string &exec = "");

      std::unique_ptr<RWebDisplayHandle> Display(const RWebDisplayArgs &args) override;

      virtual ~BrowserCreator() = default;
   };

   class ChromeCreator : public BrowserCreator {
   public:
      ChromeCreator();
      virtual ~ChromeCreator() = default;
      bool IsActive() const override { return !fProg.empty(); }
      void ProcessGeometry(std::string &, const RWebDisplayArgs &args) override;
      std::string MakeProfile(std::string &exec, bool) override;
   };

   class FirefoxCreator : public BrowserCreator {
   public:
      FirefoxCreator();
      virtual ~FirefoxCreator() = default;
      bool IsActive() const override { return !fProg.empty(); }
      std::string MakeProfile(std::string &exec, bool batch) override;
   };

   static std::map<std::string, std::unique_ptr<Creator>> &GetMap();

   static std::unique_ptr<Creator> &FindCreator(const std::string &name, const std::string &libname = "");

public:

   RWebDisplayHandle(const std::string &url) : fUrl(url) {}

   // required virtual destructor for correct cleanup at the end
   virtual ~RWebDisplayHandle() = default;

   const std::string &GetUrl() const { return fUrl; }

   void SetContent(const std::string &cont) { fContent = cont; }
   const std::string &GetContent() const { return fContent; }

   static std::unique_ptr<RWebDisplayHandle> Display(const RWebDisplayArgs &args);

   static bool DisplayUrl(const std::string &url);

   static bool ProduceImage(const std::string &fname, const std::string &json, int width = 800, int height = 600);
};

}
}

#endif
