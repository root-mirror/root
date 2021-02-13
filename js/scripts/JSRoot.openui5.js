/// @file JSRoot.openui5.js
/// Bootstraping of OpenUI5 functionality in JSROOT
/// Openui5 loaded directly in the script

JSROOT.define(['jquery', 'jquery-ui'], () => {

   "use strict";

   // very simple - openui5 was loaded before and will be used as is
   if (typeof sap == 'object')
      return sap;

   let resolveFunc, rootui5sys;

   JSROOT._.completeUI5Loading = function() {
      sap.ui.loader.config({
         paths: {
            jsroot: JSROOT.source_dir,
            rootui5: rootui5sys
         }
      });

      if (resolveFunc) {
         resolveFunc(sap);
         resolveFunc = null;
      }
   }

   function TryOpenOpenUI(sources) {

      // where to take openui5 sources
      let src = sources.shift();

      if ((src.indexOf("roothandler")==0) && (src.indexOf("://")<0)) src = src.replace(/\:\//g,"://");

      let element = document.createElement("script");
      element.setAttribute('type', "text/javascript");
      element.setAttribute('id', "sap-ui-bootstrap");
      // use nojQuery while we are already load jquery and jquery-ui, later one can use directly sap-ui-core.js

      // this is location of openui5 scripts when working with THttpServer or when scripts are installed inside JSROOT
      element.setAttribute('src', src + "resources/sap-ui-core-nojQuery.js"); // latest openui5 version

      element.setAttribute('data-sap-ui-libs', JSROOT.openui5libs || "sap.m, sap.ui.layout, sap.ui.unified, sap.ui.commons");

      element.setAttribute('data-sap-ui-theme', JSROOT.openui5theme || 'sap_belize');
      element.setAttribute('data-sap-ui-compatVersion', 'edge');
      // element.setAttribute('data-sap-ui-bindingSyntax', 'complex');

      element.setAttribute('data-sap-ui-preload', 'async'); // '' to disable Component-preload.js

      element.setAttribute('data-sap-ui-evt-oninit', "JSROOT._.completeUI5Loading()");

      element.onerror = function() {
         // remove failed element
         element.parentNode.removeChild(element);
         // and try next
         TryOpenOpenUI(sources);
      }

      element.onload = function() {
         console.log('Load openui5 from ' + src);
      }

      document.getElementsByTagName("head")[0].appendChild(element);
   }

   rootui5sys = JSROOT.source_dir.replace(/jsrootsys/g, "rootui5sys");
   if (rootui5sys == JSROOT.source_dir) {
      // if jsrootsys location not detected, try to guess it
      if (window.location.port && (window.location.pathname.indexOf("/win") >= 0) && (!JSROOT.openui5src || JSROOT.openui5src == 'nojsroot' || JSROOT.openui5src == 'jsroot'))
         rootui5sys = window.location.origin + window.location.pathname + "../rootui5sys/";
      else
         rootui5sys = undefined;
   }

   let openui5_sources = [],
       openui5_dflt = "https://openui5.hana.ondemand.com/1.82.2/",
       openui5_root = "";

   if (rootui5sys) openui5_root = rootui5sys + "distribution/";

   if (typeof JSROOT.openui5src == 'string') {
      switch (JSROOT.openui5src) {
         case "nodefault": openui5_dflt = ""; break;
         case "default": openui5_sources.push(openui5_dflt); openui5_dflt = ""; break;
         case "nojsroot": /* openui5_root = ""; */ break;
         case "jsroot": openui5_sources.push(openui5_root); openui5_root = ""; break;
         default: openui5_sources.push(JSROOT.openui5src); break;
      }
   }

   if (openui5_root && (openui5_sources.indexOf(openui5_root) < 0)) openui5_sources.push(openui5_root);
   if (openui5_dflt && (openui5_sources.indexOf(openui5_dflt) < 0)) openui5_sources.push(openui5_dflt);

   // return Promise let loader wait before dependent source will be invoked
   return new Promise(resolve => {
      resolveFunc = resolve;
      TryOpenOpenUI(openui5_sources);
   });

});
