/** @fileoverview Core methods of JavaScript ROOT
  * @namespace JSROOT */

(function( factory ) {
   if ( typeof define === "function" && define.amd ) {

      let jsroot = factory({}),
          norjs = (typeof requirejs=='undefined'),
          paths = {};

      jsroot._.amd = true; // inidcation that require will be used for loading of functionality

      for (let src in jsroot._.sources)
         if (src != 'JSRoot.core')
            paths[src] = jsroot._.get_module_src(jsroot._.sources[src]);

      if (norjs) {
         // just define locations
         require({ paths: paths });
      } else {
         let cfg_paths;
         if ((requirejs.s!==undefined) && (requirejs.s.contexts !== undefined) && ((requirejs.s.contexts._!==undefined) &&
               requirejs.s.contexts._.config!==undefined)) cfg_paths = requirejs.s.contexts._.config.paths;
         else console.warn("Require.js paths changed - please contact JSROOT developers");

         // check if modules are already loaded
         for (let p in paths)
            if (requirejs.defined(p) || (cfg_paths && (p in cfg_paths)))
               delete paths[p];

         // configure all dependencies
         requirejs.config({
            paths: paths,
            shim: {
               'jquery-ui': { deps: ['jquery'] },
               'jqueryui-mousewheel': { deps: ['jquery-ui'] },
               'jqueryui-touch-punch': { deps: ['jquery-ui'] }
            }
         });
      }

      define( jsroot );

      if (norjs || !require.specified("JSRoot.core"))
         define('JSRoot.core', [], jsroot);

      if (norjs || !require.specified("jsroot"))
         define('jsroot', [], jsroot);

      globalThis.JSROOT = jsroot;

      jsroot._.init();

   } else if (typeof exports === 'object' && typeof module !== 'undefined') {
      // processing with Node.js

      //  mark JSROOT as used with Node.js
      exports.batch_mode = exports.nodejs = (typeof global==='object') && global.process && (Object.prototype.toString.call(global.process) === '[object process]');

      factory(exports);

      globalThis.JSROOT = exports;

   } else {

      if ((typeof JSROOT != 'undefined') && !JSROOT._workaround)
         throw new Error("JSROOT is already defined", "JSRoot.core.js");

      let jsroot = {};

      factory(jsroot);

      globalThis.JSROOT = jsroot;

      if (globalThis.sap && globalThis.sap.ui) {
         console.log('Use SAP loader');

         // actiavate SAP loader
         jsroot._.sap = true;

         let rootui5sys = undefined;
         if (jsroot.source_dir.indexOf("jsrootsys") >= 0)
            rootui5sys = jsroot.source_dir.replace(/jsrootsys/g, "rootui5sys");
         sap.ui.loader.config({
            paths: {
               jsroot: jsroot.source_dir,
               rootui5: rootui5sys
            }
         });

         if (globalThis.jQuery)
            jsroot._.modules['jquery'] = { module: globalThis.jQuery };
      }

      jsroot._.init();
   }
} (function(JSROOT) {

   "use strict";

   /** @summary JSROOT version id
     * @desc For the JSROOT release the string in format "major.minor.patch" like "6.0.0"
     * For the ROOT release string is "ROOT major.minor.patch" like "ROOT 6.24.00" */
   JSROOT.version_id = "dev";

   /** @summary JSROOT version date
     * @desc Release date in format day/month/year like "14/01/2021"*/
   JSROOT.version_date = "10/02/2021";

   /** @summary JSROOT version id and date
     * @desc Produced by concatenation of {@link JSROOT.version_id} and {@link JSROOT.version_date}
     * Like "6.0.0 14/01/2021" */
   JSROOT.version = JSROOT.version_id + " " + JSROOT.version_date;

   /** @summary Location of JSROOT scripts
     * @desc Automatically detected and used to load other JSROOT scripts when required */
   JSROOT.source_dir = "";

   if (JSROOT.batch_mode === undefined)
      /** @summary Indicates if JSROOT runs in batch mode */
      JSROOT.batch_mode = false;

   if (JSROOT.nodejs === undefined)
      /** @summary Indicates if JSROOT runs inside Node.js */
      JSROOT.nodejs = false;

   //openuicfg // DO NOT DELETE, used to configure openui5 usage like JSROOT.openui5src = "nojsroot";

   /** @summary internal data
     * @memberof JSROOT
     * @private */
   let _ = {
      modules: {},            ///< list of modules
      source_min: false,      ///< is minified sources are used
      use_full_libs: false,   ///< is full libraries are used
      id_counter: 1           ///< unique id contner, starts from 1
   };

   let source_fullpath = "";
   let browser = { isOpera: false, isFirefox: true, isSafari: false, isChrome: false, isWin: false, touches: false  };

   if ((typeof document !== "undefined") && (typeof window !== "undefined")) {
      let script = document.currentScript;
      if (script && (typeof script.src == "string")) {
         const pos = script.src.indexOf("scripts/JSRoot.core.");
         if (pos >= 0) {
            source_fullpath = script.src;
            JSROOT.source_dir = source_fullpath.substr(0, pos);
            _.source_min = source_fullpath.indexOf("scripts/JSRoot.core.min.js") >= 0;
            console.log(`Set JSROOT.source_dir to ${JSROOT.source_dir}, ${JSROOT.version}`);
         }
      }

      browser.isOpera = !!window.opera || navigator.userAgent.indexOf(' OPR/') >= 0;
      browser.isFirefox = typeof InstallTrigger !== 'undefined';
      browser.isSafari = Object.prototype.toString.call(window.HTMLElement).indexOf('Constructor') > 0;
      browser.isChrome = !!window.chrome && !browser.isOpera;
      browser.isWin = navigator.platform.indexOf('Win') >= 0;
      browser.isChromeHeadless = navigator.userAgent.indexOf('HeadlessChrome') >= 0;
      browser.touches = ('ontouchend' in document); // identify if touch events are supported
   }

   _.sources = {
         'd3'                   : { src: 'd3', libs: true, extract: "d3", node: "d3" },
         'jquery'               : { src: 'jquery', libs: true,  extract: "$" },
         'jquery-ui'            : { src: 'jquery-ui', libs: true, extract: "$", dep: 'jquery' },
         'jqueryui-mousewheel'  : { src: 'jquery.mousewheel', onlymin: true, extract: "$", dep: 'jquery-ui' },
         'jqueryui-touch-punch' : { src: 'touch-punch', onlymin: true, extract: "$", dep: 'jquery-ui' },
         'rawinflate'           : { src: 'rawinflate', libs: true },
         'zstd-codec'           : { src: '../../zstd/zstd-codec.min', extract: "ZstdCodec", node: "zstd-codec" },
         'mathjax'              : { src: 'https://cdn.jsdelivr.net/npm/mathjax@3.1.2/es5/tex-svg', extract: "MathJax", node: "mathjax" },
         'dat.gui'              : { src: 'dat.gui', libs: true, extract: "dat" },
         'three'                : { src: 'three', libs: true, extract: "THREE", node: "three" },
         'threejs_jsroot'       : { src: 'three.extra', libs: true }
    };

    ['core', 'base3d','csg','geobase','geom','geoworker','gpad','hierarchy','hist','hist3d','interactive','io','menu','jq2d','latex',
      'math','more','openui5','painter','tree','v7gpad','v7hist','v7hist3d','v7more','webwindow']
         .forEach(item => _.sources[item] = { src: "JSRoot." + item });

   _.get_module_src = function(entry, fullyQualified) {
      if (entry.src.indexOf('http') == 0)
         return _.amd ? entry.src : entry.src + ".js";

      if (_.sap)
         return "jsroot/scripts/" + entry.src + ((_.source_min || entry.libs || entry.onlymin) ? ".min" : "");

      let dir = (entry.libs && _.use_full_libs && !_.source_min) ? JSROOT.source_dir + "libs/" : JSROOT.source_dir + "scripts/",
          ext = (_.source_min || (entry.libs && !_.use_full_libs) || entry.onlymin) ? ".min" : "";
      if (_.amd) return dir + entry.src + ext;
      let res = dir + entry.src + ext + ".js";

      if (fullyQualified && JSROOT.settings.NoCache) res += "?stamp=" + JSROOT.settings.NoCache;
      return res;
   }

   /** @summary Specialized JSROOT constants, used in {@link JSROOT.settings}
     * @namespace
     * @private */
   JSROOT.constants = {
      /** @summary Kind of 3D rendering, used for {@link JSROOT.settings.Render3D}
        * @namespace
        * @private */
      Render3D: {
         /** @summary Default 3D rendering, normally WebGL, if not supported - SVG*/
         Default: 0,
         /** @summary Use normal WebGL rendering and place as interactive Canvas element on HTML page */
         WebGL: 1,
         /** @summary Use WebGL rendering, but convert into svg image, not interactive */
         WebGLImage: 2,
         /** @summary Use SVG rendering, slow, inprecise and not interactive, nor recommendet */
         SVG: 3,
         fromString: function(s) {
            if ((s === "webgl") || (s == "gl")) return this.WebGL;
            if (s === "img") return this.WebGLImage;
            if (s === "svg") return this.SVG;
            return this.Default;
         }
      },
      /** @summary Way to embed 3D into SVG, used for {@link JSROOT.settings.Embed3D}
        * @namespace
        * @private */
      Embed3D: {
         /** @summary Do not embed 3D drawing, use complete space */
         NoEmbed: -1,
         /** @summary Default embeding mode, on Firefox is really ```Embed```, on all other ```Overlay``` */
         Default: 0,
         /** @summary WebGL canvas not inserted into SVG, but just overlayed The only way how Chrome browser can be used */
         Overlay: 1,
         /** @summary Really embed WebGL Canvas into SVG, only works with Firefox */
         Embed: 2,
         /** @summary Embeding, but when SVG rendering or SVG image converion is used
           * @private */
         EmbedSVG: 3,
         fromString: function(s) {
            if (s === "embed") return this.Embed;
            if (s === "overlay") return this.Overlay;
            return this.Default;
         }
      },
      /** @summary How to use latex in text drawing, used for {@link JSROOT.settings.Latex}
        * @namespace
        * @private */
      Latex: {
         /** @summary do not use Latex at all for text drawing */
         Off: 0,
         /** @summary convert only known latex symbols */
         Symbols: 1,
         /** @summary normal latex processing */
         Normal: 2,
         /** @summary use MathJax for complex cases, otherwise simple SVG text */
         MathJax: 3,
         /** @summary always use MathJax for text rendering */
         AlwaysMathJax: 4,
         fromString: function(s) {
            if (!s || (typeof s !== 'string'))
               return this.Normal;
            switch(s){
               case "off": return this.Off;
               case "symbols": return this.Symbols;
               case "MathJax":
               case "mathjax":
               case "math": return this.MathJax;
               case "AlwaysMathJax":
               case "alwaysmath":
               case "alwaysmathjax": return this.AlwaysMathJax;
            }
            let code = parseInt(s);
            return (!isNaN(code) && (code >= this.Off) && (code <= this.AlwaysMathJax)) ? code : this.Normal;
         }
      }
   };

   /** @summary Central JSROOT settings, independent from {@link JSROOT.gStyle}
     * @namespace */
   JSROOT.settings = {
      /** @summary Render of 3D drawing methods, see {@link JSROOT.constants.Render3D} for possible values */
      Render3D: JSROOT.constants.Render3D.Default,
      /** @summary Render of 3D drawing methods in batch mode, see {@link JSROOT.constants.Render3D} for possible values */
      Render3DBatch: JSROOT.constants.Render3D.Default,
      /** @summary Way to embed 3D drawing in SVG, see {@link JSROOT.constants.Embed3D} for possible values */
      Embed3D: JSROOT.constants.Embed3D.Default,
      /** @summary Enable or disable tooltips, default on */
      Tooltip: true,
      /** @summary Time in msec for appearance of tooltips, 0 - no animation */
      TooltipAnimation: 500,
      /** @summary Enables context menu usage */
      ContextMenu: true,
      /** @summary Global zooming flag, enable/disable any kind of interactive zooming */
      Zooming: true,
      /** @summary Zooming with the mouse events */
      ZoomMouse: true,
      /** @summary Zooming with mouse wheel */
      ZoomWheel: true,
      /** @summary Zooming on touch devices */
      ZoomTouch: true,
      /** @summary Enables move and resize of elements like statbox, title, pave, colz  */
      MoveResize: true,
      /** @summary enables drag and drop functionality */
      DragAndDrop: true,
      /** @summary Show progress box */
      ProgressBox: true,
      /** @summary Show additional tool buttons on the canvas, false - disabled, true - enabled, 'popup' - only toggle button */
      ToolBar: 'popup',
      /** @summary Position of toolbar 'left' left-bottom corner on canvas, 'right' - right-bottom corner on canvas, opposite on sub-pads */
      ToolBarSide: 'left',
      /** @summary display tool bar vertical (default false) */
      ToolBarVert: false,
      /** @summary if drawing inside particular div can be enlarged on full window */
      CanEnlarge: true,
      /** @summary if frame position can be adjusted to let show axis or colz labels */
      CanAdjustFrame: false,
      /** @summary calculation of text size consumes time and can be skipped to improve performance (but with side effects on text adjustments) */
      ApproxTextSize: false,
      /** @summary Histogram drawing optimization: 0 - disabled, 1 - only for large (>5000 1d bins, >50 2d bins) histograms, 2 - always */
      OptimizeDraw: 1,
      /** @summary Automatically create stats box, default on */
      AutoStat: true,
      /** @summary Default frame position in NFC */
      FrameNDC: { fX1NDC: 0.07, fY1NDC: 0.12, fX2NDC: 0.95, fY2NDC: 0.88 },
      /** @summary size of pad, where many features will be deactivated like text draw or zooming  */
      SmallPad: { width: 150, height: 100 },
      /** @summary Default color palette id  */
      Palette: 57,
      /** @summary Configures Latex usage, see {@link JSROOT.constants.Latex} for possible values */
      Latex: JSROOT.constants.Latex.Normal,
      /** @summary Grads per segment in TGeo spherical shapes like tube */
      GeoGradPerSegm: 6,
      /** @summary Enables faces compression after creation of composite shape  */
      GeoCompressComp: true,
      /** @summary if true, ignore all kind of URL options in the browser URL */
      IgnoreUrlOptions: false,
      /** @summary how many items shown on one level of hierarchy */
      HierarchyLimit: 250,
      /** @summary custom format for all X values, when not specified {@link JSROOT.gStyle.fStatFormat} is used */
      XValuesFormat : undefined,
      /** @summary custom format for all Y values, when not specified {@link JSROOT.gStyle.fStatFormat} is used */
      YValuesFormat : undefined,
      /** @summary custom format for all Z values, when not specified {@link JSROOT.gStyle.fStatFormat} is used */
      ZValuesFormat : undefined,
      /** @summary Let detect and solve problem when browser returns wrong content-length parameter
        * @desc See [jsroot#189]{@link https://github.com/root-project/jsroot/issues/189} for more info
        * Can be enabled by adding "wrong_http_response" parameter to URL when using JSROOT UI
        * @default false */
      HandleWrongHttpResponse: false,
      /** @summary Configures keybord key press handling
        * @desc Can be disabled to prevent keys heandling in complex HTML layouts
        * @default true */
      HandleKeys: true,
     /** @summary Let tweak browser caching
       * @desc When specified, extra URL parameter like ```?stamp=unique_value``` append to each JSROOT script loaded
       * In such case browser will be forced to load JSROOT functionality disregards of server cache settings
       * @default false */
      NoCache: false
   };

   /** @namespace
     * @memberof JSROOT
     * @summary Insiance of TStyle object like in ROOT
     * @desc Includes default draw styles, can be changed after loading of JSRoot.core.js
     * or can be load from the file providing style=itemname in the URL
     * See [TStyle docu]{@link https://root.cern/doc/master/classTStyle.html} "Private attributes" section for more detailed info about each value */
   let gStyle = {
      /** @summary Default log x scale */
      fOptLogx: 0,
      /** @summary Default log y scale */
      fOptLogy: 0,
      /** @summary Default log z scale */
      fOptLogz: 0,
      fOptDate: 0,
      fOptFile: 0,
      /** @summary Draw histogram title */
      fOptTitle: 1,
      fPadBottomMargin: 0.1,
      fPadTopMargin: 0.1,
      fPadLeftMargin: 0.1,
      fPadRightMargin: 0.1,
      /** @summary TPad.fGridx default value */
      fPadGridX: false,
      /** @summary TPad.fGridy default value */
      fPadGridY: false,
      fPadTickX: 0,
      fPadTickY: 0,
      /** @summary Default color of stat box */
      fStatColor: 0,
      /** @summary Default color of text in stat box */
      fStatTextColor: 1,
      fStatBorderSize: 1,
      fStatFont: 42,
      fStatFontSize: 0,
      fStatStyle: 1001,
      /** @summary Printing format for stats */
      fStatFormat: "6.4g",
      fStatX: 0.98,
      fStatY: 0.935,
      fStatW: 0.2,
      fStatH: 0.16,
      fTitleAlign: 23,
      fTitleColor: 0,
      fTitleTextColor: 1,
      fTitleBorderSize: 0,
      fTitleFont: 42,
      fTitleFontSize: 0.05,
      fTitleStyle: 0,
      /** @summary X position of top left corner of title box */
      fTitleX: 0.5,
      /** @summary Y position of top left corner of title box  */
      fTitleY: 0.995,
      /** @summary Width of title box */
      fTitleW: 0,
      /** @summary Height of title box */
      fTitleH: 0,
      /** @summary Printing format for fit parameters */
      fFitFormat: "5.4g",
      fOptStat: 1111,
      fOptFit: 0,
      fNumberContours: 20,
      fGridColor: 0,
      fGridStyle: 3,
      fGridWidth: 1,
      fFrameFillColor: 0,
      fFrameFillStyle: 1001,
      fFrameLineColor: 1,
      fFrameLineWidth: 1,
      fFrameLineStyle: 1,
      fFrameBorderSize: 1,
      fFrameBorderMode: 0,
      /** @summary size in pixels of end error for E1 draw options */
      fEndErrorSize: 2,
      /** @summary X size of the error marks for the histogram drawings */
      fErrorX: 0.5,
      /** @summary when true, BAR and LEGO drawing using base = 0  */
      fHistMinimumZero: false,
      /** @summary format for bin content */
      fPaintTextFormat: "g",
      /** @summary default time offset, UTC time at 01/01/95   */
      fTimeOffset: 788918400
   };

   /** @summary Method returns current document in use
     * @private */
   _.get_document = function() {
      if (JSROOT.nodejs)
         return _.nodejs_document;
      if (typeof document !== 'undefined')
         return document;
      if (typeof window == 'object')
         return window.document;
      return udefined;
   }

   function jsroot_require(need, factoryFunc) {

      if (!need && !factoryFunc)
         return Promise.resolve(null);

      if (typeof need == "string") need = need.split(";");

      need = need.filter(elem => !!elem);

      need.forEach((name,indx) => {
         if ((name.indexOf("load:")==0) || (name.indexOf("user:")==0))
            need[indx] = name.substr(5);
         else if (name == "jq")
            need[indx] = "jq2d";
         else if (name == "2d")
            need[indx] = "painter";
         else if (name == "v6")
            need[indx] = "gpad";
         else if (name == "v7")
            need[indx] = "v7gpad";
      });

      // loading with require.js

      if (_.amd) {
         if (!factoryFunc)
            return new Promise(resolve => {
               if (need.length > 0)
                  require(need, resolve);
               else
                  resolve();
            });

         if (need.length > 0)
            define(need, factoryFunc);
         else
            factoryFunc();
         return;
      }

      // loading inside node.js
      if (JSROOT.nodejs) {
         let arr = [];

         for (let k = 0; k < need.length; ++k) {
            let m = _.modules[need[k]];
            if (!m) {
               let src = _.sources[need[k]], modname;
               if (!src) throw Error("No module found " + need[k]);
               if (src.node)
                  modname = src.node;
               else if (src.libs)
                  modname = "./" + src.src + ".min.js";
               else
                  modname = "./" + src.src + ".js";
               let load = require(modname);
               if (load === undefined) load = 1;
               m = _.modules[need[k]] = { module: load };
            }
            arr.push(m.module);
         }

         if (factoryFunc) {
            factoryFunc(...arr);
            // TODO: how to access module.exports = res; of calling function
            return;
         }

         return Promise.resolve(arr.length == 1 ? arr[0] : arr);
      }

      // loading with sap.ui.require - but not mathjax
      if (_.sap && (need[0] != "mathjax")) {
         let req = [], reqindx = [], res = [];
         for (let k = 0; k < need.length; ++k) {
            let m = _.modules[need[k]];
            if (m) {
               res[k] = m.module;
            } else {
               let src = _.sources[need[k]];
               if (!src) {
                  req.push(need[k]);
               } else {
                  req.push(_.get_module_src(src));
               }
               reqindx.push(k);
            }
         }

         function decode_sap_results(args) {
            for (let i = 0; i < reqindx.length; ++i) {
               let k = reqindx[i]; // index in original request
               let src = _.sources[need[k]];
               if (src && src.extract) {
                  let m = _.modules[need[k]];
                  if (!m) m = _.modules[need[k]] = { module: globalThis[src.extract] };
                  res[k] = m.module;
               } else {
                  res[k] = args[i];
               }
            }
         }

         if (factoryFunc)
            return sap.ui.define(req, function() {
               decode_sap_results(arguments);
               return factoryFunc(...res);
            });
         else
            return new Promise(resolveFunc => {
               sap.ui.require(req, function() {
                  decode_sap_results(arguments);
                  resolveFunc(res.length == 1 ? res[0] : res);
               });
            });
      }

      let thisModule, thisSrc;

      if (factoryFunc) {
         if (document.currentScript) {
            thisSrc = document.currentScript.src;
            let separ = (typeof thisSrc == 'string') ? thisSrc.indexOf('?') : -1;
            if (separ > 0) thisSrc = thisSrc.substr(0, separ);
            thisModule = thisSrc;
            for (let mod in _.sources)
               if (_.get_module_src(_.sources[mod]) == thisSrc) {
                  thisModule = mod; break;
               }
         }
         if (!thisModule)
            throw Error("Cannot define module for " + (thisSrc || "uncknown script"));
      }

      // direct loading

      function finish_loading(m, res) {
         // check if promise was returned
         if (res && (typeof res == 'object') && (typeof res.then == 'function'))
            return res.then(pres => finish_loading(m, pres));

         m.module = res || 1; // just to have some value
         let waiting = m.waiting;
         delete m.loading; // clear loading flag
         delete m.waiting;

         if (waiting) waiting.forEach(func => func(true));
      }

      function handle_func(req, is_ok) {
         if (req.processed) return;
         if (!is_ok) return req.failed();
         let arr = [];
         for (let k = 0; k < req.need.length; ++k) {
            let m = _.modules[req.need[k]];
            if (!m) return req.failed();
            if (m.module === undefined) return; // not yet ready
            arr.push(m.module);
         }

         req.processed = true;

         if (req.thisModule) {

            let m = _.modules[req.thisModule], res;

            if (req.factoryFunc)
               res = req.factoryFunc(...arr);

            finish_loading(m, res);
         }

         if (req.resolve)
             req.resolve(arr.length == 1 ? arr[0] : arr);
      }

      function load_module(req, m) {
         let element = document.createElement("script");
         element.setAttribute('type', "text/javascript");
         element.setAttribute('src', m.src);
         document.getElementsByTagName("head")[0].appendChild(element);
         if (_.debug_output)
            _.debug_output.innerHTML = `Loading ${m.src} ...`;

         if (!m.jsroot || m.extract)
            element.onload = () => finish_loading(m, m.extract ? globalThis[m.extract] : 1); // mark script loaded
         element.onerror = () => { m.failure = true; req.failed(); }
      }

      function after_depend_load(req,d,m) {
         if (d.module === undefined)
            return req.failed("Dependend fail to load");

         if (!m.waiting) m.waiting = [];
         m.waiting.push(handle_func.bind(this, req));

         if (!m.loading) {
            m.loading = true;
            load_module(req,m);
         }
      }

      function analyze(resolve, reject) {
         let handler, any_dep, srcs = [],
             req = { need: need, thisModule: thisModule, factoryFunc: factoryFunc, resolve: resolve, reject: reject,
                     failed: function(msg) { this.processed = true; if (this.reject) this.reject(Error(msg || "JSROOT.require failed")); } };

         if (req.factoryFunc && req.thisModule) {

            let m = _.modules[req.thisModule];
            if (!m)
               m = _.modules[req.thisModule] = { jsroot: true, src: thisSrc, loading: true };
         }

         for (let k = 0; k < need.length; ++k) {
            let m = _.modules[need[k]];

            if (m && (m.module !== undefined)) continue;

            if (!m) {
               let jsmodule = _.sources[need[k]];

               m = _.modules[need[k]] = {};
               if (jsmodule) {
                  m.jsroot = true;
                  m.src = _.get_module_src(jsmodule, true);
                  m.extract = jsmodule.extract;
                  m.dep = jsmodule.dep; // copy dependence
              } else {
                  m.src = need[k];
               }
            }

            if (m.failure)
               // module loading failed, no nee to continue
               return req.failed(`Loading of module ${need[k]} failed`);

            if (m.dep) {
               let d = _.modules[m.dep];
               if (!d)
                  return req.failed(`Dependend module ${m.dep} not found`);
               if (d.module === undefined) {
                  any_dep = true;
                  if (!d.waiting) d.waiting = [];
                  d.waiting.push(after_depend_load.bind(this,req,d,m));
                  continue;
               }
            }

            if (!m.loading) {
               m.loading = true;
               srcs.push(m);
            }
            if (!m.waiting) m.waiting = [];
            if (!handler) handler = handle_func.bind(this, req);
            m.waiting.push(handler);
         }

         if (!handler && !any_dep)
            return handle_func(req, true);

         srcs.forEach(m => load_module(req,m));
      }

      if (factoryFunc)
         analyze();
      else
         return new Promise(function(resolve,reject) {
            analyze(resolve,reject);
         });
   }

   /** @summary Central method to load JSROOT functionality
     * @desc
     * Following components can be specified
     *    - 'io'     TFile functionality
     *    - 'tree'   TTree support
     *    - 'painter' d3.js plus basic painting functions
     *    - 'gpad'   basic 2d graphic (TCanvas/TPad/TFrame)
     *    - 'hist'   histograms 2d drawing (SVG)
     *    - 'hist3d' histograms 3d drawing (WebGL)
     *    - 'more'   extra 2d graphic (TGraph, TF1)
     *    - 'geom'   TGeo support
     *    - 'v7gpad' ROOT v7 RPad/RCanvas/RFrame
     *    - 'v7hist' ROOT v7 histograms 2d drawing (SVG)
     *    - 'v7hist3d' ROOT v7 histograms 3d drawing (WebGL)
     *    - 'v7more' ROOT v7 special classes
     *    - 'math'   some methods from TMath class
     *    - 'hierarchy' hierarchy browser
     *    - 'jq2d'   jQuery-dependent part of hierarchy
     *    - 'openui5' OpenUI5 and related functionality
     * @param {Array|string} req - list of required components (as array or string separated by semicolon)
     * @returns {Promise} with array of requirements (or single element) */
   JSROOT.require = function(req) {
      return jsroot_require(req);
   }

   /** @summary Define JSROOT module
     * @desc Should be only used for JSROOT modules
     * @param {Array|string} req - requirements, see {@link JSROOT.require} for more details
     * @param {Function} factoryFunc - called when requirements are fulfilled, with requested modules
     * @private */
   JSROOT.define = function(req, factoryFunc) {
      jsroot_require(req, factoryFunc);
   }

   /** @summary Generate mask for given bit
     * @param {number} n bit number
     * @returns {Number} produced mask
     * @private */
   JSROOT.BIT = function(n) { return 1 << n; }

   /** @summary Seed simple random generator
     * @param {number} i seed value
     * @private */
   JSROOT.seed = function(i) {
      i = Math.abs(i);
      if (i > 1e8) i = Math.abs(1e8 * Math.sin(i)); else
      if (i < 1) i*=1e8;
      this.m_w = Math.round(i);
      this.m_z = 987654321;
   }

   /** @summary Simple random generator
     * @desc Works like Math.random(), but with configurable seed - see {@link JSROOT.seed}
     * @returns {number} random value between 0 (inclusive) and 1.0 (exclusive)
     * @private */
   JSROOT.random = function() {
      if (this.m_z===undefined) return Math.random();
      this.m_z = (36969 * (this.m_z & 65535) + (this.m_z >> 16)) & 0xffffffff;
      this.m_w = (18000 * (this.m_w & 65535) + (this.m_w >> 16)) & 0xffffffff;
      let result = ((this.m_z << 16) + this.m_w) & 0xffffffff;
      result /= 4294967296;
      return result + 0.5;
   }

   /** @summary Just copy (not clone) all fields from source to the target object
     * @desc Simple replacement of jQuery.extend method
     * @memberof JSROOT
     * @private */
   let extend = (tgt, src) => {
      if ((src === null) || (typeof src !== 'object')) return tgt;
      if ((tgt === null) || (typeof tgt !== 'object')) tgt = {};

      for (let k in src)
         tgt[k] = src[k];

      return tgt;
   }

   /** @summary Adds specific methods to the object.
     * @desc JSROOT implements some basic methods for different ROOT classes.
     * @param {object} obj - object where methods are assigned
     * @param {string} [typename] - optional typename, if not specified, obj._typename will be used
     * @memberof JSROOT
     * @private */
   let addMethods = (obj, typename) => {
      extend(obj, JSROOT.getMethods(typename || obj._typename, obj));
   }

   /** @summary Should be used to parse JSON string produced with TBufferJSON class
     * @desc Replace all references inside object like { "$ref": "1" }
     * @param {object|string} json  object where references will be replaced
     * @returns {object} parsed object */
   JSROOT.parse = function(json) {

      if (!json) return null;

      let obj = (typeof json == 'string') ? JSON.parse(json) : json,
          map = [], newfmt = undefined;

      function unref_value(value) {
         if ((value===null) || (value===undefined)) return;

         if (typeof value === 'string') {
            if (newfmt || (value.length < 6) || (value.indexOf("$ref:") !== 0)) return;
            let ref = parseInt(value.substr(5));
            if (isNaN(ref) || (ref < 0) || (ref >= map.length)) return;
            newfmt = false;
            return map[ref];
         }

         if (typeof value !== 'object') return;

         let proto = Object.prototype.toString.apply(value);

         // scan array - it can contain other objects
         if ((proto.indexOf('[object')==0) && (proto.indexOf('Array]')>0)) {
             for (let i = 0; i < value.length; ++i) {
                let res = unref_value(value[i]);
                if (res!==undefined) value[i] = res;
             }
             return;
         }

         let ks = Object.keys(value), len = ks.length;

         if ((newfmt!==false) && (len===1) && (ks[0]==='$ref')) {
            const ref = parseInt(value['$ref']);
            if (isNaN(ref) || (ref < 0) || (ref >= map.length)) return;
            newfmt = true;
            return map[ref];
         }

         if ((newfmt!==false) && (len>1) && (ks[0]==='$arr') && (ks[1]==='len')) {
            // this is ROOT-coded array
            let arr = null, dflt = (value.$arr==="Bool") ? false : 0;
            switch (value.$arr) {
               case "Int8": arr = new Int8Array(value.len); break;
               case "Uint8": arr = new Uint8Array(value.len); break;
               case "Int16": arr = new Int16Array(value.len); break;
               case "Uint16": arr = new Uint16Array(value.len); break;
               case "Int32": arr = new Int32Array(value.len); break;
               case "Uint32": arr = new Uint32Array(value.len); break;
               case "Float32": arr = new Float32Array(value.len); break;
               case "Int64":
               case "Uint64":
               case "Float64": arr = new Float64Array(value.len); break;
               default: arr = new Array(value.len); break;
            }
            for (let k=0;k<value.len;++k) arr[k] = dflt;

            if (value.b !== undefined) {
               // base64 coding

               let atob_func = JSROOT.nodejs ? require('atob') : window.atob;

               let buf = atob_func(value.b);

               if (arr.buffer) {
                  let dv = new DataView(arr.buffer, value.o || 0),
                      len = Math.min(buf.length, dv.byteLength);
                  for (let k=0; k<len; ++k)
                     dv.setUint8(k, buf.charCodeAt(k));
               } else {
                  throw new Error('base64 coding supported only for native arrays with binary data');
               }
            } else {
               // compressed coding
               let nkey = 2, p = 0;
               while (nkey<len) {
                  if (ks[nkey][0]=="p") p = value[ks[nkey++]]; // position
                  if (ks[nkey][0]!=='v') throw new Error('Unexpected member ' + ks[nkey] + ' in array decoding');
                  let v = value[ks[nkey++]]; // value
                  if (typeof v === 'object') {
                     for (let k = 0; k < v.length; ++k) arr[p++] = v[k];
                  } else {
                     arr[p++] = v;
                     if ((nkey<len) && (ks[nkey][0]=='n')) {
                        let cnt = value[ks[nkey++]]; // counter
                        while (--cnt) arr[p++] = v;
                     }
                  }
               }
            }

            return arr;
         }

         if ((newfmt!==false) && (len===3) && (ks[0]==='$pair') && (ks[1]==='first') && (ks[2]==='second')) {
            newfmt = true;
            let f1 = unref_value(value.first),
                s1 = unref_value(value.second);
            if (f1!==undefined) value.first = f1;
            if (s1!==undefined) value.second = s1;
            value._typename = value['$pair'];
            delete value['$pair'];
            return; // pair object is not counted in the objects map
         }

         // add object to object map
         map.push(value);

         // add methods to all objects, where _typename is specified
         if (value._typename) addMethods(value);

         for (let k = 0; k < len; ++k) {
            const i = ks[k],
                 res = unref_value(value[i]);
            if (res !== undefined) value[i] = res;
         }
      }

      unref_value(obj);

      return obj;
   }

   /** @summary Make deep clone of the object, including all sub-objects
     * @returns {object} cloned object
     * @private */
   JSROOT.clone = function(src, map, nofunc) {
      if (!src) return null;

      if (!map) {
         map = { obj: [], clones: [], nofunc: nofunc };
      } else {
         const i = map.obj.indexOf(src);
         if (i >= 0) return map.clones[i];
      }

      let proto = Object.prototype.toString.apply(src);

      // process normal array
      if (proto === '[object Array]') {
         let tgt = [];
         map.obj.push(src);
         map.clones.push(tgt);
         for (let i = 0; i < src.length; ++i)
            if (typeof src[i] === 'object')
               tgt.push(JSROOT.clone(src[i], map));
            else
               tgt.push(src[i]);

         return tgt;
      }

      // process typed array
      if ((proto.indexOf('[object ') == 0) && (proto.indexOf('Array]') == proto.length-6)) {
         let tgt = [];
         map.obj.push(src);
         map.clones.push(tgt);
         for (let i = 0; i < src.length; ++i)
            tgt.push(src[i]);

         return tgt;
      }

      let tgt = {};
      map.obj.push(src);
      map.clones.push(tgt);

      for (let k in src) {
         if (typeof src[k] === 'object')
            tgt[k] = JSROOT.clone(src[k], map);
         else if (!map.nofunc || (typeof src[k]!=='function'))
            tgt[k] = src[k];
      }

      return tgt;
   }

   /** @summary Parse response from multi.json request
     * @desc Method should be used to parse JSON code, produced by multi.json request of THttpServer
     * @param {string} json string to parse
     * @returns {Array} array of parsed elements */
   JSROOT.parseMulti = function(json) {
      if (!json) return null;
      let arr = JSON.parse(json);
      if (arr && arr.length)
         for (let i=0;i<arr.length;++i)
            arr[i] = JSROOT.parse(arr[i]);
      return arr;
   }

   /** @summary Method converts JavaScript object into ROOT-like JSON
     * @desc Produced JSON can be used in JSROOT.parse() again
     * When performed properly, JSON can be used in [TBufferJSON::fromJSON()]{@link https://root.cern/doc/master/classTBufferJSON.html#a2ecf0daacdad801e60b8093a404c897d} method to read data back with C++
     * @param {object} obj - JavaScript object to convert
     * @param {number} [spacing] - optional line spacing in JSON
     * @returns {string} produced JSON code */
   JSROOT.toJSON = function(obj, spacing) {
      if (!obj || typeof obj !== 'object') return "";

      let map = []; // map of stored objects

      function copy_value(value) {
         if (typeof value === "function") return undefined;

         if ((value===undefined) || (value===null) || (typeof value !== 'object')) return value;

         let proto = Object.prototype.toString.apply(value);

         // typed array need to be converted into normal array, otherwise looks strange
         if ((proto.indexOf('[object ') == 0) && (proto.indexOf('Array]') == proto.length-6)) {
            let arr = new Array(value.length)
            for (let i = 0; i < value.length; ++i)
               arr[i] = copy_value(value[i]);
            return arr;
         }

         // this is how reference is code
         let refid = map.indexOf(value);
         if (refid >= 0) return { $ref: refid };

         let ks = Object.keys(value), len = ks.length, tgt = {};

         if ((len == 3) && (ks[0]==='$pair') && (ks[1]==='first') && (ks[2]==='second')) {
            // special handling of pair objects which does not included into objects map
            tgt.$pair = value.$pair;
            tgt.first = copy_value(value.first);
            tgt.second = copy_value(value.second);
            return tgt;
         }

         map.push(value);

         for (let k = 0; k < len; ++k) {
            let name = ks[k];
            tgt[name] = copy_value(value[name]);
         }

         return tgt;
      }

      let tgt = copy_value(obj);

      return JSON.stringify(tgt, null, spacing);
   }

   /** @summary decodes URL options after '?' mark
     * @desc Following options supported ?opt1&opt2=3
     * @param {string} [url] URL string with options, document.URL will be used when not specified
     * @returns {Object} with ```.has(opt)``` and ```.get(opt,dflt)``` methods
     * @example
     * let d = JSROOT.decodeUrl("any?opt1&op2=3");
     * console.log(`Has opt1 ${d.has("opt1")}`);     // true
     * console.log(`Get opt1 ${d.get("opt1")}`);     // ""
     * console.log(`Get opt2 ${d.get("opt2")}`);     // "3"
     * console.log(`Get opt3 ${d.get("opt3","-")}`); // "-" */
   JSROOT.decodeUrl = function(url) {
      let res = {
         opts: {},
         has: function(opt) { return this.opts[opt] !== undefined; },
         get: function(opt,dflt) { let v = this.opts[opt]; return v!==undefined ? v : dflt; }
      }

      if (!url || (typeof url !== 'string')) {
         if (JSROOT.settings.IgnoreUrlOptions || (typeof document === 'undefined')) return res;
         url = document.URL;
      }
      res.url = url;

      let p1 = url.indexOf("?");
      if (p1 < 0) return res;
      url = decodeURI(url.slice(p1+1));

      while (url.length > 0) {

         // try to correctly handle quotes in the URL
         let pos = 0, nq = 0, eq = -1, firstq = -1;
         while ((pos < url.length) && ((nq!==0) || ((url[pos]!=="&") && (url[pos]!=="#")))) {
            switch (url[pos]) {
               case "'": if (nq >= 0) nq = (nq+1)%2; if (firstq < 0) firstq = pos; break;
               case '"': if (nq <= 0) nq = (nq-1)%2; if (firstq < 0) firstq = pos; break;
               case '=': if ((firstq < 0) && (eq < 0)) eq = pos; break;
            }
            pos++;
         }

         if ((eq < 0) && (firstq < 0)) {
            res.opts[url.substr(0,pos)] = "";
         } if (eq > 0) {
            let val = url.slice(eq+1, pos);
            if (((val[0]==="'") || (val[0]==='"')) && (val[0]===val[val.length-1])) val = val.substr(1, val.length-2);
            res.opts[url.substr(0,eq)] = val;
         }

         if ((pos >= url.length) || (url[pos] == '#')) break;

         url = url.substr(pos+1);
      }

      return res;
   }

   /** @summary Find function with given name.
     * @desc Function name may include several namespaces like 'JSROOT.Painter.drawFrame'
     * If function starts with ., it should belong to JSROOT.Painter
     * @private */
   JSROOT.findFunction = function(name) {
      if (typeof name === 'function') return name;
      if (typeof name !== 'string') return null;
      let names, elem;
      if (name[0] == ".") { // special shortcut for JSROOT.Painter.
         names = name.substr(1).split('.');
         elem = JSROOT.Painter;
      } else {
         names = name.split('.');
         if (names[0] === 'JSROOT') {
            elem = JSROOT;
            names.shift();
         } else {
            elem = globalThis;
         }
      }

      for (let n = 0; elem && (n < names.length); ++n)
         elem = elem[names[n]];

      return (typeof elem == 'function') ? elem : null;
   }

   /** @summary Old request method, kept only for internal use
     * @private */
   JSROOT.NewHttpRequest = function(url, kind, user_accept_callback, user_reject_callback) {

      let xhr = JSROOT.nodejs ? new (require("xhr2"))() : new XMLHttpRequest();

      xhr.http_callback = (typeof user_accept_callback == 'function') ? user_accept_callback.bind(xhr) : function() {};
      xhr.error_callback = (typeof user_reject_callback == 'function') ? user_reject_callback.bind(xhr) : function(err) { console.warn(err.message); this.http_callback(null); }.bind(xhr);

      if (!kind) kind = "buf";

      let method = "GET", async = true, p = kind.indexOf(";sync");
      if (p > 0) { kind = kind.substr(0,p); async = false; }
      switch (kind) {
         case "head": method = "HEAD"; break;
         case "posttext": method = "POST"; kind = "text"; break;
         case "postbuf":  method = "POST"; kind = "buf"; break;
         case "post":
         case "multi":  method = "POST"; kind = buf; break;
      }

      xhr.kind = kind;

      if (JSROOT.settings.HandleWrongHttpResponse && (method == "GET") && (typeof xhr.addEventListener === 'function'))
         xhr.addEventListener("progress", function(oEvent) {
            if (oEvent.lengthComputable && this.expected_size && (oEvent.loaded > this.expected_size)) {
               this.did_abort = true;
               this.abort();
               this.error_callback(Error('Server sends more bytes ' + oEvent.loaded + ' than expected ' + this.expected_size + '. Abort I/O operation'), 598);
            }
         }.bind(xhr));

      xhr.onreadystatechange = function() {

         if (this.did_abort) return;

         if ((this.readyState === 2) && this.expected_size) {
            let len = parseInt(this.getResponseHeader("Content-Length"));
            if (!isNaN(len) && (len > this.expected_size) && !JSROOT.settings.HandleWrongHttpResponse) {
               this.did_abort = true;
               this.abort();
               return this.error_callback(Error('Server response size ' + len + ' larger than expected ' + this.expected_size + '. Abort I/O operation'), 599);
            }
         }

         if (this.readyState != 4) return;

         if ((this.status != 200) && (this.status != 206) && !browser.qt5 &&
             // in these special cases browsers not always set status
             !((this.status == 0) && ((url.indexOf("file://")==0) || (url.indexOf("blob:")==0)))) {
               return this.error_callback(Error('Fail to load url ' + url), this.status);
         }

         if (this.nodejs_checkzip && (this.getResponseHeader("content-encoding") == "gzip")) {
            // special handling of gzipped JSON objects in Node.js
            let zlib = require('zlib'),
                res = zlib.unzipSync(Buffer.from(this.response)),
                obj = JSON.parse(res); // zlib returns Buffer, use JSON to parse it
            return this.http_callback(JSROOT.parse(obj));
         }

         switch(this.kind) {
            case "xml": return this.http_callback(this.responseXML);
            case "text": return this.http_callback(this.responseText);
            case "object": return this.http_callback(JSROOT.parse(this.responseText));
            case "multi": return this.http_callback(JSROOT.parseMulti(this.responseText));
            case "head": return this.http_callback(this);
         }

         // if no response type is supported, return as text (most probably, will fail)
         if (this.responseType === undefined)
            return this.http_callback(this.responseText);

         if ((this.kind == "bin") && ('byteLength' in this.response)) {
            // if string representation in requested - provide it

            let filecontent = "", u8Arr = new Uint8Array(this.response);
            for (let i = 0; i < u8Arr.length; ++i)
               filecontent += String.fromCharCode(u8Arr[i]);

            return this.http_callback(filecontent);
         }

         this.http_callback(this.response);
      }

      xhr.open(method, url, async);

      if ((kind == "bin") || (kind == "buf")) xhr.responseType = 'arraybuffer';

      if (JSROOT.nodejs && (method == "GET") && (kind === "object") && (url.indexOf('.json.gz')>0)) {
         xhr.nodejs_checkzip = true;
         xhr.responseType = 'arraybuffer';
      }

      return xhr;
   }

   /** @summary Submit asynchronoues http request
     * @desc Following requests kind can be specified:
     *    - "bin" - abstract binary data, result as string
     *    - "buf" - abstract binary data, result as ArrayBuffer (default)
     *    - "text" - returns req.responseText
     *    - "object" - returns JSROOT.parse(req.responseText)
     *    - "multi" - returns correctly parsed multi.json request
     *    - "xml" - returns req.responseXML
     *    - "head" - returns request itself, uses "HEAD" request method
     *    - "post" - creates post request, submits req.send(post_data)
     *    - "postbuf" - creates post request, expectes binary data as response
     * @param {string} url - URL for the request
     * @param {string} kind - kind of requested data
     * @param {string} [post_data] - data submitted with post kind of request
     * @returns {Promise} Promise for requested data, result type depends from the kind
     * @example
     * JSROOT.httpRequest("https://root.cern/js/files/thstack.json.gz", "object")
     *       .then(obj => console.log(`Get object of type ${obj._typename}`))
     *       .catch(err => console.error(err.message)); */
   JSROOT.httpRequest = function(url, kind, post_data) {
      return new Promise(function(accept, reject) {
         JSROOT.NewHttpRequest(url, kind, accept, reject).send(post_data || null);
      });
   }

   /** @summary Load script or CSS file into the browser
     * @desc Normal JSROOT functionality should be loaded via {@link JSROOT.require} method
     * @param {String} url - script or css file URL (or array, in this case they all loaded secuentially)
     * @returns {Promise} */
   JSROOT.loadScript = function(url) {
      if (!url)
         return Promise.resolve(true);

      if (typeof url != 'string') {
         let scripts = url, loadNext = () => {
            if (!scripts.length) return Promise.resolve(true);
            return JSROOT.loadScript(scripts.shift()).then(loadNext, loadNext);
         }
         return loadNext();
      }

      if (url.indexOf("$$$")===0) {
         url = url.slice(3);
         if ((url.indexOf("style/")==0) && (url.indexOf('.css') < 0))
            url += _.source_min ? '.min.css' : ".css";
         url = JSROOT.source_dir + url;
      }

      let element, isstyle = url.indexOf(".css") > 0;

      if (JSROOT.nodejs) {
         let res = null;
         if (!isstyle) {
            if ((url.indexOf("http:") == 0) || (url.indexOf("https:") == 0))
               return JSROOT.httpRequest(url,"text").then(txt => eval(txt));
            res = require(url);
         }

         return Promise.resolve(res);
      }

      function match_url(src) {
         if (src == url) return true;
         let indx = src.indexOf(url);
         return (indx > 0) && (indx + url.length == src.length) && (src[indx-1] == "/");
      }

      if (isstyle) {
         let styles = document.getElementsByTagName('link');
         for (let n = 0; n < styles.length; ++n) {
            if (!styles[n].href || (styles[n].type !== 'text/css') || (styles[n].rel !== 'stylesheet')) continue;
            if (match_url(styles[n].href))
               return Promise.resolve();
         }

      } else {
         let scripts = document.getElementsByTagName('script');
         for (let n = 0; n < scripts.length; ++n)
            if (match_url(scripts[n].src))
               return Promise.resolve();
      }

      if (isstyle) {
         element = document.createElement("link");
         element.setAttribute("rel", "stylesheet");
         element.setAttribute("type", "text/css");
         element.setAttribute("href", url);
      } else {
         element = document.createElement("script");
         element.setAttribute('type', "text/javascript");
         element.setAttribute('src', url);
      }

      return new Promise((resolve, reject) => {
         element.onload = () => resolve(true);
         element.onerror = () => reject(Error(`Fail to load ${url}`));
         document.getElementsByTagName("head")[0].appendChild(element);
      });
   }

   // Open ROOT file, defined in JSRoot.io.js
   JSROOT.openFile = (filename, callback) => {
      return jsroot_require("io").then(() => JSROOT.openFile(filename, callback));
   }

   // Draw object, defined in JSRoot.painter.js
   JSROOT.draw = (divid, obj, opt) => {
      return jsroot_require("painter").then(() => JSROOT.draw(divid, obj, opt));
   }

   // Redaraw object, defined in JSRoot.painter.js
   JSROOT.redraw = (divid, obj, opt) => {
      return jsroot_require("painter").then(() => JSROOT.redraw(divid, obj, opt));
   }

   // Dummy, when painter is not yet loaded, should happens nothing
   JSROOT.cleanup = () => {}

   // Create SVG, defined in JSRoot.painter.js
   JSROOT.makeSVG = args => {
      return jsroot_require("painter").then(() => JSROOT.makeSVG(args));
   }

   /** @summary Method to build main JSROOT GUI
     * @param {string|object} gui_element - where gui should be started, element id or just HTML element
     * @param {string} [gui_kind = "gui"] - kind of the gui: "gui", "online", "draw"
     * @returns {Promise} when ready
     * @private */
   JSROOT.buildGUI = (gui_element, gui_kind) => {
      let d = JSROOT.decodeUrl(),
          nobrowser = d.has('nobrowser'),
          requirements = ["hierarchy"];

      if (typeof gui_element == 'string')
         gui_element = document.getElementById(gui_element);

      if (!gui_element) {
         console.log('Fail to find element for GUI drawing');
         return Promise.resolve(false);
      }

      if (d.has('libs')) _.use_full_libs = true;

      if (gui_kind == "nobrowser") {
          gui_kind = "gui"; nobrowser = true;
      } else if (gui_kind == "draw") {
         nobrowser = true;
      } else if (gui_kind != "online") {
         gui_kind = "gui";
      }

      if (!nobrowser) requirements.push("jq2d");

      let user_scripts = d.get("autoload") || d.get("load");

      if (user_scripts) requirements.push("painter");

      _.debug_output = gui_element;

      return jsroot_require(requirements).then(() => JSROOT.loadScript(user_scripts)).then(() => {
         gui_element.innerHTML = "";
         delete _.debug_output;
         return nobrowser ? JSROOT.buildNobrowserGUI(gui_element, gui_kind) : JSROOT.buildGUI(gui_element, gui_kind);
      });
   }

   /** @summary Create some ROOT classes
     * @desc Supported classes: "TObject", "TNamed", "TList", "TAxis", "TLine", "TText", "TLatex", "TPad", "TCanvas"
     * @param {string} typename - ROOT class name
     * @memberof JSROOT
     * @example
     * let obj = JSROOT.create("TNamed");
     * obj.fName = "name";
     * obj.fTitle = "title"; */
   let create = (typename, target) => {
      let obj = target || {};

      switch (typename) {
         case 'TObject':
             extend(obj, { fUniqueID: 0, fBits: 0 });
             break;
         case 'TNamed':
            extend(obj, { fUniqueID: 0, fBits: 0, fName: "", fTitle: "" });
            break;
         case 'TList':
         case 'THashList':
            extend(obj, { name: typename, arr : [], opt : [] });
            break;
         case 'TAttAxis':
            extend(obj, { fNdivisions: 510, fAxisColor: 1,
                                 fLabelColor: 1, fLabelFont: 42, fLabelOffset: 0.005, fLabelSize: 0.035, fTickLength: 0.03,
                                 fTitleOffset: 1, fTitleSize: 0.035, fTitleColor: 1, fTitleFont : 42 });
            break;
         case 'TAxis':
            create("TNamed", obj);
            create("TAttAxis", obj);
            extend(obj, { fNbins: 1, fXmin: 0, fXmax: 1, fXbins : [], fFirst: 0, fLast: 0,
                                 fBits2: 0, fTimeDisplay: false, fTimeFormat: "", fLabels: null, fModLabs: null });
            break;
         case 'TAttLine':
            extend(obj, { fLineColor: 1, fLineStyle: 1, fLineWidth: 1 });
            break;
         case 'TAttFill':
            extend(obj, { fFillColor: 0, fFillStyle: 0 } );
            break;
         case 'TAttMarker':
            extend(obj, { fMarkerColor: 1, fMarkerStyle: 1, fMarkerSize: 1. });
            break;
         case 'TLine':
            create("TObject", obj);
            create("TAttLine", obj);
            extend(obj, { fX1: 0, fX2: 1, fY1: 0, fY2: 1 });
            break;
         case 'TBox':
            create("TObject", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            extend(obj, { fX1: 0, fX2: 1, fY1: 0, fY2: 1 });
            break;
         case 'TPave':
            create("TBox", obj);
            extend(obj, { fX1NDC : 0., fY1NDC: 0, fX2NDC: 1, fY2NDC: 1,
                                 fBorderSize: 0, fInit: 1, fShadowColor: 1,
                                 fCornerRadius: 0, fOption: "blNDC", fName: "title" });
            break;
         case 'TAttText':
            extend(obj, { fTextAngle: 0, fTextSize: 0, fTextAlign: 22, fTextColor: 1, fTextFont: 42});
            break;
         case 'TPaveText':
            create("TPave", obj);
            create("TAttText", obj);
            extend(obj, { fLabel: "", fLongest: 27, fMargin: 0.05, fLines: create("TList") });
            break;
         case 'TPaveStats':
            create("TPaveText", obj);
            extend(obj, { fOptFit: 0, fOptStat: 0, fFitFormat: "", fStatFormat: "", fParent: null });
            break;
         case 'TLegend':
            create("TPave", obj);
            create("TAttText", obj);
            extend(obj, { fColumnSeparation: 0, fEntrySeparation: 0.1, fMargin: 0.25, fNColumns: 1, fPrimitives: create("TList") });
            break;
         case 'TLegendEntry':
            create("TObject", obj);
            create("TAttText", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            create("TAttMarker", obj);
            extend(obj, { fLabel: "", fObject: null, fOption: "" });
            break;
         case 'TText':
            create("TNamed", obj);
            create("TAttText", obj);
            extend(obj, { fLimitFactorSize: 3, fOriginSize: 0.04 });
            break;
         case 'TLatex':
            create("TText", obj);
            create("TAttLine", obj);
            extend(obj, { fX: 0, fY: 0 });
            break;
         case 'TObjString':
            create("TObject", obj);
            extend(obj, { fString: "" });
            break;
         case 'TH1':
            create("TNamed", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            create("TAttMarker", obj);

            extend(obj, {
               fBits: 8,
               fNcells: 0,
               fXaxis: create("TAxis"),
               fYaxis: create("TAxis"),
               fZaxis: create("TAxis"),
               fBarOffset: 0, fBarWidth: 1000, fEntries: 0.,
               fTsumw: 0., fTsumw2: 0., fTsumwx: 0., fTsumwx2: 0.,
               fMaximum: -1111., fMinimum: -1111, fNormFactor: 0., fContour: [],
               fSumw2: [], fOption: "",
               fFunctions: create("TList"),
               fBufferSize: 0, fBuffer: [], fBinStatErrOpt: 0, fStatOverflows: 2 });
            break;
         case 'TH1I':
         case 'TH1F':
         case 'TH1D':
         case 'TH1S':
         case 'TH1C':
            create("TH1", obj);
            obj.fArray = [];
            break;
         case 'TH2':
            create("TH1", obj);
            extend(obj, { fScalefactor: 1., fTsumwy: 0.,  fTsumwy2: 0, fTsumwxy: 0 });
            break;
         case 'TH2I':
         case 'TH2F':
         case 'TH2D':
         case 'TH2S':
         case 'TH2C':
            create("TH2", obj);
            obj.fArray = [];
            break;
         case 'TH3':
            create("TH1", obj);
            extend(obj, { fTsumwy: 0.,  fTsumwy2: 0, fTsumwz: 0.,  fTsumwz2: 0, fTsumwxy: 0, fTsumwxz: 0, fTsumwyz: 0 });
            break;
         case 'TH3I':
         case 'TH3F':
         case 'TH3D':
         case 'TH3S':
         case 'TH3C':
            create("TH3", obj);
            obj.fArray = [];
            break;
         case 'THStack':
            create("TNamed", obj);
            extend(obj, { fHists: create("TList"), fHistogram: null, fMaximum: -1111, fMinimum: -1111 });
            break;
         case 'TGraph':
            create("TNamed", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            create("TAttMarker", obj);
            extend(obj, { fFunctions: create("TList"), fHistogram: null,
                                 fMaxSize: 0, fMaximum: -1111, fMinimum: -1111, fNpoints: 0, fX: [], fY: [] });
            break;
         case 'TGraphAsymmErrors':
            create("TGraph", obj);
            extend(obj, { fEXlow: [], fEXhigh: [], fEYlow: [], fEYhigh: []});
            break;
         case 'TMultiGraph':
            create("TNamed", obj);
            extend(obj, { fFunctions: create("TList"), fGraphs: create("TList"),
                                 fHistogram: null, fMaximum: -1111, fMinimum: -1111 });
            break;
         case 'TGraphPolargram':
            create("TNamed", obj);
            create("TAttText", obj);
            create("TAttLine", obj);
            extend(obj, { fRadian: true, fDegree: false, fGrad: false, fPolarLabelColor: 1, fRadialLabelColor: 1,
                                 fAxisAngle: 0, fPolarOffset: 0.04, fPolarTextSize: 0.04, fRadialOffset: 0.025, fRadialTextSize: 0.035,
                                 fRwrmin: 0, fRwrmax: 1, fRwtmin: 0, fRwtmax: 2*Math.PI, fTickpolarSize: 0.02,
                                 fPolarLabelFont: 62, fRadialLabelFont: 62, fCutRadial: 0, fNdivRad: 508, fNdivPol: 508 });
            break;
         case 'TPolyLine':
            create("TObject", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            extend(obj, { fLastPoint: -1, fN: 0, fOption: "", fX: null, fY: null });
            break;
         case 'TGaxis':
            create("TLine", obj);
            create("TAttText", obj);
            extend(obj, { fChopt: "", fFunctionName: "", fGridLength: 0,
                                 fLabelColor: 1, fLabelFont: 42, fLabelOffset: 0.005, fLabelSize: 0.035,
                                 fName: "", fNdiv: 12, fTickSize: 0.02, fTimeFormat: "",
                                 fTitle: "", fTitleOffset: 1, fTitleSize: 0.035,
                                 fWmax: 100, fWmin: 0 });
            break;
         case 'TAttPad':
            extend(obj, { fLeftMargin: gStyle.fPadLeftMargin,
                                 fRightMargin: gStyle.fPadRightMargin,
                                 fBottomMargin: gStyle.fPadBottomMargin,
                                 fTopMargin: gStyle.fPadTopMargin,
                                 fXfile: 2, fYfile: 2, fAfile: 1, fXstat: 0.99, fYstat: 0.99, fAstat: 2,
                                 fFrameFillColor: gStyle.fFrameFillColor,
                                 fFrameFillStyle: gStyle.fFrameFillStyle,
                                 fFrameLineColor: gStyle.fFrameLineColor,
                                 fFrameLineWidth: gStyle.fFrameLineWidth,
                                 fFrameLineStyle: gStyle.fFrameLineStyle,
                                 fFrameBorderSize: gStyle.fFrameBorderSize,
                                 fFrameBorderMode: gStyle.fFrameBorderMode });
            break;
         case 'TPad':
            create("TObject", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            create("TAttPad", obj);
            extend(obj, { fX1: 0, fY1: 0, fX2: 1, fY2: 1, fXtoAbsPixelk: 1, fXtoPixelk: 1,
                                 fXtoPixel: 1, fYtoAbsPixelk: 1, fYtoPixelk: 1, fYtoPixel: 1,
                                 fUtoAbsPixelk: 1, fUtoPixelk: 1, fUtoPixel: 1, fVtoAbsPixelk: 1,
                                 fVtoPixelk: 1, fVtoPixel: 1, fAbsPixeltoXk: 1, fPixeltoXk: 1,
                                 fPixeltoX: 1, fAbsPixeltoYk: 1, fPixeltoYk: 1, fPixeltoY: 1,
                                 fXlowNDC: 0, fYlowNDC: 0, fXUpNDC: 0, fYUpNDC: 0, fWNDC: 1, fHNDC: 1,
                                 fAbsXlowNDC: 0, fAbsYlowNDC: 0, fAbsWNDC: 1, fAbsHNDC: 1,
                                 fUxmin: 0, fUymin: 0, fUxmax: 0, fUymax: 0, fTheta: 30, fPhi: 30, fAspectRatio: 0,
                                 fNumber: 0, fLogx: gStyle.fOptLogx, fLogy: gStyle.fOptLogy, fLogz: gStyle.fOptLogz,
                                 fTickx: gStyle.fPadTickX,
                                 fTicky: gStyle.fPadTickY,
                                 fPadPaint: 0, fCrosshair: 0, fCrosshairPos: 0, fBorderSize: 2,
                                 fBorderMode: 0, fModified: false,
                                 fGridx: gStyle.fPadGridX,
                                 fGridy: gStyle.fPadGridY,
                                 fAbsCoord: false, fEditable: true, fFixedAspectRatio: false,
                                 fPrimitives: create("TList"), fExecs: null,
                                 fName: "pad", fTitle: "canvas" });

            break;
         case 'TAttCanvas':
            extend(obj, { fXBetween: 2, fYBetween: 2, fTitleFromTop: 1.2,
                                 fXdate: 0.2, fYdate: 0.3, fAdate: 1 });
            break;
         case 'TCanvas':
            create("TPad", obj);
            extend(obj, { fNumPaletteColor: 0, fNextPaletteColor: 0, fDISPLAY: "$DISPLAY",
                                 fDoubleBuffer: 0, fRetained: true, fXsizeUser: 0,
                                 fYsizeUser: 0, fXsizeReal: 20, fYsizeReal: 10,
                                 fWindowTopX: 0, fWindowTopY: 0, fWindowWidth: 0, fWindowHeight: 0,
                                 fCw: 500, fCh: 300, fCatt: create("TAttCanvas"),
                                 kMoveOpaque: true, kResizeOpaque: true, fHighLightColor: 5,
                                 fBatch: true, kShowEventStatus: false, kAutoExec: true, kMenuBar: true });
            break;
         case 'TGeoVolume':
            create("TNamed", obj);
            create("TAttLine", obj);
            create("TAttFill", obj);
            extend(obj, { fGeoAtt:0, fFinder: null, fMedium: null, fNodes: null, fNtotal: 0, fNumber: 0, fRefCount: 0, fShape: null, fVoxels: null });
            break;
         case 'TGeoNode':
            create("TNamed", obj);
            extend(obj, { fGeoAtt:0, fMother: null, fNovlp: 0, fNumber: 0, fOverlaps: null, fVolume: null });
            break;
         case 'TGeoNodeMatrix':
            create("TGeoNode", obj);
            extend(obj, { fMatrix: null });
            break;
         case 'TGeoTrack':
            create("TObject", obj);
            create("TAttLine", obj);
            create("TAttMarker", obj);
            extend(obj, { fGeoAtt: 0, fNpoints: 0, fPoints: [] });
            break;
      }

      obj._typename = typename;
      addMethods(obj, typename);
      return obj;
   }

   /** @summary Create histogram object of specified type
     * @param {string} typename - histogram typename like 'TH1I' or 'TH2F'
     * @param {number} nbinsx - number of bins on X-axis
     * @param {number} [nbinsy] - number of bins on Y-axis (for 2D/3D histograms)
     * @param {number} [nbinsz] - number of bins on Z-axis (for 3D histograms)
     * @returns {Object} created histogram object
     * @example
     * let h1 = JSROOT.createHistogram("TH1I", 20);
     * h1.fName = "Hist1";
     * h1.fTitle = "Histogram title";
     * h1.fXaxis.fTitle = "xaxis";
     * h1.fYaxis.fTitle = "yaxis";
     * h1.fXaxis.fLabelSize = 0.02; */
   JSROOT.createHistogram = (typename, nbinsx, nbinsy, nbinsz) => {
      let histo = create(typename);
      if (!histo.fXaxis || !histo.fYaxis || !histo.fZaxis) return null;
      histo.fName = "hist"; histo.fTitle = "title";
      if (nbinsx) extend(histo.fXaxis, { fNbins: nbinsx, fXmin: 0, fXmax: nbinsx });
      if (nbinsy) extend(histo.fYaxis, { fNbins: nbinsy, fXmin: 0, fXmax: nbinsy });
      if (nbinsz) extend(histo.fZaxis, { fNbins: nbinsz, fXmin: 0, fXmax: nbinsz });
      switch (parseInt(typename[2])) {
         case 1: if (nbinsx) histo.fNcells = nbinsx+2; break;
         case 2: if (nbinsx && nbinsy) histo.fNcells = (nbinsx+2) * (nbinsy+2); break;
         case 3: if (nbinsx && nbinsy && nbinsz) histo.fNcells = (nbinsx+2) * (nbinsy+2) * (nbinsz+2); break;
      }
      if (histo.fNcells > 0) {
         switch (typename[3]) {
            case "C": histo.fArray = new Int8Array(histo.fNcells); break;
            case "S": histo.fArray = new Int16Array(histo.fNcells); break;
            case "I": histo.fArray = new Int32Array(histo.fNcells); break;
            case "F": histo.fArray = new Float32Array(histo.fNcells); break;
            case "L": histo.fArray = new Float64Array(histo.fNcells); break;
            case "D": histo.fArray = new Float64Array(histo.fNcells); break;
            default: histo.fArray = new Array(histo.fNcells); break;
         }
         for (let i=0;i<histo.fNcells;++i) histo.fArray[i] = 0;
      }
      return histo;
   }

   /** @summary Creates TPolyLine object
     * @param {number} npoints - number of points
     * @param {boolean} [use_int32] - use Int32Array type for points, default is Float32Array */
   JSROOT.createTPolyLine = (npoints, use_int32) => {
      let poly = create("TPolyLine");
      if (npoints) {
         poly.fN = npoints;
         if (use_int32) {
            poly.fX = new Int32Array(npoints);
            poly.fY = new Int32Array(npoints);
         } else {
            poly.fX = new Float32Array(npoints);
            poly.fY = new Float32Array(npoints);
         }
      }
      return poly;
   }

   /** @summary Creates TGraph object
     * @param {number} npoints - number of points in TGraph
     * @param {array} [xpts] - array with X coordinates
     * @param {array} [ypts] - array with Y coordinates */
   JSROOT.createTGraph = (npoints, xpts, ypts) => {
      let graph = extend(create("TGraph"), { fBits: 0x408, fName: "graph", fTitle: "title" });

      if (npoints>0) {
         graph.fMaxSize = graph.fNpoints = npoints;

         const usex = (typeof xpts == 'object') && (xpts.length === npoints);
         const usey = (typeof ypts == 'object') && (ypts.length === npoints);

         for (let i = 0; i < npoints; ++i) {
            graph.fX.push(usex ? xpts[i] : i/npoints);
            graph.fY.push(usey ? ypts[i] : i/npoints);
         }
      }

      return graph;
   }

   /** @summary Creates THStack object
     * @desc As arguments one could specify any number of histograms objects
     * @example
     * let nbinsx = 20;
     * let h1 = JSROOT.createHistogram("TH1F", nbinsx);
     * let h2 = JSROOT.createHistogram("TH1F", nbinsx);
     * let h3 = JSROOT.createHistogram("TH1F", nbinsx);
     * let stack = JSROOT.createTHStack(h1, h2, h3); */
   JSROOT.createTHStack = function() {
      let stack = create("THStack");
      for(let i=0; i<arguments.length; ++i)
         stack.fHists.Add(arguments[i], "");
      return stack;
   }

   /** @summary Creates TMultiGraph object
     * @desc As arguments one could specify any number of TGraph objects
     * @example
     * let gr1 = JSROOT.createTGraph(100);
     * let gr2 = JSROOT.createTGraph(100);
     * let gr3 = JSROOT.createTGraph(100);
     * let mgr = JSROOT.createTMultiGraph(gr1, gr2, gr3); */
   JSROOT.createTMultiGraph = function() {
      let mgraph = create("TMultiGraph");
      for(let i=0; i<arguments.length; ++i)
          mgraph.fGraphs.Add(arguments[i], "");
      return mgraph;
   }

   let methodsCache = {}; // variable used to keep methods for known classes

   /** @summary Returns methods for given typename
     * @private */
   JSROOT.getMethods = function(typename, obj) {

      let m = methodsCache[typename],
          has_methods = (m!==undefined);

      if (!has_methods) m = {};

      // Due to binary I/O such TObject methods may not be set for derived classes
      // Therefore when methods requested for given object, check also that basic methods are there
      if ((typename=="TObject") || (typename=="TNamed") || (obj && (obj.fBits!==undefined)))
         if (m.TestBit === undefined) {
            m.TestBit = function (f) { return (this.fBits & f) != 0; };
            m.InvertBit = function (f) { this.fBits = this.fBits ^ (f & 0xffffff); };
         }

      if (has_methods) return m;

      if ((typename === 'TList') || (typename === 'THashList')) {
         m.Clear = function() {
            this.arr = [];
            this.opt = [];
         }
         m.Add = function(obj,opt) {
            this.arr.push(obj);
            this.opt.push((opt && typeof opt=='string') ? opt : "");
         }
         m.AddFirst = function(obj,opt) {
            this.arr.unshift(obj);
            this.opt.unshift((opt && typeof opt=='string') ? opt : "");
         }
         m.RemoveAt = function(indx) {
            this.arr.splice(indx, 1);
            this.opt.splice(indx, 1);
         }
      }

      if ((typename === "TPaveText") || (typename === "TPaveStats")) {
         m.AddText = function(txt) {
            let line = create("TLatex");
            line.fTitle = txt;
            this.fLines.Add(line);
         }
         m.Clear = function() {
            this.fLines.Clear();
         }
      }

      if ((typename.indexOf("TF1") == 0) || (typename === "TF2")) {
         m.addFormula = function(obj) {
            if (!obj) return;
            if (this.formulas === undefined) this.formulas = [];
            this.formulas.push(obj);
         }

         m.evalPar = function(x, y) {
            if (! ('_func' in this) || (this._title !== this.fTitle)) {

              let _func = this.fTitle, isformula = false, pprefix = "[";
              if (_func === "gaus") _func = "gaus(0)";
              if (this.fFormula && typeof this.fFormula.fFormula == "string") {
                 if (this.fFormula.fFormula.indexOf("[](double*x,double*p)")==0) {
                    isformula = true; pprefix = "p[";
                    _func = this.fFormula.fFormula.substr(21);
                 } else {
                    _func = this.fFormula.fFormula;
                    pprefix = "[p";
                 }
                 if (this.fFormula.fClingParameters && this.fFormula.fParams) {
                    for (let i=0;i<this.fFormula.fParams.length;++i) {
                       let regex = new RegExp('(\\[' + this.fFormula.fParams[i].first + '\\])', 'g'),
                           parvalue = this.fFormula.fClingParameters[this.fFormula.fParams[i].second];
                       _func = _func.replace(regex, (parvalue < 0) ? "(" + parvalue + ")" : parvalue);
                    }
                 }
              }

              if ('formulas' in this)
                 for (let i=0;i<this.formulas.length;++i)
                    while (_func.indexOf(this.formulas[i].fName) >= 0)
                       _func = _func.replace(this.formulas[i].fName, this.formulas[i].fTitle);
              _func = _func.replace(/\b(abs)\b/g, 'TMath::Abs')
                           .replace(/TMath::Exp\(/g, 'Math.exp(')
                           .replace(/TMath::Abs\(/g, 'Math.abs(');
              if (typeof JSROOT.Math == 'object') {
                 this._math = JSROOT.Math;
                 _func = _func.replace(/TMath::Prob\(/g, 'this._math.Prob(')
                              .replace(/TMath::Gaus\(/g, 'this._math.Gaus(')
                              .replace(/TMath::BreitWigner\(/g, 'this._math.BreitWigner(')
                              .replace(/xygaus\(/g, 'this._math.gausxy(this, x, y, ')
                              .replace(/gaus\(/g, 'this._math.gaus(this, x, ')
                              .replace(/gausn\(/g, 'this._math.gausn(this, x, ')
                              .replace(/expo\(/g, 'this._math.expo(this, x, ')
                              .replace(/landau\(/g, 'this._math.landau(this, x, ')
                              .replace(/landaun\(/g, 'this._math.landaun(this, x, ')
                              .replace(/ROOT::Math::/g, 'this._math.');
              }
              for (let i=0;i<this.fNpar;++i) {
                 let parname = pprefix + i + "]";
                 while(_func.indexOf(parname) != -1)
                    _func = _func.replace(parname, '('+this.GetParValue(i)+')');
              }
              _func = _func.replace(/\b(sin)\b/gi, 'Math.sin')
                           .replace(/\b(cos)\b/gi, 'Math.cos')
                           .replace(/\b(tan)\b/gi, 'Math.tan')
                           .replace(/\b(exp)\b/gi, 'Math.exp')
                           .replace(/\b(pow)\b/gi, 'Math.pow')
                           .replace(/pi/g, 'Math.PI');
              for (let n=2;n<10;++n)
                 _func = _func.replace('x^'+n, 'Math.pow(x,'+n+')');

              if (isformula) {
                 _func = _func.replace(/x\[0\]/g,"x");
                 if (this._typename==="TF2") {
                    _func = _func.replace(/x\[1\]/g,"y");
                    this._func = new Function("x", "y", _func).bind(this);
                 } else {
                    this._func = new Function("x", _func).bind(this);
                 }
              } else
              if (this._typename==="TF2")
                 this._func = new Function("x", "y", "return " + _func).bind(this);
              else
                 this._func = new Function("x", "return " + _func).bind(this);

              this._title = this.fTitle;
            }

            return this._func(x, y);
         }
         m.GetParName = function(n) {
            if (this.fParams && this.fParams.fParNames) return this.fParams.fParNames[n];
            if (this.fFormula && this.fFormula.fParams) {
               for (let k=0;k<this.fFormula.fParams.length;++k)
                  if(this.fFormula.fParams[k].second == n)
                     return this.fFormula.fParams[k].first;
            }
            if (this.fNames && this.fNames[n]) return this.fNames[n];
            return "p"+n;
         }
         m.GetParValue = function(n) {
            if (this.fParams && this.fParams.fParameters) return this.fParams.fParameters[n];
            if (this.fFormula && this.fFormula.fClingParameters) return this.fFormula.fClingParameters[n];
            if (this.fParams) return this.fParams[n];
            return undefined;
         }
         m.GetParError = function(n) {
            return this.fParErrors ? this.fParErrors[n] : undefined;
         }
         m.GetNumPars = function() {
            return this.fNpar;
         }
      }

      if (((typename.indexOf("TGraph") == 0) || (typename == "TCutG")) && (typename != "TGraphPolargram") && (typename != "TGraphTime")) {
         // check if point inside figure specified by the TGraph
         m.IsInside = function(xp,yp) {
            let i, j = this.fNpoints - 1, x = this.fX, y = this.fY, oddNodes = false;

            for (i=0; i<this.fNpoints; ++i) {
               if ((y[i]<yp && y[j]>=yp) || (y[j]<yp && y[i]>=yp)) {
                  if (x[i]+(yp-y[i])/(y[j]-y[i])*(x[j]-x[i])<xp) {
                     oddNodes = !oddNodes;
                  }
               }
               j=i;
            }

            return oddNodes;
         }
      }

      if (typename.indexOf("TH1") == 0 ||
          typename.indexOf("TH2") == 0 ||
          typename.indexOf("TH3") == 0) {
         m.getBinError = function(bin) {
            //   -*-*-*-*-*Return value of error associated to bin number bin*-*-*-*-*
            //    if the sum of squares of weights has been defined (via Sumw2),
            //    this function returns the sqrt(sum of w2).
            //    otherwise it returns the sqrt(contents) for this bin.
            if (bin >= this.fNcells) bin = this.fNcells - 1;
            if (bin < 0) bin = 0;
            if (bin < this.fSumw2.length)
               return Math.sqrt(this.fSumw2[bin]);
            return Math.sqrt(Math.abs(this.fArray[bin]));
         }
         m.setBinContent = function(bin, content) {
            // Set bin content - only trivial case, without expansion
            this.fEntries++;
            this.fTsumw = 0;
            if ((bin>=0) && (bin<this.fArray.length))
               this.fArray[bin] = content;
         }
      }

      if (typename.indexOf("TH1") == 0) {
         m.getBin = function(x) { return x; }
         m.getBinContent = function(bin) { return this.fArray[bin]; }
         m.Fill = function(x, weight) {
            let axis = this.fXaxis,
                bin = 1 + Math.floor((x - axis.fXmin) / (axis.fXmax - axis.fXmin) * axis.fNbins);
            if (bin < 0) bin = 0; else
            if (bin > axis.fNbins + 1) bin = axis.fNbins + 1;
            this.fArray[bin] += ((weight===undefined) ? 1 : weight);
            this.fEntries++;
         }
      }

      if (typename.indexOf("TH2") == 0) {
         m.getBin = function(x, y) { return (x + (this.fXaxis.fNbins+2) * y); }
         m.getBinContent = function(x, y) { return this.fArray[this.getBin(x, y)]; }
         m.Fill = function(x, y, weight) {
            let axis1 = this.fXaxis, axis2 = this.fYaxis,
                bin1 = 1 + Math.floor((x - axis1.fXmin) / (axis1.fXmax - axis1.fXmin) * axis1.fNbins),
                bin2 = 1 + Math.floor((y - axis2.fXmin) / (axis2.fXmax - axis2.fXmin) * axis2.fNbins);
            if (bin1 < 0) bin1 = 0; else
            if (bin1 > axis1.fNbins + 1) bin1 = axis1.fNbins + 1;
            if (bin2 < 0) bin2 = 0; else
            if (bin2 > axis2.fNbins + 1) bin2 = axis2.fNbins + 1;
            this.fArray[bin1 + (axis1.fNbins+2)*bin2] += ((weight===undefined) ? 1 : weight);
            this.fEntries++;
         }
      }

      if (typename.indexOf("TH3") == 0) {
         m.getBin = function(x, y, z) { return (x + (this.fXaxis.fNbins+2) * (y + (this.fYaxis.fNbins+2) * z)); }
         m.getBinContent = function(x, y, z) { return this.fArray[this.getBin(x, y, z)]; }
         m.Fill = function(x, y, z, weight) {
            let axis1 = this.fXaxis, axis2 = this.fYaxis, axis3 = this.fZaxis,
                bin1 = 1 + Math.floor((x - axis1.fXmin) / (axis1.fXmax - axis1.fXmin) * axis1.fNbins),
                bin2 = 1 + Math.floor((y - axis2.fXmin) / (axis2.fXmax - axis2.fXmin) * axis2.fNbins),
                bin3 = 1 + Math.floor((z - axis3.fXmin) / (axis3.fXmax - axis3.fXmin) * axis3.fNbins);
            if (bin1 < 0) bin1 = 0; else
            if (bin1 > axis1.fNbins + 1) bin1 = axis1.fNbins + 1;
            if (bin2 < 0) bin2 = 0; else
            if (bin2 > axis2.fNbins + 1) bin2 = axis2.fNbins + 1;
            if (bin3 < 0) bin3 = 0; else
            if (bin3 > axis3.fNbins + 1) bin3 = axis3.fNbins + 1;
            this.fArray[bin1 + (axis1.fNbins+2)* (bin2+(axis2.fNbins+2)*bin3)] += ((weight===undefined) ? 1 : weight);
            this.fEntries++;
         }
      }

      if (typename.indexOf("TProfile") == 0) {
         if (typename.indexOf("TProfile2D") == 0) {
            m.getBin = function(x, y) { return (x + (this.fXaxis.fNbins+2) * y); }
            m.getBinContent = function(x, y) {
               let bin = this.getBin(x, y);
               if (bin < 0 || bin >= this.fNcells) return 0;
               if (this.fBinEntries[bin] < 1e-300) return 0;
               if (!this.fArray) return 0;
               return this.fArray[bin]/this.fBinEntries[bin];
            }
            m.getBinEntries = function(x, y) {
               let bin = this.getBin(x, y);
               if (bin < 0 || bin >= this.fNcells) return 0;
               return this.fBinEntries[bin];
            }
         } else {
            m.getBin = function(x) { return x; }
            m.getBinContent = function(bin) {
               if (bin < 0 || bin >= this.fNcells) return 0;
               if (this.fBinEntries[bin] < 1e-300) return 0;
               if (!this.fArray) return 0;
               return this.fArray[bin]/this.fBinEntries[bin];
            }
         }
         m.getBinEffectiveEntries = function(bin) {
            if (bin < 0 || bin >= this.fNcells) return 0;
            let sumOfWeights = this.fBinEntries[bin];
            if ( !this.fBinSumw2 || this.fBinSumw2.length != this.fNcells) {
               // this can happen  when reading an old file
               return sumOfWeights;
            }
            let sumOfWeightsSquare = this.fBinSumw2[bin];
            return (sumOfWeightsSquare > 0) ? sumOfWeights * sumOfWeights / sumOfWeightsSquare : 0;
         }
         m.getBinError = function(bin) {
            if (bin < 0 || bin >= this.fNcells) return 0;
            let cont = this.fArray[bin],               // sum of bin w *y
                sum  = this.fBinEntries[bin],          // sum of bin weights
                err2 = this.fSumw2[bin],               // sum of bin w * y^2
                neff = this.getBinEffectiveEntries(bin);  // (sum of w)^2 / (sum of w^2)
            if (sum < 1e-300) return 0;                  // for empty bins
            const EErrorType = { kERRORMEAN: 0, kERRORSPREAD: 1, kERRORSPREADI: 2, kERRORSPREADG: 3 };
            // case the values y are gaussian distributed y +/- sigma and w = 1/sigma^2
            if (this.fErrorMode === EErrorType.kERRORSPREADG)
               return 1.0/Math.sqrt(sum);
            // compute variance in y (eprim2) and standard deviation in y (eprim)
            let contsum = cont/sum, eprim = Math.sqrt(Math.abs(err2/sum - contsum*contsum));
            if (this.fErrorMode === EErrorType.kERRORSPREADI) {
               if (eprim != 0) return eprim/Math.sqrt(neff);
               // in case content y is an integer (so each my has an error +/- 1/sqrt(12)
               // when the std(y) is zero
               return 1.0/Math.sqrt(12*neff);
            }
            // if approximate compute the sums (of w, wy and wy2) using all the bins
            //  when the variance in y is zero
            // case option "S" return standard deviation in y
            if (this.fErrorMode === EErrorType.kERRORSPREAD) return eprim;
            // default case : fErrorMode = kERRORMEAN
            // return standard error on the mean of y
            return eprim/Math.sqrt(neff);
         }
      }

      if (typename == "TAxis") {
         m.GetBinLowEdge = function(bin) {
            if (this.fNbins <= 0) return 0;
            if ((this.fXbins.length > 0) && (bin > 0) && (bin <= this.fNbins)) return this.fXbins[bin-1];
            return this.fXmin + (bin-1) * (this.fXmax - this.fXmin) / this.fNbins;
         }
         m.GetBinCenter = function(bin) {
            if (this.fNbins <= 0) return 0;
            if ((this.fXbins.length > 0) && (bin > 0) && (bin < this.fNbins)) return (this.fXbins[bin-1] + this.fXbins[bin])/2;
            return this.fXmin + (bin-0.5) * (this.fXmax - this.fXmin) / this.fNbins;
         }
      }

      if (typeof JSROOT.getMoreMethods == "function")
         JSROOT.getMoreMethods(m, typename, obj);

      methodsCache[typename] = m;
      return m;
   }

   /** @summary Add methods for specified type.
     * @desc Will be automatically applied when decoding JSON string
     * @private */
   JSROOT.registerMethods = function(typename, m) {
      methodsCache[typename] = m;
   }

   /** @summary Returns true if object represents basic ROOT collections
     * @desc Checks if type is TList or TObjArray or TClonesArray or TMap or THashList
     * @param {object} lst - object to check
     * @param {string} [typename] - or just typename to check
     * @private */
   JSROOT.isRootCollection = function(lst, typename) {
      if (lst && (typeof lst === 'object')) {
         if ((lst.$kind === "TList") || (lst.$kind === "TObjArray")) return true;
         if (!typename) typename = lst._typename;
      }
      if (!typename) return false;
      return (typename === 'TList') || (typename === 'THashList') || (typename === 'TMap') ||
             (typename === 'TObjArray') || (typename === 'TClonesArray');
   }

   // Connects web window
   JSROOT.connectWebWindow = function(arg) {
      if (typeof arg == 'function') arg = { callback: arg };

      if (arg.openui5src) JSROOT.openui5src = arg.openui5src;
      if (arg.openui5libs) JSROOT.openui5libs = arg.openui5libs;
      if (arg.openui5theme) JSROOT.openui5theme = arg.openui5theme;
      if (!arg.ignoreUrl) {
         let url = JSROOT.decodeUrl();
         if (url.has('nogl')) JSROOT.settings.Render3D = JSROOT.constants.Render3D.SVG;
         if (url.has('libs')) JSROOT._.use_full_libs = true;
      }

      let prereq = "webwindow;";
      if (arg && arg.prereq) prereq += arg.prereq;

      return jsroot_require(prereq).then(() => {
         if (arg && arg.prereq_logdiv && document) {
            let elem = document.getElementById(arg.prereq_logdiv);
            if (elem) elem.innerHTML = '';
            delete arg.prereq_logdiv;
         }
         if (arg && arg.prereq) delete arg.prereq;
         return JSROOT.connectWebWindow(arg);
      });
   }

   /** @summary Initialize JSROOT
     * @desc Called when main JSRoot.core.js script is loaded.
     * @private */
   _.init = function() {

      if (!source_fullpath) return this;

      let d = JSROOT.decodeUrl(source_fullpath);

      if (d.has('nocache')) JSROOT.settings.NoCache = (new Date).getTime(); // use timestamp to overcome cache limitation
      if (d.has('wrong_http_response') || JSROOT.decodeUrl().has('wrong_http_response'))
         JSROOT.settings.HandleWrongHttpResponse = true; // server may send wrong content length by partial requests, use other method to control this
      if (d.has('nosap')) _.sap = undefined; // let ignore sap loader even with openui5 loaded

      return this;
   }

   JSROOT._ = _;
   JSROOT.browser = browser;
   JSROOT.gStyle = gStyle;
   JSROOT.create = create;
   JSROOT.extend = extend;
   JSROOT.addMethods = addMethods;

   return JSROOT;

}));
