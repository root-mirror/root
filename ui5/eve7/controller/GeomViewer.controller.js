sap.ui.define(['sap/ui/core/Component',
               'sap/ui/core/mvc/Controller',
               'sap/ui/core/Control',
               'sap/m/Text',
               'sap/m/CheckBox',
               "sap/ui/layout/HorizontalLayout",
               "sap/ui/table/Column",
               "sap/ui/model/json/JSONModel",
               "rootui5/browser/model/BrowserModel"
],function(Component, Controller, CoreControl, mText, mCheckBox,
           HorizontalLayout, tableColumn, JSONModel, BrowserModel) {

   "use strict";

   let geo = null; // naming convention for JSROOT.GEO

   var geomColorBox = CoreControl.extend("rootui5.eve7.controller.ColorBox", { // call the new Control type "my.ColorBox" and let it inherit from sap.ui.core.Control

      // the control API:
      metadata : {
         properties : {           // setter and getter are created behind the scenes, incl. data binding and type validation
            "color" : {type: "sap.ui.core.CSSColor", defaultValue: "#fff"} // you can give a default value and more
         }
      },

      // the part creating the HTML:
      renderer : function(oRm, oControl) { // static function, so use the given "oControl" instance instead of "this" in the renderer function
         // if (!oControl.getVisible()) return;

         oRm.write("<div");
         oRm.writeControlData(oControl);  // writes the Control ID and enables event handling - important!
         oRm.addStyle("background-color", oControl.getColor());  // write the color property; UI5 has validated it to be a valid CSS color
         oRm.writeStyles();
         oRm.addClass("geomColorBox");      // add a CSS class for styles common to all control instances
         oRm.writeClasses();              // this call writes the above class plus enables support for Square.addStyleClass(...)
         oRm.write(">");
         oRm.write("</div>"); // no text content to render; close the tag
      }

/*
      // an event handler:
      onclick : function(evt) {   // is called when the Control's area is clicked - no further event registration required
         sap.ui.require([
            'sap/ui/unified/ColorPickerPopover'
         ], function (ColorPickerPopover) {
            if (!this.oColorPickerPopover) {
               this.oColorPickerPopover = new ColorPickerPopover({
                  change: this.handleChange.bind(this)
               });
            }
            this.oColorPickerPopover.setColorString(this.getColor());
            this.oColorPickerPopover.openBy(this);
         }.bind(this));
      },

      handleChange: function (oEvent) {
         var newColor = oEvent.getParameter("colorString");
         this.setColor(newColor);
         // TODO: fire a "change" event, in case the application needs to react explicitly when the color has changed
         // but when the color is bound via data binding, it will be updated also without this event
      }
*/

   });

   /** Central geometry viewer contoller
    * All TGeo functionality is loaded after main ui5 rendering is performed,
    * To start drawing, following stages should be completed:
    *    - ui5 element is rendered (onAfterRendering is called)
    *    - TGeo-related JSROOT functionality is loaded
    *    - REveGeomDrawing object delivered from the server
    * Only after all this stages are completed, one could start to analyze  */

   return Controller.extend("rootui5.geom.controller.GeomViewer", {
      onInit: function () {

         let viewData = this.getView().getViewData();

         this.websocket = viewData.conn_handle;
         this._embeded = viewData.embeded;

         // this is code for the Components.js
         // this.websocket = Component.getOwnerComponentFor(this.getView()).getComponentData().conn_handle;

         this.websocket.setReceiver(this);
         this.websocket.connect(viewData.conn_href);

         this.queue = []; // received draw messages

         // if true, most operations are performed locally without involving server
         this.standalone = this.websocket.kind == "file";

         this.cfg = { standalone: this.websocket.kind == "file" };
         this.cfg_model = new JSONModel(this.cfg);
         this.getView().setModel(this.cfg_model);

         let nobrowser = this.websocket.getUserArgs('nobrowser') || JSROOT.decodeUrl().has('nobrowser');

         if (nobrowser) {
            // remove main area - plain geometry drawing
            // if master activated - immediately show control
            let app = this.byId("geomViewerApp");
            app.setMode(sap.m.SplitAppMode.HideMode);
            app.setInitialMaster(this.createId("geomControl"));
            app.removeMasterPage(this.byId("geomHierarchy"));
            this.byId("geomControl").setShowNavButton(false);
         } else {

            // create model only for browser - no need for anybody else
            this.model = new BrowserModel();

            this.model.useIndexSuffix = false;

            let t = this.byId("treeTable");

            t.setModel(this.model);

            // let vis_selected_handler = this.visibilitySelected.bind(this);

            this.model.assignTreeTable(t);

            t.addColumn(new tableColumn({
               label: "Description",
               template: new HorizontalLayout({
                  content: [
                     //new mCheckBox({ enabled: true, visible: true, selected: "{node_visible}", select: vis_selected_handler }),
                     //new geomColorBox({color:"{color}", visible: "{color_visible}"}),
                     new mText({text:"{name}", wrapping: false })
                  ]
               })
            }));

            // catch re-rendering of the table to assign handlers
            t.addEventDelegate({
               onAfterRendering: function() { this.assignRowHandlers(); }
            }, this);
         }

         JSROOT.require("geom").then(_geo => {
            geo = _geo;
            sap.ui.define(['rootui5/eve7/lib/EveElements'], EveElements => {
               this.creator = new EveElements();
               this.creator.useIndexAsIs = JSROOT.decodeUrl().has('useindx');
               this.checkSendRequest();
            });
         });
      },

      /** @summary invoked when visibility checkbox clicked */
      visibilitySelected: function(oEvent) {
         let nodeid = this.getRowNodeId(oEvent.getSource());
         if (nodeid < 0) {
            console.error('Fail to identify nodeid');
            return;
         }

         let msg = "SETVI" + (oEvent.getParameter("selected") ? "1:" : "0:") + JSON.stringify(nodeid);

         // send info message to client to change visibility
         this.websocket.send(msg);
      },

      assignRowHandlers: function() {
         let rows = this.byId("treeTable").getRows();
         for (let k = 0; k < rows.length; ++k) {
            rows[k].$().hover(this.onRowHover.bind(this, rows[k], true), this.onRowHover.bind(this, rows[k], false));
         }
      },

      /** @summary Send REveGeomRequest data to geometry viewer */
      sendViewerRequest: function(_oper, args) {
         let req = { oper: _oper, path: [], stack: [] };
         JSROOT.extend(req, args);
         this.websocket.send("GVREQ:" + JSON.stringify(req));
      },

      /** Process reply on REveGeomRequest */
      processViewerReply: function(repl) {
         if (!repl || (typeof repl != "object") || !repl.oper)
            return false;

         if (repl.oper == "HOVER") {

            this._hover_stack = repl.stack || null;
            if (this.geo_painter)
               this.geo_painter.highlightMesh(null, 0x00ff00, null, undefined, this._hover_stack, true);

         } else if (repl.oper == "HIGHL") {

            this.highlighRowWithPath(repl.path);

         }
      },

      /** @brief Handler for mouse-hover event
       * @desc Used to highlight correspondent volume on geometry drawing */
      onRowHover: function(row, is_enter) {

         // ignore hover event when drawing not exists
         if (!this.isDrawPageActive()) return;

         // property of current entry, not used now
         let ctxt = row.getBindingContext(),
             prop = ctxt ? ctxt.getProperty(ctxt.getPath()) : null;

         if (!this.standalone) {
            let req = is_enter && prop && prop.path && prop.isLeaf ? prop.path : [ "OFF" ];
            // avoid multiple time submitting same request
            if (this.comparePaths(this._last_hover_req, req) === 1000) return;

            this._last_hover_req = req;
            return this.sendViewerRequest("HOVER", { path: req });
         }

         if (this.geo_painter && this.geo_clones) {
            let strpath = "";

            if (prop && prop.path && is_enter)
               strpath = prop.path.join("/");

            // remember current element with hover stack
            this._hover_stack = strpath ? this.geo_clones.findStackByName(strpath) : null;

            this.geo_painter.highlightMesh(null, 0x00ff00, null, undefined, this._hover_stack, true);
         }
      },

      /** @summary Return nodeid for the row */
      getRowNodeId: function(row) {
         let ctxt = row.getBindingContext(),
             ttt = ctxt ? ctxt.getProperty(ctxt.getPath()) : null;
         return ttt && (ttt.id !== undefined) ? ttt.id : -1;
      },

      /** @summary Return arrys of ids for this row  */
      getRowIds: function(row) {
         let ctxt = row.getBindingContext();
         if (!ctxt) return null;

         let path = ctxt.getPath(), lastpos = 0, ids = [];

         while (lastpos>=0) {
            lastpos = path.indexOf("/childs", lastpos+1);

            let ttt = ctxt.getProperty(path.substr(0,lastpos));

            if (!ttt || (ttt.id===undefined)) {
               // it is not an error - sometime TableTree does not have displayed items
               // console.error('Fail to extract node id for path ' + path.substr(0,lastpos) + ' full path ' + ctxt.getPath());
               return null;
            }

            ids.push(ttt.id);
         }
         return ids;
      },

      /** @summary try to produce stack out of row path */
      getRowStack: function(row) {
         let ids = this.getRowIds(row);
         return ids ? this.geo_clones.buildStackByIds(ids) : null;
      },

      /** @summary Callback from geo painter when mesh object is highlighted. Use for update of TreeTable */
      highlightMesh: function(active_mesh, color, geo_object, geo_index, geo_stack) {
         if (!this.standalone) {
            let req = geo_stack ? geo_stack : [];
            // avoid multiple time submitting same request
            if (geo.isSameStack(this._last_highlight_req, req)) return;
            this._last_highlight_req = req;
            return this.sendViewerRequest("HIGHL", { stack: req });
         }

         let hpath = "";

         if (this.geo_clones && geo_stack) {
            let info = this.geo_clones.resolveStack(geo_stack);
            if (info && info.name) hpath = info.name;
         }

         this.highlighRowWithPath(hpath.split("/"));
      },

      /** @summary compare two paths to verify that both are the same
        * @returns 1000 if both are equivalent or maximal match length */
      comparePaths: function(path1, path2) {
         if (!path1) path1 = [];
         if (!path2) path2 = [];
         let len = Math.min(path1.length, path2.length);
         for (let i = 0; i < len; i++)
            if (path1[i] != path2[i])
               return i-1;

         return path1.length == path2.length ? 1000 : len;
      },

      /** @summary Highlights row with specified path */
      highlighRowWithPath: function(path) {
         let rows = this.byId("treeTable").getRows(), best_cmp = 0, best_indx = 0;

         for (let i=0;i<rows.length;++i) {
            rows[i].$().css("background-color", "");
            if (path && (path[0] != "OFF")) {
               let ctxt = rows[i].getBindingContext(),
                   prop = ctxt ? ctxt.getProperty(ctxt.getPath()) : null;

               if (prop && prop.path) {
                  let cmp = this.comparePaths(prop.path, path);
                  if (cmp > best_cmp) { best_cmp = cmp; best_indx = i; }
               }
            }
         }

         if (best_cmp > 0)
            rows[best_indx].$().css("background-color", best_cmp == 1000 ? "yellow" : "lightgrey");
      },

      createGeoPainter: function(drawopt) {

         if (this.geo_painter) {
            this.geo_painter.clearDrawings();
         } else {
            let geomDrawing = this.byId("geomDrawing");
            this.geo_painter = JSROOT.Painter.createGeoPainter(geomDrawing.getDomRef(), null, drawopt);
            this.geo_painter.setMouseTmout(0);
            // this.geo_painter.setDepthMethod("dflt");
            this.geo_painter.ctrl.notoolbar = true;
            // this.geo_painter.showControlOptions = this.showControl.bind(this);

            this.geo_painter.setCompleteHandler(this.completeGeoDrawing.bind(this));

            this.geom_model = new JSONModel(this.geo_painter.ctrl);
            this.geo_painter.ctrl.cfg = {}; // dummy config until real config is received
            this.byId("geomControl").setModel(this.geom_model);
            geomDrawing.setGeomPainter(this.geo_painter);
            this.geo_painter.addHighlightHandler(this);
         }

         this.geo_painter.activateInBrowser = this.activateInTreeTable.bind(this);

         this.geo_painter.assignClones(this.geo_clones);
      },

      /** @summary Extract shapes from binary data using appropriate draw message
       * @desc Draw message is vector of REveGeomVisible objects, including info where shape is in raw data */
      extractRawShapes: function(draw_msg, recreate) {

         let nodes = null, old_gradpersegm = 0;

         // array for descriptors for each node
         // if array too large (>1M), use JS object while only ~1K nodes are expected to be used
         if (recreate) {
            if (draw_msg.kind !== "draw") return false;
            nodes = (draw_msg.numnodes > 1e6) ? { length: draw_msg.numnodes } : new Array(draw_msg.numnodes); // array for all nodes
         }

         for (let cnt=0; cnt < draw_msg.nodes.length; ++cnt) {
            let node = draw_msg.nodes[cnt];
            this.formatNodeElement(node);
            if (nodes)
               nodes[node.id] = node;
            else
               this.geo_clones.updateNode(node);
         }

         if (recreate) {
            this.geo_clones = new geo.ClonedNodes(null, nodes);
            this.geo_clones.name_prefix = this.geo_clones.getNodeName(0);
            // normally only need when making selection, not used in geo viewer
            // this.geo_clones.setMaxVisNodes(draw_msg.maxvisnodes);
            // this.geo_clones.setVisLevel(draw_msg.vislevel);
            // parameter need for visualization with transparency
            // TODO: provide from server
            this.geo_clones.maxdepth = 20;
         }

         let nsegm = 0;
         if (draw_msg.cfg)
            nsegm = draw_msg.cfg.nsegm;
         else if (this.geom_model)
            nsegm = this.geom_model.getProperty("/cfg/nsegm");

         if (nsegm) {
            old_gradpersegm = geo.GradPerSegm;
            geo.GradPerSegm = 360 / Math.max(nsegm,6);
         }

         for (let cnt = 0; cnt < draw_msg.visibles.length; ++cnt) {
            let item = draw_msg.visibles[cnt], rd = item.ri;

            // entry may be provided without shape - it is ok
            if (!rd) continue;

            item.server_shape = rd.server_shape =
               this.createServerShape(rd, nsegm);
         }

         if (old_gradpersegm)
            geo.GradPerSegm = old_gradpersegm;

         return true;
      },

      /** @summary Create single shape from provided raw data. If nsegm changed, shape will be recreated */
      createServerShape: function(rd, nsegm) {

         if (rd.server_shape && ((rd.nsegm===nsegm) || !rd.shape))
            return rd.server_shape;

         rd.nsegm = nsegm;

         let g = null, off = 0;

         if (rd.shape) {
            // case when TGeoShape provided as is
            g = geo.createGeometry(rd.shape);
         } else {

            if (!rd.raw || (rd.raw.length==0)) {
               console.error('No raw data at all');
               return null;
            }

            if (!rd.raw.buffer) {
               console.error('No raw buffer');
               return null;
            }

            if (rd.sz[0]) {
               rd.vtxBuff = new Float32Array(rd.raw.buffer, off, rd.sz[0]);
               off += rd.sz[0]*4;
            }

            if (rd.sz[1]) {
               rd.nrmBuff = new Float32Array(rd.raw.buffer, off, rd.sz[1]);
               off += rd.sz[1]*4;
            }

            if (rd.sz[2]) {
               rd.idxBuff = new Uint32Array(rd.raw.buffer, off, rd.sz[2]);
               off += rd.sz[2]*4;
            }

            g = this.creator.makeEveGeometry(rd);
         }

         // shape handle is similar to created in JSROOT.GeoPainter
         return {
            _typename: "$$Shape$$", // indicate that shape can be used as is
            ready: true,
            geom: g,
            nfaces: geo.numGeometryFaces(g)
         };
      },

      /** @summary function to accumulate and process all drawings messages
       * @desc if not all scripts are loaded, messages are quied and processed later */
      checkDrawMsg: function(kind, msg) {
         console.log('Get message kind ', kind);
         if (kind) {
            if (!msg)
               return console.error("No message is provided for " + kind);

            msg.kind = kind;

            this.queue.push(msg);
         }


         if (!this.creator ||            // complete JSROOT/EVE7 TGeo functionality is loaded
            !this.queue.length ||        // drawing messages are created
            !this.renderingDone) return; // UI5 rendering is performed

         // only from here we can start to analyze messages and create TGeo painter, clones objects and so on

         msg = this.queue.shift();

         console.log('Process message kind ', msg.kind);

         switch (msg.kind) {
            case "draw":

               // keep for history
               this.last_draw_msg = msg;

               // here we should decode render data
               this.extractRawShapes(msg, true);

               // after clones are existing - ensure geo painter is there
               this.createGeoPainter(msg.cfg ? msg.cfg.drawopt : "");

               // assign configuration to the control
               if (msg.cfg) {
                  this.geo_painter.ctrl.cfg = msg.cfg;
                  this.geo_painter.ctrl.show_config = true;
                  this.geom_model.refresh();
               }

               console.log('start drawing');

               this.geo_painter.prepareObjectDraw(msg.visibles, "__geom_viewer_selection__");

               // TODO: handle completion of geo drawing

               // this is just start drawing, main work will be done asynchronous
               break;

            case "found":
               // only extend nodes and decode shapes
               if (this.extractRawShapes(msg))
                  this.paintFoundNodes(msg.visibles, true);
               break;

            case "append":
               this.extractRawShapes(msg);
               this.appendNodes(msg.visibles);
               break;
         }
      },

      completeGeoDrawing: function() {
         if (this.geom_model)
            this.geom_model.refresh();
      },

      onWebsocketOpened: function(handle) {
         this.isConnected = true;

         if (this.model)
            this.model.sendFirstRequest(this.websocket);

         // when connection established, checked if we can submit request
         this.checkSendRequest();
      },

      onWebsocketClosed: function() {
         // when connection closed, close panel as well
         console.log('CLOSE WINDOW WHEN CONNECTION CLOSED');

         if (window && !this._embeded) window.close();

         this.isConnected = false;
      },

      /** Entry point for all data from server */
      onWebsocketMsg: function(handle, msg, offset) {

         // binary data can be send only as addition to draw message
         // here data can be placed in the queue and processed when all other prerequicities are done
         if (typeof msg != "string")
            return console.error("Geom viewer do not uses binary messages len = " + mgs.byteLength);

         let mhdr = msg.substr(0,6);
         msg = msg.substr(6);

         console.log(mhdr, msg.length, msg.substr(0,70), "...");

         switch (mhdr) {
         case "DESCR:":  // browser hierarchy
            this.parseDescription(msg, true);
            break;
         case "FESCR:":  // searching hierarchy
            this.parseDescription(msg, false);
            break;
         case "BREPL:":   // browser reply
            if (this.model)
               this.model.processResponse(JSON.parse(msg));
            break;
         case "FOUND:":  // text message for found query
            this.showTextInBrowser(msg);
            this.paintFoundNodes(null); // nothing can be shown
            break;
         case "MODIF:":
            this.modifyDescription(msg);
            break;
         case "GDRAW:":   // normal drawing of geometry
            this.checkDrawMsg("draw", JSROOT.parse(msg)); // use JSROOT.parse while refs are used
            break;
         case "FDRAW:":   // drawing of found nodes
            this.checkDrawMsg("found", JSROOT.parse(msg));
            break;
         case "APPND:":
            this.checkDrawMsg("append", JSROOT.parse(msg));
            break;
         case "GVRPL:":
            this.processViewerReply(JSROOT.parse(msg));
            break;
         case "NINFO:":
            this.provideNodeInfo(JSROOT.parse(msg));
            break;
         case "RELOAD":
            this.paintFoundNodes(null);
            this.doReload(true);
            break;
         case "DROPT:":
            this.applyDrawOptions(msg);
            break;
         case "IMAGE:":
            this.produceImage(msg);
            break;
         default:
            console.error('Non recognized msg ' + mhdr + ' len=' + msg.length);
         }
      },

      /** @summary Format REveGeomNode data to be able use it in list of clones */
      formatNodeElement: function(elem) {
         elem.kind = 2; // special element for geom viewer, used in TGeoPainter
         elem.vis = 2; // visibility is alwys on
         let m = elem.matr;
         delete elem.matr;
         if (!m || !m.length) return;

         if (m.length == 16) {
            elem.matrix = m;
         } else {
            let nm = elem.matrix = new Array(16);
            for (let k = 0; k < 16; ++k) nm[k] = 0;
            nm[0] = nm[5] = nm[10] = nm[15] = 1;

            if (m.length == 3) {
               // translation martix
               nm[12] = m[0]; nm[13] = m[1]; nm[14] = m[2];
            } else if (m.length == 4) {
               // scale matrix
               nm[0] = m[0]; nm[5] = m[1]; nm[10] = m[2]; nm[15] = m[3];
            } else if (m.length == 9) {
               // rotation matrix
               nm[0] = m[0]; nm[4] = m[1]; nm[8]  = m[2];
               nm[1] = m[3]; nm[5] = m[4]; nm[9]  = m[5];
               nm[2] = m[6]; nm[6] = m[7]; nm[10] = m[8];
            } else {
               console.error('wrong number of elements in the matrix ' + m.length);
            }
         }
      },

      /** @summary Parse compact geometry description
       * @desc Used only to initialize hierarchy browser with full Tree,
       * later should be done differently */
      parseDescription: function(msg, is_original) {

         if (!this.model) return;

         let descr = JSON.parse(msg), br = this.byId("treeTable");

         br.setNoData("");
         br.setShowNoData(false);

         let topnode = this.buildTreeNode(descr, [], 0, is_original ? 1 : 999);
         if (this.standalone)
            this.fullModel = topnode;

         this.model.setFullModel(topnode);
      },

      /** TO BE CHANGED !!! When single node element is modified on the server side */
      modifyDescription: function(msg) {
         let arr = JSON.parse(msg), can_refresh = true;

         if (!arr || !this.geo_clones) return;

         console.error('modifyDescription should be changed');

         return;

         for (let k=0;k<arr.length;++k) {
            let moditem = arr[k];

            this.formatNodeElement(moditem);

            let item = this.geo_clones.nodes[moditem.id];

            if (!item)
               return console.error('Fail to find item ' + moditem.id);

            item.vis = moditem.vis;
            item.matrix = moditem.matrix;

            let dnode = this.originalCache ? this.originalCache[moditem.id] : null;

            if (dnode) {
               // here we can modify only node which was changed
               dnode.title = moditem.name;
               dnode.color_visible = false;
               dnode.node_visible = moditem.vis != 0;
            } else {
               can_refresh = false;
            }

            if (!moditem.vis && this.geo_painter)
               this.geo_painter.removeDrawnNode(moditem.id);
         }

         if (can_refresh) {
            this.model.refresh();
         } else {
            // rebuild complete tree for TreeBrowser
         }

      },

      buildTreeNode: function(nodes, cache, indx, expand_lvl) {
         let tnode = cache[indx];
         if (tnode) return tnode;
         if (!expand_lvl) expand_lvl = 0;

         let node = nodes[indx];

         cache[indx] = tnode = { name: node.name, id: indx, color_visible: false, node_visible: node.vis != 0 };

         if (expand_lvl > 0) tnode.expanded = true;

         if (node.color) {
            tnode.color = "rgb(" + node.color + ")";
            tnode.color_visisble = true;
         }

         if (node.chlds && (node.chlds.length>0)) {
            tnode.childs = [];
            tnode.nchilds = node.chlds.length;
            for (let k = 0; k < tnode.nchilds; ++k)
               tnode.childs.push(this.buildTreeNode(nodes, cache, node.chlds[k], expand_lvl-1));
         } else {
            tnode.end_node = true; // TODO: no need for such flag
         }

         return tnode;
      },

      /** @summary search main drawn nodes for matches */
      findMatchesFromDraw: function(func) {
         let matches = [];

         if (this.last_draw_msg && this.last_draw_msg.visibles && this.geo_clones)
            for (let k = 0; k < this.last_draw_msg.visibles.length; ++k) {
               let item = this.last_draw_msg.visibles[k];
               let res = this.geo_clones.resolveStack(item.stack);
               if (func(res.node))
                  matches.push({ stack: item.stack, color: item.color });
            }

         return matches;
      },

      /** @summary Show special message insted of nodes hierarchy */
      showTextInBrowser: function(text) {
         let br = this.byId("treeTable");
         br.collapseAll();
         if (!text || (text === "RESET")) {
            br.setNoData("");
            br.setShowNoData(false);

            this.model.setNoData(false);
            this.model.refresh();

         } else {
            br.setNoData(text);
            br.setShowNoData(true);
            this.model.setNoData(true);
            this.model.refresh();
         }
      },

      /** @summary Show found nodes in the browser, used for offline */
      showFoundNodes: function(matches) {

         let br = this.byId("treeTable");

         let nodes = [];
         for (let k = 0; k < matches.length; ++k)
             this.appendStackToTree(nodes, matches[k].stack, matches[k].color);

         br.setNoData("");
         br.setShowNoData(false);
         this.model.setFullModel(nodes[0]);
      },

      /** @summary Here one tries to append only given stack to the tree
        * @desc used to build partial tree with visible objects
        * Used only in standalone mode */
      appendStackToTree: function(tnodes, stack, color) {
         let prnt = null, node = null, path = [];
         for (let i = -1; i < stack.length; ++i) {
            let indx = (i < 0) ? 0 : node.chlds[stack[i]];
            node = this.geo_clones.nodes[indx];
            path.push(node.name);
            let tnode = tnodes[indx];
            if (!tnode)
               tnodes[indx] = tnode = { name: node.name, path: path.slice(), id: indx, color_visible: false, node_visible: true };

            if (prnt) {
               if (!prnt.childs) prnt.childs = [];
               if (prnt.childs.indexOf(tnode) < 0)
                  prnt.childs.push(tnode);
               prnt.nchilds = prnt.childs.length;
               prnt.expanded = true;
            }
            prnt = tnode;
         }

         prnt.end_node = true;
         prnt.color = color ? "rgb(" + color + ")" : "";
         prnt.color_visible = prnt.color.length > 0;
      },

      /** @summary Paint extra node - or remove them from painting */
      paintFoundNodes: function(visibles, append_more) {
         if (!this.geo_painter) return;

         if (append_more)
            this.geo_painter.appendMoreNodes(visibles || null);

         if (visibles && visibles.length && (visibles.length < 100)) {
            let dflt = Math.max(this.geo_painter.ctrl.transparency, 0.98);
            this.geo_painter.changedGlobalTransparency(function(node) {
               if (node.stack)
                  for (let n=0;n<visibles.length;++n)
                     if (geo.isSameStack(node.stack, visibles[n].stack))
                        return 0;
               return dflt;
            });

         } else {
            this.geo_painter.changedGlobalTransparency();
         }
      },

      appendNodes: function(nodes) {
         if (this.geo_painter) this.geo_painter.prepareObjectDraw(nodes, "__geom_viewer_append__");
      },

      showMoreNodes: function(matches) {
         if (!this.geo_painter) return;
         this.geo_painter.appendMoreNodes(matches);
         if (this._hover_stack)
            this.geo_painter.highlightMesh(null, 0x00ff00, null, undefined, this._hover_stack, true);
      },

      omBeforeRendering: function() {
         this.renderingDone = false;
      },

      onAfterRendering: function() {
         this.renderingDone = true;

         this.checkSendRequest();
      },

      onAfterMasterOpen: function() {
      },

      checkSendRequest: function(force) {
         if (force) this.ask_getdraw = false;

         if (this.isConnected && this.renderingDone) {

            if (this.creator && !this.ask_getdraw) {
               this.websocket.send("GETDRAW");
               this.ask_getdraw = true;
            }
         }
      },

      /** @summary method called from geom painter when specific node need to be activated in the browser
       * @desc Due to complex indexing in TreeTable it is not trivial to select special node */
      activateInTreeTable: function(itemnames, force) {

         if (!force || !itemnames || !this.model) return;

         let index = this.model.expandNodeByPath(itemnames[0]),
             tt = this.byId("treeTable");

         if ((index > 0) && tt)
            tt.setFirstVisibleRow(Math.max(0, index - Math.round(tt.getVisibleRowCount()/2)));
      },

      /** @summary Submit node search query to server, ignore in offline case */
      submitSearchQuery: function(query, from_handler) {

         if (!from_handler) {
            // do not submit immediately, but after very short timeout
            // if user types very fast - only last selection will be shown
            if (this.search_handler) clearTimeout(this.search_handler);
            this.search_handler = setTimeout(this.submitSearchQuery.bind(this, query, true), 1000);
            return;
         }

         delete this.search_handler;

         this.websocket.send("SEARCH:" + (query || ""));
      },

      /** when new draw options send from server */
      applyDrawOptions: function(opt) {
         if (!this.geo_painter) return;

         this.geo_painter.setAxesDraw(opt.indexOf("axis") >= 0);

         this.geo_painter.setAutoRotate(opt.indexOf("rotate") >= 0);
      },

      /** @summary when new query entered in the seach field */
      onSearch : function(oEvt) {
         let query = oEvt.getSource().getValue();
         if (!query) {
            this.paintFoundNodes(null); // remove all search results
            this.doReload(false);
         } else if (!this.standalone) {
            this.submitSearchQuery(query);
         } else {
            let lst = this.findMatchesFromDraw(function(node) {
               return node.name.indexOf(query)==0;
            });

            if (lst && lst.length) {
               this.showFoundNodes(lst);
               this.paintFoundNodes(lst);
            } else {
               this.showTextInBrowser("No found nodes");
               this.paintFoundNodes(null);
            }
         }
      },

      onCellClick: function(oEvent) {

         let tt = this.byId("treeTable"),
             first = tt.getFirstVisibleRow() || 0,
             rowindx = oEvent.getParameters().rowIndex - first,
             row = (rowindx >=0) ? tt.getRows()[rowindx] : null,
             ctxt = row ? row.getBindingContext() : null,
             prop = ctxt ? ctxt.getProperty(ctxt.getPath()) : null;

         if(prop && this.isInfoPageActive())
            if (this.standalone) {
               this.processInfoOffline(prop.path, prop.id);
            } else {
               this.sendViewerRequest("INFO", { path: prop.path });
            }
      },

      /** Try to provide as much info as possible offline */
      processInfoOffline: function(path, id) {
         let model = new JSONModel({ path: path, strpath: path.join("/")  });

         this.byId("geomInfo").setModel(model);

         if (this.geo_clones && path) {
            let stack = this.geo_clones.findStackByName(path.join("/"));

            let info = stack ? this.geo_clones.resolveStack(stack) : null;

            let build_shape = null;

            // this can be moved into GeoPainter later
            if (info && (info.id !== undefined) && this.geo_painter && this.geo_painter._draw_nodes) {
               for (let k = 0; k < this.geo_painter._draw_nodes.length; ++k) {
                  let item = this.geo_painter._draw_nodes[k];
                  if ((item.nodeid == info.id) && item.server_shape) {
                     build_shape = item.server_shape;
                     break;
                  }
               }
            }

            this.drawNodeShape(build_shape, true);
         }
      },

      /** @summary This is reply on INFO request */
      provideNodeInfo: function(info) {

         info.strpath = info.path.join("/"); // only for display

         let model = new JSONModel(info);

         this.byId("geomInfo").setModel(model);

         let server_shape = null;

         if (info.ri)
            server_shape = this.createServerShape(info.ri, 0);

         this.drawNodeShape(server_shape, false);
      },

      drawNodeShape: function(server_shape, skip_cleanup) {

         let nodeDrawing = this.byId("nodeDrawing");

         nodeDrawing.setGeomPainter(null);

         this.node_painter_active = false;
         if (this.node_painter) {
            this.node_painter.clearDrawings();
            delete this.node_painter;
            delete this.node_model;
         }

         if (!server_shape) {
            this.byId("geomControl").setModel(this.geom_model);
            return;
         }

         this.node_painter = JSROOT.Painter.createGeoPainter(nodeDrawing.getDomRef(), server_shape, "");
         this.node_painter.setMouseTmout(0);
         this.node_painter.ctrl.notoolbar = true;
         this.node_painter_active = true;

         // this.node_painter.setDepthMethod("dflt");
         this.node_model = new JSONModel(this.node_painter.ctrl);

         nodeDrawing.setGeomPainter(this.node_painter, skip_cleanup);
         this.byId("geomControl").setModel(this.node_model);
         this.node_painter.prepareObjectDraw(server_shape, "");
      },

      /** Save as png image */
      pressSaveButton: function() {
         this.produceImage("");
      },

      produceImage: function(name) {
         let painter = (this.node_painter_active && this.node_painter) ? this.node_painter : this.geo_painter;
         if (!painter) return;

         let dataUrl = painter.createSnapshot(this.standalone ? "geometry.png" : "asis");
         if (!dataUrl) return;
         let separ = dataUrl.indexOf("base64,");
         if ((separ>=0) && this.websocket && !this.standalone)
            this.websocket.send("IMAGE:" + name + "::" + dataUrl.substr(separ+7));
      },

      /** @summary Reload geometry description and base drawing, normally not required */
      onRealoadPress: function () {
         this.doReload(true);
      },

      doReload: function(force) {
         if (this.standalone) {
            this.showTextInBrowser();
            this.paintFoundNodes(null);
            if (this.model)
               this.model.setFullModel(this.fullModel);
         } else {
            this.checkSendRequest(force);

            if (this.model) {
               this.model.clearFullModel();
               this.model.reloadMainModel(force);
            }
         }
      },

      isDrawPageActive: function() {
         let app = this.byId("geomViewerApp");
         let curr = app ? app.getCurrentDetailPage() : null;
         return curr ? curr.getId() == this.createId("geomDraw") : false;
      },

      isInfoPageActive: function() {
         let app = this.byId("geomViewerApp");
         let curr = app ? app.getCurrentDetailPage() : null;
         return curr ? curr.getId() == this.createId("geomInfo") : false;
      },

      onInfoPress: function() {
         let app = this.byId("geomViewerApp"), ctrlmodel;

         if (this.isInfoPageActive()) {
            app.toDetail(this.createId("geomDraw"));
            this.node_painter_active = false;
            ctrlmodel = this.geom_model;
         } else {
            app.toDetail(this.createId("geomInfo"));
            this.node_painter_active = true;
            ctrlmodel = this.node_model;
         }

         if (ctrlmodel) {
            this.byId("geomControl").setModel(ctrlmodel);
            ctrlmodel.refresh();
         }

      },

      /** Quit ROOT session */
      onQuitRootPress: function() {
         if (!this.standalone)
            this.websocket.send("QUIT_ROOT");
      },

      onPressMasterBack: function() {
         this.byId("geomViewerApp").backMaster();
      },

      onPressDetailBack: function() {
         this.byId("geomViewerApp").backDetail();
      },

      showControl: function() {
         this.byId("geomViewerApp").toMaster(this.createId("geomControl"));
      },

      sendConfig: function() {
         if (!this.standalone && this.geo_painter && this.geo_painter.ctrl.cfg) {
            let cfg = this.geo_painter.ctrl.cfg;
            cfg.build_shapes = parseInt(cfg.build_shapes);
            this.websocket.send("CFG:" + JSROOT.toJSON(cfg));
         }
      },

      /** @summary configuration handler changes,
        * @desc after short timeout send updated config to server  */
      configChanged: function() {
         if (this.config_tmout)
            clearTimeout(this.config_tmout);

         this.config_tmout = setTimeout(this.sendConfig.bind(this), 500);
      },

      processPainterChange: function(func, arg) {
         let painter = (this.node_painter_active && this.node_painter) ? this.node_painter : this.geo_painter;

         if (painter && (typeof painter[func] == 'function'))
            painter[func](arg);
      },

      lightChanged: function() {
         this.processPainterChange('changedLight');
      },

      sliderXchange: function() {
         this.processPainterChange('changedClipping', 0);
      },

      sliderYchange: function() {
         this.processPainterChange('changedClipping', 1);
      },

      sliderZchange: function() {
         this.processPainterChange('changedClipping', 2);
      },

      clipChanged: function() {
         this.processPainterChange('changedClipping', -1);
      },

      hightlightChanged: function() {
         this.processPainterChange('changedHighlight');
      },

      transparencyChange: function() {
         this.processPainterChange('changedGlobalTransparency');
      },

      wireframeChanged: function() {
         this.processPainterChange('changedWireFrame');
      },

      backgroundChanged: function(oEvent) {
         this.processPainterChange('changedBackground', oEvent.getParameter('value'));
      },

      axesChanged: function() {
         this.processPainterChange('changedAxes');
      },

      autorotateChanged: function() {
         this.processPainterChange('changedAutoRotate');
      },

      cameraReset: function() {
         this.processPainterChange('focusCamera');
      },

      depthTestChanged: function() {
         this.processPainterChange('changedDepthTest');
      },

      depthMethodChanged: function() {
         this.processPainterChange('changedDepthMethod');
      },

      sliderTransChange: function() {
         this.processPainterChange('changedTransformation');
      },

      pressTransReset: function() {
         this.processPainterChange('changedTransformation', 'reset');
      },

      pressReset: function() {
         this.processPainterChange('resetAdvanced');
         this.byId("geomControl").getModel().refresh();
      },

      ssaoChanged: function() {
         this.processPainterChange('changedSSAO');
      }
   });

});
