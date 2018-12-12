sap.ui.define([
   'sap/ui/core/mvc/Controller',
   'sap/ui/model/json/JSONModel',
   "sap/ui/core/ResizeHandler"
], function (Controller, JSONModel, ResizeHandler) {

   "use strict";

   return Controller.extend("eve.GL", {

      onInit : function() {
         var id = this.getView().getId();
         console.log("eve.GL.onInit id = ", id);

         var data = this.getView().getViewData();
         // console.log("VIEW DATA", data);

         if (data.standalone && data.conn_handle) {
            this.mgr = new JSROOT.EVE.EveManager();
            this.mgr.UseConnection(data.conn_handle);
            this.standalone = data.standalone;
            this.mgr.RegisterUpdate(this, "onManagerUpdate");
            
         } else {
            this.mgr = data.mgr;
            this.elementid = data.elementid;
            this.kind = data.kind;
         }

         ResizeHandler.register(this.getView(), this.onResize.bind(this));
         this.fast_event = [];

         this._load_scripts = false;
         this._render_html = false;
         this.geo_painter = null;
         this.painter_ready = false;

         this.mgr.RegisterHighlight(this, "onElementHighlight");

         JSROOT.AssertPrerequisites("geom;user:evedir/EveElements.js", this.onLoadScripts.bind(this));

         // this.checkScences();
      },

      onLoadScripts: function() {
         this._load_scripts = true;
         // only when scripts loaded, one could create objects
         this.creator = new JSROOT.EVE.EveElements();
         this.checkScences();
      },
      
      onManagerUpdate: function() {
         // called when manager was updated, need only in standalone modes to detect own element id
         if (!this.standalone || this.elementid) return;
         
         var viewers = this.mgr.FindViewers();

         // first check number of views to create
         var found = null;
         for (var n=0;n<viewers.length;++n) {
            if (viewers[n].fName.indexOf(this.standalone) == 0) { found = viewers[n]; break; }
         }
         if (!found) return;
         
         this.elementid = found.fElementId;
         this.kind = (found.fName == "Default Viewer") ? "3D" : "2D";
         this.checkScences();
         
      },

      // function called from GuiPanelController
      onExit: function() {
         if (this.mgr) this.mgr.Unregister(this);
      },

      onElementChanged: function(id, element) {
         // console.log("!!!GL DETECT CHANGED", id);

         this.checkScences();
      },

      onAfterRendering: function() {

         this._render_html = true;

         // TODO: should be specified somehow in XML file
         this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%").parent().css("overflow", "hidden");

         this.checkScences();
      },

      checkScences: function() {
         
         if (!this._load_scripts || !this._render_html || !this.elementid) return;
         
         if (!this.register_for_change) {
            this.register_for_change = true;
            
            // only when rendering completed - register for modify events
            var element = this.mgr.GetElement(this.elementid);

            // loop over scene and add dependency
            for (var k=0;k<element.childs.length;++k) {
               var scene = element.childs[k];
               this.mgr.Register(scene.fSceneId, this, "onElementChanged");
            }
         }

         // start drawing only when all scenes has childs
         // this is configured view
         var element = this.mgr.GetElement(this.elementid), allok = true, anyextras = false;

         // loop over scene and add dependency
         for (var k=0;k<element.childs.length;++k) {
            var scene = element.childs[k];
            if (!scene) { allok = false; break; }
            var realscene = this.mgr.GetElement(scene.fSceneId);
            if (this.hasExtras(realscene)) anyextras = true;
         }

         if (allok && anyextras) this.drawGeometry();
      },

      drawGeometry: function() {

         // console.log("start geometry drawing", this.getView().getId());

 /*           var shape = {
      _typename: "TGeoBBox",
      fUniqueID: 0, fBits: 0x3000000, fName: "BOX", fTitle: "",
      fShapeId: 256, fShapeBits: 1024, fDX: 200, fDY: 300, fDZ: 400, fOrigin: [0,0,0]
      };

      var geom_obj = JSROOT.extend(JSROOT.Create("TEveGeoShapeExtract"),
                      { fTrans: null, fShape: shape, fRGBA: [0, 1, 0, 0.2], fElements: null, fRnrSelf: true });
 */
         var options = "", geom_obj = null;
         if (this.kind != "3D") options = "ortho_camera";

         if (this.geo_painter) {

            // when geo painter alreay exists - clear all our additional objects
            this.geo_painter.clearExtras();

            this.geo_painter.ResetReady();

         } else {

            // TODO: should be specified somehow in XML file
            this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%");

            this.geo_painter = JSROOT.Painter.CreateGeoPainter(this.getView().getDomRef(), null, options);
         }

         this.painter_ready = false;
         // assign callback function - when needed
         this.geo_painter.WhenReady(this.onGeomertyDrawn.bind(this));

         // now loop over all  scene and create three.js objects

         // top scene element
         var element = this.mgr.GetElement(this.elementid);

         // loop over scene and add dependency
         for (var k=0;k<element.childs.length;++k) {
            var scene = element.childs[k];
            if (!scene) continue;
            var realscene = this.mgr.GetElement(scene.fSceneId);

            if (realscene && realscene.childs)
               this.createExtras(realscene.childs);
         }

         // if geometry detected in the scenes, it will be used to display

         this.geo_painter.AssignObject(geom_obj);

         this.geo_painter.prepareObjectDraw(geom_obj); // and now start everything

         // AMT temporary here, should be set in camera instantiation time
         if (this.geo_painter._camera.type == "OrthographicCamera") {
            this.geo_painter._camera.left = -this.getView().$().width();
            this.geo_painter._camera.right = this.getView().$().width();
            this.geo_painter._camera.top = this.getView().$().height();
            this.geo_painter._camera.bottom = -this.getView().$().height();
            this.geo_painter._camera.updateProjectionMatrix();
            this.geo_painter.Render3D();
         }
         // JSROOT.draw(this.getView().getDomRef(), obj, "", this.onGeomertyDrawn.bind(this));
      },

      onGeomertyDrawn: function(painter) {
         this.painter_ready = true;
         this.geo_painter._highlight_handlers = [ this ]; // register ourself for highlight handling
         this.last_highlight = null;
      },

      /// function called by GeoPainter, when mesh is highlighted
      /// forward message to the EveManager to redistribute event to all other drawings
      HighlightMesh: function(mesh, color, geo_object) {
         if (this.last_highlight === geo_object) return;
         this.last_highlight = geo_object;
         this.mgr.ProcessHighlight(this, geo_object, geo_object ? 0 : 100);
      },
      
      /// invoked from the manager
      
      onElementHighlight: function(masterid) {
         if (!this.painter_ready || !this.geo_painter) return;

         // masterid used as identifier, no any recursions
         this.geo_painter.HighlightMesh(null, null, masterid, null, true);
      },
      
      makeGLRepresentation: function(elem) {
         var fname = elem.render_data.rnr_func;
         var obj3d = this.creator[fname](elem, elem.render_data);
         if (obj3d) {
            obj3d._typename = "THREE.Mesh";

            // SL: this is just identifier for highlight, required to show items on other places, set in creator
            // obj3d.geo_object = elem.fMasterId || elem.fElementId;
            // obj3d.geo_name = elem.fName; // used for highlight

            //AMT: reference needed in MIR callback
            obj3d.eveId = elem.fElementId;

            if (elem.render_data.matrix) {
               obj3d.matrixAutoUpdate = false;
               obj3d.matrix.fromArray( elem.render_data.matrix );
               obj3d.updateMatrixWorld(true);
            }
            return obj3d;
         }
      },
      
      createExtras: function(arr, toplevel) {
         if (!arr) return;
         for (var k=0;k<arr.length;++k) {
            var elem = arr[k];
            if (elem.render_data) {
               var fname = elem.render_data.rnr_func, obj3d = null;
               if (!this.creator[fname]) {
                  console.error("Function " +fname + " missing in creator");
               } else {
                  var obj3d =  this.makeGLRepresentation(elem);
                  if (obj3d)
                     this.geo_painter.addExtra(obj3d);
               }
            }
            this.createExtras(elem.childs);
         }
      },
      
      hasExtras: function(elem) {
         if (!elem) return false;
         if (elem.render_data) return true;
         if (elem.childs)
            for (var k=0;k<elem.childs.length;++k) 
               if (this.hasExtras(elem.childs[k])) return true;
         return false;
      },

      onResize: function(event) {
         // use timeout
         // console.log("resize painter")
         if (this.resize_tmout) clearTimeout(this.resize_tmout);
         this.resize_tmout = setTimeout(this.onResizeTimeout.bind(this), 300); // minimal latency
      },

      onResizeTimeout: function() {
         delete this.resize_tmout;

         // TODO: should be specified somehow in XML file
         this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%");
         if (this.geo_painter)
            this.geo_painter.CheckResize();
      },

      colorChanged: function(el) {
         this.replaceElement(el);
      },

      replaceElement: function(el) {
         if (!this.geo_painter) return;
         
         var mesh = this.getMesh(el.fElementId);
         
         var ex = this.geo_painter._extraObjects;
         for (var i=0; i < ex.arr.length; ++i) {
            if (ex.arr[i].geo_object == el.fElementId) {
               ex.RemoveAt(i);
            }
         }
         this.geo_painter.getExtrasContainer().remove(mesh);
         var obj3d =  this.makeGLRepresentation(el);
         if (obj3d) {
            this.geo_painter.addExtra(obj3d);
            this.geo_painter.getExtrasContainer().add(obj3d);
         }
         this.geo_painter.Render3D(-1);
      },

      getMesh : function(elementId) {
         if (this.geo_painter && this.geo_painter._extraObjects ) {
            var ex = this.geo_painter._extraObjects;
            for (var i=0; i < ex.arr.length; ++i) {
               if (ex.arr[i].eveId == elementId)
                  return ex.arr[i];
            }
         }
         return null;
      },

      visibilityChanged: function(el, msg) {
         console.log("visibilit element changed ", this, msg);
         var mesh = this.getMesh(el.fElementId);
         console.log("MESH visibility change ", mesh , el.fRnrSelf);
         mesh.visible = el.fRnrSelf;
         this.geo_painter.Render3D(-1);
      },

      visibilityChildrenChanged: function(el, msg) {
         console.log("visibilit children changed ", this.mgr, el);
         el.fRnrChildren = msg.fRnrChildren;
         if (el.childs) {
            for ( var i = 0; i < el.childs.length; ++i)
            {
               // console.log("visChildren", el.fElementId, "loop child ", el.childs[i] )
               var mesh = this.getMesh(el.childs[i].fElementId);
               el.childs[i].fRnrSelf =  msg.fRnrChildren;;
               mesh.visible = msg.fRnrChildren;
               this.visibilityChildrenChanged(el.childs[i], msg);
               this.geo_painter.Render3D(-1);
            }
         }
      },
      elementAdded: function(em) {
      },
      elementRemoved: function() {
      },
      beginChanges: function() {
      },
      endChanges: function() {
      }
   });

});
