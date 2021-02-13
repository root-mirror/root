sap.ui.define([
   'sap/ui/core/Component',
   'sap/ui/core/UIComponent',
   'sap/ui/core/mvc/Controller',
   "sap/ui/core/ResizeHandler",
   'rootui5/eve7/lib/EveManager',
   'rootui5/eve7/lib/EveScene'
], function (Component, UIComponent, Controller, ResizeHandler, EveManager, EveScene) {

   "use strict";

   return Controller.extend("rootui5.eve7.Controller.GL", {

      //==============================================================================
      // Initialization, bootstrap, destruction & cleanup.
      //==============================================================================

      onInit : function()
      {
         // var id = this.getView().getId();

         let viewData = this.getView().getViewData();
         if (viewData)
         {
            this.setupManagerAndViewType(viewData);
         }
         else
         {
            UIComponent.getRouterFor(this).getRoute("View").attachPatternMatched(this.onViewObjectMatched, this);
         }

         this._load_scripts = false;
         this._render_html  = false;
         this.htimeout = 250;

         ResizeHandler.register(this.getView(), this.onResize.bind(this));

         JSROOT.require("geom").then(() => this.onLoadScripts());
      },

      onLoadScripts: function()
      {
         this._load_scripts = true;

         this.checkViewReady();
      },

      onViewObjectMatched: function(oEvent)
      {
         let args = oEvent.getParameter("arguments");

         console.log('ON MATCHED', args.viewName);

         console.log('MORE DATA', JSROOT.$eve7tmp);

         console.log('COMPONENT DATA', Component.getOwnerComponentFor(this.getView()).getComponentData());

         this.setupManagerAndViewType(Component.getOwnerComponentFor(this.getView()).getComponentData(),
                                      args.viewName, JSROOT.$eve7tmp);

         delete JSROOT.$eve7tmp;

         this.checkViewReady();
      },

      // Initialization that can be done immediately onInit or later through UI5 bootstrap callbacks.
      setupManagerAndViewType: function(data, viewName, moredata)
      {
         if (viewName)
         {
            data.standalone = viewName;
            data.kind       = viewName;
         }

         // console.log("VIEW DATA", data);

         if (moredata && moredata.mgr)
         {
            this.mgr        = moredata.mgr;
            this.eveViewerId  = moredata.eveViewerId;
            this.kind       = moredata.kind;
            this.standalone = viewName;

            this.checkViewReady();
         }
         else if (data.standalone && data.conn_handle)
         {
            this.mgr        = new EveManager();
            this.standalone = data.standalone;
            this.mgr.UseConnection(data.conn_handle);
         }
         else
         {
            this.mgr       = data.mgr;
            this.eveViewerId = data.eveViewerId;
            this.kind      = data.kind;
         }

         this.mgr.RegisterController(this);
         this.mgr.RegisterGlController(this);
      },

      // Called when HTML parent/container rendering is complete.
      onAfterRendering: function()
      {
         this._render_html = true;

         // TODO: should be specified somehow in XML file
         this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%").parent().css("overflow", "hidden");

         this.checkViewReady();
      },

      OnEveManagerInit: function()
      {
         // called when manager was updated, need only in standalone modes to detect own element id
         if (!this.standalone) return;

         let viewers = this.mgr.FindViewers();

         // first check number of views to create
         let found = null;
         for (let n = 0; n < viewers.length; ++n)
         {
            if (viewers[n].fName.indexOf(this.standalone) == 0)
            {
               found = viewers[n];
               break;
            }
         }
         if (!found) return;

         this.eveViewerId = found.fElementId;
         this.kind      = (found.fName == "Default Viewer") ? "3D" : "2D";

         this.checkViewReady();
      },

      // Function called from GuiPanelController.
      onExit: function()
      {
         // QQQQ EveManager does not have Unregister ... nor UnregisterController
         if (this.mgr) this.mgr.Unregister(this);
         // QQQQ plus, we should unregister this as gl-controller, too
      },

      // Checks if all initialization is performed and startup renderer.
      checkViewReady: function()
      {
         if (!this.mgr || !this._load_scripts || !this._render_html || !this.eveViewerId || this.viewer_class) return;

         this.viewer_class = this.mgr.handle.getUserArgs("GLViewer");
         if ((this.viewer_class != "JSRoot") && (this.viewer_class != "Three") && (this.viewer_class != "RCore"))
            this.viewer_class = "Three";

         this.htimeout = this.mgr.handle.getUserArgs("HTimeout");
         if (this.htimeout === undefined) this.htimeout = 250;

         // when "Reset" - reset camera position
         this.dblclick_action = this.mgr.handle.getUserArgs("DblClick");

         sap.ui.require(['rootui5/eve7/lib/GlViewer' + this.viewer_class],
               function(GlViewer) {
                  this.viewer = new GlViewer(this.viewer_class);
                  this.viewer.init(this);
               }.bind(this));
      },

      //==============================================================================
      // Common functions between THREE and GeoPainter
      //==============================================================================

      /** returns container for 3d objects */
      getSceneContainer: function(scene_name)
      {
         let parent = this.viewer.get_top_scene();

         for (let k = 0; k < parent.children.length; ++k)
         {
            if (parent.children[k]._eve_name === scene_name)
            {
               return parent.children[k];
            }
         }

         let obj3d = this.viewer.make_object();
         obj3d._eve_name = scene_name;
         parent.add(obj3d);
         return obj3d;
      },

      createScenes: function()
      {
         if (this.created_scenes !== undefined) return;
         this.created_scenes = [];

         // only when rendering completed - register for modify events
         let element = this.mgr.GetElement(this.eveViewerId);

         // loop over scene and add dependency
         for (let scene of element.childs)
         {
            this.created_scenes.push(new EveScene(this.mgr, scene, this));
         }
      },

      redrawScenes: function()
      {
         for (let s of this.created_scenes)
         {
            s.redrawScene();
         }
      },

      /// invoked from ResizeHandler
      onResize: function(event)
      {
         if (this.resize_tmout) clearTimeout(this.resize_tmout);
         this.resize_tmout = setTimeout(this.onResizeTimeout.bind(this), 250); // small latency
      },

      onResizeTimeout: function()
      {
         delete this.resize_tmout;

         // console.log("onResizeTimeout", this.camera);

         // TODO: should be specified somehow in XML file
         this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%");

         if (this.viewer)
            this.viewer.onResizeTimeout();
      },

      /** Called from JSROOT context menu when object selected for browsing */
      invokeBrowseOf: function(obj_id) {
         this.mgr.SendMIR("BrowseElement(" + obj_id + ")", 0, "ROOT::Experimental::REveManager");
      }

   });

});
