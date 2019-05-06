sap.ui.define([
   'sap/ui/core/Component',
   'sap/ui/core/UIComponent',
   'sap/ui/core/mvc/Controller',
   'sap/ui/model/json/JSONModel',
   "sap/ui/core/ResizeHandler",
   'rootui5/eve7/lib/EveManager'
], function (Component, UIComponent, Controller, JSONModel, ResizeHandler, EveManager) {

   "use strict";

   // for debug purposes - do not create geometry painter, just three.js renderer
   var direct_threejs = false;

   var EveScene = null;

   return Controller.extend("rootui5.eve7.controller.GL", {

      onInit : function()
      {
         var id = this.getView().getId();

         var viewData = this.getView().getViewData();
         if (viewData) {
            this.createXXX(viewData);
         } else {
            var oRouter = UIComponent.getRouterFor(this);
            oRouter.getRoute("View").attachPatternMatched(this._onObjectMatched, this);
         }

         ResizeHandler.register(this.getView(), this.onResize.bind(this));
         this.fast_event = [];

         this._load_scripts = false;
         this._render_html = false;
         this.geo_painter = null;

         JSROOT.AssertPrerequisites("geom", this.onLoadScripts.bind(this));
      },

      onLoadScripts: function()
      {
         var pthis = this;

         // one only can load EveScene after geometry painter
         sap.ui.define(['rootui5/eve7/lib/EveScene', 'rootui5/eve7/lib/OutlinePass', 'rootui5/eve7/lib/FXAAShader'], function (_EveScene) {
            EveScene = _EveScene;
            pthis._load_scripts = true;
            pthis.checkViewReady();
         });
      },

      _onObjectMatched: function(oEvent) {
         var args = oEvent.getParameter("arguments");

         this.createXXX(Component.getOwnerComponentFor(this.getView()).getComponentData(), args.viewName, JSROOT.$eve7tmp);

         delete JSROOT.$eve7tmp;
      },

      createXXX: function(data, viewName, moredata) {

         if (viewName) {
            data.standalone = viewName;
            data.kind = viewName;
         }
         //var data = this.getView().getViewData();
         // console.log("VIEW DATA", data);

         if (moredata && moredata.mgr) {
            this.mgr = moredata.mgr;
            this.elementid = moredata.elementid;
            this.kind = moredata.kind;
            this.standalone = viewName;

            this.checkViewReady();

         } else if (data.standalone && data.conn_handle) {
            this.mgr = new EveManager();
            this.mgr.UseConnection(data.conn_handle);
            this.standalone = data.standalone;
            this.mgr.RegisterUpdate(this, "onManagerUpdate");
         } else {
            this.mgr = data.mgr;
            this.elementid = data.elementid;
            this.kind = data.kind;
         }

      },

      // MT-HAKA
      createThreejsRenderer: function()
      {
         if (!direct_threejs || this.renderer) return;

         this.scene      = new THREE.Scene();
         this.camera     = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, 100000 );
         this.rot_center = new THREE.Vector3(0,0,0);

         // this.controls = new THREE.OrbitControls( this.camera );
         //var controls = new THREE.FirstPersonControls( camera );

         this.renderer = new THREE.WebGLRenderer();
         this.renderer.setPixelRatio( window.devicePixelRatio );
         this.renderer.setSize( window.innerWidth, window.innerHeight );

         // this.scene.fog = new THREE.FogExp2( 0xaaaaaa, 0.05 );
         this.renderer.setClearColor( 0xffffff, 1 );

         this.dom_registered = false;
         //document.body.appendChild( this.renderer.domElement );

         //document.getElementById('EveViewer9').appendChild( this.renderer.domElement );

         // this.getView().getDomRef().appendChild( this.renderer.domElement );

         // -------

         // var sphere = new THREE.SphereGeometry( 0.1, 8, 8 );
         // var lamp = new THREE.DirectionalLight( 0xff5050, 0.5 );
         var lampR = new THREE.PointLight( 0xff5050, 0.7 );
         // lampR.add(new THREE.Mesh( sphere, new THREE.MeshBasicMaterial( { color: lampR.color } ) ));
         lampR.position.set(2, 2, -2);
         this.scene.add( lampR );

         var lampG = new THREE.PointLight( 0x50ff50, 0.7 );
         // lampG.add(new THREE.Mesh( sphere, new THREE.MeshBasicMaterial( { color: lampG.color } ) ));
         lampG.position.set(-2, 2, 2);
         this.scene.add( lampG );

         var lampB = new THREE.PointLight( 0x5050ff, 0.7 );
         // lampB.add(new THREE.Mesh( sphere, new THREE.MeshBasicMaterial( { color: lampB.color } ) ));
         lampB.position.set(2, 2, 2);
         this.scene.add( lampB );

         //var plane = new THREE.GridHelper(20, 20, 0x80d080, 0x8080d0);
         //this.scene.add(plane);
      },

      /** returns container for 3d objects */
      getThreejsContainer: function(name)
      {
         var prnt = null;

         if (!direct_threejs)
            prnt = this.geo_painter.getExtrasContainer();
         else
            prnt = this.scene;

         for (var k=0;k<prnt.children.length;++k)
            if (prnt.children[k]._eve_name === name)
               return prnt.children[k];

         var obj3d = new THREE.Object3D();
         obj3d._eve_name = name;
         prnt.add(obj3d);
         return obj3d;
      },

      // MT-HAKA
      render: function()
      {
         if (!direct_threejs) {
            if (this.geo_painter) {
               if (!this.first_time_render) {
                  this.first_time_render = true;
                  this.geo_painter.adjustCameraPosition(true);
               }
               this.geo_painter.Render3D();
            }
            return;
         }

         if ( ! this.dom_registered)
         {
            this.getView().getDomRef().appendChild( this.renderer.domElement );

            //this.controls = new THREE.OrbitControls( this.camera);
            this.controls = new THREE.OrbitControls( this.camera, this.getView().getDomRef() );

            this.controls.addEventListener( 'change', this.render.bind(this) );

            this.dom_registered = true;

            // Setup camera
            var sbbox = new THREE.Box3();
            sbbox.setFromObject( this.scene );

            //var center = boundingBox.getCenter();
            this.controls.target = this.rot_center;

            var maxV = new THREE.Vector3; maxV.subVectors(sbbox.max, this.rot_center);
            var minV = new THREE.Vector3; minV.subVectors(sbbox.min, this.rot_center);

            var posV = new THREE.Vector3; posV = maxV.multiplyScalar(2);

            this.camera.position.set( posV.x, posV.y, posV.z );
            this.camera.lookAt(this.rot_center);

            console.log("scene bbox ", sbbox, ", camera_pos ", posV, ", look_at ", this.rot_center);
         }

         // console.log(this.camera);

         //console.log(this.controls);
         //console.log(this.getView().getDomRef());
         //console.log(this.renderer.domElement);

         // requestAnimationFrame( this.render.bind(this) );

         // this.controls.update( );

         this.renderer.render( this.scene, this.camera );
      },

      onManagerUpdate: function()
      {
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
         this.checkViewReady();

      },

      // function called from GuiPanelController
      onExit: function()
      {
         if (this.mgr) this.mgr.Unregister(this);
      },

      onAfterRendering: function()
      {

         this._render_html = true;

         // TODO: should be specified somehow in XML file
         this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%").parent().css("overflow", "hidden");

         this.checkViewReady();
      },

      createScenes: function()
      {
         if (this.created_scenes !== undefined) return;
         this.created_scenes = [];

         // only when rendering completed - register for modify events
         var element = this.mgr.GetElement(this.elementid);

         // loop over scene and add dependency
         for (var k=0;k<element.childs.length;++k)
         {
            var scene = element.childs[k];

            var handler = new EveScene(this.mgr, scene, this);

            this.created_scenes.push(handler);
            this.mgr.addSceneHandler(handler);
         }
      },

      redrawScenes: function() {
         for (var k=0;k<this.created_scenes.length;++k)
            this.created_scenes[k].redrawScene();
      },


      /** checks if all initialization is performed */
      checkViewReady: function()
      {
         if (!this._load_scripts || !this._render_html || !this.elementid) return;

         if (direct_threejs) {
            this.createThreejsRenderer();
            this.createScenes();
            this.redrawScenes();
            return;
         }

         if (this.geo_painter) {
            this.redrawScenes();
            return;
         }


         var options = "outline";
         // options += " black, ";
         if (this.kind != "3D") options += ", ortho_camera";


         // TODO: should be specified somehow in XML file
         this.getView().$().css("overflow", "hidden").css("width", "100%").css("height", "100%");

         this.geo_painter = JSROOT.Painter.CreateGeoPainter(this.getView().getDomRef(), null, options);

         // function used by TGeoPainter to create OutlineShader - for the moment remove from JSROOT
         this.geo_painter.createOutline = function(w,h) {
            this._outlinePass = new THREE.OutlinePass( new THREE.Vector2( w, h ), this._scene, this._camera );
            this._outlinePass.edgeStrength = 5.5;
            this._outlinePass.edgeGlow = 0.7;
            this._outlinePass.edgeThickness = 1.5;
            this._outlinePass.usePatternTexture = false;
            this._outlinePass.downSampleRatio = 1;
            this._outlinePass.glowDownSampleRatio = 3;

            // const sh = THREE.OutlinePass.selection_enum["select"]; // doesnt stand for spherical harmonics :P
            // THREE.OutlinePass.selection_atts[sh].visibleEdgeColor.set('#dd1111');
            // THREE.OutlinePass.selection_atts[sh].hiddenEdgeColor.set('#1111dd');

            this._effectComposer.addPass( this._outlinePass );

            this._effectFXAA = new THREE.ShaderPass( THREE.FXAAShader );
            this._effectFXAA.uniforms[ 'resolution' ].value.set( 1 / w, 1 / h );
            this._effectFXAA.renderToScreen = true;
            this._effectComposer.addPass( this._effectFXAA );
         }

         // assign callback function - when needed
         this.geo_painter.WhenReady(this.onGeoPainterReady.bind(this));

         this.geo_painter.AssignObject(null);

         this.geo_painter.prepareObjectDraw(null); // and now start everything
      },

      onGeoPainterReady: function(painter) {

         // AMT temporary here, should be set in camera instantiation time
         if (this.geo_painter._camera.type == "OrthographicCamera") {
            this.geo_painter._camera.left = -this.getView().$().width();
            this.geo_painter._camera.right = this.getView().$().width();
            this.geo_painter._camera.top = this.getView().$().height();
            this.geo_painter._camera.bottom = -this.getView().$().height();
            this.geo_painter._camera.updateProjectionMatrix();
         }

         painter.eveGLcontroller = this;
         painter._controls.ProcessMouseMove = function(intersects) {
            var active_mesh = null, tooltip = null, resolve = null, names = [], geo_object, geo_index;

            // try to find mesh from intersections
            for (var k=0;k<intersects.length;++k) {
               var obj = intersects[k].object, info = null;
               if (!obj) continue;
               if (obj.geo_object) info = obj.geo_name; else
                  if (obj.stack) info = painter.GetStackFullName(obj.stack);
               if (info===null) continue;

               if (info.indexOf("<prnt>")==0)
                  info = painter.GetItemName() + info.substr(6);

               names.push(info);

               if (!active_mesh) {
                  active_mesh = obj;
                  tooltip = info;
                  geo_object = obj.geo_object;
                  if (obj.get_ctrl) {
                     geo_index = obj.get_ctrl().extractIndex(intersects[k]);
                     if ((geo_index !== undefined) && (typeof tooltip == "string")) tooltip += " indx:" + JSON.stringify(geo_index);
                  }
                  if (active_mesh.stack) resolve = painter.ResolveStack(active_mesh.stack);
               }
            }

            // painter.HighlightMesh(active_mesh, undefined, geo_object, geo_index); AMT override
            if (active_mesh && active_mesh.get_ctrl()){
               active_mesh.get_ctrl().setHighlight( 0xffaa33, geo_index);
            }
            else {
               var sl = painter.eveGLcontroller.created_scenes;
               for (var k=0; k < sl.length; ++k)
                   sl[k].clearHighlight();
            }


            if (painter.options.update_browser) {
               if (painter.options.highlight && tooltip) names = [ tooltip ];
               painter.ActivateInBrowser(names);
            }

            if (!resolve || !resolve.obj) return tooltip;

            var lines = JSROOT.GEO.provideInfo(resolve.obj);
            lines.unshift(tooltip);

            return { name: resolve.obj.fName, title: resolve.obj.fTitle || resolve.obj._typename, lines: lines };
         }

         // this.geo_painter._highlight_handlers = [ this ]; // register ourself for highlight handling
         this.last_highlight = null;

         // outlinePass passthrough
         this.outlinePass = this.geo_painter._outlinePass;

         // this is try to remap standard THREE.js OutlinePass with specialized
         /* if (this.outlinePass) {
            this.outlinePass.id2obj_map = {}; // FIXME:!!!!!!!!

            this.outlinePass.oldRender = this.outlinePass.render;

            this.outlinePass.render = function() {
               this._selectedObjects = Object.values(this.id2obj_map).flat();
               this.oldRender.apply(this, arguments);
            }
         }
         */
         var sz = this.geo_painter.size_for_3d();
         this.geo_painter._effectComposer.setSize( sz.width, sz.height);
         this.geo_painter._effectFXAA.uniforms[ 'resolution' ].value.set( 1 / sz.width, 1 / sz.height );

         // create only when geo painter is ready
         this.createScenes();
         this.redrawScenes();
      },

      /// invoked from the manager
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

         if (this.geo_painter){
            this.geo_painter.CheckResize();
            if (this.geo_painter._effectFXAA)
               this.geo_painter._effectFXAA.uniforms[ 'resolution' ].value.set( 1 / this.geo_painter._scene_width, 1 / this.geo_painter._scene_height );
         }
      },
   });

});
