sap.ui.define([
   'rootui5/eve7/lib/GlViewer',
   'rootui5/eve7/lib/EveElements',
   'rootui5/eve7/lib/OutlinePass',
   'rootui5/eve7/lib/FXAAShader'
], function(GlViewer, EveElements) {

   "use strict";

   function GlViewerJSRoot(viewer_class)
   {
      GlViewer.call(this, viewer_class);
   };

   GlViewerJSRoot.prototype = Object.assign(Object.create(GlViewer.prototype), {

      constructor: GlViewerJSRoot,

      init: function(controller)
      {
         GlViewer.prototype.init.call(this, controller);
         //super.init(controller);

         this.creator = new EveElements(controller);
         this.creator.useIndexAsIs = JSROOT.decodeUrl().has('useindx');

         this.createGeoPainter();
      },

      //==============================================================================

      make_object: function(name)
      {
         return new THREE.Object3D;
      },

      get_top_scene: function()
      {
         return this.geo_painter.getExtrasContainer();
      },

      //==============================================================================

      createGeoPainter: function()
      {
         let options = "outline";
         options += ", mouse_click"; // process mouse click events
         // options += " black, ";
         if (this.controller.kind != "3D") options += ", ortho_camera";

         // TODO: should be specified somehow in XML file
         this.get_view().$().css("overflow", "hidden").css("width", "100%").css("height", "100%");

         this.geo_painter = JSROOT.Painter.createGeoPainter(this.get_view().getDomRef(), null, options);

         this.geo_painter._geom_viewer = true; // disable several JSROOT features

         // function used by TGeoPainter to create OutlineShader - for the moment remove from JSROOT
         this.geo_painter.createOutline = function(w,h)
         {
            // this here will be TGeoPainter!

            this.outline_pass = new THREE.OutlinePass( new THREE.Vector2( w, h ), this._scene, this._camera );
            this.outline_pass.edgeStrength = 5.5;
            this.outline_pass.edgeGlow = 0.7;
            this.outline_pass.edgeThickness = 1.5;
            this.outline_pass.usePatternTexture = false;
            this.outline_pass.downSampleRatio = 1;
            this.outline_pass.glowDownSampleRatio = 3;

            // const sh = THREE.OutlinePass.selection_enum["select"]; // doesnt stand for spherical harmonics :P
            // THREE.OutlinePass.selection_atts[sh].visibleEdgeColor.set('#dd1111');
            // THREE.OutlinePass.selection_atts[sh].hiddenEdgeColor.set('#1111dd');

            this._effectComposer.addPass( this.outline_pass );

            this.fxaa_pass = new THREE.ShaderPass( THREE.FXAAShader );
            this.fxaa_pass.uniforms[ 'resolution' ].value.set( 1 / w, 1 / h );
            this.fxaa_pass.renderToScreen = true;
            this._effectComposer.addPass( this.fxaa_pass );
         }

         this.geo_painter.setMouseTmout(this.controller.htimeout);

         this.geo_painter.assignObject(null);

         this.geo_painter.prepareObjectDraw(null) // and now start everything
             .then(() => this.onGeoPainterReady(this.geo_painter));
      },

      onGeoPainterReady: function(painter)
      {
         // AMT temporary here, should be set in camera instantiation time
         if (this.geo_painter._camera.type == "OrthographicCamera")
         {
            this.geo_painter._camera.left   = -this.get_width();
            this.geo_painter._camera.right  =  this.get_width();
            this.geo_painter._camera.top    =  this.get_height();
            this.geo_painter._camera.bottom = -this.get_height();
            this.geo_painter._camera.updateProjectionMatrix();
         }

         painter.eveGLcontroller = this.controller;

         /** Handler for single mouse click, provided by basic control, used in GeoPainter */
         painter._controls.ProcessSingleClick = function(intersects)
         {
            if (!intersects) return;
            var intersect = null;
            for (var k=0;k<intersects.length;++k) {
               if (intersects[k].object.get_ctrl) {
                  intersect = intersects[k];
                  break;
               }
            }
            if (intersect) {
               var c = intersect.object.get_ctrl();
               c.elementSelected(c.extractIndex(intersect));
            }
         }

         /** Handler of mouse double click - either ignore or reset camera position */
         if (this.controller.dblclick_action != "Reset")
            painter._controls.processDblClick = function() { }

         painter._controls.ProcessMouseMove = function(intersects)
         {
            var active_mesh = null, tooltip = null, resolve = null, names = [], geo_object, geo_index;

            // try to find mesh from intersections
            for (var k=0;k<intersects.length;++k)
            {
               var obj = intersects[k].object, info = null;
               if (!obj) continue;
               if (obj.geo_object) info = obj.geo_name; else
                  if (obj.stack) info = painter.getStackFullName(obj.stack);
               if (info===null) continue;

               if (info.indexOf("<prnt>")==0)
                  info = painter.getItemName() + info.substr(6);

               names.push(info);

               if (!active_mesh) {
                  active_mesh = obj;
                  tooltip = info;
                  geo_object = obj.geo_object;
                  if (obj.get_ctrl) {
                     geo_index = obj.get_ctrl().extractIndex(intersects[k]);
                     if ((geo_index !== undefined) && (typeof tooltip == "string")) tooltip += " indx:" + JSON.stringify(geo_index);
                  }
                  if (active_mesh.stack) resolve = painter.resolveStack(active_mesh.stack);
               }
            }

            // painter.highlightMesh(active_mesh, undefined, geo_object, geo_index); AMT override
            if (active_mesh && active_mesh.get_ctrl())
            {
               active_mesh.get_ctrl().elementHighlighted(geo_index);
            }
            else
            {
               var sl = painter.eveGLcontroller.created_scenes;
               for (var k=0; k < sl.length; ++k)
                  sl[k].clearHighlight();
            }

            if (painter.options.update_browser) {
               if (painter.options.highlight && tooltip) names = [ tooltip ];
               painter.activateInBrowser(names);
            }

            if (!resolve || !resolve.obj) return tooltip;

            var lines = JSROOT.GEO.provideInfo(resolve.obj);
            lines.unshift(tooltip);

            return { name: resolve.obj.fName, title: resolve.obj.fTitle || resolve.obj._typename, lines: lines };
         }

         // this.geo_painter._highlight_handlers = [ this ]; // register ourself for highlight handling
         this.last_highlight = null;

         // outline_pass passthrough
         this.outline_pass = this.geo_painter.outline_pass;

         var sz = this.geo_painter.getSizeFor3d();
         this.geo_painter._effectComposer.setSize( sz.width, sz.height);
         this.geo_painter.fxaa_pass.uniforms[ 'resolution' ].value.set( 1 / sz.width, 1 / sz.height );

         this.geo_painter._controls.contextMenu = this.jsrootOrbitContext.bind(this);

         // create only when geo painter is ready
         this.controller.createScenes();
         this.controller.redrawScenes();

         this.geo_painter.adjustCameraPosition(true);
         this.render();
      },

      /** @summary Used together with the geo painter for processing context menu */
      jsrootOrbitContext: function(evnt, intersects) {

         var browseHandler = this.controller.invokeBrowseOf.bind(this.controller);

         JSROOT.Painter.createMenu(evnt, this.geo_painter).then(menu => {
            var numitems = 0;
            if (intersects)
               for (var n=0;n<intersects.length;++n)
                  if (intersects[n].object.geo_name) numitems++;

            if (numitems === 0) {
               // default JSROOT context menu
               menu.painter.fillContextMenu(menu);
            } else {
               var many = numitems > 1;

               if (many) menu.add("header: Items");

               for (var n=0;n<intersects.length;++n) {
                  var obj = intersects[n].object;
                  if (!obj.geo_name) continue;

                  menu.add((many ? "sub:" : "header:") + obj.geo_name, obj.geo_object, browseHandler);

                  menu.add("Browse", obj.geo_object, browseHandler);

                  var wireframe = menu.painter.accessObjectWireFrame(obj);

                  if (wireframe!==undefined)
                     menu.addchk(wireframe, "Wireframe", n, function(indx) {
                        var m = intersects[indx].object.material;
                        m.wireframe = !m.wireframe;
                        this.render3D();
                     });


                  // not yet working
                  // menu.add("Focus", n, function(indx) { this.focusCamera(intersects[indx].object); });

                  if (many) menu.add("endsub:");
               }
            }

            // show menu
            menu.show();
         });
      },

      //==============================================================================
      remoteToolTip: function ()
      {
         // to be implemented
      },

      //==============================================================================

      render: function()
      {
         this.geo_painter.render3D();
      },

      //==============================================================================

      onResizeTimeout: function()
      {
         this.geo_painter.checkResize();
         if (this.geo_painter.fxaa_pass)
            this.geo_painter.fxaa_pass.uniforms[ 'resolution' ].value.set( 1 / this.geo_painter._scene_width, 1 / this.geo_painter._scene_height );
      },
   });

   return GlViewerJSRoot;
});
