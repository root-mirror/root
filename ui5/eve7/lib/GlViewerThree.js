sap.ui.define([
   'rootui5/eve7/lib/GlViewer',
   'rootui5/eve7/lib/EveElements',
   'rootui5/eve7/lib/OrbitControlsEve',
   'rootui5/eve7/lib/OutlinePass',
   'rootui5/eve7/lib/FXAAShader'
], function(GlViewer, EveElements) {

   "use strict";

   function GlViewerThree(viewer_class)
   {
      GlViewer.call(this, viewer_class);
   }

   GlViewerThree.prototype = Object.assign(Object.create(GlViewer.prototype), {

      constructor: GlViewerThree,

      g_highlight_update: function(mgr)
      {
         let sa = THREE.OutlinePass.selection_atts;
         let gs = mgr.GetElement(mgr.global_selection_id);
         let gh = mgr.GetElement(mgr.global_highlight_id);

         if (gs && gh) {
            sa[0].visibleEdgeColor.setStyle(JSROOT.Painter.getColor(gs.fVisibleEdgeColor));
            sa[0].hiddenEdgeColor .setStyle(JSROOT.Painter.getColor(gs.fHiddenEdgeColor));
            sa[1].visibleEdgeColor.setStyle(JSROOT.Painter.getColor(gh.fVisibleEdgeColor));
            sa[1].hiddenEdgeColor .setStyle(JSROOT.Painter.getColor(gh.fHiddenEdgeColor));
         }
      },

      init: function(controller)
      {
         GlViewer.prototype.init.call(this, controller);
         //super.init(controller);

         this.creator = new EveElements(controller);
         this.creator.useIndexAsIs = JSROOT.decodeUrl().has('useindx');

         if(!GlViewerThree.g_global_init_done)
         {
            GlViewerThree.g_global_init_done = true;

            this.controller.mgr.RegisterSelectionChangeFoo(this.g_highlight_update.bind(this));
            this.g_highlight_update(this.controller.mgr);
         }

         this.createThreejsRenderer();
         this.controller.createScenes();
         this.controller.redrawScenes();
         this.setupThreejsDomAndEventHandlers();
      },

      //==============================================================================

      make_object: function(name)
      {
         return new THREE.Object3D;
      },

      get_top_scene: function()
      {
         return this.scene;
      },

      //==============================================================================
      // THREE renderer creation, DOM/event handler setup, reset
      //==============================================================================

      createThreejsRenderer: function()
      {
         var w = this.get_width();
         var h = this.get_height();

         // console.log("createThreejsRenderer", this.kind, "w=", w, "h=", h);

         this.scene = new THREE.Scene();
         // this.scene.fog = new THREE.FogExp2( 0xaaaaaa, 0.05 );

         if (this.controller.kind === "3D")
         {
            this.camera = new THREE.PerspectiveCamera(75, w / h, 1, 5000);
         }
         else
         {
            this.camera = new THREE.OrthographicCamera(-w/2, w/2, -h/2, h/2, 0, 2000);
         }
         this.scene.add(this.camera);

         this.rot_center = new THREE.Vector3(0,0,0);

         this.renderer = new THREE.WebGLRenderer();
         this.renderer.setPixelRatio( window.devicePixelRatio );
         this.renderer.setSize( w, h );

         this.renderer.setClearColor( 0xffffff, 1 );

         // -------- Raycaster, lights, composer & FXAA and Outline passes.

         this.raycaster = new THREE.Raycaster();
         this.raycaster.params.Points.threshold = 4;   // ???
         this.raycaster.params.Line.threshold = 2.5;   // new three.js r121

         // Lights are positioned in resetRenderer

         this.point_lights = new THREE.Object3D;
         this.point_lights.add( new THREE.PointLight( 0xffffff, 0.7 )); // R
         this.point_lights.add( new THREE.PointLight( 0xffffff, 0.7 )); // G
         this.point_lights.add( new THREE.PointLight( 0xffffff, 0.7 )); // B
         this.scene.add(this.point_lights);

         // var plane = new THREE.GridHelper(20, 20, 0x80d080, 0x8080d0);
         // this.scene.add(plane);

         this.composer = new THREE.EffectComposer(this.renderer);
         this.composer.addPass(new THREE.RenderPass(this.scene, this.camera));

         this.outline_pass = new THREE.OutlinePass(new THREE.Vector2(w, h), this.scene, this.camera);
         this.outline_pass.edgeStrength        = 5.5;
         this.outline_pass.edgeGlow            = 0.7;
         this.outline_pass.edgeThickness       = 1.5;
         this.outline_pass.usePatternTexture   = false;
         this.outline_pass.downSampleRatio     = 1;
         this.outline_pass.glowDownSampleRatio = 3;

         // This does not work ... seems it is not standard pass?
         // this.outline_pass.renderToScreen = true;
         // Tried hacking with this, but would apparently need to load it somehow, sigh.
         // var copyPass = new ShaderPass( CopyShader );
         // this.composer.addPass( new THREE.ShaderPass(CopyShader) );

         this.composer.addPass(this.outline_pass);

         this.fxaa_pass = new THREE.ShaderPass( THREE.FXAAShader );
         this.fxaa_pass.uniforms.resolution.value.set(0.5 / w, 0.5 / h);
         this.fxaa_pass.renderToScreen = true;

         this.composer.addPass( this.fxaa_pass );
      },

      setupThreejsDomAndEventHandlers: function()
      {
         this.get_view().getDomRef().appendChild( this.renderer.domElement );

         // Setup tooltip
         this.ttip = document.createElement('div');
         this.ttip.setAttribute('class', 'eve_tooltip');
         this.ttip_text = document.createElement('div');
         this.ttip.appendChild(this.ttip_text);
         this.get_view().getDomRef().appendChild(this.ttip);

         // Setup controls
         this.controls = new THREE.OrbitControlsEve( this.camera, this.get_view().getDomRef() );

         this.controls.addEventListener('change', this.render.bind(this));

         // Setup some event pre-handlers
         var glc = this;

         this.renderer.domElement.addEventListener('mousemove', function(event) {

            if (event.movementX == 0 && event.movementY == 0)
               return;

            glc.removeMouseupListener();

            if (event.buttons === 0) {
               glc.removeMouseMoveTimeout();
               glc.mousemove_timeout = setTimeout(glc.onMouseMoveTimeout.bind(glc, event.offsetX, event.offsetY), glc.controller.htimeout);
            } else {
               glc.clearHighlight();
            }
         });

         this.renderer.domElement.addEventListener('mouseleave', function(event) {

            glc.removeMouseMoveTimeout();
            glc.clearHighlight();
            glc.removeMouseupListener();
         });

         this.renderer.domElement.addEventListener('mousedown', function(event) {

            glc.removeMouseMoveTimeout();
            if (event.buttons != 1 && event.buttons != 2)  glc.clearHighlight();
            glc.removeMouseupListener();

            // console.log("GLC::mousedown", this, glc, event, event.offsetX, event.offsetY);

            glc.mouseup_listener = function(event2)
            {
               this.removeEventListener('mouseup', glc.mouseup_listener);

               if (event.buttons == 1) // Selection on mouseup without move
               {
                  glc.handleMouseSelect(event2);
               }
               else if (event.buttons == 2) // Context menu on delay without move
               {
                  // Was needed for "on press with timeout"
                  // glc.controls.resetMouseDown(event);

               }
            }

            this.addEventListener('mouseup', glc.mouseup_listener);
         });

         this.renderer.domElement.addEventListener('dblclick', function(event) {
            if (glc.controller.dblclick_action == "Reset")
               glc.resetThreejsRenderer();
         });

         // Key-handlers go on window ...

         window.addEventListener('keydown', function(event) {

            // console.log("GLC::keydown", event.key, event.code, event);

            let handled = true;

            if (event.key == "t")
            {
               glc.scene.traverse( function( node ) {

                  if ( node.material && node.material.linewidth )
                  {
                     if ( ! node.material.linewidth_orig) node.material.linewidth_orig = node.material.linewidth;

                     node.material.linewidth *= 1.2;
                  }
               });
            }
            else if (event.key == "e")
            {
               glc.scene.traverse( function( node ) {

                  if ( node.material && node.material.linewidth )
                  {
                     if ( ! node.material.linewidth_orig) node.material.linewidth_orig = node.material.linewidth;

                     node.material.linewidth *= 0.8;
                  }
               });
            }
            else if (event.key == "r")
            {
               glc.scene.traverse( function( node ) {

                  if ( node.material && node.material.linewidth && node.material.linewidth_orig )
                  {
                     node.material.linewidth = node.material.linewidth_orig;
                  }
               });
            }
            else
            {
               handled = false;
            }

            if (handled)
            {
               // // // event.stopPropagation();
               // event.preventDefault();
               // event.stopImmediatePropagation();

               glc.render();
            }
         });

         // This will also call render().
         this.resetThreejsRenderer();
      },

      /** Reset camera, lights based on scene bounding box. */
      resetThreejsRenderer: function()
      {
         let sbbox = new THREE.Box3();
         sbbox.setFromObject( this.scene );

         let posV = new THREE.Vector3; posV.subVectors(sbbox.max, this.rot_center);
         let negV = new THREE.Vector3; negV.subVectors(sbbox.min, this.rot_center);

         let extV = new THREE.Vector3; extV = negV; extV.negate(); extV.max(posV);
         let extR = extV.length();

         let lc = this.point_lights.children;
         lc[0].position.set( extR, extR, -extR);
         lc[1].position.set(-extR, extR,  extR);
         lc[2].position.set( extR, extR,  extR);

         let s = 1.02;
         if (this.camera.isPerspectiveCamera)
         {
            let posC = new THREE.Vector3(-s * extR, s * extR, -s * extR);

            this.camera.position.copy(posC);

            this.controls.screenSpacePanning = true;

            // console.log("resetThreejsRenderer 3D scene bbox ", sbbox, ", camera_pos ", posC, ", look_at ", this.rot_center);
         }
         else
         {
            let posC = new THREE.Vector3(0, 0, 1000);

            this.camera.position.copy(posC);
            let ex, ey;
            if (extV.x > extV.y)
            {
               ex = extV.x;
               ey = ex / this.get_width() * this.get_height();
               if (ey < extV.y)
                  s *= extV.y/ey;
            }
            else {
               ey = extV.y;
               ex = ey / this.get_height() * this.get_width();
               if (ex < extV.x)
                  s *= extV.x/ex;
            }

            ex *= s;
            ey *= s;
            this.camera.left   = -ex;
            this.camera.right  =  ex;
            this.camera.top    =  ey;
            this.camera.bottom = -ey;

            this.controls.resetOrthoPanZoom();

            this.controls.screenSpacePanning = true;
            this.controls.enableRotate = false;

            // console.log("resetThreejsRenderer 2D scene bbox ex ey", sbbox, ex, ey, ", camera_pos ", posC, ", look_at ", this.rot_center);
         }
         this.controls.target.copy( this.rot_center );

         this.composer.reset();

         this.controls.update();
      },


      //==============================================================================

      render: function()
      {
         // Render through composer:
         this.composer.render( this.scene, this.camera );

         // or directly through renderer:
         // this.renderer.render( this.scene, this.camera );
      },

      //==============================================================================

      onResizeTimeout: function()
      {
         let w = this.get_width();
         let h = this.get_height();

         if (this.camera.isPerspectiveCamera)
         {
            this.camera.aspect = w / h;
         }
         else
         {
            this.camera.left  =  this.camera.bottom / h * w;
            this.camera.right = -this.camera.left;
            this.camera.updateProjectionMatrix();
         }
         this.camera.updateProjectionMatrix();

         this.renderer.setSize(w, h);
         this.outline_pass.setSize(w, h);
         this.fxaa_pass.uniforms.resolution.value.set(0.5 / w, 0.5 / h);

         this.composer.reset();
         this.controls.update();
         this.render();
      },


      //==============================================================================
      // THREE renderer event handlers etc.
      //==============================================================================

      //------------------------------------------------------------------------------
      // Highlight & Mouse move timeout handling
      //------------------------------------------------------------------------------

      clearHighlight: function()
      {
         if (this.highlighted_scene)
         {
            this.highlighted_scene.clearHighlight(); // XXXX should go through manager
            this.highlighted_scene = 0;

            this.ttip.style.display = "none";
         }
      },

      removeMouseMoveTimeout: function()
      {
         if (this.mousemove_timeout)
         {
            clearTimeout(this.mousemove_timeout);
            delete this.mousemove_timeout;
         }
      },

      /** Get three.js intersect object at specified mouse position */
      getIntersectAt: function(x, y) {
         let w = this.get_width();
         let h = this.get_height();

         // console.log("GLC::onMouseMoveTimeout", this, event, x, y);

         let mouse = new THREE.Vector2( ((x + 0.5) / w) * 2 - 1, -((y + 0.5) / h) * 2 + 1 );

         this.raycaster.setFromCamera(mouse, this.camera);

         let intersects = this.raycaster.intersectObjects(this.scene.children, true);

         let o = null, c = null;

         for (let i = 0; i < intersects.length; ++i)
         {
            if (intersects[i].object.get_ctrl)
            {
               intersects[i].mouse = mouse;
               intersects[i].w = w;
               intersects[i].h = h;
               return intersects[i];
            }
         }
      },

      onMouseMoveTimeout: function(x, y)
      {
         delete this.mousemove_timeout;

         var intersect = this.getIntersectAt(x,y);

         if (!intersect)
            return this.clearHighlight();

         var c = intersect.object.get_ctrl();

         var mouse = intersect.mouse;

         c.elementHighlighted(c.extractIndex(intersect));

         this.highlighted_scene = c.obj3d.scene;

         if (c.obj3d && c.obj3d.eve_el)
            this.ttip_text.innerHTML = c.getTooltipText(intersect);
         else
            this.ttip_text.innerHTML = "";

         let del  = this.controller.getView().getDomRef();
         let offs = (mouse.x > 0 || mouse.y < 0) ? this.getRelativeOffsets(del) : null;

         if (mouse.x <= 0) {
            this.ttip.style.left  = (x + del.offsetLeft + 10) + "px";
            this.ttip.style.right = null;
         } else {
            this.ttip.style.right = (intersect.w - x + offs.right + 10) + "px";
            this.ttip.style.left  = null;
         }
         if (mouse.y >= 0) {
            this.ttip.style.top    = (y + del.offsetTop + 10) + "px";
            this.ttip.style.bottom = null;
         } else {
            this.ttip.style.bottom = (intersect.h - y + offs.bottom + 10) + "px";
            this.ttip.style.top = null;
         }

         // show tooltip only in static mode
         if (this.handle.kind != "file") this.ttip.style.display= "block";
      },
      remoteToolTip: function (msg)
      {
         this.ttip_text.innerHTML = msg;
         if (this.highlighted_scene) this.ttip.style.display= "block";
      },
      getRelativeOffsets: function(elem)
      {
         // Based on:
         // https://stackoverflow.com/questions/3000887/need-to-calculate-offsetright-in-javascript

         let r = { left: 0, right: 0, top:0, bottom: 0 };

         let parent = elem.offsetParent;

         while (parent && getComputedStyle(parent).position === 'relative')
         {
            r.top    += elem.offsetTop;
            r.left   += elem.offsetLeft;
            r.right  += parent.offsetWidth  - (elem.offsetLeft + elem.offsetWidth);
            r.bottom += parent.offsetHeight - (elem.offsetTop  + elem.offsetHeight);

            elem   = parent;
            parent = parent.offsetParent;
         }

         return r;
      },

      //------------------------------------------------------------------------------
      // Mouse button handlers, selection, context menu
      //------------------------------------------------------------------------------

      removeMouseupListener: function()
      {
         if (this.mouseup_listener)
         {
            this.renderer.domElement.removeEventListener('mouseup', this.mouseup_listener);
            this.mouseup_listener = 0;
         }
      },

      showContextMenu: function(event, menu)
      {
         // console.log("GLC::showContextMenu", this, menu)

         // See js/scripts/JSRootPainter.jquery.js JSROOT.Painter.createMenu(), menu.add()


         var intersect = this.getIntersectAt(event.offsetX, event.offsetY);

         menu.add("header:Context Menu");

         if (intersect) {
            if (intersect.object.eve_el)
               menu.add("Browse to " + (intersect.object.eve_el.fName || "element"), intersect.object.eve_el.fElementId, this.controller.invokeBrowseOf.bind(this.controller));
         }

         menu.add("Reset camera", this.resetThreejsRenderer);

         menu.add("separator");

         let fff = this.defaultContextMenuAction;
         menu.add("sub:Sub Test");
         menu.add("Foo",     'foo', fff);
         menu.add("Bar",     'bar', fff);
         menu.add("Baz",     'baz', fff);
         menu.add("endsub:");

         menu.show(event);
      },

      defaultContextMenuAction: function(arg)
      {
         console.log("GLC::defaultContextMenuAction", this, arg);
      },

      handleMouseSelect: function(event)
      {
         var intersect = this.getIntersectAt(event.offsetX, event.offsetY);

         if (intersect) {
            var c = intersect.object.get_ctrl();
            c.event = event;
            c.elementSelected(c.extractIndex(intersect));
            this.highlighted_scene = intersect.object.scene;
         } else {
            // XXXX HACK - handlersMIR senders should really be in the mgr

            this.controller.created_scenes[0].processElementSelected(null, [], event);
         }
      },

   });

   return GlViewerThree;

});
