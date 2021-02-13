sap.ui.define([
   'rootui5/eve7/lib/GlViewer',
   'rootui5/eve7/lib/EveElementsRCore'
], function(GlViewer, EveElements) {

   "use strict";

   function GlViewerRCore(viewer_class)
   {
      GlViewer.call(this, viewer_class);
   }

   GlViewerRCore.prototype = Object.assign(Object.create(GlViewer.prototype), {

      constructor: GlViewerRCore,

      init: function(controller)
      {
         GlViewer.prototype.init.call(this, controller);
         // super.init(controller);

         var pthis = this;
         import("../../eve7/rnr_core/RenderCore.js").then((module) => {
            console.log("GLC onInit RenderCore loaded");
           // alert("Step 1: controller says: RnrCore loaded")
            pthis.RCore = module;

            pthis.creator = new EveElements(controller);
            pthis.creator.useIndexAsIs = JSROOT.decodeUrl().has('useindx');

            pthis.createRCoreRenderer();
            pthis.controller.createScenes();
            pthis.controller.redrawScenes();
            pthis.setupRCoreDomAndEventHandlers();
         });

      },

      //==============================================================================

      RC: function()
      {
         return this.RCore;
      },

      make_object: function(name)
      {
         let RC = this.RC();

         // return new RC.Object3D();
         return new RC.Group();
      },

      get_top_scene: function()
      {
         return this.scene;
      },

      //==============================================================================

      createRCoreRenderer: function()
      {
         let RC = this.RC();

         var w = this.get_width();
         var h = this.get_height();

         this.canvas = document.createElement('canvas');
         this.get_view().getDomRef().appendChild(this.canvas);

         this.renderer = new RC.MeshRenderer(this.canvas, RC.WEBGL2);
         this.renderer.clearColor = "#FFFFFFFF";
         this.renderer.addShaderLoaderUrls("rootui5sys/eve7/rnr_core/shaders");

         this.camera = new RC.PerspectiveCamera(120, w / h, 10, 10000);
         this.camera.position = new RC.Vector3(-500, 0, 0);
         this.camera.lookAt(new RC.Vector3(0, 0, 0), new RC.Vector3(0, 1, 0));

         this.scene = new RC.Scene();
      },

      setupRCoreDomAndEventHandlers: function()
      {
         // this.get_view().getDomRef().appendChild( this.renderer.domElement );

         // This will also call render().
         this.resetRCoreRenderer();
      },

      resetRCoreRenderer: function()
      {

         this.onResizeTimeout();

         // this.render();
      },

      //==============================================================================

      render: function()
      {
         this.renderer.render( this.scene, this.camera );
      },

      //==============================================================================

      onResizeTimeout: function()
      {
         let w = this.get_width();
         let h = this.get_height();

         this.canvas.width  = w;
         this.canvas.height = h;

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
         // this.camera.updateProjectionMatrix();

         // this.renderer.setSize(w, h);
         this.renderer.updateViewport(w, h);

         //this.outline_pass.setSize(w, h);
         //this.fxaa_pass.uniforms.resolution.value.set(0.5 / w, 0.5 / h);

         //this.composer.reset();
         //this.controls.update();

         this.render();
      },

   });

   return GlViewerRCore;
});
