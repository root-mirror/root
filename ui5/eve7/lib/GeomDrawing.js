sap.ui.define(['sap/ui/core/Control',
               "sap/ui/core/ResizeHandler"
],function(CoreControl, ResizeHandler) {

   "use strict";

   var GeomDrawing = CoreControl.extend("rootui5.eve7.lib.GeomDrawing", {

      metadata : {
         properties : {           // setter and getter are created behind the scenes, incl. data binding and type validation
            "color" : {type: "sap.ui.core.CSSColor", defaultValue: "#fff"} // you can give a default value and more
         }
      },

      onBeforeRendering: function() {
         // remove Canvas and painter from DOM
         if (this.geom_painter) {
            this.geom_painter.clear3dCanvas();
            this.geom_painter.clearTopPainter();
         }
      },

      // the part creating the HTML:
      renderer : function(oRm, oControl) { // static function, so use the given "oControl" instance instead of "this" in the renderer function
         oRm.write("<div");
         oRm.writeControlData(oControl);  // writes the Control ID and enables event handling - important!
         // oRm.addStyle("background-color", oControl.getColor());  // write the color property; UI5 has validated it to be a valid CSS color
         oRm.addStyle("width", "100%");
         oRm.addStyle("height", "100%");
         oRm.addStyle("overflow", "hidden");
         oRm.writeStyles();
         oRm.writeClasses();              // this call writes the above class plus enables support for Square.addStyleClass(...)
         oRm.write(">");
         oRm.write("</div>"); // no text content to render; close the tag
      },

      onAfterRendering: function() {
         ResizeHandler.register(this, this.onResize.bind(this));

         if (this.geom_painter) {
            // this should be moved to GeomPainter itself !!!

            // work with and without canvas
            this.geom_painter.setDom(this.getDomRef());

            var size = this.geom_painter.getSizeFor3d();

            this.geom_painter.add3dCanvas(size, this.geom_painter._renderer.domElement);

            // set as main canvas painter or just top painter
            this.geom_painter.setAsMainPainter();

            this.geom_painter.render3D();
         }
      },

      setGeomPainter: function(painter, skip_cleanup) {

         if (this.geom_painter) {
            if (this.geom_skip_cleanup) {
               // workaround, done for nodes drawing to avoid deletion of 3D objects
               this.geom_painter._clones = null;
               this.geom_painter._clones_owner = false;
            }

            this.geom_painter.cleanup();
            delete this.geom_painter;
            delete this.geom_skip_cleanup;
         }

         if (painter) {
            this.geom_painter = painter;
            this.geom_skip_cleanup = skip_cleanup;
         }
      },

      onResize: function() {
         if (this.resize_tmout) clearTimeout(this.resize_tmout);
         this.resize_tmout = setTimeout(this.onResizeTimeout.bind(this), 100); // minimal latency
      },

      onResizeTimeout: function() {
         delete this.resize_tmout;
         if (this.geom_painter)
            this.geom_painter.checkResize();
      }
   });

   return GeomDrawing;

});

