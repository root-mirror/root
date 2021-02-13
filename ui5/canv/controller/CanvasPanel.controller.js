sap.ui.define([
   'sap/ui/core/mvc/Controller',
   'sap/ui/core/ResizeHandler'
], function (Controller, ResizeHandler) {
   "use strict";

   return Controller.extend("rootui5.canv.controller.CanvasPanel", {

      preserveCanvasContent: function () {
         // workaround, openui5 does not preserve DOM elements when calling onBeforeRendering
         let dom = this.getView().getDomRef();
         if (this.canvas_painter && dom && dom.children.length && !this._mainChild) {
            this._mainChild = dom.children[0];
            dom.removeChild(this._mainChild);
         }
      },

      onBeforeRendering: function() {
         this.preserveCanvasContent();
         this._has_after_rendering = false;
      },

      setPainter: function(painter) {
         this.canvas_painter = painter;
      },

      getPainter: function() {
         return this.canvas_painter;
      },

      onAfterRendering: function() {
         // workaround for openui5 problem - called before actual dimension of HTML element is assigned
         // using timeout and resize event to handle it correctly
         this._has_after_rendering = true;
         this.invokeResizeTimeout(10);
      },

      hasValidSize: function() {
         return (this.getView().$().width() > 0) && (this.getView().$().height() > 0);
      },

      invokeResizeTimeout: function(tmout) {
        if (this.resize_tmout) {
            clearTimeout(this.resize_tmout);
            delete this.resize_tmout;
        }

        if (this.hasValidSize() && this._has_after_rendering)
           this.onResizeTimeout();
        else
           this.resize_tmout = setTimeout(this.onResizeTimeout.bind(this), tmout);
      },

      /** @summary Handle openui5 resize glitch
        * @desc onAfterRendering method does not provide valid dimension of the HTML element
        * One should wait either resize event or timeout and check if valid size is there
        * Only then normal rendering can be started
        * Method also used to check resize events  */
      onResizeTimeout: function() {
         delete this.resize_tmout;

         if (!this.hasValidSize())
             return this.invokeResizeTimeout(5000); // very rare check if something changed

         let check_resize = true;

         if (this._has_after_rendering) {

            this._has_after_rendering = false;
            check_resize = false;

            let dom = this.getView().getDomRef();

            if (dom && this._mainChild) {
               dom.appendChild(this._mainChild)
               delete this._mainChild;
               check_resize = true;
            }

            if (dom && !dom.children.length) {
               let d = document.createElement("div");
               d.style = "position:relative;left:0;right:0;top:0;bottom:0;height:100%;width:100%";
               dom.appendChild(d);
            }

            if (this.canvas_painter) {
               this.canvas_painter.setDom(dom.lastChild);
               this.canvas_painter.setPadName("");
            }

            if (this.canvas_painter && this.canvas_painter._window_handle) {
               this.canvas_painter.useWebsocket(this.canvas_painter._window_handle);
               delete this.canvas_painter._window_handle;
            }
         }

         if (this.canvas_painter && check_resize)
            this.canvas_painter.checkCanvasResize();
      },

      onInit: function() {
         ResizeHandler.register(this.getView(), this.invokeResizeTimeout.bind(this, 200));
      },

      onExit: function() {
         if (this.canvas_painter) {
            this.canvas_painter.cleanup();
            delete this.canvas_painter;
         }
      }
   });

});
