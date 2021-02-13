sap.ui.define([
   'sap/ui/core/mvc/Controller',
   'sap/ui/core/ResizeHandler'
], function (Controller, ResizeHandler) {
   "use strict";

   return Controller.extend("rootui5.canv.controller.Panel", {

      onBeforeRendering: function() {
         console.log("Cleanup Panel", this.getView().getId());
         if (this.object_painter) {
            this.object_painter.cleanup();
            delete this.object_painter;
         }
         this.rendering_perfromed = false;
      },

      drawObject: function(obj, options) {

         return new Promise(resolveFunc => {
            this.draw_callbacks.push(resolveFunc);

            if (!this.rendering_perfromed) {
               this.panel_data = { object: obj, opt: options };
               return;
            }

            this.object = obj;
            d3.select(this.getView().getDomRef()).style('overflow','hidden');

            JSROOT.draw(this.getView().getDomRef(), obj, options).then(painter => {
               console.log("object painting finished");
               this.object_painter = painter;
               let arr = this.draw_callbacks;
               this.draw_callbacks = [];
               arr.forEach(cb => cb(painter));
            });
         });
      },

      drawModel: function(model) {
         if (!model) return;
         if (!this.rendering_perfromed) {
            this.panel_data = model;
            return;
         }

         if (model.object) {
            this.drawObject(model.object, model.opt);
         } else if (model.jsonfilename) {
            JSROOT.httpRequest(model.jsonfilename, 'object')
                  .then(obj => this.drawObject(obj, model.opt));
         } else if (model.filename) {
            JSROOT.openFile(model.filename)
                  .then(file => file.readObject(model.itemname))
                  .then(obj => this.drawObject(obj, model.opt));
         }
      },

      /** method to access object painter
          Promise is used while painting may be not finished when painter is requested
            let panel = sap.ui.getCore().byId("YourPanelId");
            let object_painter = null;
            panel.getController().getPainter().then(painter => {
               object_painter = painter;
            });
      */
      getPainter: function() {
         return new Promise(resolveFunc => {
            if (this.object_painter)
               resolveFunc(this.object_painter);
            else
               this.draw_callbacks.push(resolveFunc);
         });
      },

      getRenderPromise: function() {
         if (this.rendering_perfromed)
            return Promise.resolve(true);

         return new Promise(resolve => {
            if (!this.funcs) this.funcs = [];
            this.funcs.push(resolve);
         });
      },

      onAfterRendering: function() {
         this.rendering_perfromed = true;

         if (this.funcs) {
            let arr = this.funcs;
            delete this.funcs;
            arr.forEach(func => func(true));
         }

         if (this.panel_data) this.drawModel(this.panel_data);
      },

      onResize: function(event) {
         // use timeout
         if (this.resize_tmout) clearTimeout(this.resize_tmout);
         this.resize_tmout = setTimeout(this.onResizeTimeout.bind(this), 300); // minimal latency
      },

      onResizeTimeout: function() {
         delete this.resize_tmout;
         if (this.object_painter)
            this.object_painter.checkResize();
      },

      onInit: function() {

         this.draw_callbacks = []; // list of callbacks

         this.rendering_perfromed = false;

         let id = this.getView().getId();

         console.log("Initialization of JSROOT Panel", id);

         let oModel = sap.ui.getCore().getModel(id);
         if (!oModel && (id.indexOf("__xmlview0--")==0)) oModel = sap.ui.getCore().getModel(id.substr(12));

         if (oModel)
            this.panel_data = oModel.getData();

         ResizeHandler.register(this.getView(), this.onResize.bind(this));
      },

      onExit: function() {
         console.log("Exit from JSROOT Panel", this.getView().getId());

         if (this.object_painter) {
            this.object_painter.cleanup();
            delete this.object_painter;
         }
      }
   });

});
