sap.ui.define([], function() {

   "use strict";

   function GlViewer(viewer_class)
   {
      this.viewer_class = viewer_class;

      // console.log(this.get_name() + " - constructor");
   };

   GlViewer.prototype = {

      init: function(controller)
      {
         // console.log(this.get_name() + ".init()");

         if (this.controller) throw new Error(this.get_name() + "already initialized.");

         this.controller = controller;
      },

      //==============================================================================

      get_name:   function() { return "EVE.GlViewer" + this.viewer_class; },
      get_view:   function() { return this.controller.getView(); },
      get_width:  function() { return this.controller.getView().$().width(); },
      get_height: function() { return this.controller.getView().$().height(); },

      //==============================================================================

      make_object:   function(name) { return null; },
      get_top_scene: function()     { return null; },

   };

   return GlViewer;

});
