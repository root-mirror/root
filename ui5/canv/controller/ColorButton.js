sap.ui.define([
   'sap/m/Button',
   'sap/m/ButtonRenderer',
   'sap/m/Dialog',
   'sap/ui/commons/ColorPicker'
], function (Button, ButtonRenderer, Dialog, ColorPicker) {

   "use strict";

   var ColorButton = Button.extend("rootui5.canv.controller.ColorButton", {
      metadata: {
         properties: {
            attrcolor : { type : "string", group : "Misc", defaultValue : null }
         }
      },
      renderer: ButtonRenderer.render,
      init: function() {
         // svg images are always loaded without @2
         this.addEventDelegate({
            onAfterRendering: function() { this._setColor(); }
         }, this);
      }
   });

   ColorButton.prototype._setColor = function() {
      this.$().children().css('background-color', this.getProperty("attrcolor"));
   }

   ColorButton.prototype.firePress = function(args) {
      // if (Button.prototype.firePress)
      //   Button.prototype.firePress.call(this, args);

      var that = this;

      if (!this.colorPicker)
         this.colorPicker = new ColorPicker();

      if (!this.colorDialog) {
         this.colorDialog = new Dialog({
            title: 'Select color',
            content: this.colorPicker,
            beginButton: new Button({
               text: 'Apply',
               press: function () {
                  if (that.colorPicker) {
                     var col = that.colorPicker.getColorString();
                      that.setProperty("attrcolor", col);
                  }
                  that.colorDialog.close();
               }
            }),
            endButton: new Button({
               text: 'Cancel',
               press: function () {
                  that.colorDialog.close();
               }
            })
         });
      }

      var col = this.getProperty("attrcolor");

      this.colorPicker.setColorString(col);
      this.colorDialog.open();
   }

   return ColorButton;

});
