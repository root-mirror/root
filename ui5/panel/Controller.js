sap.ui.define([
   'sap/ui/core/mvc/Controller',
   'sap/ui/model/json/JSONModel'
], function (Controller, JSONModel) {
   "use strict";

   return Controller.extend("rootui5.panel.Controller", {

      // one could define in derived classes followin methods
      //   onPanelInit: function() {},
      //   onPanelReceive: function(msg, offset) {},
      //   onPanelExit: function(),


      onInit: function() {
         var data = this.getView().getViewData();

         if (data && data.handle) {
            this.websocket = data.handle;
            this.websocket.setReceiver(this); // redirect websocket handling on controller itself
            this.websocket.send("PANEL_READY"); // confirm panel creation, only then GUI can send commands
         }

         // TODO: use more specific API between Canvas and Panel
         if (data && data.masterPanel) {
            this.masterPanel = data.masterPanel;
         }

         if (typeof this.onPanelInit == "function")
            this.onPanelInit();
      },

      onWebsocketOpened: function(handle) {
         console.log('Connection established - should never happen');
      },

      onWebsocketMsg: function(handle, msg, offset) {
         if (typeof this.onPanelReceive == 'function')
            this.onPanelReceive(msg, offset);
         else
            console.log('GuiPanel Get message ' + msg);
      },

      onWebsocketClosed: function(handle) {
          console.log('GuiPanel closed');
          if (window) window.open('','_self').close(); // window.close();
          delete this.websocket; // remove reference on websocket
      },

      setPanelTitle: function(title) {
         if (!this.masterPanel && document)
            document.title = title;
      },

      panelSend: function(msg) {
         if (this.websocket)
            this.websocket.send(msg);
         else
            console.error('No connection available to send message');
      },

      onExit: function() {
         if (typeof this.onPanelExit == 'function')
            this.onPanelExit();
         if (this.websocket) {
            this.websocket.close();
            delete this.websocket;
         }
      },

      /** Method should be used to close panel
       * Depending from used window manager different functionality can be used here */
      closePanel: function() {
         if (this.masterPanel) {
            if (this.masterPanel.showLeftArea) this.masterPanel.showLeftArea("");
         } else {
            if (window) window.open('','_self').close(); // window.close();
         }
      }

   });

});
