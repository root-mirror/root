sap.ui.define([
   'rootui5/panel/Controller',
   'sap/ui/model/json/JSONModel',
   'sap/ui/model/Filter',
   'sap/ui/model/FilterOperator'
], function (GuiPanelController, JSONModel, Filter, FilterOperator) {

   "use strict";

   return GuiPanelController.extend("rootui5.fitpanel.controller.FitPanel", {

      //function called from GuiPanelController
      onPanelInit : function() {

         // WORKAROUND, need to be FIXED IN THE FUTURE
         if (window && window.location && window.location.hostname && window.location.hostname.indexOf("github.io")>=0)
            JSROOT.loadScript('../rootui5/fitpanel/style/style.css');
         else
            JSROOT.loadScript('rootui5sys/fitpanel/style/style.css');

         // need dummy ranges to initialize RangeSliders,
         var data = {
               fDim: 1,
               fMinRangeX: -1,
               fMaxRangeX: 1,
               fStepX: 0.1,
               fRangeX: [-1,1],
               fMinRangeY: -1,
               fMaxRangeY: 1,
               fStepY: 0.1,
               fRangeY: [-1,1]

         };
         this.getView().setModel(new JSONModel(data));
      },

      // returns actual model object of class RFitPanelModel
      data: function() {
        return this.getView().getModel().getData();
      },

      // cause refresh of complete fit panel
      refresh: function() {
         this.doing_refresh = true;
         this.getView().getModel().refresh();
         this.doing_refresh = false;
      },

      sendModel: function(prefix) {
         if (!prefix || (typeof prefix!="string")) {
            // this is protection against infinite loop
            // may happen if by refresh of model any callbacks are activated and trying update server side
            // this should be prevented
            if (this.doing_refresh) return;
            prefix = "UPDATE:";
         }

         if (this.websocket)
            this.websocket.send(prefix + this.getView().getModel().getJSON());
      },

      // Assign the new JSONModel to data
      onWebsocketMsg: function(handle, msg) {

         if(msg.startsWith("MODEL:")) {
            var data = JSROOT.parse(msg.substr(6));

            if(data) {
               this.getView().setModel(new JSONModel(data));

               this.verifySelectedMethodMin(data);

               this.refresh();
            }
         } else if (msg.startsWith("PARS:")) {

            this.data().fFuncPars = JSROOT.parse(msg.substr(5));

            this.refresh();
         }
      },

      // Update Button
      doUpdate: function() {
         if (this.websocket)
            this.websocket.send("RELOAD");
      },

      // Fit Button
      doFit: function() {
         this.sendModel("DOFIT:");
      },

      // Draw Button
      doDraw: function() {
         this.sendModel("DODRAW:");
      },

      onPanelExit: function(){
      },

      // when selected data is changing - cause update of complete model
      onSelectedDataChange: function() {
         this.sendModel();
      },

      // when change function many elements may be changed - resync model
      onSelectedFuncChange: function() {
         this.sendModel();
      },

      // approve current fSelectMethodMin value - and change if require
      verifySelectedMethodMin: function(data) {

         this.getView().byId("MethodMin").getBinding("items").filter(new Filter("lib", FilterOperator.EQ, data.fLibrary));

         var best = 0, selected = null;

         // first find selected item
         for (var k=0;k<data.fMethodMinAll.length;++k) {
            if (data.fMethodMinAll[k].id == data.fSelectMethodMin) {
               selected = data.fMethodMinAll[k];
               break;
            }
         }

         for (var k=0;k<data.fMethodMinAll.length;++k) {
            var item = data.fMethodMinAll[k];
            if (item.lib != data.fLibrary) continue;
            if (item.id === data.fSelectMethodMin) return;
            if (!best) best = item.id;
            // if name matches with previous selection
            if (selected && (selected.text == item.text)) best = item.id;
         }

         data.fSelectMethodMin = best;
      },

      //change the combo box in Minimization Tab --- Method depending on Radio Buttons values
      selectMinimizationLibrary: function() {
         this.verifySelectedMethodMin(this.data());

         // refresh all UI elements
         this.refresh();
      },

      onContourPar1Change: function() {
         var data = this.data();
         if (data.fContourPar1Id == data.fContourPar2Id) {
            var par2 = parseInt(data.fContourPar2Id);
            if (par2 > 0) par2--; else par2 = 1;
            data.fContourPar2Id = par2.toString();
            this.refresh();
         }
      },

      onContourPar2Change: function() {
         var data = this.data();
         if (data.fContourPar1Id == data.fContourPar2Id) {
            var par1 = parseInt(data.fContourPar1Id);
            if (par1 > 0) par1--; else par1 = 1;
            data.fContourPar1Id = par1.toString();
            this.refresh();
         }
      },

      pressApplyPars: function() {
         var json = JSROOT.toJSON(this.data().fFuncPars);

         if (this.websocket)
            this.websocket.send("SETPARS:" + json);
      }

   });

   return
});
