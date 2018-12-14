sap.ui.define(['sap/ui/core/mvc/Controller',
               'sap/ui/model/json/JSONModel',
               'sap/ui/layout/Splitter',
               'sap/ui/layout/SplitterLayoutData',
               "sap/ui/core/ResizeHandler",
               "eve/GeomDraw"
],function(Controller, JSONModel, Splitter, SplitterLayoutData, ResizeHandler, GeomDraw) {
   "use strict";

   return Controller.extend("eve.Geom", {
      onInit: function () {
         
         console.log('GEOM CONTROLLER INIT');
         
         console.log('USE CONNECTION', this.getView().getViewData().conn_handle);
         
         var data = {
               Nodes: [
                {
                  title: "1",
                  childs: [ { title: "1.1" } , { title: "1.2" },  { title: "1.3" } ]
                },
                {
                   title: "2",
                   childs: [ { title: "2.1" } , { title: "2.2" },  { title: "2.3" } ]
                 },
                 {
                    title: "3",
                    childs: [ { title: "3.1" } , { title: "3.2" },  { title: "3.3" } ]
                 }
              ]
         };
         
         
         var model = new JSONModel(data);
         this.getView().setModel(model);
         
         // PART 2: instantiate Control and place it onto the page

         var myControl = new GeomDraw({color:"#f00"});
         //myControl.placeAt("content");

         // ok, add another instance...:
         //new my.ColorBox({color:"green"}).placeAt("content");
         
         
         this.getView().byId("mainSplitter").addContentArea(myControl);
      },

      onAfterRendering: function(){
      },

      onToolsMenuAction : function (oEvent) {

         var item = oEvent.getParameter("item");

         switch (item.getText()) {
            case "GED Editor": this.getView().byId("Summary").getController().toggleEditor(); break;
         }
      },
      
      showHelp : function(oEvent) {
         alert("User support: root-webgui@cern.ch");
      },
      
      showUserURL : function(oEvent) {
         sap.m.URLHelper.redirect("https://github.com/alja/jsroot/blob/dev/eve7.md", true);
      }
   });
});
