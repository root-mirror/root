sap.ui.define(['sap/ui/core/Component',
               'sap/ui/core/UIComponent',
               'sap/ui/core/mvc/Controller',
               'sap/ui/layout/Splitter',
               'sap/ui/layout/SplitterLayoutData',
               'sap/m/library',
               'sap/m/Button',
               'sap/m/MenuItem',
               'rootui5/eve7/lib/EveManager'
], function(Component, UIComponent, Controller, Splitter, SplitterLayoutData, MobileLibrary, mButton, mMenuItem, EveManager) {

   "use strict";

   return Controller.extend("rootui5.eve7.controller.Main", {
      onInit: function () {

         console.log('MAIN CONTROLLER INIT');

         this.mgr = new EveManager();

         var conn_handle = Component.getOwnerComponentFor(this.getView()).getComponentData().conn_handle;
         this.mgr.UseConnection(conn_handle);
         // this.mgr.UseConnection(this.getView().getViewData().conn_handle);

         // method to found summary controller by ID and set manager to it
         var elem = this.byId("Summary");
         var ctrl = elem.getController();
         ctrl.SetMgr(this.mgr);

         this.mgr.RegisterUpdate(this, "onManagerUpdate");
      },

      getHandle: function () {
         return this.handle;
      },

      UpdateCommandsButtons: function(cmds) {
          if (!cmds || this.commands) return;

          var toolbar = this.byId("otb1");

          this.commands = cmds;
          for (var k=cmds.length-1;k>=0;--k) {
             var btn = new mButton({
                // text: "ButtonNew",
                icon: cmds[k].icon,
                tooltip: cmds[k].name,
                press: this.mgr.executeCommand.bind(this.mgr, cmds[k])
              });
             toolbar.insertContent(btn, 0);
          }
      },
      
      viewItemPressed: function (oEvent) {
         var item = oEvent.getSource();
         console.log('item pressed', item.getText());
         
         var name = item.getText();
         if (name.indexOf(" ")>0) name = name.substr(0, name.indexOf(" "));
         
         var oRouter = UIComponent.getRouterFor(this);
         oRouter.navTo("View", { viewName: name });
      },
      
      updateViewers: function(loading_done) {
         var viewers = this.mgr.FindViewers();

         // first check number of views to create
         var staged = [];
         for (var n=0;n<viewers.length;++n) {
            var el = viewers[n];
            if (!el.$view_created && el.fRnrSelf) staged.push(el);
         }
         if (staged.length == 0) return;
         
         console.log("FOUND viewers", viewers.length, "not yet exists", staged.length);
         
         if (staged.length > 1) {
            var vMenu = this.getView().byId("menuViewId");
            var item = new mMenuItem({text:"Browse to"});
            vMenu.addItem(item);
            for (var n=0;n<staged.length;++n) 
               item.addItem(new mMenuItem({text: staged[n].fName, press: this.viewItemPressed.bind(this) }));
         }

         var main = this, vv = null, count = 0, sv = this.getView().byId("MainAreaSplitter");

         for (var n=0;n<staged.length;++n) {
            var elem = staged[n];
            var viewid = "EveViewer" + elem.fElementId;

            // create missing view
            elem.$view_created = true;
            console.log("Creating view", viewid);
            count++;

            var oLd = undefined;
            if ((count == 1) && (staged.length > 1))
               oLd = new SplitterLayoutData({ resizable: true, size: "50%" });

            var vtype = "rootui5.eve7.view.GL";
            if (elem.fName === "Table") vtype = "rootui5.eve7.view.EveTable"; // AMT temporary solution
            
            var oOwnerComponent = Component.getOwnerComponentFor(this.getView());
            var view = oOwnerComponent.runAsOwner(function() {
               return new JSROOT.sap.ui.xmlview({
                  id: viewid,
                  viewName: vtype,
                  viewData: { mgr: main.mgr, elementid: elem.fElementId, kind: (count==1) ? "3D" : "2D" },
                  layoutData: oLd
               });
            });

            if (count == 1) {
               sv.addContentArea(view);
               continue;
            }

            if (!vv) {
               vv = new Splitter("SecondaryViewSplitter", { orientation : "Vertical" });
               sv.addContentArea(vv);
            }

            vv.addContentArea(view);
         }
      },

      onManagerUpdate: function() {
         console.log("manager updated");
         this.UpdateCommandsButtons(this.mgr.commands);
         this.updateViewers();
      },

      /*
       * processWaitingMsg: function() { for ( var i = 0; i <
       * msgToWait.length; ++i ) {
       * this.OnWebsocketMsg(handleToWait, msgToWait[i]); }
       * handleToWait = 0; msgToWait = []; },
       */
      event: function() {
         // this._event = lst;

      },

      setMainVerticalSplitterHeight: function(){
         var mainViewHeight = document.body.clientHeight;
         var mainToolbarHeight = 49;
         var height = mainViewHeight - mainToolbarHeight;
         var splitter =  this.getView().byId("MainAreaSplitter");
         if (splitter) {
            // console.log("set splitter height >>> " , height);
            splitter.setHeight(height + "px");
         }
      },

      onAfterRendering: function(){
         var me = this;
         setTimeout(
            function(){
               $(window).on("resize", function() {
                  me.setMainVerticalSplitterHeight();
               });
               me.setMainVerticalSplitterHeight();
            }, 100);
      },

      onToolsMenuAction : function (oEvent) {

         var item = oEvent.getParameter("item");

         switch (item.getText()) {
            case "GED Editor": this.getView().byId("Summary").getController().toggleEditor(); break;
          //  case "Event statusbar": this.toggleShowStatus(); break;
          //  case "Toolbar": this.toggleToolBar(); break;
          //  case "Tooltip info": this.toggleToolTip(); break;
         }
      },

      showHelp : function(oEvent) {
         alert("User support: root-webgui@cern.ch");
      },
      showUserURL : function(oEvent) {
         MobileLibrary.URLHelper.redirect("https://github.com/alja/jsroot/blob/dev/eve7.md", true);
      }
   });
});
