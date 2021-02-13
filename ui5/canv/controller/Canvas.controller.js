sap.ui.define([
   'sap/ui/core/mvc/Controller',
   'sap/ui/core/Component',
   'sap/ui/model/json/JSONModel',
   'sap/ui/core/mvc/XMLView',
   'sap/m/MessageToast',
   'sap/m/Dialog',
   'sap/m/List',
   'sap/m/InputListItem',
   'sap/m/Input',
   'sap/m/Button',
   'sap/ui/layout/Splitter',
   'sap/ui/layout/SplitterLayoutData'
], function (Controller, Component, JSONModel, XMLView, MessageToast, Dialog, List, InputListItem, Input, Button, Splitter, SplitterLayoutData) {
   "use strict";

   let CController = Controller.extend("rootui5.canv.controller.Canvas", {
      onInit : function() {
         this._Page = this.getView().byId("CanvasMainPage");

         let id = this.getView().getId();
         console.log("Initialization CANVAS id = " + id);

         this.bottomVisible = false;

         let model = new JSONModel({ GedIcon: "", StatusIcon: "", ToolbarIcon: "", TooltipIcon: "sap-icon://accept",
                                     StatusLbl1:"", StatusLbl2:"", StatusLbl3:"", StatusLbl4:"" });
         this.getView().setModel(model);

         let vd = this.getView().getViewData();
         let cp = vd ? vd.canvas_painter : null;

         if (!cp) cp = Component.getOwnerComponentFor(this.getView()).getComponentData().canvas_painter;

         if (cp) {
            this.getView().byId("MainPanel").getController().setPainter(cp);

            cp.executeObjectMethod = this.executeObjectMethod.bind(this);

            // overwriting method of canvas with standalone handling of GED
            cp.activateGed = this.openuiActivateGed.bind(this);
            cp.removeGed = this.cleanupIfGed.bind(this);
            cp.hasGed = this.isGedEditor.bind(this);

            cp.hasEventStatus = this.isStatusShown.bind(this);
            cp.activateStatusBar = this.toggleShowStatus.bind(this);
            cp.showCanvasStatus = this.showCanvasStatus.bind(this); // used only for UI5, otherwise global func
            cp.showMessage = this.showMessage.bind(this);
            cp.showSection = this.showSection.bind(this);

            cp.showUI5ProjectionArea = this.showProjectionArea.bind(this);
            cp.drawInUI5ProjectionArea = this.drawInProjectionArea.bind(this);

            cp.showUI5Panel = this.showPanelInLeftArea.bind(this);
         }

         // this.toggleGedEditor();
      },

      executeObjectMethod: function(painter, method, menu_obj_id) {

         if (method.fArgs!==undefined) {
            this.showMethodsDialog(painter, method, menu_obj_id);
            return true;
         }

         if (method.fName == "Inspect") {
            painter.showInspector();
            return true;
         }

         if (method.fName == "FitPanel") {
            this.showLeftArea("FitPanel");
            return true;
         }

         if (method.fName == "Editor") {
            this.openuiActivateGed(painter);
            return true;
         }

         return false; // not processed

      },

      /** @summary function used to activate GED in full canvas */
      openuiActivateGed: function(painter, kind, mode) {

         let canvp = this.getCanvasPainter();

         return this.showGeEditor(true).then(() => {
            canvp.selectObjectPainter(painter);

            if (typeof canvp.processChanges == 'function')
               canvp.processChanges("sbits", canvp);

            return true;
         });
      },

      getCanvasPainter : function(also_without_websocket) {
         let elem = this.getView().byId("MainPanel");

         let p = elem ? elem.getController().getPainter() : null;

         return (p && (p._websocket || also_without_websocket)) ? p : null;
      },

      closeMethodDialog : function(painter, method, menu_obj_id) {

         let args = "";

         if (method) {
            let cont = this.methodDialog.getContent();

            let items = cont[0].getItems();

            if (method.fArgs.length !== items.length)
               alert('Mismatch between method description' + method.fArgs.length + ' and args list in dialog ' + items.length);

            // console.log('ITEMS', method.fArgs.length, items.length);

            for (let k=0;k<method.fArgs.length;++k) {
               let arg = method.fArgs[k];
               let value = items[k].getContent()[0].getValue();

               if (value==="") value = arg.fDefault;

               if ((arg.fTitle=="Option_t*") || (arg.fTitle=="const char*")) {
                  // check quotes,
                  // TODO: need to make more precise checking of escape characters
                  if (!value) value = '""';
                  if (value[0]!='"') value = '"' + value;
                  if (value[value.length-1] != '"') value += '"';
               }

               args += (k>0 ? "," : "") + value;
            }
         }

         this.methodDialog.close();
         this.methodDialog.destroy();

         if (painter && method && args) {

            if (painter.executeMenuCommand(method, args)) return;
            let exec = method.fExec;
            if (args) exec = exec.substr(0,exec.length-1) + args + ')';
            // invoked only when user press Ok button
            console.log('execute method for object ' + menu_obj_id + ' exec= ' + exec);

            let canvp = this.getCanvasPainter();

            if (canvp)
               canvp.sendWebsocket('OBJEXEC:' + menu_obj_id + ":" + exec);
         }
      },

      showMethodsDialog : function(painter, method, menu_obj_id) {

         // TODO: deliver class name together with menu items
         method.fClassName = painter.getClassName();
         if ((menu_obj_id.indexOf("#x")>0) || (menu_obj_id.indexOf("#y")>0) || (menu_obj_id.indexOf("#z")>0)) method.fClassName = "TAxis";

         let items = [];

         for (let n=0;n<method.fArgs.length;++n) {
            let arg = method.fArgs[n];
            arg.fValue = arg.fDefault;
            if (arg.fValue == '\"\"') arg.fValue = "";
            let item = new InputListItem({
               label: arg.fName + " (" +arg.fTitle + ")",
               content: new Input({ placeholder: arg.fName, value: arg.fValue })
            });
            items.push(item);
         }

         this.methodDialog = new Dialog({
            title: method.fClassName + '::' + method.fName,
            content: new List({
                items: items
             }),
             beginButton: new Button({
               text: 'Cancel',
               press: this.closeMethodDialog.bind(this)
             }),
             endButton: new Button({
               text: 'Ok',
               press: this.closeMethodDialog.bind(this, painter, method, menu_obj_id)
             })
         });

         // this.getView().getModel().setProperty("/Method", method);
         //to get access to the global model
         // this.getView().addDependent(this.methodDialog);

         this.methodDialog.addStyleClass("sapUiSizeCompact");

         this.methodDialog.open();
      },

      onFileMenuAction : function (oEvent) {
         //let oItem = oEvent.getParameter("item"),
         //    sItemPath = "";
         //while (oItem instanceof sap.m.MenuItem) {
         //   sItemPath = oItem.getText() + " > " + sItemPath;
         //   oItem = oItem.getParent();
         //}
         //sItemPath = sItemPath.substr(0, sItemPath.lastIndexOf(" > "));

         let p = this.getCanvasPainter();
         if (!p) return;

         let name = oEvent.getParameter("item").getText();

         switch (name) {
            case "Close canvas": p.onWebsocketClosed(); p.closeWebsocket(true); break;
            case "Interrupt": p.sendWebsocket("INTERRUPT"); break;
            case "Quit ROOT": p.sendWebsocket("QUIT"); break;
            case "Canvas.png":
            case "Canvas.jpeg":
            case "Canvas.svg":
               p.saveCanvasAsFile(name);
               break;
            case "Canvas.root":
            case "Canvas.pdf":
            case "Canvas.ps":
            case "Canvas.C":
               p.sendSaveCommand(name);
               break;
         }

         MessageToast.show("Action triggered on item: " + name);
      },

      onCloseCanvasPress : function() {
         let p = this.getCanvasPainter();
         if (p) {
            p.onWebsocketClosed();
            p.closeWebsocket(true);
         }
      },

      onInterruptPress : function() {
         let p = this.getCanvasPainter();
         if (p) p.sendWebsocket("INTERRUPT");
      },

      onQuitRootPress : function() {
         let p = this.getCanvasPainter();
         if (p) p.sendWebsocket("QUIT");
      },

      onReloadPress : function() {
         let p = this.getCanvasPainter();
         if (p) p.sendWebsocket("RELOAD");
      },

      isGedEditor : function() {
         return this.getView().getModel().getProperty("/LeftArea") == "Ged";
      },

      showGeEditor : function(new_state) {
         return this.showLeftArea(new_state ? "Ged" : "");
      },

      cleanupIfGed: function() {
         let ged = this.getLeftController("Ged"),
             p = this.getCanvasPainter();
         if (p) p.registerForPadEvents(null);
         if (ged) ged.cleanupGed();
         if (p && p.processChanges) p.processChanges("sbits", p);
      },

      getLeftController: function(name) {
         if (this.getView().getModel().getProperty("/LeftArea") != name) return null;
         let split = this.getView().byId("MainAreaSplitter");
         return split ? split.getContentAreas()[0].getController() : null;
      },

      toggleGedEditor : function() {
         this.showGeEditor(!this.isGedEditor());
      },

      showPanelInLeftArea: function(panel_name, panel_handle) {

         let split = this.getView().byId("MainAreaSplitter");
         let curr = this.getView().getModel().getProperty("/LeftArea");
         if (!split || (curr === panel_name))
            return Promise.resolve(false);

         // first need to remove existing
         if (curr) {
            console.log('REMOVE CURRENT AREA', curr);
            this.cleanupIfGed();
            split.removeContentArea(split.getContentAreas()[0]);
         }

         this.getView().getModel().setProperty("/LeftArea", panel_name);
         this.getView().getModel().setProperty("/GedIcon", (panel_name=="Ged") ? "sap-icon://accept" : "");

         if (!panel_handle || !panel_name)
            return Promise.resolve(false);

         let oLd = new SplitterLayoutData({
            resizable : true,
            size      : "250px"
         });

         let viewName = panel_name;
         if (viewName.indexOf(".") < 0) viewName = "rootui5.canv.view." + panel_name;

         let can_elem = this.getView().byId("MainPanel");

         return XMLView.create({
            viewName: viewName,
            viewData: { handle: panel_handle, masterPanel: this },
            layoutData: oLd,
            height: (panel_name == "Panel") ? "100%" : undefined
         }).then(oView => {
            // workaround, while CanvasPanel.onBeforeRendering called too late
            can_elem.getController().preserveCanvasContent();
            split.insertContentArea(oView, 0);
            return true;
         });

      },

      // TODO: sync with showPanelInLeftArea, it is more or less same
      showLeftArea: function(panel_name) {
         let split = this.getView().byId("MainAreaSplitter");
         let curr = this.getView().getModel().getProperty("/LeftArea");
         if (!split || (curr === panel_name))
            return Promise.resolve(null);

         // first need to remove existing
         if (curr) {
            this.cleanupIfGed();
            split.removeContentArea(split.getContentAreas()[0]);
         }

         this.getView().getModel().setProperty("/LeftArea", panel_name);
         this.getView().getModel().setProperty("/GedIcon", (panel_name=="Ged") ? "sap-icon://accept" : "");

         if (!panel_name) return Promise.resolve(null);

         let oLd = new SplitterLayoutData({
            resizable: true,
            size: "250px"
         });

         let canvp = this.getCanvasPainter();

         let viewName = "rootui5.canv.view." + panel_name;
         if (panel_name == "FitPanel") viewName = "rootui5.fitpanel.view.FitPanel";

         let can_elem = this.getView().byId("MainPanel");

         return XMLView.create({
            viewName: viewName,
            viewData: { masterPanel: this },
            layoutData: oLd,
            height: (panel_name == "Panel") ? "100%" : undefined
         }).then(oView => {

            // workaround, while CanvasPanel.onBeforeRendering called too late
            can_elem.getController().preserveCanvasContent();

            split.insertContentArea(oView, 0);

            if (panel_name === "Ged") {
               let ged = oView.getController();
               if (canvp && ged && (typeof canvp.registerForPadEvents == "function")) {
                  canvp.registerForPadEvents(ged.padEventsReceiver.bind(ged));
                  canvp.selectObjectPainter(canvp);
               }
            }

            return oView.getController();
         });
      },

      getBottomController : function() {
         if (!this.bottomVisible) return null;
         let split = this.getView().byId("MainAreaSplitter"),
             cont = split.getContentAreas(),
             vsplit = cont[cont.length-1],
             vcont = vsplit.getContentAreas(),
             bottom = vcont[vcont.length-1];
         return bottom ? bottom.getController() : null;
      },

      drawInProjectionArea : function(can, opt) {
         let ctrl = this.getBottomController();
         if (!ctrl) ctrl = this.getLeftController("Panel");

         if (ctrl && ctrl.drawObject)
            return ctrl.drawObject(can, opt);

         return Promise.resolve(null);
      },

      showProjectionArea : function(kind) {
         let bottom = null;
         return this.showBottomArea(kind == "X")
             .then(area => { bottom = area; return this.showLeftArea(kind == "Y" ? "Panel" : ""); })
             .then(left => {

               let ctrl = bottom || left;

               if (!ctrl || ctrl.getView().getDomRef())
                  return Promise.resolve(!!ctrl);

               return ctrl.getRenderPromise();
            });
      },

      showBottomArea : function(is_on) {

         if (this.bottomVisible == is_on)
            return Promise.resolve(this.getBottomController());

         let split = this.getView().byId("MainAreaSplitter");
         if (!split) return Promise.resolve(null);

         let cont = split.getContentAreas();

         this.bottomVisible = !this.bottomVisible;

         if (!this.bottomVisible) {
            // vertical splitter exists - toggle it

            let vsplit = cont[cont.length-1],
                main = vsplit.removeContentArea(0);

            vsplit.destroyContentAreas();
            split.removeContentArea(vsplit);
            split.addContentArea(main);
            return Promise.resolve(null);
         }

         // remove panel with normal drawing
         split.removeContentArea(cont[cont.length-1]);
         let vsplit = new Splitter({orientation: "Vertical"});

         split.addContentArea(vsplit);

         vsplit.addContentArea(cont[cont.length-1]);

         let oLd = new SplitterLayoutData({
            resizable : true,
            size      : "200px"
         });

         return XMLView.create({
            viewName : "rootui5.canv.view.Panel",
            layoutData: oLd,
            height: "100%"
         }).then(oView => {
            vsplit.addContentArea(oView);
            return oView.getController();
         });

      },

      showCanvasStatus : function (text1,text2,text3,text4) {
         let model = this.getView().getModel();
         model.setProperty("/StatusLbl1", text1);
         model.setProperty("/StatusLbl2", text2);
         model.setProperty("/StatusLbl3", text3);
         model.setProperty("/StatusLbl4", text4);
      },

      isStatusShown : function() {
         return this._Page.getShowFooter();
      },

      toggleShowStatus : function(new_state) {
         if (new_state === undefined) new_state = !this.isStatusShown();

         this._Page.setShowFooter(new_state);
         this.getView().getModel().setProperty("/StatusIcon", new_state ? "sap-icon://accept" : "");

         let canvp = this.getCanvasPainter();
         if (canvp) canvp.processChanges("sbits", canvp);
      },

      toggleToolBar : function(new_state) {
         if (new_state === undefined) new_state = !this.getView().getModel().getProperty("/ToolbarIcon");

         this._Page.setShowSubHeader(new_state);

         this.getView().getModel().setProperty("/ToolbarIcon", new_state ? "sap-icon://accept" : "");
      },

      toggleToolTip : function(new_state) {
         if (new_state === undefined) new_state = !this.getView().getModel().getProperty("/TooltipIcon");

         this.getView().getModel().setProperty("/TooltipIcon", new_state ? "sap-icon://accept" : "");

         let p = this.getCanvasPainter(true);
         if (p) p.setTooltipAllowed(new_state);
      },

      setShowMenu: function(new_state) {
         this._Page.setShowHeader(new_state);
      },

      onViewMenuAction: function (oEvent) {

         let item = oEvent.getParameter("item");

         switch (item.getText()) {
            case "Editor": this.toggleGedEditor(); break;
            case "Event statusbar": this.toggleShowStatus(); break;
            case "Toolbar": this.toggleToolBar(); break;
            case "Tooltip info": this.toggleToolTip(); break;
         }
      },

      onToolsMenuAction : function(oEvent) {
         let item = oEvent.getParameter("item"),
             name = item.getText();

         if (name != "Fit panel") return;

         let curr = this.getView().getModel().getProperty("/LeftArea");

         this.showLeftArea(curr == "FitPanel" ? "" : "FitPanel");
      },

      showMessage : function(msg) {
         MessageToast.show(msg);
      },

      showSection : function(that, on) {
         // this function call when section state changed from server side
         switch(that) {
            case "Menu": this.setShowMenu(on); break;
            case "StatusBar": this.toggleShowStatus(on); break;
            case "Editor": return this.showGeEditor(on);
            case "ToolBar": this.toggleToolBar(on); break;
            case "ToolTips": this.toggleToolTip(on); break;
         }
         return Promise.resolve(true);
      }
   });

   return CController;

});
