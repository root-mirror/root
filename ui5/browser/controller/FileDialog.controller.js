sap.ui.define(['rootui5/panel/Controller',
               "sap/ui/model/json/JSONModel",
               "sap/ui/core/Fragment",
               "sap/m/Link",
               "sap/m/Text",
               "sap/m/Button",
               "sap/m/ButtonType",
               "sap/m/Dialog"
],function(GuiPanelController, JSONModel, Fragment, Link, Text, Button, ButtonType, Dialog) {

   "use strict";

   /** FileDialog controller */

   let FileDialog = GuiPanelController.extend("rootui5.browser.controller.FileDialog", {

      //function called from GuiPanelController
      onPanelInit : function() {

         // TODO: provide functionality via BrowserModel - once we know that exactly we need

         // create model only for browser - no need for anybody else
         // this.model = new BrowserModel();

         // copy extra attributes from element to node in the browser
         // later can be done automatically
         /* this.model.addNodeAttributes = function(node, elem) {
            node.icon = elem.icon;
            node.fsize = elem.fsize;
            node.mtime = elem.mtime;
            node.ftype = elem.ftype;
            node.fuid = elem.fuid;
            node.fgid = elem.fgid;
            node.className = elem.className
         };
         */

         this.kind = "None"; // not yet known
         this.oModel = new JSONModel({ canEditName: (this.kind == "SaveAs") || (this.kind == "NewFile"),
                                       dialogTitle: "Dialog Title",
                                       fileName: "",
                                       filesList: [{name:"first.txt", counter: 11}, {name:"second.txt", counter: 22}, {name:"third.xml", counter: 33}],
                                       fileExt: "AllFiles",
                                       showFileExt: false,
                                       fileExtList: [] });
                                      // fileExtList: [{ id: "AllFiles", text: "All files (*.*)" }, { id: "png", text: "png files (*.png)"}, { id: "cxx", text: "CXX files (*.cxx)"}] });
         this.getView().setModel(this.oModel);

         this.own_window = true;

         Fragment.load({
            name: "rootui5.browser.view.filedialog",
            controller: this,
            id: "FileDialogFragment"
         }).then(oFragment => {
            this.fragment = oFragment;
            this.getView().addDependent(this.fragment);
            this.getView().byId("dialogPage").addContent(this.fragment);
            this.fragment.setModel(this.oModel);

            if (this._init_msg) {
               this.processInitMsg(this._init_msg);
               delete this._init_msg;
            }
         });

      },

      /** @summary Set path to the Breadcrumb element */
      updateBReadcrumbs: function(split) {
         this._currentPath = split;

         let oBreadcrumbs = sap.ui.core.Fragment.byId("FileDialogFragment", "breadcrumbs");
         if (!oBreadcrumbs)
            return;

         oBreadcrumbs.removeAllLinks();
         oBreadcrumbs.setCurrentLocationText("");
         for (let i=-1; i<split.length; i++) {
            let txt = i<0 ? "/": split[i];
            if (i === split.length-1) {
               oBreadcrumbs.setCurrentLocationText(txt);
            } else {
               let link = new Link({text: txt});
               link.attachPress(this, this.onBreadcrumbsPress, this);
               oBreadcrumbs.addLink(link);
            }
         }
      },

      /** @summary Returns coded in Breadcrumb path
       * @desc If selectedId specified, return path up to that element id */
      getBreadcrumbPath: function(selectedId) {
         let oBreadcrumbs = sap.ui.core.Fragment.byId("FileDialogFragment", "breadcrumbs");
         if (selectedId && oBreadcrumbs) {
            let oLinks = oBreadcrumbs.getLinks(), path = [];
            for (let i = 0; i < oLinks.length; i++) {
               if (i>0) path.push(oLinks[i].getText());
               if (oLinks[i].getId() === selectedId) return path;
            }
         }

         return this._currentPath.slice(); // make copy of original array
      },

      // returns full file name as array
      getFullFileName: function() {
         let path = this.getBreadcrumbPath();
         path.push(this.oModel.getProperty("/fileName"));
         return path;
      },

      /** Handler for Breadcrumbs press event */
      onBreadcrumbsPress: function(oEvent) {
         let path = this.getBreadcrumbPath(oEvent.getSource().sId);
         if (this.kind == "OpenFile")
            this.oModel.setProperty("/fileName", "");
         this.websocket.send("CHPATH:" + JSON.stringify(path));
      },

      /** @summary Handler for List item press event */
      onItemPress: function(event) {
         let item = event.getParameters().listItem;
         if (!item) return;

         let tm = new Date().getTime();

         if (this._last_item === item) {
            let diff = tm - this._last_tm;
            if ((diff < 300) && (this.kind == "OpenFile"))
               return this.onOkPress();
         }

         this._last_item = item;
         this._last_tm = tm;

         let ctxt = item.getBindingContext(),
             prop = ctxt ? ctxt.getProperty(ctxt.getPath()) : null;

         if (prop && (prop.icon == "sap-icon://folder-blank")) {
            if (this.kind == "OpenFile")
               this.oModel.setProperty("/fileName", "");

            let path = this.getBreadcrumbPath();
            path.push(item.getTitle());
            return this.websocket.send("CHPATH:" + JSON.stringify(path));
         }

         this.oModel.setProperty("/fileName", item.getTitle());
      },

      /** @summary When selected file extenstion changed */
      onFileExtChanged: function() {
         let extName = this.oModel.getProperty("/fileExt");
         this.websocket.send("CHEXT:" + extName);
      },

      processInitMsg: function(msg) {
         let cfg = JSON.parse(msg);

         if (!this.dialog) {
            // when not an embedded dialog, update configuration from server
            this.kind = cfg.kind;
            this.oModel.setProperty("/dialogTitle", cfg.title);
            this.oModel.setProperty("/canEditName", (this.kind == "SaveAs") || (this.kind == "NewFile"));
         }

         if (this.own_window && document)
             document.title = cfg.title || this.kind + " dialog";

         this.updateBReadcrumbs(cfg.path);

         this.oModel.setProperty("/fileName", cfg.fname);

         if (cfg.filters && cfg.filters.length) {
            let arr = [];
            for (let k = 0; k < cfg.filters.length; ++k) {
               let fname = cfg.filters[k];
               let p = fname.indexOf("(");
               if (p>0) fname = fname.substr(0,p);
               arr.push({id: fname.trim(), text: cfg.filters[k]});
            }
            this.oModel.setProperty("/fileExt", cfg.filter);
            this.oModel.setProperty("/fileExtList", arr);
            this.oModel.setProperty("/showFileExt", true);
         } else {
            this.oModel.setProperty("/fileExt", "All files");
            this.oModel.setProperty("/showFileExt", false);
         }

         this.oModel.setProperty("/filesList", cfg.brepl.nodes);
      },

      processChangePathMsg: function(msg) {
         let cfg = JSON.parse(msg);
         this.updateBReadcrumbs(cfg.path);
         this.oModel.setProperty("/filesList", cfg.brepl.nodes);
      },

      /** @summary Close file dialog */
      closeFileDialog: function(fname) {
         // add more logic when FileDialog embed into main window
         if (this.did_close) return;

         this.did_close = true;

         if (fname !== undefined)
            this._completeDialog(fname ? "onOk" : "onCancel", fname);

         if (this.websocket) {
            this.websocket.close();
            delete this.websocket;
         }

         if (this.dialog) {
            this.dialog.close();
            this.dialog.destroy();
            delete this.dialog;
         }

         if (this.own_window && window)
            window.open('','_self').close();
      },

      onWebsocketOpened: function(handle) {
         if (this.model)
            this.model.sendFirstRequest(this.websocket);
      },

      onWebsocketClosed: function() {
         // when connection closed, close panel as well
         this.closeFileDialog();
      },

      /** @summary Entry point for all data from server */
      onWebsocketMsg: function(handle, msg) {
         if (typeof msg != "string")
            return console.error("Browser do not uses binary messages len = " + mgs.byteLength);

         let mhdr = msg.split(":")[0];
         msg = msg.substr(mhdr.length+1);

         switch (mhdr) {
         case "INMSG":
            if (!this.fragment)
               this._init_msg = msg;
            else
               this.processInitMsg(msg);
            break;
         case "CHMSG":
            this.processChangePathMsg(msg);
            break;
         case "SELECT_CONFIRMED": // when selected file can be used for SaveAs operation
            this.closeFileDialog(msg);
            break;
         case "NEED_CONFIRM": // need confirmation warning
            this.showWarningDialog();
            break;
         case "NOSELECT_CONFIRMED":
            this.closeFileDialog("");
            break;

         default:
            console.error('Non recognized msg ' + mhdr + ' len=' + msg.length);
         }
      },

      /** Shown when warning message about overwritten file should appear */
      showWarningDialog: function() {
         let oWarnDlg = new Dialog({
            title: "Warning",
            type: "Message",
            state: "Warning",
            content: new Text({ text: "File already exists. Overwrite? "}),
            beginButton: new Button({
               text: 'Cancel',
               press: ()  => oWarnDlg.close()
            }),
            endButton: new Button({
               text: 'Ok',
               type: ButtonType.Emphasized,
               press: () => {
                  oWarnDlg.close();
                  this.websocket.send("DLG_CONFIRM_SELECT");
               }
            }),
            afterClose: () => oWarnDlg.destroy()
         });

         oWarnDlg.open();
      },

      /** @summary method used to complete dialog */
      _completeDialog: function(funcname, arg) {
         if (!this.dialog_args)
            return;

         if (!funcname) funcname = "onFailure";

         if (typeof this.dialog_args[funcname] == "function")
            this.dialog_args[funcname](arg);

         delete this.dialog_args;
      },

      /** @summary Start SaveAs dialog
        * @private */
      _initDialog: async function(kind, args) {

         if (!args || typeof args != "object")
            return null;

         this.dialog_args = args || {};
         if (!args.websocket)
            return this._completeDialog("onFailure");

         let fname = args.filename || "";
         let p = Math.max(fname.lastIndexOf("/"), fname.lastIndexOf("\\"));
         if (p>0) fname = fname.substr(p+1);

         this.kind = kind;
         this.oModel = new JSONModel({ canEditName: (this.kind == "SaveAs") || (this.kind == "NewFile"),
                                       dialogTitle: args.title || "Title",
                                       fileName: fname, // will be returned from the server, just for initialization
                                       filesList: [{name:"first.txt", counter: 11}, {name:"second.txt", counter: 22}, {name:"third.xml", counter: 33}],
                                       fileExt: "All files",
                                       showFileExt: false,
                                       fileExtList: [] });
                                       // fileExtList: [{ id: "AllFiles", text: "All files (*.*)" }, { id: "png", text: "png files (*.png)"}, { id: "cxx", text: "CXX files (*.cxx)"}] });

         // create extra channel for the FileDialog
         this.websocket = args.websocket.createChannel();
         // assign ourself as receiver of all
         this.websocket.setReceiver(this);

         this.fragment = await Fragment.load({
            name: "rootui5.browser.view.filedialog",
            controller: this,
            id: "FileDialogFragment"
         });

         this.fragment.setModel(this.oModel);

         this.dialog = new Dialog({
            title: "{/dialogTitle}",
            contentWidth: args.width || "70%",
            contentHeight: args.height || "50%",
            resizable: (args.resizable === undefined) ? true : args.resizable,
            draggable: (args.draggable === undefined) ? true : args.draggable,
            content: this.fragment,
            beginButton: new Button({
               text: 'Cancel',
               press: () => this.onCancelPress()
            }),
            endButton: new Button({
               text: 'Ok',
               enabled: "{= ${/fileName} !== '' }",
               press: () => this.onOkPress()
            })
         });

         this.dialog.addStyleClass("sapUiSizeCompact");

         this.dialog.setModel(this.oModel);

         this.dialog.open();

         if (this._init_msg) {
            // probably never happens here, but keep it for completnece
            this.processInitMsg(this._init_msg);
            delete this._init_msg;
         }

         let server_args = [ this.kind, args.filename || "", this.websocket.getChannelId().toString() ];

         // add at the end filter and filter array
         if (args.filter || args.filters) {
            server_args.push(args.filter || "Any files");
            server_args = server_args.concat(args.filters || ["Any files (*)"]);
         }

         args.websocket.send("FILEDIALOG:" + JSON.stringify(server_args));

         return this;
      },

      /** @summary Press Ok button id Dialog,
        * @desc send selected file name and wait if confirmation required */
      onOkPress: function() {
         let fullname = this.getFullFileName();

         if (this.websocket)
            this.websocket.send("DLGSELECT:" + JSON.stringify(fullname));
         else
            this.closeFileDialog();
      },

      /** Press Cancel button Dialog */
      onCancelPress: function() {
         if (this.websocket)
            this.websocket.send("DLGNOSELECT");
         else
            this.closeFileDialog();
      },

      onBeforeRendering: function() {
         this.renderingDone = false;
      },

      onAfterRendering: function() {
         this.renderingDone = true;
      }
   });

   /** Function to initiate SaveAs dialog from client side
    * Following arguments has to be specified:
    * args.websocket - current available connection, used to send "FILEDIALOG:" request
    * args.filename - initial file name in the dialog
    * args.title - title used for the dialog
    * args.onOk - handler when file is selected and "Ok" button is pressed
    * args.onCancel - handler when "Cancel" button is pressed
    * args.onFailure - handler when any failure appears, dialog will be closed afterwards
    */
   FileDialog.SaveAs = function(args) {
      let controller = new FileDialog();

      return controller._initDialog("SaveAs", args);
   }

   FileDialog.NewFile = function(args) {
      let controller = new FileDialog();

      return controller._initDialog("NewFile", args);
   }

   FileDialog.OpenFile = function(args) {
      let controller = new FileDialog();

      return controller._initDialog("OpenFile", args);
   }


   return FileDialog;

});
