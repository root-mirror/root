sap.ui.define([
   'sap/ui/core/mvc/Controller',
   "sap/ui/model/json/JSONModel",
   "sap/ui/core/mvc/XMLView",
   "sap/m/CustomTreeItem",
   "sap/m/FlexBox",
   "sap/m/CheckBox",
   "sap/m/Text",
   "sap/m/Button",
   "sap/ui/layout/SplitterLayoutData",
   "rootui5/eve7/controller/Ged.controller"
], function(Controller, JSONModel, XMLView, CustomTreeItem,
            FlexBox, mCheckBox, mText, mButton, SplitterLayoutData, GedController) {

   "use strict";

   var EveSummaryCustomItem = CustomTreeItem.extend('rootui5.eve7.lib.EveSummaryCustomItem', {
      renderer: {},

      metadata: {
         properties: {
            elementId: 'string',
            background: 'string',
            mainColor: 'string',
            showCheckbox: 'boolean',
            showRnrChildren: 'boolean'
         }
      },

      onAfterRendering: function()
      {
         var btn = this.getContent()[0].getItems()[1];
         btn.$().css('background-color', this.getMainColor());
      }

   });

   return Controller.extend("rootui5.eve7.controller.Summary", {

      onInit: function () {

         var data = [{ fName: "Event" }];

         this.summaryElements = {}; // object with all elements, used for fast access to elements by id

         var oTree = this.getView().byId("tree");
         this.expandLevel = 2;

         var oModel = new JSONModel();
         oModel.setData([]);
         oModel.setSizeLimit(10000);
         oModel.setDefaultBindingMode("OneWay");
         this.getView().setModel(oModel, "treeModel");

         var oItemTemplate = new EveSummaryCustomItem({
            content: [
                new FlexBox({
                   width: "100%",
                   alignItems: "Start",
                   justifyContent: "SpaceBetween",
                   items: [
                     new FlexBox({
                        alignItems: "Start",
                        items: [
                           new mCheckBox({ visible: "{treeModel>fShowCheckbox}", selected: "{treeModel>fSelected}", select: this.clickItemSelected.bind(this) }),
                           new mText({text:" {treeModel>fName}", tooltip: "{treeModel>fTitle}" , renderWhitespace: true, wrapping: false })
                         ]
                      }),
                      new mButton({ id: "detailBtn", visible: "{treeModel>fShowButton}", icon: "sap-icon://edit", type: "Transparent", tooltip: "Actiavte GED", press: this.pressGedButton.bind(this) })
                    ]
                })
            ],

            elementId: "{treeModel>fElementId}",
            mainColor: "{treeModel>fMainColor}",
            showCheckbox: "{treeModel>fShowCheckbox}",
            showRnrChildren: "{treeModel>fShowRnrChildren}"
         });

         oItemTemplate.addStyleClass("eveSummaryItem");
         oItemTemplate.attachBrowserEvent("mouseenter", this.onMouseEnter, this);
         oItemTemplate.attachBrowserEvent("mouseleave", this.onMouseLeave, this);
         oTree.bindItems("treeModel>/", oItemTemplate);

         this.template = oItemTemplate;

         this.rebuild = false;
      },

      clickItemSelected: function(oEvent) {
         var item = oEvent.getSource().getParent().getParent().getParent();

         var selected = oEvent.getSource().getSelected();

         var id = item.getElementId();

         if (id === undefined)
            return console.log('clickItemSelected: fail to extract element id');

         var elem = this.mgr.GetElement(id);

         if (!elem)
            return console.log('clickItemSelected: fail to find element', id);

         var method = item.getShowRnrChildren() ? "SetRnrChildren" : GedController.GetRnrSelfMethod(elem._typename)

         this.mgr.SendMIR(method + "(" + selected + ")", elem.fElementId, elem._typename);
      },

      pressGedButton: function(oEvent) {
         var item = oEvent.getSource().getParent().getParent();

         // var path = item.getBindingContext("treeModel").getPath(),
         //    ttt = item.getBindingContext("treeModel").getProperty(path);

         this.showGedEditor(item.getElementId());
      },

      SetMgr: function(mgr) {
         this.mgr = mgr;
         this.mgr.RegisterController(this);
         this.selected = {}; // container of selected objects
      },

      OnEveManagerInit: function() {
         var model = this.getView().getModel("treeModel");
         model.setData(this.createModel());
         model.refresh();

         var oTree = this.getView().byId("tree");
         oTree.expandToLevel(this.expandLevel);

         // hide editor
         if (this.ged)
            this.ged.getController().closeGedEditor();

         var scenes = this.mgr.getSceneElements();
         for (var i = 0; i < scenes.length; ++i) {
            this.mgr.RegisterSceneReceiver(scenes[i].fElementId, this);
         }
      },

      onToggleOpenState: function(oEvent) {
      },

      processHighlight: function(kind, evid, force) {

         if (!force) {
            if (this._trigger_timer)
               clearTimeout(this._trigger_timer);

            this._trigger_timer = setTimeout(this.processHighlight.bind(this,kind,evid,true), 200);
            return;
         }

         delete this._trigger_timer;

         var objid = 0;

         if (kind != "leave") {
            var tree = this.getView().byId("tree"),
                items = tree.getItems(true), item = null;
            for (var n = 0; n < items.length; ++n)
               if (items[n].getId() == evid) {
                  item = items[n]; break;
               }

            if (item) objid = item.getElementId() || 0;
         }

         // FIXME: provide more generic code which should
         this.mgr.SendMIR("NewElementPicked(" + objid + ",false,false)",
                          this.mgr.global_highlight_id, "ROOT::Experimental::REveSelection");
      },

      onMouseEnter: function(oEvent) {
         this.processHighlight("enter", oEvent.target.id);
      },

      onMouseLeave: function(oEvent) {
         this.processHighlight("leave");
      },

      GetSelectionColor: function(selection_obj) {
         return selection_obj.fName == "Global Highlight" ? "rgb(230, 230, 230)" : "rgb(66, 124, 172)";
      },

      FindTreeItemForEveElement:function(element_id) {
         var items = this.getView().byId("tree").getItems();
         for (var n = 0; n<items.length;++n)
            if (items[n].getElementId() == element_id)
               return items[n];
         return null;
      },

      SelectElement: function(selection_obj, element_id, sec_idcs) {
         var item = this.FindTreeItemForEveElement(element_id);
         if (item) {
            var color = this.GetSelectionColor(selection_obj);
            item.$().css("background-color", color);
            if (this.ged && (this.mgr.global_selection_id == selection_obj.fElementId)) this.ged.getController().updateSecondarySelectionGED(element_id, sec_idcs);
         }

      },

      UnselectElement: function (selection_obj, element_id) {
         var item = this.FindTreeItemForEveElement(element_id);
         if (item) {
            var color = this.GetSelectionColor(selection_obj);
            var cc = item.$().css("background-color");
            if (cc == color)
               item.$().css("background-color", "");
             if (this.ged && (this.mgr.global_selection_id == selection_obj.fElementId)) this.ged.getController().updateSecondarySelectionGED();
         }
      },

      toggleEditor: function() {
         var showid = -1;
         if (!this.ged || !this.ged.getController().isGedVisible())
            for (var id in this.summaryElements) {
               if (this.summaryElements[id].fShowButton) { showid = id; break; }
            }

         if (showid > 0)
            this.showGedEditor(showid);
         else if (this.ged)
            this.ged.getController().closeGedEditor();
      },

      showGedEditor: function(elementId) {

         var sumSplitter = this.byId("sumSplitter");

         if (!this.ged) {
            var pthis = this;

            XMLView.create({
               viewName: "rootui5.eve7.view.Ged",
               layoutData: new SplitterLayoutData("sld", {size: "30%"}),
               height: "100%"
            }).then(function(oView) {
               pthis.ged = oView;
               pthis.ged.getController().setManager(pthis.mgr);
               pthis.ged.getController().showGedEditor(sumSplitter, elementId);

            });
         } else {
            this.ged.getController().showGedEditor(sumSplitter, elementId);
         }
      },

      anyVisible: function(arr) {
         if (!arr) return false;
         for (var k=0;k<arr.length;++k) {
            if (arr[k].fName) return true;
         }
         return false;
      },

      /** Set summary element attributes from original element */
      setElementsAttributes: function(newelem, elem) {
         newelem.fShowCheckbox = false;
         newelem.fShowRnrChildren = false;
         newelem.fElementId = elem.fElementId;
         newelem.fShowButton = false;
         newelem.fMainColor = "";

         if (GedController.canEditClass(elem._typename)) {
            newelem.fShowButton = true;

            if (!elem.childs) {
               newelem.fShowCheckbox = true;
               newelem.fSelected = elem.fRnrSelf;
            }  else if (elem.fRnrChildren !== undefined) {
               newelem.fShowCheckbox = true;
               newelem.fShowRnrChildren = true;
               newelem.fSelected = elem.fRnrChildren;
            }

            if (elem.fMainColor) {
               newelem.fMainColor = JSROOT.Painter.getColor(elem.fMainColor);
            }
         }
      },

      createModel: function() {
         this.summaryElements = {};

         /*
         var src = this.mgr.childs[0].childs[2].childs;
         for (var i = 0; i < src.length; i++) {
            if (src[i].fName == "Collections")
               src = src[i].childs;
         }
         */
         let src = this.mgr.childs;
         return this.createSummaryModel([], src, "/");
      },

      createSummaryModel: function(tgt, src, path) {
         for (var n=0;n<src.length;++n) {
            var elem = src[n];

            var newelem = { fName: elem.fName, fTitle: elem.fTitle || elem.fName, id: elem.fElementId, fHighlight: "None", fBackground: "", fMainColor: "", fSelected: false };

            this.setElementsAttributes(newelem, elem);

            newelem.path = path + n;
            newelem.masterid = elem.fMasterId || elem.fElementId;

            tgt.push(newelem);

            this.summaryElements[newelem.id] = newelem;

            if ((elem.childs !== undefined) && this.anyVisible(elem.childs))
               newelem.childs = this.createSummaryModel([], elem.childs, newelem.path + "/childs/");
         }

         return tgt;
      },

      beginChanges: function() {
        // this.rebuild=false;
      },

      elementsRemoved: function(ids) {
         this.rebuild = true;
      },

      sceneElementChange: function(msg) {

         if (this.ged)
            this.ged.getController().updateGED(msg.fElementId);

         var newelem = this.summaryElements[msg.fElementId];

         var elem = this.mgr.GetElement(msg.fElementId);

         if (newelem && elem)
            this.setElementsAttributes(newelem, elem);

         // console.log('SUMMURY: detect changed', elem.id, elem.path);

         if (msg.rnr_self_changed)
            elem.fSelected = msg.fRnrSelf;

         this.any_changed = true;

      },

      endChanges: function() {
         if (this.rebuild) {
            var oTree = this.getView().byId("tree");
            oTree.unbindItems();

            var model = this.getView().getModel("treeModel");
            model.setData(this.createModel());
            model.refresh();

            this.getView().setModel(model, "treeModel");
            oTree.bindItems("treeModel>/", this.template);
            oTree.setModel(model, "treeModel");

            oTree.expandToLevel(this.expandLevel);

            if (this.ged)
               this.ged.getController().closeGedEditor();

            this.rebuild = false;
         } else if (this.any_changed) {
            var model = this.getView().getModel("treeModel");
            model.refresh();

            this.any_changed = false;
         }
      },

      /** Invoked via EveManager when specified element should be focused */
      BrowseElement: function(elid) {
         var summaryElement = this.summaryElements[elid];
         if (!summaryElement) return;

         var oTree = this.getView().byId("tree"),
             element_path = summaryElement.path,
             bestindx = 1, bestlen = 0;

         while (bestindx >= 0) {

            bestindx = -1;

            var items = oTree.getItems();

            for (var k = 0; k < items.length; ++k) {
               var item = items[k],
                   model = item.getBindingContext("treeModel"),
                   path = model.getPath();

               if (element_path == path) {
                  var dom = item.$()[0];
                  if (dom) dom.scrollIntoView();
                  return;
               }

               if ((element_path.substr(0, path.length) == path) && (path.length > bestlen)) {
                  bestindx = k;
                  bestlen = path.length;
               }
            }

            if (bestindx >= 0) oTree.expand(bestindx);
         }

      }
   });
});
