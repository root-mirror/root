sap.ui.define([
   "sap/ui/core/mvc/Controller",
   "sap/ui/model/json/JSONModel",
   "sap/m/Button",
   "sap/m/Input",
   "sap/m/StepInput",
   "sap/m/CheckBox",
   "sap/m/Text",
   "sap/m/ColorPalettePopover",
   "sap/ui/layout/HorizontalLayout"
], function (Controller, JSONModel, Button, mInput, mStepInput, mCheckBox, mText, ColorPalettePopover, HorizontalLayout) {
   "use strict";

   var UI5PopupColors = {
         aliceblue: 'f0f8ff',
         antiquewhite: 'faebd7',
         aqua: '00ffff',
         aquamarine: '7fffd4',
         azure: 'f0ffff',
         beige: 'f5f5dc',
         bisque: 'ffe4c4',
         black: '000000',
         blanchedalmond: 'ffebcd',
         blue: '0000ff',
         blueviolet: '8a2be2',
         brown: 'a52a2a',
         burlywood: 'deb887',
         cadetblue: '5f9ea0',
         chartreuse: '7fff00',
         chocolate: 'd2691e',
         coral: 'ff7f50',
         cornflowerblue: '6495ed',
         cornsilk: 'fff8dc',
         crimson: 'dc143c',
         cyan: '00ffff',
         darkblue: '00008b',
         darkcyan: '008b8b',
         darkgoldenrod: 'b8860b',
         darkgray: 'a9a9a9',
         darkgrey: 'a9a9a9',
         darkgreen: '006400',
         darkkhaki: 'bdb76b',
         darkmagenta: '8b008b',
         darkolivegreen: '556b2f',
         darkorange: 'ff8c00',
         darkorchid: '9932cc',
         darkred: '8b0000',
         darksalmon: 'e9967a',
         darkseagreen: '8fbc8f',
         darkslateblue: '483d8b',
         darkslategray: '2f4f4f',
         darkslategrey: '2f4f4f',
         darkturquoise: '00ced1',
         darkviolet: '9400d3',
         deeppink: 'ff1493',
         deepskyblue: '00bfff',
         dimgray: '696969',
         dimgrey: '696969',
         dodgerblue: '1e90ff',
         firebrick: 'b22222',
         floralwhite: 'fffaf0',
         forestgreen: '228b22',
         fuchsia: 'ff00ff',
         gainsboro: 'dcdcdc',
         ghostwhite: 'f8f8ff',
         gold: 'ffd700',
         goldenrod: 'daa520',
         gray: '808080',
         grey: '808080',
         green: '008000',
         greenyellow: 'adff2f',
         honeydew: 'f0fff0',
         hotpink: 'ff69b4',
         indianred: 'cd5c5c',
         indigo: '4b0082',
         ivory: 'fffff0',
         khaki: 'f0e68c',
         lavender: 'e6e6fa',
         lavenderblush: 'fff0f5',
         lawngreen: '7cfc00',
         lemonchiffon: 'fffacd',
         lightblue: 'add8e6',
         lightcoral: 'f08080',
         lightcyan: 'e0ffff',
         lightgoldenrodyellow: 'fafad2',
         lightgray: 'd3d3d3',
         lightgrey: 'd3d3d3',
         lightgreen: '90ee90',
         lightpink: 'ffb6c1',
         lightsalmon: 'ffa07a',
         lightseagreen: '20b2aa',
         lightskyblue: '87cefa',
         lightslategray: '778899',
         lightslategrey: '778899',
         lightsteelblue: 'b0c4de',
         lightyellow: 'ffffe0',
         lime: '00ff00',
         limegreen: '32cd32',
         linen: 'faf0e6',
         magenta: 'ff00ff',
         maroon: '800000',
         mediumaquamarine: '66cdaa',
         mediumblue: '0000cd',
         mediumorchid: 'ba55d3',
         mediumpurple: '9370db',
         mediumseagreen: '3cb371',
         mediumslateblue: '7b68ee',
         mediumspringgreen: '00fa9a',
         mediumturquoise: '48d1cc',
         mediumvioletred: 'c71585',
         midnightblue: '191970',
         mintcream: 'f5fffa',
         mistyrose: 'ffe4e1',
         moccasin: 'ffe4b5',
         navajowhite: 'ffdead',
         navy: '000080',
         oldlace: 'fdf5e6',
         olive: '808000',
         olivedrab: '6b8e23',
         orange: 'ffa500',
         orangered: 'ff4500',
         orchid: 'da70d6',
         palegoldenrod: 'eee8aa',
         palegreen: '98fb98',
         paleturquoise: 'afeeee',
         palevioletred: 'db7093',
         papayawhip: 'ffefd5',
         peachpuff: 'ffdab9',
         peru: 'cd853f',
         pink: 'ffc0cb',
         plum: 'dda0dd',
         powderblue: 'b0e0e6',
         purple: '800080',
         red: 'ff0000',
         rosybrown: 'bc8f8f',
         royalblue: '4169e1',
         saddlebrown: '8b4513',
         salmon: 'fa8072',
         sandybrown: 'f4a460',
         seagreen: '2e8b57',
         seashell: 'fff5ee',
         sienna: 'a0522d',
         silver: 'c0c0c0',
         skyblue: '87ceeb',
         slateblue: '6a5acd',
         slategray: '708090',
         slategrey: '708090',
         snow: 'fffafa',
         springgreen: '00ff7f',
         steelblue: '4682b4',
         tan: 'd2b48c',
         teal: '008080',
         thistle: 'd8bfd8',
         tomato: 'ff6347',
         turquoise: '40e0d0',
         violet: 'ee82ee',
         wheat: 'f5deb3',
         white: 'ffffff',
         whitesmoke: 'f5f5f5',
         yellow: 'ffff00',
         yellowgreen: '9acd32',
         transparent: '00000000'
   };// colorButton colors


   // TODO: move to separate file
   var EVEColorButton = Button.extend("rootui5.eve7.controller.EVEColorButton", {
      // when default value not specified - openui tries to load custom
      renderer: {}, // ButtonRenderer.render,

      metadata: {
         properties: {
            background: 'string'
         }
      },

      onAfterRendering: function() {
         this.$().children().css("background-color", this.getBackground());
      }

   });

    var EVEColorPopup = ColorPalettePopover.extend("rootui5.eve7.controller.EVEColorPopup", {
        // when default value not specified - openui tries to load custom
        defaultColors : ['gold','darkorange', 'indianred','rgb(102,51,0)', 'cyan',// 'magenta'
                             'blue', 'lime', 'gray','slategray','rgb(204, 198, 170)',
                             'white', 'black','red' , 'rgb(102,154,51)', 'rgb(200, 0, 200)'],

        parseRGB : function(val) {
            let rgb, regex = /rgb\((\d+)\,\s?(\d+)\,\s?(\d+)\)/,
                found = val.match(regex);
            if (found) {
                console.log("match color ", found);
                rgb = { r: found[1], g: found[2], b: found[3] };
            } else {
                let hex = UI5PopupColors[val];

                // Expand shorthand form (e.g. "03F") to full form (e.g. "0033FF")
                let shorthandRegex = /^#?([a-f\d])([a-f\d])([a-f\d])$/i;

                hex = hex.replace(shorthandRegex, function(m, r, g, b) {
                    return r + r + g + g + b + b;
                });

                rgb = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);

                rgb = rgb ? { r: parseInt(rgb[1], 16), g: parseInt(rgb[2], 16), b: parseInt(rgb[3], 16) } : null;
            }
            return rgb;
        }

    });

   var GedController = Controller.extend("rootui5.eve7.controller.Ged", {

      onInit : function() {
         this.oModel = new JSONModel({ title: "GED title", "widgetlist" : [] });
         this.getView().setModel(this.oModel, "ged");

         this.ged_visible = false;
         this.ged_id = -1;
      },

      onExit : function() {
      },

      setManager: function(mgr) {
         this.mgr = mgr;
      },

      isGedVisible : function() {
         return this.ged_visible;
      },

      closeGedEditor: function() {
         if (this.ged_visible) {
            var prnt = this.getView().getParent();
            if (prnt) prnt.removeContentArea(this.getView());
            this.ged_visible = false;
         }

         this.ged_id = -1;
         this.editorElement = null;
      },

      showGedEditor: function(sumSplitter, elementId) {

         if (this.ged_visible && (elementId == this.ged_id))
            return this.closeGedEditor();

         var editorElement = this.mgr ? this.mgr.GetElement(elementId) : null;
         if (!editorElement)
            return this.closeGedEditor();

         if (!this.ged_visible)
            sumSplitter.addContentArea(this.getView());

         this.ged_id = elementId;
         this.ged_visible = true;

         this.editorElement = editorElement;

         var title = this.editorElement.fName + " (" +  this.editorElement._typename.substring(20) + " )" ;
         this.oModel.setProperty("/title", title);
         this.buildEditor();
      },

      buildEditor: function() {
         let gedFrame =  this.getView().byId("GED");
         gedFrame.unbindElement();
         gedFrame.destroyContent();
         this.secSelectList = 0;

         let t = this.editorElement._typename;
         if (t.indexOf("ROOT::Experimental::")==0) t = t.substring(20);
         let fn = "build" + t + "Setter";
         if (typeof this[fn] === "function")
            this[fn](this.editorElement);
         else
            this.buildREveElementSetter(this.editorElement);
      },

      buildREveElementSetter : function(el)
      {
         this.makeBoolSetter(el.fRnrSelf, "RnrSelf");
         this.makeBoolSetter(el.fRnrChildren, "RnrChildren");
         this.makeColorSetter(el.fMainColor, "MainColor");
      },

      buildREveSelectionSetter : function(el)
      {
         this.makeColorSetter(el.fVisibleEdgeColor, "VisibleEdgeColor");
         this.makeColorSetter(el.fHiddenEdgeColor, "HiddenEdgeColor");
      },

      buildREveJetConeSetter : function(el)
      {
         this.makeBoolSetter(el.fRnrSelf, "RnrSelf");
         this.makeBoolSetter(el.fRnrChildren, "RnrChildren");
         this.makeColorSetter(el.fMainColor, "MainColor");
         this.makeNumberSetter(el.fNDiv, "NDiv");
      },

      buildREveTrackSetter : function(el)
      {
         this.buildREveElementSetter(el);
         this.makeNumberSetter(el.fLineWidth, "LineWidth");
      },

      buildREveDataCollectionSetter : function(el)
      {
         this.makeBoolSetter(el.fRnrSelf, "RnrSelf");
         this.makeColorSetter(el.fMainColor, "MainColor");
         this.makeStringSetter(el.fFilterExpr, "FilterExpr");
      },
      buildREveDataItemListSetter : function(el)
      {
         let pthis = this;
         let gedFrame =  this.getView().byId("GED");
         let list = new sap.m.List({});
         this.secSelectList = list;

         list.addStyleClass("eveSummaryItem");
	 list.setMode("MultiSelect");
	 list.setIncludeItemInSelection(true);
	 list.addStyleClass("eveNoSelectionCheckBox");
         let citems = el.items;

         // CAUTION! This state is only valid for the last click event.
	 // If the next itemPress is triggered by a keyboard or touch event, it will still
	 // read this outdated ctrlKeyPressed information!!
	 // So ALL events causing itemPress must clear/set ctrlKeyPressed
	 // or ctrlKeyPressed must be reset to false after a short timeout.
	 //
	 // Also, it is not tested whether for all types of events, the direct browser
	 // event is coming BEFORE the itemPress event handler invocation!
         var ctrlKeyPressed = false;
         list.attachBrowserEvent("click", function(e) {
	    ctrlKeyPressed = e.ctrlKey;
	 });

         let lastLabel = "item_"+ (citems.length -1)
         let makeItem = function(i) {
            let iid = "item_"+ i;
            let fout = citems[i].fFiltered;
	    var item  = new sap.m.CustomListItem( iid, {type:sap.m.ListType.Active});
	    item.addStyleClass("sapUiTinyMargin");

            // item info
	    let label = new sap.m.Label({text: iid});
            label.addStyleClass("sapUiTinyMarginBeginEnd");

            // rnr self
	    let rb = new mCheckBox({
               selected: citems[i].fRnrSelf,
               text: "RnrSelf",
               select: function(oEvent)
               {
                  let value = oEvent.getSource().getSelected();
                  let mir =  "SetItemVisible( " + i + ", " + value + " )";
                  pthis.mgr.SendMIR(mir, el.fElementId, el._typename );
               }
            });

            rb.addStyleClass("sapUiTinyMarginEnd");

            let col_widget = new EVEColorButton( {
               icon : "sap-icon://palette",
               background: JSROOT.Painter.getColor(citems[i].fColor),
               press: function () {
                  let oCPPop = new EVEColorPopup( {
                     colorSelect: function(event) {
                        let rgb = this.parseRGB(event.getParameters().value);
                        let mir = "SetItemColorRGB(" + i + ", " + rgb.r + ", " + rgb.g +  ", " + rgb.b + ")";
                        pthis.mgr.SendMIR(mir, el.fElementId, el._typename );
                        // console.log("color mir -  .... ", mir);
                     }
                  });
                  oCPPop.openBy(this);
               }
            });
            col_widget.addStyleClass("sapUiTinyMarginBeginEnd");
            if (fout){
               label.addStyleClass("eveTableCellUnfiltered");
               rb.setEnabled(false);
               col_widget.setEnabled(false);
            }

            let box = new sap.m.HBox({
               items : [ label, rb, col_widget ]
            });

            item.addContent(box);
            list.addItem(item);

         };

         for (let i = 0; i < citems.length; ++i ) {
            if (!citems[i].fFiltered)
               makeItem(i);
         }
         for (let i = 0; i < citems.length; ++i ) {
            if (citems[i].fFiltered)
               makeItem(i);
         }
         list.attachItemPress(function(oEvent) {
	    let p = oEvent.getParameters("item");
	    let idx = p.listItem.sId.substring(5);
            let secIdcs = [idx];
            let is_multi = false;

            if(!ctrlKeyPressed)
	    {
	       let selected = list.getSelectedItems();
	       console.log("selected items ", selected, "idx = ",  idx);
	       for (let s = 0; s < selected.length; s++) {
		  if (selected[s].sId !=  p.listItem.sId)
		     list.setSelectedItem(selected[s], false);
	       }
	    }
            let fcall = "ProcessSelection(" + pthis.mgr.global_selection_id + `, ${is_multi}, true`;
                              fcall += ", { " + secIdcs.join(", ")  + " }";
                              fcall += ")";
                              pthis.mgr.SendMIR(fcall, el.fElementId, el._typename);

	 });
         gedFrame.addContent(list);
      },

      buildREveCaloDataHistSetter : function(el)
      {
         let si = el.sliceInfos;

         for (let i = 0; i < si.length; i++)
         {
            let pthis = this;
            let col_widget = new EVEColorButton( {
               background: JSROOT.Painter.getColor(si[i].color),
               press: function () {
                  let oCPPop = new EVEColorPopup( {
                     colorSelect: function(event) {
                        let rgb = this.parseRGB(event.getParameters().value);
                        let mir = "SetSliceColor(" + i + ", TColor::GetColor(" + rgb.r + ", " + rgb.g +  ", " + rgb.b + "))";
                        console.log("color mir -  .... ", mir);
                        pthis.mgr.SendMIR(mir, pthis.editorElement.fElementId, pthis.editorElement._typename);
                     }
                  });
                  oCPPop.openBy(this);
               }
            });

            let name = si[i].name;
            let label = new mText({ text:name });
            label.addStyleClass("sapUiTinyMargin");

            let cx = new mText({ text:"Color:"});
            cx.addStyleClass("sapUiTinyMargin");

            let fx = new mText({ text:"Threshold:"});
            fx.addStyleClass("sapUiTinyMargin");

            let in_widget = new mStepInput({
               displayValuePrecision: 3,
               min : 0,
               value:  si[i].threshold,
               step : 0.1,
               change: function (event)
               {
                  let mir =  "SetSliceThreshold( " + i + ", " + event.getParameter("value") + " )";
                  pthis.mgr.SendMIR(mir, pthis.editorElement.fElementId, pthis.editorElement._typename );
               }
            });
            in_widget.setWidth("100px");
            let frame = new HorizontalLayout({
               content : [label, cx, col_widget, fx, in_widget ]
            });
            let gedFrame =  this.getView().byId("GED");
            gedFrame.addContent(frame);
         }
      },

      makeBoolSetter : function(val, labelName, funcName, gedFrame)
      {
         if (!gedFrame)
            gedFrame =  this.getView().byId("GED");


         if (!funcName)
            funcName = "Set" + labelName;

         let gcm = this;
         let widget = new mCheckBox({
            selected: val,

            select: function(oEvent)
            {
               console.log("Bool setter select event", oEvent.getSource());
               let value = oEvent.getSource().getSelected();
               let mir =  funcName + "( " + value + " )";
               gcm.mgr.SendMIR(mir, gcm.editorElement.fElementId, gcm.editorElement._typename );
            }
         });

         let label = new mText({ text: labelName });
         label.addStyleClass("sapUiTinyMargin");

         let frame = new HorizontalLayout({
            content : [widget, label]
         });

         gedFrame.addContent(frame);
      },

      makeColorSetter : function(val, labelName, funcName, gedFrame)
      {
         if (!gedFrame)
            gedFrame =  this.getView().byId("GED");

         if (!funcName)
            funcName = "Set" + labelName + "RGB";


         let pthis = this;
         let widget = new EVEColorButton( {
            icon: "sap-icon://palette",
            background: JSROOT.Painter.getColor(val),
            press: function () {
               let oCPPop = new EVEColorPopup( {
                  colors: this.defaultColors,
                  colorSelect: function(event) {
                     let rgb = this.parseRGB(event.getParameters().value);
                     let mir =  funcName + "((UChar_t)" + rgb.r + ", (UChar_t)" + rgb.g +  ", (UChar_t)" + rgb.b + ")";
                     pthis.mgr.SendMIR(mir, pthis.editorElement.fElementId, pthis.editorElement._typename);
                  }
               });
               oCPPop.openBy(this);
            }
         });

         let label = new mText({ text: labelName });
         label.addStyleClass("sapUiTinyMargin");

         let frame = new HorizontalLayout({
            content : [widget, label]
         });
         gedFrame.addContent(frame);
      },

      makeNumberSetter : function(val, labelName, funcName, gedFrame)
      {
         if (!gedFrame)
            gedFrame =  this.getView().byId("GED");

         if (!funcName)
            funcName = "Set" + labelName;

         let gcm = this;
         let widget = new mInput({
            value: val,
            change: function (event)
            {
               let value = event.getParameter("value");
               let mir =  funcName + "( " + value + " )";
               gcm.mgr.SendMIR(mir, gcm.editorElement.fElementId, gcm.editorElement._typename );
            }
         });
         widget.setType(sap.m.InputType.Number);
         let label = new mText({ text: labelName });
         label.addStyleClass("sapUiTinyMargin");

         let frame = new HorizontalLayout({
            content : [widget, label]
         });
         gedFrame.addContent(frame);
      },

      makeStringSetter : function(val, labelName, funcName, gedFrame)
      {
         if (!gedFrame)
            gedFrame =  this.getView().byId("GED");

         if (!funcName)
            funcName = "Set" + labelName;

         let gcm = this;
         let widget = new mInput({
            value: val,
            change: function (event)
            {
               let value = event.getParameter("value");
               let mir =  funcName + "( \"" + value + "\" )";
               gcm.mgr.SendMIR(mir, gcm.editorElement.fElementId, gcm.editorElement._typename );
            }
         });
         widget.setType(sap.m.InputType.String);
         widget.setWidth("250px"); // AMT this should be handled differently

         let label = new mText({ text: labelName });
         label.addStyleClass("sapUiTinyMargin");

         let frame = new HorizontalLayout({
            content : [widget, label]
         });
         gedFrame.addContent(frame);
      },

      updateGED: function(elementId) {
         if (this.ged_visible && this.editorElement && (this.editorElement.fElementId == elementId)) {
            this.buildEditor();
         }
      },

      updateSecondarySelectionGED:function(elementId, sec_idcs) {
         if (this.secSelectList)
         {
            if (this.editorElement.fElementId == elementId){
               let selected = this.secSelectList.getSelectedItems();
               for (let s = 0; s < selected.length; s++)
                  this.secSelectList.setSelectedItem(selected[s], false);


               for (let i =0; i < sec_idcs.length; ++i) {
                  let sid = "item_"+sec_idcs[i];
                  this.secSelectList.setSelectedItemById(sid, true);
               }
            }
            else
               this.secSelectList.removeSelections();
         }
   }

   });
   GedController.canEditClass = function(typename) {
      return true;
   };

   /** Return method to toggle rendering self */
   GedController.GetRnrSelfMethod = function(typename) {
      var desc = this.canEditClass(typename);
      if (desc)
         for (var k=0;k<desc.length;++k)
            if ((desc[k].member == "fRnrSelf") && desc[k].name)
               return "Set" + desc[k].name;

      return "SetRnrSelf";
   }


   return GedController;

});
