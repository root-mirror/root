sap.ui.define([
   'sap/ui/core/mvc/Controller',
   'sap/ui/model/json/JSONModel',
   'sap/ui/core/Fragment',
   'sap/ui/core/HTML'
], function (Controller, JSONModel, Fragment, HTML) {
   "use strict";

   return Controller.extend("rootui5.canv.controller.Ged", {

      currentPainter: null,

      gedFragments : [],

      onInit : function() {
         console.log('init GED editor');
         var model = new JSONModel({ SelectedClass: "none" });
         this.getView().setModel(model);
      },

      onExit : function() {
         console.log('exit GED editor');
         this.cleanupGed();

      },

      cleanupGed : function() {
         console.log('Clenup GED editor');

         // empty fragments
         this.getView().byId("ged_page").removeAllContent();

         // set dummy model
         this.getView().setModel(new JSONModel({ SelectedClass: "none" }));

         // remove references
         this.currentPainter = null;
         this.currentPadPainter = null;
         this.currentPlace = undefined;

         // TODO: deregsiter for all events

      },

      addFragment : function(page, kind, model) {
         let fragm = this.gedFragments[kind];

         if (!fragm)
            return Fragment.load({
               name: "rootui5.canv.view." + kind,
               type: "XML",
               controller: this
            }).then(function(_page, _kind, _model, oFragm) {
               this.gedFragments[_kind] = oFragm;
               this.addFragment(_page, _kind, _model);
            }.bind(this, page, kind, model));

         fragm.ged_fragment = true; // mark as ged fragment

         let html = new HTML();
         html.setContent("<hr>");
         html.setTooltip(kind);
         page.addContent(html);

         fragm.setModel(model);
         page.addContent(fragm);
      },

     /** @summary Produce exec string for WebCanas to set color value
       * @desc Color can be id or string, but should belong to list of known colors
       * For higher color numbers TColor::GetColor(r,g,b) will be invoked to ensure color will be created on server side
       * @private */
      getColorExec: function(painter, col, method) {
         let id = -1, arr = JSROOT.Painter.root_colors;
         if (typeof col == "string") {
            if (!col || (col == "none")) id = 0; else
               for (let k = 1; k < arr.length; ++k)
                  if (arr[k] == col) { id = k; break; }
            if ((id < 0) && (col.indexOf("rgb") == 0)) id = 9999;
         } else if (!isNaN(col) && arr[col]) {
            id = col;
            col = arr[id];
         }

         if (id < 0) return "";

         if (id >= 50) {
            // for higher color numbers ensure that such color exists
            let c = d3.color(col);
            id = "TColor::GetColor(" + c.r + "," + c.g + "," + c.b + ")";
         }

         return "exec:" + method + "(" + id + ")";
      },

      /// function called when user changes model property
      /// data object includes _kind, _painter and _handle (optionally)
      modelPropertyChange : function(evnt, data) {
         let pars = evnt.getParameters();
         console.log('Model property changes', pars.path, pars.value, data._kind);

         if (data._handle) {
            //var subname = pars.path.substr(1);
            //if (subname in data._handle) data._handle[subname] = pars.value;

            if (typeof data._handle.verifyDirectChange === 'function')
                data._handle.verifyDirectChange(data._painter);
            data._handle.changed = true;
         }

         let exec = "", item = pars.path.substr(1), obj;

         if (data._painter) {
            obj = data._painter.snapid ? data._painter.getObject() : null;
         }

         if (obj) {
            if ((data._kind === "TAttLine") && (obj.fLineColor!==undefined) && (obj.fLineStyle!==undefined) && (obj.fLineWidth!==undefined)) {
               if (item == "attline/width")
                  exec = "exec:SetLineWidth(" + pars.value + ")";
               else if (item == "attline/style")
                  exec = "exec:SetLineStyle(" + pars.value + ")";
               else if (item == "attline/color")
                  exec = this.getColorExec(data._painter, pars.value, "SetLineColor");
            } else if ((data._kind === "TAttFill") && (obj.fFillColor!==undefined) && (obj.fFillStyle!==undefined))  {
               if (item == "attfill/pattern")
                  exec = "exec:SetFillStyle(" + pars.value + ")";
               else if (item == "attfill/color")
                  exec = this.getColorExec(data._painter, pars.value, "SetFillColor");
            }
         }

         if (data._painter)
            data._painter.interactiveRedraw("pad", exec); // TODO: some objects can readraw directly, no need to redraw pad
         else if (this.currentPadPainter)
            this.currentPadPainter.redraw();
      },


      getAxisHandle: function() {
         if (this.currentPainter)
            switch (this.currentPlace) {
               case "xaxis": return this.currentPainter.x_handle;
               case "yaxis": return this.currentPainter.y_handle;
               case "zaxis": return this.currentPainter.z_handle;
            }
         return null
      },

      setAxisModel : function(model) {
         let obj =  this.currentPainter.getObject(this.currentPlace),
             painter = this.getAxisHandle();

         let data = {
             specialRefresh: 'setAxisModel',
             axis: obj,
             axiscolor: painter.lineatt.color,
             color_label: this.currentPadPainter.getColor(obj.fLabelColor),
             center_label: obj.TestBit(JSROOT.EAxisBits.kCenterLabels),
             vert_label: obj.TestBit(JSROOT.EAxisBits.kLabelsVert),
             color_title: this.currentPadPainter.getColor(obj.fTitleColor),
             center_title: obj.TestBit(JSROOT.EAxisBits.kCenterTitle),
             rotate_title: obj.TestBit(JSROOT.EAxisBits.kRotateTitle),
         };

         model.setData(data);
      },

      processAxisModelChange: function(evnt, data) {
         let pars = evnt.getParameters(),
             item = pars.path.substr(1),
             exec = "",
             painter = this.currentPainter,
             kind = this.currentPlace,
             axis = painter.getObject(kind);

         // while axis painter is temporary object, we should not try change it attributes

         if (!this.currentPadPainter || !axis) return;

         if ((typeof kind == 'string') && (kind.indexOf("axis")==1)) kind = kind.substr(0,1);

         console.log(`Change axis ${kind} item ${item} value ${pars.value}`);

         switch(item) {
            case "axis/fTitle":
               exec = `exec:SetTitle("${pars.value}")`;
               break;
            case "axiscolor":
               axis.fAxisColor = this.currentPadPainter.addColor(pars.value);
               exec = this.getColorExec(painter, pars.value, "SetAxisColor");
               break;
            case "color_label":
               axis.fLabelColor = this.currentPadPainter.addColor(pars.value);
               exec = this.getColorExec(painter, pars.value, "SetLabelColor");
               break;
            case "center_label":
               axis.InvertBit(JSROOT.EAxisBits.kCenterLabels);
               exec = `exec:CenterLabels(${pars.value ? true : false})`;
               break;
            case "vert_label":
               axis.InvertBit(JSROOT.EAxisBits.kLabelsVert);
               exec = `exec:SetBit(TAxis::kLabelsVert,${pars.value ? true : false})`;
               break;
            case "axis/fLabelOffset":
               exec = `exec:SetLabelOffset(${pars.value})`;
               break;
            case "axis/fLabelSize":
               exec = `exec:SetLabelOffset(${pars.value})`;
               break;
            case "color_title":
               axis.fLabelColor = this.currentPadPainter.addColor(pars.value);
               exec = this.getColorExec(painter, pars.value, "SetTitleColor");
               break;
            case "center_title":
               axis.InvertBit(JSROOT.EAxisBits.kCenterTitle);
               exec = `exec:CenterTitle(${pars.value ? true : false})`;
               break;
            case "rotate_title":
               axis.InvertBit(JSROOT.EAxisBits.kRotateTitle);
               exec = `exec:RotateTitle(${pars.value ? true : false})`;
               break;
            case "axis/fTickLength":
               exec = `exec:SetTickLength(${pars.value})`;
               break;
            case "axis/fTitleOffset":
               exec = `exec:SetTitleOffset(${pars.value})`;
               break;
            case "axis/fTitleSize":
               exec = `exec:SetTitleSize(${pars.value})`;
               break;
         }

         // TAxis belongs to main painter like TH1, therefore submit commands there
         let main = this.currentPainter.getMainPainter(true);

         if (main && main.snapid) {
            console.log('Invoke interactive redraw ', main.snapid, kind)
            main.interactiveRedraw("pad", exec, kind);
         } else {
            this.currentPadPainter.redraw();
         }
      },

      setRAxisModel: function(model) {
         let handle = this.getAxisHandle();
         if (!handle) return;

         let data = {
            logbase: handle.logbase || 0,
            handle: handle,
            ticks_size: handle.ticksSize/handle.scaling_size,
            labels_offset: handle.labelsOffset/handle.scaling_size,
            labels_rotate: handle.labelsFont.angle != 0,
            title_offset: handle.titleOffset/handle.scaling_size,
            title_rotate: handle.isTitleRotated(),
            specialRefresh: 'setRAxisModel'
         };

         if (Math.abs(data.logbase - 2.7) < 0.1) data.logbase = 3; // for ln

         model.setData(data);
      },

      processRAxisModelChange: function(evnt, data) {
         let pars = evnt.getParameters(),
             item = pars.path.substr(1);

         let handle = this.getAxisHandle();
         if (!handle) return;

         switch(item) {
            case "logbase": handle.changeAxisLog((pars.value == 3) ? Math.exp(1) : pars.value); break;
            case "handle/ticksColor": handle.changeAxisAttr(1, "ticks_color_name", pars.value); break;
            case "ticks_size": handle.changeAxisAttr(1, "ticks_size", pars.value); break;
            case "handle/ticksSide": handle.changeAxisAttr(1, "ticks_side", pars.value); break;
            case "labels_offset": handle.changeAxisAttr(1, "labels_offset", pars.value); break;
            case "labels_rotate": handle.changeAxisAttr(1, "labels_angle", pars.value ? 180 : 0); break;
            case "handle/fTitle": handle.changeAxisAttr(1, "title", pars.value); break;
            case "title_offset": handle.changeAxisAttr(1, "title_offset", pars.value); break;
            case "handle/titlePos": handle.changeAxisAttr(1, "title_position", pars.value); break;
            case "title_rotate": handle.changeAxisAttr(1, "title_angle", pars.value ? 180 : 0); break;
         }
      },

      processHistModelChange : function(evnt, data) {
         // let pars = evnt.getParameters();
         let opts = data.options;

         opts.Mode3D = opts.Mode3Dindx > 0;

         opts.Lego = parseInt(opts.Lego);

         let cl = this.getView().getModel().getProperty("/SelectedClass");

         if ((typeof cl == "string") && opts.Mode3D && (cl.indexOf("ROOT::Experimental::RHist") == 0))
            opts.Lego = 12;

         opts.Contor = parseInt(opts.Contor);
         opts.ErrorKind = parseInt(opts.ErrorKind);

         if (this.currentPainter)
            this.currentPainter.interactiveRedraw("pad","drawopt");
      },


      onObjectSelect : function(padpainter, painter, place) {

         if ((this.currentPainter === painter) && (place === this.currentPlace)) return;

         this.currentPadPainter = padpainter;
         this.currentPainter = painter;
         this.currentPlace = place;

         let obj = painter.getObject(place), selectedClass = "";

         if (place == "xaxis" && painter.x_handle) {
            painter = painter.x_handle;
            selectedClass = painter.getAxisType();
         } else if (place == "yaxis" && painter.y_handle) {
            painter = painter.y_handle;
            selectedClass = painter.getAxisType();
         } else if (place == "zaxis" && painter.z_handle) {
            painter = painter.z_handle;
            selectedClass = painter.getAxisType();
         } else {
            selectedClass = obj ? obj._typename : painter.getObjectHint();
         }

         this.getView().getModel().setProperty("/SelectedClass", selectedClass);

         let oPage = this.getView().byId("ged_page");
         oPage.removeAllContent();


         if (selectedClass == "RAttrAxis") {
            let model = new JSONModel({});
            this.setRAxisModel(model);
            this.addFragment(oPage, "RAxis", model);
            model.attachPropertyChange({ _kind: "RAttrAxis" }, this.processRAxisModelChange, this);
            return;
         }

         if (painter.lineatt && painter.lineatt.used && !painter.lineatt.not_standard) {
            let model = new JSONModel( { attline: painter.lineatt } );
            model.attachPropertyChange({ _kind: "TAttLine", _painter: painter, _handle: painter.lineatt }, this.modelPropertyChange, this);

            this.addFragment(oPage, "TAttLine", model);
         }

         if (painter.fillatt && painter.fillatt.used) {
            let model = new JSONModel( { attfill: painter.fillatt } );
            model.attachPropertyChange({ _kind: "TAttFill", _painter: painter, _handle: painter.fillatt }, this.modelPropertyChange, this);

            this.addFragment(oPage, "TAttFill", model);
         }

         if (painter.markeratt && painter.markeratt.used) {
            let model = new JSONModel( { attmark: painter.markeratt } );
            model.attachPropertyChange({ _kind: "TAttMarker", _painter: painter, _handle: painter.markeratt }, this.modelPropertyChange, this);

            this.addFragment(oPage, "TAttMarker", model);
         }

         if (typeof painter.processTitleChange == 'function') {
            let tobj = painter.processTitleChange("check");
            if (tobj) {
               let model = new JSONModel({ tnamed: tobj });
               model.attachPropertyChange( {}, painter.processTitleChange, painter );
               this.addFragment(oPage, "TNamed", model);
            }
         }

         if (selectedClass == "TAxis") {
            let model = new JSONModel({});
            this.setAxisModel(model);
            this.addFragment(oPage, "Axis", model);
            model.attachPropertyChange({ _kind: "TAxis" }, this.processAxisModelChange, this);
         }

         if (typeof painter.getHisto == 'function') {

            painter.options.Mode3Dindx = painter.options.Mode3D ? 1 : 0;

            let model = new JSONModel({ opts : painter.options });

            // model.attachPropertyChange({}, painter.processTitleChange, painter);
            this.addFragment(oPage, "Hist", model);

            model.attachPropertyChange({ options: painter.options }, this.processHistModelChange, this);
         }
      },

      onObjectRedraw : function(padpainter, painter) {
         if ((this.currentPadPainter !== padpainter) || (this.currentPainter !== painter)) return;

         // console.log('GED sees selected object redraw');

         let page = this.getView().byId("ged_page");
         let cont = page.getContent();

         for (let n = 0; n < cont.length; ++n)
            if (cont[n] && cont[n].ged_fragment) {
               let model = cont[n].getModel();

               let func = model.getProperty("/specialRefresh");
               if (func)
                  this[func](model);
               else
                  model.refresh();
            }
      },

      onPadRedraw : function(padpainter) {
         if (this.currentPadPainter === padpainter)
            this.onObjectRedraw(this.currentPadPainter, this.currentPainter);
      },

      padEventsReceiver : function(evnt) {
         if (!evnt) return;

         if (evnt.what == "select")
            this.onObjectSelect(evnt.padpainter, evnt.painter, evnt.place);
         else if (evnt.what == "redraw")
            this.onObjectRedraw(evnt.padpainter, evnt.painter);
         else if (evnt.what == "padredraw")
            this.onPadRedraw(evnt.padpainter);
      }

   });

});
