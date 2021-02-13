/// @file EveScene.js

// TODO: add dependency from JSROOT components

sap.ui.define(['rootui5/eve7/lib/EveManager'], function(EveManager) {

   "use strict";

   /// constructor, handle for REveScene class

   function EveScene(mgr, scene, glctrl)
   {
      this.mgr     = mgr;
      this.scene   = scene;
      this.id      = scene.fSceneId;
      this.glctrl  = glctrl;
      this.creator = glctrl.viewer.creator;
      this.id2obj_map  = new Map; // base on element id
      this.mid2obj_map = new Map; // base on master id

      this.first_time = true;

      // register ourself for scene events
      this.mgr.RegisterSceneReceiver(scene.fSceneId, this);

      // AMT temporary solution ... resolve with callSceneReceivers in EveManager.js
      scene.eve_scene = this;
   }

   //==============================================================================
   // Render object creation / management
   //==============================================================================

   EveScene.prototype.makeGLRepresentation = function(elem)
   {
      if ( ! elem.render_data) return null;

      let fname = elem.render_data.rnr_func;
      let obj3d = this.creator[fname](elem, elem.render_data);

      if (obj3d)
      {
         // MT ??? why?, it can really be anything, even just container Object3D
         obj3d._typename = "THREE.Mesh";

         // add reference to a streamed eve element to obj3d
         obj3d.eve_el = elem;

         // SL: this is just identifier for highlight, required to show items on other places, set in creator
         obj3d.geo_object = elem.fMasterId || elem.fElementId;
         obj3d.geo_name   = elem.fName; // used for highlight
         obj3d.scene  = this; // required for get changes when highlight/selection is changed

         if (elem.render_data.matrix)
         {
            obj3d.matrixAutoUpdate = false;
            obj3d.matrix.fromArray( elem.render_data.matrix );
            obj3d.updateMatrixWorld(true);
         }

         return obj3d;
      }
   }

   EveScene.prototype.getObj3D = function(elementId, is_master)
   {
      let map = is_master ? this.mid2obj_map : this.id2obj_map;
      return map.get(elementId);
   }

   EveScene.prototype.create3DObjects = function(all_ancestor_children_visible, prnt, res3d)
   {
      if (prnt === undefined) {
         prnt = this.mgr.GetElement(this.id);
         res3d = [];
      }

      if (!prnt || !prnt.childs) return res3d;

      for (let k = 0; k < prnt.childs.length; ++k)
      {
         let elem = prnt.childs[k];
         if (elem.render_data)
         {
            let fname = elem.render_data.rnr_func, obj3d = null;
            if ( ! this.creator[fname])
            {
               console.error("Function " + fname + " missing in creator");
            }
            else
            {
               let obj3d = this.makeGLRepresentation(elem);
               if (obj3d)
               {
                  // MT - should maintain hierarchy ????
                  // Easier to remove ... but might need sub-class of
                  // Object3D to separate "graphical" children and structural children.

                  res3d.push(obj3d);

                  this.id2obj_map.set(elem.fElementId, obj3d);
                  if (elem.fMasterId) this.mid2obj_map.set(elem.fMasterId, obj3d);

                  obj3d.visible = elem.fRnrSelf && all_ancestor_children_visible;
                  obj3d.all_ancestor_children_visible = all_ancestor_children_visible;
               }
            }
         }

         this.create3DObjects(elem.fRnrChildren && all_ancestor_children_visible, elem, res3d);
      }

      return res3d;
   }

   //==============================================================================

   //==============================================================================

   /** method insert all objects into three.js container */
   EveScene.prototype.redrawScene = function()
   {
      if ( ! this.glctrl) return;

      let res3d = this.create3DObjects(true);
      if ( ! res3d.length && this.first_time) return;

      let cont = this.glctrl.getSceneContainer("scene" + this.id);
      while (cont.children.length > 0)
         cont.remove(cont.children[0]);

      for (let k = 0; k < res3d.length; ++k)
         cont.add(res3d[k]);

      this.applySelectionOnSceneCreate(this.mgr.global_selection_id);
      this.applySelectionOnSceneCreate(this.mgr.global_highlight_id);

      this.first_time = false;
   }

   EveScene.prototype.update3DObjectsVisibility = function(arr, all_ancestor_children_visible)
   {
      if (!arr) return;

      for (let k = 0; k < arr.length; ++k)
      {
         let elem = arr[k];
         if (elem.render_data)
         {
            let obj3d = this.getObj3D(elem.fElementId);
            if (obj3d)
            {
               obj3d.visible = elem.fRnrSelf && all_ancestor_children_visible;
               obj3d.all_ancestor_children_visible = all_ancestor_children_visible;
            }
         }

         this.update3DObjectsVisibility(elem.childs, elem.fRnrChildren && all_ancestor_children_visible);
      }
   }

   EveScene.prototype.onSceneCreate = function(id)
   {
      this.redrawScene();
   }

   //==============================================================================
   // Scene changes processing
   //==============================================================================

   EveScene.prototype.beginChanges = function()
   {
   }

   EveScene.prototype.endChanges = function()
   {
      if (this.glctrl)
         this.glctrl.viewer.render();
   }

   EveScene.prototype.elementAdded = function(el)
   {
      if ( ! this.glctrl) return;

      let obj3d =  this.makeGLRepresentation(el);
      if ( ! obj3d) return;

      // AMT this is an overkill, temporary solution
      let scene = this.mgr.GetElement(el.fSceneId);
      this.update3DObjectsVisibility(scene.childs, true);

      let container = this.glctrl.getSceneContainer("scene" + this.id);

      container.add(obj3d);

      this.id2obj_map.set(el.fElementId, obj3d);
      if (el.fMasterId) this.mid2obj_map.set(el.fMasterId, obj3d);
   }

   EveScene.prototype.replaceElement = function(el)
   {
      if ( ! this.glctrl) return;

      let obj3d = this.getObj3D(el.fElementId);
      let all_ancestor_children_visible = obj3d.all_ancestor_children_visible;
      let visible = obj3d.visible;

      let container = this.glctrl.getSceneContainer("scene" + this.id);

      container.remove(obj3d);

      obj3d = this.makeGLRepresentation(el);
      obj3d.all_ancestor_children_visible = all_ancestor_children_visible;
      obj3d.visible = visible;
      container.add(obj3d);


      this.id2obj_map.set(el.fElementId, obj3d);
      if (el.fMasterId) this.mid2obj_map.set(el.fMasterId, obj3d);

      this.glctrl.viewer.render();
   }

   EveScene.prototype.elementRemoved = function()
   {
      // XXXXX how is this empty? not called?
   }

   EveScene.prototype.elementsRemoved = function(ids)
   {
      for (let i = 0; i < ids.length; i++)
      {
         let elId  = ids[i];
         let obj3d = this.getObj3D(elId);
         if (!obj3d) {
            let el = this.mgr.GetElement(elId);
            if (el && el.render_data) {
               console.log("ERROR EveScene.prototype.elementsRemoved can't find obj3d ",this.mgr.GetElement(el));
            }
            continue;
         }

         let container = this.glctrl.getSceneContainer("scene" + this.id);
         container.remove(obj3d);

         this.id2obj_map.delete(elId);
      }
   }

   EveScene.prototype.sceneElementChange = function(msg)
   {
      let el = this.mgr.GetElement(msg.fElementId);

      // visibility
      if (msg.changeBit & this.mgr.EChangeBits.kCBVisibility) {
         // self
         if (msg.rnr_self_changed)
         {
            let obj3d = this.getObj3D( el.fElementId );
            if (obj3d)
            {
               obj3d.visible = obj3d.all_ancestor_children_visible && el.fRnrSelf;
            }
         }
         // children
         if (msg.rnr_children_changed && el.childs)
         {
            let scene = this.mgr.GetElement(el.fSceneId);
            this.update3DObjectsVisibility(scene.childs, true);
         }
      }

      // other change bits
      if (el.render_data) {
         if ((el.changeBit & this.mgr.EChangeBits.kCBObjProps) || (el.changeBit & this.mgr.EChangeBits.kCBColorSelection))
         {
            this.replaceElement(el);
         }
      }
   }

   //==============================================================================
   // Selection handling
   //==============================================================================

   /** interactive handler. Calculates selection state, apply to element and distribute to other scene */
   EveScene.prototype.processElementSelected = function(obj3d, indx, event)
   {
      // MT BEGIN
      // console.log("EveScene.prototype.processElementSelected", obj3d, col, indx, evnt);

      let is_multi  = event && event.ctrlKey ? true : false;
      let is_secsel = indx !== undefined;

      let fcall = "NewElementPicked(" + (obj3d ? obj3d.eve_el.fElementId : 0) + `, ${is_multi}, ${is_secsel}`;
      if (is_secsel)
      {
         fcall += ", { " + (Array.isArray(indx) ? indx.join(", ") : indx) + " }";
      }
      fcall += ")";

      this.mgr.SendMIR(fcall, this.mgr.global_selection_id, "ROOT::Experimental::REveSelection");

      return true;
   }

   /** interactive handler */
   EveScene.prototype.processElementHighlighted = function(obj3d, indx, evnt)
   {
      if (this.mgr.MatchSelection(this.mgr.global_selection_id, obj3d.eve_el, indx))
         return true;


      // Need check for duplicates before call server, else server will un-higlight highlighted element
      // console.log("EveScene.prototype.processElementHighlighted", obj3d.eve_el.fElementId, indx, evnt);
      if (this.mgr.MatchSelection(this.mgr.global_highlight_id, obj3d.eve_el, indx))
         return true;

      // when send queue below threshold, ignre highlight
      if (this.mgr.CheckSendThreshold())
         return true;

      let is_multi  = false;
      let is_secsel = indx !== undefined;

      let fcall = "NewElementPicked(" + obj3d.eve_el.fElementId + `, ${is_multi}, ${is_secsel}`;
      if (is_secsel)
         fcall += ", { " + (Array.isArray(indx) ? indx.join(", ") : indx) + " }";
      fcall += ")";

      this.mgr.SendMIR(fcall, this.mgr.global_highlight_id, "ROOT::Experimental::REveSelection");

      return true;
   }


   EveScene.prototype.clearHighlight = function()
   {
      // QQQQ This will have to change for multi client support.
      // Highlight will always be multi and we will have to track
      // which highlight is due to our connection.

      // when send queue below threshold, ignre highlight
      if (this.mgr.CheckSendThreshold())
         return true;

      let is_multi  = false;
      let is_secsel = false;

      let so = this.mgr.GetElement(this.mgr.global_highlight_id);

      if (so && so.prev_sel_list && so.prev_sel_list.length)
      {
         let fcall = "NewElementPicked(" + 0 + `, ${is_multi}, ${is_secsel}` + ")";

         this.mgr.SendMIR(fcall, this.mgr.global_highlight_id, "ROOT::Experimental::REveSelection");
      }

      return true;
   }

   EveScene.prototype.applySelectionOnSceneCreate = function(selection_id)
   {
      let selection_obj = this.mgr.GetElement(selection_id);
      if ( ! selection_obj || ! selection_obj.prev_sel_list) return;

      var pthis = this;
      selection_obj.prev_sel_list.forEach(function(rec) {

         let prl = pthis.mgr.GetElement(rec.primary);
         if (prl && prl.fSceneId == pthis.id)
         {
            pthis.SelectElement(selection_obj, rec.primary, rec.sec_idcs, rec.extra );
         }
         else // XXXXX why else ... should we not process all of them?!!!!
         {
            for (let impId of rec.implied)
            {
               let eli = pthis.mgr.GetElement(impId);
               if (eli && eli.fSceneId == pthis.id)
               {
                  // console.log("CHECK select IMPLIED", pthis);
                  pthis.SelectElement(selection_obj, impId, rec.sec_idcs, rec.extra);
               }
            }
         }
      });
   }

   EveScene.prototype.SelectElement = function(selection_obj, element_id, sec_idcs, extra)
   {
      let obj3d = this.getObj3D( element_id );
      if (!obj3d) return;

      let opass = this.glctrl.viewer.outline_pass;
      opass.id2obj_map[element_id] = opass.id2obj_map[element_id] || [];

      if (opass.id2obj_map[element_id][selection_obj.fElementId] !== undefined)
      {
         return;
      }

      let stype  = selection_obj.fName.endsWith("Selection") ? "select" : "highlight";
      let estype = THREE.OutlinePass.selection_enum[stype];
      let oe = this.mgr.GetElement(element_id);
      // console.log("EveScene.SelectElement ", selection_obj.fName, oe.fName, selection_obj.fElementId, this.glctrl.viewer.outline_pass.id2obj_map);

      let res = {
         "sel_type" : estype,
         "sec_sel"  : (oe.fSecondarySelect && sec_idcs.length > 0) ? true: false,
         "geom"     : []
      };

      // exit if you try to highlight an object that has already been selected
      if (estype == THREE.OutlinePass.selection_enum["highlight"] &&
          opass.id2obj_map[element_id][this.mgr.global_selection_id] !== undefined)
      {
         if (!res.sec_sel)
         return;
      }

      if (!res.sec_sel) opass.id2obj_map[element_id] = [];

      if (obj3d.get_ctrl)
      {
         let ctrl = obj3d.get_ctrl();
         ctrl.DrawForSelection(sec_idcs, res, extra);
         opass.id2obj_map[element_id][selection_obj.fElementId] = res;

         if (stype == "highlight" && selection_obj.sel_list) {
            this.glctrl.viewer.remoteToolTip(selection_obj.sel_list[0].tooltip);
         }
      }
   }

   EveScene.prototype.UnselectElement = function(selection_obj, element_id)
   {
      let opass = this.glctrl.viewer.outline_pass;
      // console.log("EveScene.UnselectElement ", selection_obj.fName, element_id, selection_obj.fElementId, this.glctrl.viewer.outline_pass.id2obj_map);
      if (opass.id2obj_map[element_id] !== undefined)
      {
         delete opass.id2obj_map[element_id][selection_obj.fElementId];
      }
   }

   return EveScene;
});
