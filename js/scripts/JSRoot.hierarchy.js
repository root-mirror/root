/// @file JSRoot.hierarchy.js
/// Hierarchy display functionality

JSROOT.define(['d3', 'painter'], (d3, jsrp) => {

   "use strict";

   // ===========================================================================================

   /** @summary draw list content
     * @desc used to draw all items from TList or TObjArray inserted into the TCanvas list of primitives
     * @memberof JSROOT.Painter
     * @private */
   function drawList(divid, lst, opt) {
      if (!lst || !lst.arr) return Promise.resolve(null);

      let obj = {
        divid: divid,
        lst: lst,
        opt: opt,
        indx: -1,
        painter: null,
        draw_next: function() {
           while (++this.indx < this.lst.arr.length) {
              let item = this.lst.arr[this.indx],
                  opt = (this.lst.opt && this.lst.opt[this.indx]) ? this.lst.opt[this.indx] : this.opt;
              if (!item) continue;
              return JSROOT.draw(this.getDom(), item, opt).then(p => {
                 if (p && !this.painter) this.painter = p;
                 return this.draw_next(); // reenter loop
              });
           }

           return Promise.resolve(this.painter);
        }
      }

      return obj.draw_next();
   }

   // ===================== hierarchy scanning functions ==================================

   /** @summary Create hierarchy elements for TFolder object
     * @memberof JSROOT.Painter
     * @private */
   function folderHierarchy(item, obj) {

      if (!obj || !('fFolders' in obj) || (obj.fFolders===null)) return false;

      if (obj.fFolders.arr.length===0) { item._more = false; return true; }

      item._childs = [];

      for ( let i = 0; i < obj.fFolders.arr.length; ++i) {
         let chld = obj.fFolders.arr[i];
         item._childs.push( {
            _name : chld.fName,
            _kind : "ROOT." + chld._typename,
            _obj : chld
         });
      }
      return true;
   }

   /** @summary Create hierarchy elements for TTask object
     * @memberof JSROOT.Painter
     * @private */
   function taskHierarchy(item, obj) {
      // function can be used for different derived classes
      // we show not only child tasks, but all complex data members

      if (!obj || !('fTasks' in obj) || (obj.fTasks === null)) return false;

      objectHierarchy(item, obj, { exclude: ['fTasks', 'fName'] } );

      if ((obj.fTasks.arr.length===0) && (item._childs.length==0)) { item._more = false; return true; }

      // item._childs = [];

      for ( let i = 0; i < obj.fTasks.arr.length; ++i) {
         let chld = obj.fTasks.arr[i];
         item._childs.push( {
            _name : chld.fName,
            _kind : "ROOT." + chld._typename,
            _obj : chld
         });
      }
      return true;
   }

   /** @summary Create hierarchy elements for TList object
     * @memberof JSROOT.Painter
     * @private */
   function listHierarchy(folder, lst) {
      if (!JSROOT.isRootCollection(lst)) return false;

      if ((lst.arr === undefined) || (lst.arr.length === 0)) {
         folder._more = false;
         return true;
      }

      let do_context = false, prnt = folder;
      while (prnt) {
         if (prnt._do_context) do_context = true;
         prnt = prnt._parent;
      }

      // if list has objects with similar names, create cycle number for them
      let ismap = (lst._typename == 'TMap'), names = [], cnt = [], cycle = [];

      for (let i = 0; i < lst.arr.length; ++i) {
         let obj = ismap ? lst.arr[i].first : lst.arr[i];
         if (!obj) continue; // for such objects index will be used as name
         let objname = obj.fName || obj.name;
         if (!objname) continue;
         let indx = names.indexOf(objname);
         if (indx>=0) {
            cnt[indx]++;
         } else {
            cnt[names.length] = cycle[names.length] = 1;
            names.push(objname);
         }
      }

      folder._childs = [];
      for ( let i = 0; i < lst.arr.length; ++i) {
         let obj = ismap ? lst.arr[i].first : lst.arr[i];

         let item;

         if (!obj || !obj._typename) {
            item = {
               _name: i.toString(),
               _kind: "ROOT.NULL",
               _title: "NULL",
               _value: "null",
               _obj: null
            }
         } else {
           item = {
             _name: obj.fName || obj.name,
             _kind: "ROOT." + obj._typename,
             _title: (obj.fTitle || "") + " type:"  +  obj._typename,
             _obj: obj
           };

           switch(obj._typename) {
              case 'TColor': item._value = jsrp.getRGBfromTColor(obj); break;
              case 'TText': item._value = obj.fTitle; break;
              case 'TLatex': item._value = obj.fTitle; break;
              case 'TObjString': item._value = obj.fString; break;
              default: if (lst.opt && lst.opt[i] && lst.opt[i].length) item._value = lst.opt[i];
           }

           if (do_context && jsrp.canDraw(obj._typename)) item._direct_context = true;

           // if name is integer value, it should match array index
           if (!item._name || (!isNaN(parseInt(item._name)) && (parseInt(item._name)!==i))
               || (lst.arr.indexOf(obj)<i)) {
              item._name = i.toString();
           } else {
              // if there are several such names, add cycle number to the item name
              let indx = names.indexOf(obj.fName);
              if ((indx>=0) && (cnt[indx]>1)) {
                 item._cycle = cycle[indx]++;
                 item._keyname = item._name;
                 item._name = item._keyname + ";" + item._cycle;
              }
           }
         }

         folder._childs.push(item);
      }
      return true;
   }

   /** @summary Create hierarchy of TKey lists in file or sub-directory
     * @memberof JSROOT.Painter
     * @private */
   function keysHierarchy(folder, keys, file, dirname) {

      if (keys === undefined) return false;

      folder._childs = [];

      for (let i = 0; i < keys.length; ++i) {
         let key = keys[i];

         let item = {
            _name : key.fName + ";" + key.fCycle,
            _cycle : key.fCycle,
            _kind : "ROOT." + key.fClassName,
            _title : key.fTitle,
            _keyname : key.fName,
            _readobj : null,
            _parent : folder
         };

         if (key.fObjlen > 1e5) item._title += ' (size: ' + (key.fObjlen/1e6).toFixed(1) + 'MB)';

         if ('fRealName' in key)
            item._realname = key.fRealName + ";" + key.fCycle;

         if (key.fClassName == 'TDirectory' || key.fClassName == 'TDirectoryFile') {
            let dir = null;
            if (dirname && file) dir = file.getDir(dirname + key.fName);
            if (!dir) {
               item._more = true;
               item._expand = function(node, obj) {
                  // one can get expand call from child objects - ignore them
                  return keysHierarchy(node, obj.fKeys);
               }
            } else {
               // remove cycle number - we have already directory
               item._name = key.fName;
               keysHierarchy(item, dir.fKeys, file, dirname + key.fName + "/");
            }
         } else
         if ((key.fClassName == 'TList') && (key.fName == 'StreamerInfo')) {
            item._name = 'StreamerInfo';
            item._kind = "ROOT.TStreamerInfoList";
            item._title = "List of streamer infos for binary I/O";
            item._readobj = file.fStreamerInfos;
         }

         folder._childs.push(item);
      }

      return true;
   }

   /** @summary Create hierarchy for arbitrary object
     * @memberof JSROOT.Painter
     * @private */
   function objectHierarchy(top, obj, args) {
      if (!top || (obj===null)) return false;

      top._childs = [];

      let proto = Object.prototype.toString.apply(obj);

      if (proto === '[object DataView]') {

         let item = {
             _parent: top,
             _name: 'size',
             _value: obj.byteLength.toString(),
             _vclass: 'h_value_num'
         };

         top._childs.push(item);
         let namelen = (obj.byteLength < 10) ? 1 : Math.log10(obj.byteLength);

         for (let k=0;k<obj.byteLength;++k) {
            if (k % 16 === 0) {
               item = {
                 _parent: top,
                 _name: k.toString(),
                 _value: "",
                 _vclass: 'h_value_num'
               };
               while (item._name.length < namelen) item._name = "0" + item._name;
               top._childs.push(item);
            }

            let val = obj.getUint8(k).toString(16);
            while (val.length<2) val = "0"+val;
            if (item._value.length>0)
               item._value += (k%4===0) ? " | " : " ";

            item._value += val;
         }
         return true;
      }

      // check nosimple property in all parents
      let nosimple = true, do_context = false, prnt = top;
      while (prnt) {
         if (prnt._do_context) do_context = true;
         if ('_nosimple' in prnt) { nosimple = prnt._nosimple; break; }
         prnt = prnt._parent;
      }

      let isarray = (proto.lastIndexOf('Array]') == proto.length-6) && (proto.indexOf('[object')==0) && !isNaN(obj.length),
          compress = isarray && (obj.length > JSROOT.settings.HierarchyLimit),  arrcompress = false;

      if (isarray && (top._name==="Object") && !top._parent) top._name = "Array";

      if (compress) {
         arrcompress = true;
         for (let k=0;k<obj.length;++k) {
            let typ = typeof obj[k];
            if ((typ === 'number') || (typ === 'boolean') || (typ=='string' && (obj[k].length<16))) continue;
            arrcompress = false; break;
         }
      }

      if (!('_obj' in top))
         top._obj = obj;
      else
      if (top._obj !== obj) alert('object missmatch');

      if (!top._title) {
         if (obj._typename)
            top._title = "ROOT." + obj._typename;
         else
         if (isarray) top._title = "Array len: " + obj.length;
      }

      if (arrcompress) {
         for (let k=0;k<obj.length;) {

            let nextk = Math.min(k+10,obj.length), allsame = true, prevk = k;

            while (allsame) {
               allsame = true;
               for (let d=prevk;d<nextk;++d)
                  if (obj[k]!==obj[d]) allsame = false;

               if (allsame) {
                  if (nextk===obj.length) break;
                  prevk = nextk;
                  nextk = Math.min(nextk+10,obj.length);
               } else
               if (prevk !== k) {
                  // last block with similar
                  nextk = prevk;
                  allsame = true;
                  break;
               }
            }

            let item = { _parent: top, _name: k+".."+(nextk-1), _vclass: 'h_value_num' };

            if (allsame) {
               item._value = obj[k].toString();
            } else {
               item._value = "";
               for (let d=k;d<nextk;++d)
                  item._value += ((d===k) ? "[ " : ", ") + obj[d].toString();
               item._value += " ]";
            }

            top._childs.push(item);

            k = nextk;
         }
         return true;
      }

      let lastitem, lastkey, lastfield, cnt;

      for (let key in obj) {
         if ((key == '_typename') || (key[0]=='$')) continue;
         let fld = obj[key];
         if (typeof fld == 'function') continue;
         if (args && args.exclude && (args.exclude.indexOf(key)>=0)) continue;

         if (compress && lastitem) {
            if (lastfield===fld) { ++cnt; lastkey = key; continue; }
            if (cnt>0) lastitem._name += ".." + lastkey;
         }

         let item = { _parent: top, _name: key };

         if (compress) { lastitem = item;  lastkey = key; lastfield = fld; cnt = 0; }

         if (fld === null) {
            item._value = item._title = "null";
            if (!nosimple) top._childs.push(item);
            continue;
         }

         let simple = false;

         if (typeof fld == 'object') {

            proto = Object.prototype.toString.apply(fld);

            if ((proto.lastIndexOf('Array]') == proto.length-6) && (proto.indexOf('[object')==0)) {
               item._title = "array len=" + fld.length;
               simple = (proto != '[object Array]');
               if (fld.length === 0) {
                  item._value = "[ ]";
                  item._more = false; // hpainter will not try to expand again
               } else {
                  item._value = "[...]";
                  item._more = true;
                  item._expand = objectHierarchy;
                  item._obj = fld;
               }
            } else if (proto === "[object DataView]") {
               item._title = 'DataView len=' + fld.byteLength;
               item._value = "[...]";
               item._more = true;
               item._expand = objectHierarchy;
               item._obj = fld;
            } else if (proto === "[object Date]") {
               item._more = false;
               item._title = 'Date';
               item._value = fld.toString();
               item._vclass = 'h_value_num';
            } else {

               if (fld.$kind || fld._typename)
                  item._kind = item._title = "ROOT." + (fld.$kind || fld._typename);

               if (fld._typename) {
                  item._title = fld._typename;
                  if (do_context && jsrp.canDraw(fld._typename)) item._direct_context = true;
               }

               // check if object already shown in hierarchy (circular dependency)
               let curr = top, inparent = false;
               while (curr && !inparent) {
                  inparent = (curr._obj === fld);
                  curr = curr._parent;
               }

               if (inparent) {
                  item._value = "{ prnt }";
                  item._vclass = 'h_value_num';
                  item._more = false;
                  simple = true;
               } else {
                  item._obj = fld;
                  item._more = false;

                  switch(fld._typename) {
                     case 'TColor': item._value = jsrp.getRGBfromTColor(fld); break;
                     case 'TText': item._value = fld.fTitle; break;
                     case 'TLatex': item._value = fld.fTitle; break;
                     case 'TObjString': item._value = fld.fString; break;
                     default:
                        if (JSROOT.isRootCollection(fld) && (typeof fld.arr === "object")) {
                           item._value = fld.arr.length ? "[...]" : "[]";
                           item._title += ", size:"  + fld.arr.length;
                           if (fld.arr.length>0) item._more = true;
                        } else {
                           item._more = true;
                           item._value = "{ }";
                        }
                  }
               }
            }
         } else if ((typeof fld === 'number') || (typeof fld === 'boolean')) {
            simple = true;
            if (key == 'fBits')
               item._value = "0x" + fld.toString(16);
            else
               item._value = fld.toString();
            item._vclass = 'h_value_num';
         } else if (typeof fld === 'string') {
            simple = true;
            item._value = '&quot;' + fld.replace(/\&/g, '&amp;').replace(/\"/g, '&quot;').replace(/</g, '&lt;').replace(/>/g, '&gt;') + '&quot;';
            item._vclass = 'h_value_str';
         } else if (typeof fld === 'undefined') {
            simple = true;
            item._value = "undefined";
            item._vclass = 'h_value_num';
         } else {
            simple = true;
            alert('miss ' + key + '  ' + typeof fld);
         }

         if (!simple || !nosimple)
            top._childs.push(item);
      }

      if (compress && lastitem && (cnt>0)) lastitem._name += ".." + lastkey;

      return true;
   }

   /** @summary Create hierarchy for streamer info object
     * @memberof JSROOT.Painter
     * @private */
   function createStreamerInfoContent(lst) {
      let h = { _name : "StreamerInfo", _childs : [] };

      for (let i = 0; i < lst.arr.length; ++i) {
         let entry = lst.arr[i]

         if (entry._typename == "TList") continue;

         if (typeof entry.fName == 'undefined') {
            console.warn(`strange element in StreamerInfo with type ${entry._typename}`);
            continue;
         }

         let item = {
            _name : entry.fName + ";" + entry.fClassVersion,
            _kind : "class " + entry.fName,
            _title : "class:" + entry.fName + ' version:' + entry.fClassVersion + ' checksum:' + entry.fCheckSum,
            _icon: "img_class",
            _childs : []
         };

         if (entry.fTitle != '') item._title += '  ' + entry.fTitle;

         h._childs.push(item);

         if (typeof entry.fElements == 'undefined') continue;
         for ( let l = 0; l < entry.fElements.arr.length; ++l) {
            let elem = entry.fElements.arr[l];
            if (!elem || !elem.fName) continue;
            let info = elem.fTypeName + " " + elem.fName,
                title = elem.fTypeName + " type:" + elem.fType;
            if (elem.fArrayDim===1)
               info += "[" + elem.fArrayLength + "]";
            else
               for (let dim=0;dim<elem.fArrayDim;++dim)
                  info+="[" + elem.fMaxIndex[dim] + "]";
            if (elem.fBaseVersion===4294967295) info += ":-1"; else
            if (elem.fBaseVersion!==undefined) info += ":" + elem.fBaseVersion;
            info += ";";
            if (elem.fTitle != '') info += " // " + elem.fTitle;

            item._childs.push({ _name : info, _title: title, _kind: elem.fTypeName, _icon: (elem.fTypeName == 'BASE') ? "img_class" : "img_member" });
         }
         if (item._childs.length == 0) delete item._childs;
      }

      return h;
   }

   /** @summary Create hierarchy for object inspector
     * @memberof JSROOT.Painter
     * @private */
   function createInspectorContent(obj) {
      let h = { _name: "Object", _title: "", _click_action: "expand", _nosimple: false, _do_context: true };

      if ((typeof obj.fName === 'string') && (obj.fName.length > 0))
         h._name = obj.fName;

      if ((typeof obj.fTitle === 'string') && (obj.fTitle.length > 0))
         h._title = obj.fTitle;

      if (obj._typename)
         h._title += "  type:" + obj._typename;

      if (JSROOT.isRootCollection(obj)) {
         h._name = obj.name || obj._typename;
         listHierarchy(h, obj);
      } else {
         objectHierarchy(h, obj);
      }

      return h;
   }


   /** @summary Parse string value as array.
    * @desc It could be just simple string:  "value" or
    * array with or without string quotes:  [element], ['elem1',elem2]
    * @private */
   let parseAsArray = val => {

      let res = [];

      if (typeof val != 'string') return res;

      val = val.trim();
      if (val=="") return res;

      // return as array with single element
      if ((val.length<2) || (val[0]!='[') || (val[val.length-1]!=']')) {
         res.push(val); return res;
      }

      // try to split ourself, checking quotes and brackets
      let nbr = 0, nquotes = 0, ndouble = 0, last = 1;

      for (let indx = 1; indx < val.length; ++indx) {
         if (nquotes > 0) {
            if (val[indx]==="'") nquotes--;
            continue;
         }
         if (ndouble > 0) {
            if (val[indx]==='"') ndouble--;
            continue;
         }
         switch (val[indx]) {
            case "'": nquotes++; break;
            case '"': ndouble++; break;
            case "[": nbr++; break;
            case "]": if (indx < val.length - 1) { nbr--; break; }
            case ",":
               if (nbr === 0) {
                  let sub =  val.substring(last, indx).trim();
                  if ((sub.length>1) && (sub[0]==sub[sub.length-1]) && ((sub[0]=='"') || (sub[0]=="'")))
                     sub = sub.substr(1, sub.length-2);
                  res.push(sub);
                  last = indx+1;
               }
               break;
         }
      }

      if (res.length === 0)
         res.push(val.substr(1, val.length-2).trim());

      return res;
   }

   // =================================================================================================

   /**
     * @summary special layout with three different areas for browser (left), status line (bottom) and central drawing
     * Main application is normal browser in JSROOT, but later one should be able to use it in ROOT6 canvas
     *
     * @class
     * @memberof JSROOT
     * @private
     */

   function BrowserLayout(id, hpainter, objpainter) {
      this.gui_div = id;
      this.hpainter = hpainter; // painter for brwoser area (if any)
      this.objpainter = objpainter; // painter for object area (if any)
      this.browser_kind = null; // should be 'float' or 'fix'
   }

   /** @summary Selects main element */
   BrowserLayout.prototype.main = function() {
      return d3.select("#" + this.gui_div);
   }

   /** @summary Returns drawing divid */
   BrowserLayout.prototype.drawing_divid = function() {
      return this.gui_div + "_drawing";
   }

   /** @summary Check resize action */
   BrowserLayout.prototype.checkResize = function() {
      if (this.hpainter && (typeof this.hpainter.checkResize == 'function'))
         this.hpainter.checkResize();
      else if (this.objpainter && (typeof this.objpainter.checkResize == 'function')) {
         this.objpainter.checkResize(true);
      }
   }

   /** @summary method used to create basic elements
     * @desc should be called only once */
   BrowserLayout.prototype.create = function(with_browser) {
      let main = this.main();

      main.append("div").attr("id", this.drawing_divid())
                        .classed("jsroot_draw_area", true)
                        .style('position',"absolute").style('left',0).style('top',0).style('bottom',0).style('right',0);

      if (with_browser) main.append("div").classed("jsroot_browser", true);
   }

   /** @summary Create buttons in the layout */
   BrowserLayout.prototype.createBrowserBtns = function() {
      let br = this.main().select(".jsroot_browser");
      if (br.empty()) return;
      let btns = br.append("div").classed("jsroot_browser_btns", true).classed("jsroot", true);
      btns.style('position',"absolute").style("left","7px").style("top","7px");
      if (JSROOT.browser.touches) btns.style('opacity','0.2'); // on touch devices should be always visible
      return btns;
   }

   /** @summary Remove browser buttons */
   BrowserLayout.prototype.removeBrowserBtns = function() {
      this.main().select(".jsroot_browser").select(".jsroot_browser_btns").remove();
   }

   /** @summary Set browser content */
   BrowserLayout.prototype.setBrowserContent = function(guiCode) {
      let main = d3.select("#" + this.gui_div + " .jsroot_browser");
      if (main.empty()) return;

      main.insert('div', ".jsroot_browser_btns").classed('jsroot_browser_area', true)
          .style('position',"absolute").style('left',0).style('top',0).style('bottom',0).style('width','250px')
          .style('padding-left','5px')
          .style('display','flex').style('flex-direction', 'column')   /* use the flex model */
          .html("<p class='jsroot_browser_title'>title</p>" +  guiCode);
   }

   /** @summary Check if there is browser content */
   BrowserLayout.prototype.hasContent = function() {
      let main = d3.select("#" + this.gui_div + " .jsroot_browser");
      if (main.empty()) return false;
      return !main.select(".jsroot_browser_area").empty();
   }

   /** @summary Delete browser content */
   BrowserLayout.prototype.deleteContent = function() {
      let main = d3.select("#" + this.gui_div + " .jsroot_browser");
      if (main.empty()) return;

      main.selectAll("*").remove();
      delete this.browser_visible;
   }

   /** @summary Returns true when status line exists */
   BrowserLayout.prototype.hasStatus = function() {
      let main = d3.select("#"+this.gui_div+" .jsroot_browser");
      if (main.empty()) return false;

      let id = this.gui_div + "_status",
          line = d3.select("#"+id);

      return !line.empty();
   }

   /** @summary Create status line */
   BrowserLayout.prototype.createStatusLine = function(height, mode) {
      if (!this.gui_div) return Promise.resolve('');
      return JSROOT.require('jq2d').then(() => this.createStatusLine(height, mode));
   }

   // ==============================================================================


   /** @summary central function for expand of all online items
     * @private */
   function onlineHierarchy(node, obj) {
      if (obj && node && ('_childs' in obj)) {

         for (let n=0;n<obj._childs.length;++n)
            if (obj._childs[n]._more || obj._childs[n]._childs)
               obj._childs[n]._expand = onlineHierarchy;

         node._childs = obj._childs;
         obj._childs = null;
         return true;
      }

      return false;
   }

   // ==============================================================

   /** @summary Current hierarchy painter
     * @desc Instance of {@link JSROOT.HierarchyPainter} object
     * @private */
   JSROOT.hpainter = null;

   /**
     * @summary Painter of hierarchical structures
     *
     * @class
     * @memberof JSROOT
     * @param {string} name - symbolic name
     * @param {string} frameid - element id where hierarchy is drawn
     * @param {string} [backgr] - background color
     * @example
     * // create hierarchy painter in "myTreeDiv"
     * let h = new JSROOT.HierarchyPainter("example", "myTreeDiv");
     * // configure 'simple' layout in "myMainDiv"
     * // one also can specify "grid2x2" or "flex" or "tabs"
     * h.setDisplay("simple", "myMainDiv");
     * // open file and display element
     * h.openRootFile("https://root.cern/js/files/hsimple.root").then(() => h.display("hpxpy;1","colz")); */
   function HierarchyPainter(name, frameid, backgr) {
      JSROOT.BasePainter.call(this, frameid);
      this.name = name;
      this.h = null; // hierarchy
      this.with_icons = true;
      this.background = backgr;
      this.files_monitoring = !frameid; // by default files monitored when nobrowser option specified
      this.nobrowser = (frameid === null);

      // remember only very first instance
      if (!JSROOT.hpainter)
         JSROOT.hpainter = this;
   }

   HierarchyPainter.prototype = Object.create(JSROOT.BasePainter.prototype);

   /** @summary Cleanup hierarchy painter
     * @desc clear drawing and browser */
   HierarchyPainter.prototype.cleanup = function() {
      this.clearHierarchy(true);

      JSROOT.BasePainter.prototype.cleanup.call(this);

      if (JSROOT.hpainter === this)
         JSROOT.hpainter = null;
   }

   /** @summary Create file hierarchy
     * @private */
   HierarchyPainter.prototype.fileHierarchy = function(file) {
      let painter = this;

      let folder = {
         _name : file.fFileName,
         _title : (file.fTitle ? (file.fTitle + ", path ") : "")  + file.fFullURL,
         _kind : "ROOT.TFile",
         _file : file,
         _fullurl : file.fFullURL,
         _localfile : file.fLocalFile,
         _had_direct_read : false,
         // this is central get method, item or itemname can be used, returns promise
         _get : function(item, itemname) {

            let fff = this; // file item

            if (item && item._readobj)
               return Promise.resolve(item._readobj);

            if (item) itemname = painter.itemFullName(item, fff);

            function ReadFileObject(file) {
               if (!fff._file) fff._file = file;

               if (!file) return Promise.resolve(null);

               return file.readObject(itemname).then(obj => {

                  // if object was read even when item did not exist try to reconstruct new hierarchy
                  if (!item && obj) {
                     // first try to found last read directory
                     let d = painter.findItem({name:itemname, top:fff, last_exists:true, check_keys:true });
                     if ((d!=null) && ('last' in d) && (d.last!=fff)) {
                        // reconstruct only subdir hierarchy
                        let dir = file.getDir(painter.itemFullName(d.last, fff));
                        if (dir) {
                           d.last._name = d.last._keyname;
                           let dirname = painter.itemFullName(d.last, fff);
                           keysHierarchy(d.last, dir.fKeys, file, dirname + "/");
                        }
                     } else {
                        // reconstruct full file hierarchy
                        keysHierarchy(fff, file.fKeys, file, "");
                     }
                     item = painter.findItem({name:itemname, top: fff});
                  }

                  if (item) {
                     item._readobj = obj;
                     // remove cycle number for objects supporting expand
                     if ('_expand' in item) item._name = item._keyname;
                  }

                  return obj;
               });
            }

            if (fff._file) return ReadFileObject(fff._file);
            if (fff._localfile) return JSROOT.openFile(fff._localfile).then(f => ReadFileObject(f));
            if (fff._fullurl) return JSROOT.openFile(fff._fullurl).then(f => ReadFileObject(f));
            return Promise.resolve(null);
         }
      };

      keysHierarchy(folder, file.fKeys, file, "");

      return folder;
   }

   /** @summary Iterate over all items in hierarchy
     * @param {function} func - function called for every item
     * @param {object} [top] - top item to start from
     * @private */
   HierarchyPainter.prototype.forEachItem = function(func, top) {
      function each_item(item, prnt) {
         if (!item) return;
         if (prnt) item._parent = prnt;
         func(item);
         if ('_childs' in item)
            for (let n = 0; n < item._childs.length; ++n)
               each_item(item._childs[n], item);
      }

      if (typeof func == 'function')
         each_item(top || this.h);
   }

   /** @summary Search item in the hierarchy
     * @param {object|string} arg - item name or object with arguments
     * @param {string} arg.name -  item to search
     * @param {boolean} [arg.force] - specified elements will be created when not exists
     * @param {boolean} [arg.last_exists] -  when specified last parent element will be returned
     * @param {boolean} [arg.check_keys] - check TFile keys with cycle suffix
     * @param {boolean} [arg.allow_index] - let use sub-item indexes instead of name
     * @param {object} [arg.top] - element to start search from
     * @private */
   HierarchyPainter.prototype.findItem = function(arg) {

      function find_in_hierarchy(top, fullname) {

         if (!fullname || (fullname.length == 0) || !top) return top;

         let pos = fullname.length;

         if (!top._parent && (top._kind !== 'TopFolder') && (fullname.indexOf(top._name)===0)) {
            // it is allowed to provide item name, which includes top-parent like file.root/folder/item
            // but one could skip top-item name, if there are no other items
            if (fullname === top._name) return top;

            let len = top._name.length;
            if (fullname[len] == "/") {
               fullname = fullname.substr(len+1);
               pos = fullname.length;
            }
         }

         function process_child(child, ignore_prnt) {
            // set parent pointer when searching child
            if (!ignore_prnt) child._parent = top;

            if ((pos >= fullname.length-1) || (pos < 0)) return child;

            return find_in_hierarchy(child, fullname.substr(pos + 1));
         }

         while (pos > 0) {
            // we try to find element with slashes inside - start from full name
            let localname = (pos >= fullname.length) ? fullname : fullname.substr(0, pos);

            if (top._childs) {
               // first try to find direct matched item
               for (let i = 0; i < top._childs.length; ++i)
                  if (top._childs[i]._name == localname)
                     return process_child(top._childs[i]);

               // if first child online, check its elements
               if ((top._kind === 'TopFolder') && (top._childs[0]._online!==undefined))
                  for (let i = 0; i < top._childs[0]._childs.length; ++i)
                     if (top._childs[0]._childs[i]._name == localname)
                        return process_child(top._childs[0]._childs[i], true);

               // if allowed, try to found item with key
               if (arg.check_keys) {
                  let newest = null;
                  for (let i = 0; i < top._childs.length; ++i) {
                    if (top._childs[i]._keyname === localname) {
                       if (!newest || (newest._cycle < top._childs[i]._cycle)) newest = top._childs[i];
                    }
                  }
                  if (newest) return process_child(newest);
               }

               let allow_index = arg.allow_index;
               if ((localname[0] === '[') && (localname[localname.length-1] === ']') &&
                    /^\d+$/.test(localname.substr(1,localname.length-2))) {
                  allow_index = true;
                  localname = localname.substr(1,localname.length-2);
               }

               // when search for the elements it could be allowed to check index
               if (allow_index && /^\d+$/.test(localname)) {
                  let indx = parseInt(localname);
                  if (!isNaN(indx) && (indx >= 0) && (indx < top._childs.length))
                     return process_child(top._childs[indx]);
               }
            }

            pos = fullname.lastIndexOf("/", pos - 1);
         }

         if (arg.force) {
             // if did not found element with given name we just generate it
             if (top._childs === undefined) top._childs = [];
             pos = fullname.indexOf("/");
             let child = { _name: ((pos < 0) ? fullname : fullname.substr(0, pos)) };
             top._childs.push(child);
             return process_child(child);
         }

         return (arg.last_exists && top) ? { last: top, rest: fullname } : null;
      }

      let top = this.h, itemname = "";

      if (arg === null) return null; else
      if (typeof arg == 'string') { itemname = arg; arg = {}; } else
      if (typeof arg == 'object') { itemname = arg.name; if ('top' in arg) top = arg.top; } else
         return null;

      if (itemname === "__top_folder__") return top;

      if ((typeof itemname == 'string') && (itemname.indexOf("img:")==0)) return null;

      return find_in_hierarchy(top, itemname);
   }

   /** @summary Produce full string name for item
     * @param {Object} node - item element
     * @param {Object} [uptoparent] - up to which parent to continue
     * @param {boolean} [compact] - if specified, top parent is not included
     * @returns {string} produced name
     * @private */
   HierarchyPainter.prototype.itemFullName = function(node, uptoparent, compact) {

      if (node && node._kind ==='TopFolder') return "__top_folder__";

      let res = "";

      while (node) {
         // online items never includes top-level folder
         if ((node._online!==undefined) && !uptoparent) return res;

         if ((node === uptoparent) || (node._kind==='TopFolder')) break;
         if (compact && !node._parent) break; // in compact form top-parent is not included
         if (res.length > 0) res = "/" + res;
         res = node._name + res;
         node = node._parent;
      }

      return res;
   }

    /** @summary Executes item marked as 'Command'
      * @desc If command requires additional arguments, they could be specified as extra arguments arg1, arg2, ...
      * @param {String} itemname - name of command item
      * @param {Object} [elem] - HTML element for command execution
      * @param [arg1] - first optional argument
      * @param [arg2] - second optional argument and so on
      * @returns {Promise} with command result */
   HierarchyPainter.prototype.executeCommand = function(itemname, elem) {

      let hitem = this.findItem(itemname),
          url = this.getOnlineItemUrl(hitem) + "/cmd.json",
          d3node = d3.select(elem);

      if ('_numargs' in hitem)
         for (let n = 0; n < hitem._numargs; ++n) {
            let argname = "arg" + (n+1), argvalue = undefined;
            if (n+2 < arguments.length) argvalue = arguments[n+2];
            if ((argvalue===undefined) && elem)
               argvalue = prompt("Input argument " + argname + " for command " + hitem._name, "");
            if (argvalue===null) return Promise.resolve(false);
            url += ((n==0) ? "?" : "&") + argname + "=" + argvalue;
         }

      if (!d3node.empty()) {
         d3node.style('background','yellow');
         if (hitem && hitem._title) d3node.attr('title', "Executing " + hitem._title);
      }

      return JSROOT.httpRequest(url, 'text').then(res => {
         if (!d3node.empty()) {
            let col = ((res != null) && (res != 'false')) ? 'green' : 'red';
            if (hitem && hitem._title) d3node.attr('title', hitem._title + " lastres=" + res);
            d3node.style('background', col);
            setTimeout(() => d3node.style('background', ''), 2000);
            if ((col == 'green') && ('_hreload' in hitem))
               this.reload();
            if ((col == 'green') && ('_update_item' in hitem))
               this.updateItems(hitem._update_item.split(";"));
         }
         return res;
      });
   }

   /** @summary Refresh HTML for hierachy painter
     * @returns {Promise} when completed
     * @private */
   HierarchyPainter.prototype.refreshHtml = function() {
      if (!this.getDom() || JSROOT.batch_mode)
         return Promise.resolve(this);
      return JSROOT.require('jq2d').then(() => this.refreshHtml());
   }

   /** @summary Get object item with specified name
     * @desc depending from provided option, same item can generate different object types
     * @param {Object} arg - item name or config object
     * @param {string} arg.name - item name
     * @param {Object} arg.item - or item itself
     * @param {string} options - supposed draw options
     * @returns {Promise} with object like { item, obj, itemname }
     * @private */
   HierarchyPainter.prototype.getObject = function(arg, options) {

      let itemname, item, result = { item: null, obj: null };

      if (arg===null) return Promise.resolve(result);

      if (typeof arg === 'string') {
         itemname = arg;
      } else if (typeof arg === 'object') {
         if ((arg._parent!==undefined) && (arg._name!==undefined) && (arg._kind!==undefined)) item = arg; else
         if (arg.name!==undefined) itemname = arg.name; else
         if (arg.arg!==undefined) itemname = arg.arg; else
         if (arg.item!==undefined) item = arg.item;
      }

      if ((typeof itemname == 'string') && (itemname.indexOf("img:")==0)) {
         // artificial class, can be created by users
         result.obj = {_typename: "TJSImage", fName: itemname.substr(4)};
         return Promise.resolve(result);
      }

      if (item) itemname = this.itemFullName(item);
           else item = this.findItem( { name: itemname, allow_index: true, check_keys: true } );

      // if item not found, try to find nearest parent which could allow us to get inside

      let d = item ? null : this.findItem({ name: itemname, last_exists: true, check_keys: true, allow_index: true });

      // if item not found, try to expand hierarchy central function
      // implements not process get in central method of hierarchy item (if exists)
      // if last_parent found, try to expand it
      if ((d !== null) && ('last' in d) && (d.last !== null)) {
         let parentname = this.itemFullName(d.last);

         // this is indication that expand does not give us better path to searched item
         if ((typeof arg == 'object') && ('rest' in arg))
            if ((arg.rest == d.rest) || (arg.rest.length <= d.rest.length))
               return Promise.resolve(result);

         return this.expandItem(parentname, null, true).then(res => {
            if (!res) return result;
            let newparentname = this.itemFullName(d.last);
            if (newparentname.length > 0) newparentname += "/";
            return this.getObject( { name: newparentname + d.rest, rest: d.rest }, options);
         });
      }

      result.item = item;

      if ((item !== null) && (typeof item._obj == 'object')) {
         result.obj = item._obj;
         return Promise.resolve(result);
      }

      // normally search _get method in the parent items
      let curr = item;
      while (curr) {
         if (typeof curr._get == 'function')
            return curr._get(item, null, options).then(obj => { result.obj = obj; return result; });
         curr = ('_parent' in curr) ? curr._parent : null;
      }

      return Promise.resolve(result);
   }

   /** @summary Starts player for specified item
     * @desc Same as "Player" context menu
     * @param {string} itemname - item name for which player should be started
     * @param {string} [option] - extra options for the player
     * @returns {Promise} when ready*/
   HierarchyPainter.prototype.player = function(itemname, option) {
      let item = this.findItem(itemname);

      if (!item || !item._player) return Promise.resolve(null);

      return JSROOT.require(item._prereq || '').then(() => {

         let player_func = JSROOT.findFunction(item._player);
         if (!player_func) return null;

         return this.createDisplay().then(mdi => {
            let res = mdi ? player_func(this, itemname, option) : null;
            return res;
         });
      });
   }

   /** @summary Checks if item can be displayed with given draw option
     * @private */
   HierarchyPainter.prototype.canDisplay = function(item, drawopt) {
      if (!item) return false;
      if (item._player) return true;
      if (item._can_draw !== undefined) return item._can_draw;
      if (drawopt == 'inspect') return true;
      let handle = jsrp.getDrawHandle(item._kind, drawopt);
      return handle && (('func' in handle) || ('draw_field' in handle));
   }

   /** @summary Returns true if given item displayed
     * @param {string} itemname - item name */
   HierarchyPainter.prototype.isItemDisplayed = function(itemname) {
      let mdi = this.getDisplay();
      return mdi ? mdi.findFrame(itemname) !== null : false;
   }

   /** @summary Display specified item
     * @param {string} itemname - item name
     * @param {string} [drawopt] - draw option for the item
     * @returns {Promise} with created painter object */
   HierarchyPainter.prototype.display = function(itemname, drawopt) {
      let h = this,
          painter = null,
          updating = false,
          item = null,
          display_itemname = itemname,
          frame_name = itemname,
          marker = "::_display_on_frame_::",
          p = drawopt ? drawopt.indexOf(marker) : -1;

      if (p>=0) {
         frame_name = drawopt.substr(p + marker.length);
         drawopt = drawopt.substr(0, p);
      }

      function complete(respainter, err) {
         if (err) console.log('When display ', itemname, err);

         if (updating && item) delete item._doing_update;
         if (!updating) jsrp.showProgress();
         if (respainter && (typeof respainter === 'object') && (typeof respainter.setItemName === 'function')) {
            respainter.setItemName(display_itemname, updating ? null : drawopt, h); // mark painter as created from hierarchy
            if (item && !item._painter) item._painter = respainter;
         }

         return respainter || painter;
      }

      return h.createDisplay().then(mdi => {

         if (!mdi) return complete();

         item = h.findItem(display_itemname);

         if (item && ('_player' in item))
            return h.player(display_itemname, drawopt).then(res => complete(res));

         updating = (typeof drawopt == 'string') && (drawopt.indexOf("update:")==0);

         if (updating) {
            drawopt = drawopt.substr(7);
            if (!item || item._doing_update) return complete();
            item._doing_update = true;
         }

         if (item && !h.canDisplay(item, drawopt)) return complete();

         let divid = "";
         if ((typeof drawopt == 'string') && (drawopt.indexOf("divid:") >= 0)) {
            let pos = drawopt.indexOf("divid:");
            divid = drawopt.slice(pos+6);
            drawopt = drawopt.slice(0, pos);
         }

         if (!updating) jsrp.showProgress("Loading " + display_itemname);

         return h.getObject(display_itemname, drawopt).then(result => {
            if (!updating) jsrp.showProgress();

            if (!item) item = result.item;
            let obj = result.obj;

            if (!obj) return complete();

            if (!updating) jsrp.showProgress("Drawing " + display_itemname);

            if (divid.length > 0) {
               let func = updating ? JSROOT.redraw : JSROOT.draw;
               return func(divid, obj, drawopt).then(p => complete(p)).catch(err => complete(null, err));
            }

            mdi.forEachPainter((p, frame) => {
               if (p.getItemName() != display_itemname) return;
               // verify that object was drawn with same option as specified now (if any)
               if (!updating && drawopt && (p.getItemDrawOpt() != drawopt)) return;
               mdi.activateFrame(frame);

               let handle = null;
               if (obj._typename) handle = jsrp.getDrawHandle("ROOT." + obj._typename);
               if (handle && handle.draw_field && obj[handle.draw_field])
                  obj = obj[handle.draw_field];

               if ((typeof p.redrawObject == 'function') && p.redrawObject(obj, drawopt)) painter = p;
            });

            if (painter) return complete();

            if (updating) {
               console.warn(`something went wrong - did not found painter when doing update of ${display_itemname}`);
               return complete();
            }

            let frame = mdi.findFrame(frame_name, true);
            d3.select(frame).html("");
            mdi.activateFrame(frame);

            return JSROOT.draw(frame, obj, drawopt).then(p => {
               if (JSROOT.settings.DragAndDrop)
                  h.enableDrop(frame, display_itemname);
               return complete(p);
            }).catch(err => complete(null, err));

         });
      });
   }

   HierarchyPainter.prototype.enableDrag = function(/*element, itemname*/) {
      // here is not defined - implemented with jquery
   }

   HierarchyPainter.prototype.enableDrop = function(/*frame, itemname*/) {
      // here is not defined - implemented with jquery
   }

  /** @summary Drop item on specified element for drawing
    * @returns {Promise} when completed
    * @private */
   HierarchyPainter.prototype.dropItem = function(itemname, divid, opt) {

      if (opt && typeof opt === 'function') { call_back = opt; opt = ""; }
      if (opt===undefined) opt = "";

      let drop_complete = drop_painter => {
         if (drop_painter && (typeof drop_painter === 'object') && (typeof drop_painter.setItemName == 'function'))
            drop_painter.setItemName(itemname, null, this);
         return drop_painter;
      }

      if (itemname == "$legend")
         return JSROOT.require("hist").then(() => {
            let legend_painter = jsrp.produceLegend(divid, opt);
            return drop_complete(legend_painter);
         });

      return this.getObject(itemname).then(res => {

         if (!res.obj) return null;

         let main_painter = jsrp.getElementMainPainter(divid);

         if (main_painter && (typeof main_painter.performDrop === 'function'))
            return main_painter.performDrop(res.obj, itemname, res.item, opt).then(p => drop_complete(p));

         if (main_painter && main_painter.accept_drops)
            return JSROOT.draw(divid, res.obj, "same " + opt).then(p => drop_complete(p));

         this.cleanupFrame(divid);
         return JSROOT.draw(divid, res.obj, opt).then(p => drop_complete(p));
      });
   }

   /** @summary Update specified items
     * @desc Method can be used to fetch new objects and update all existing drawings
     * @param {string|array|boolean} arg - either item name or array of items names to update or true if only automatic items will be updated
     * @returns {Promise} when ready */
   HierarchyPainter.prototype.updateItems = function(arg) {

      if (!this.disp)
         return Promise.resolve(false);

      let allitems = [], options = [], only_auto_items = false, want_update_all = false;

      if (typeof arg == "string")
         arg = [ arg ];
      else if (typeof arg != 'object') {
         if (arg === undefined) arg = !this.isMonitoring();
         want_update_all = true;
         only_auto_items = !!arg;
      }

      // first collect items
      this.disp.forEachPainter(p => {
         let itemname = p.getItemName();

         if ((typeof itemname != 'string') || (allitems.indexOf(itemname) >= 0)) return;

         if (want_update_all) {
            let item = this.findItem(itemname);
            if (!item || ('_not_monitor' in item) || ('_player' in item)) return;
            if (!('_always_monitor' in item)) {
               let forced = false, handle = jsrp.getDrawHandle(item._kind);
               if (handle && ('monitor' in handle)) {
                  if ((handle.monitor === false) || (handle.monitor == 'never')) return;
                  if (handle.monitor == 'always') forced = true;
               }
               if (!forced && only_auto_items) return;
            }
         } else {
            if (arg.indexOf(itemname) < 0) return;
         }

         allitems.push(itemname);
         options.push("update:" + p.getItemDrawOpt());
      }, true); // only visible panels are considered

      // force all files to read again (normally in non-browser mode)
      if (this.files_monitoring && !only_auto_items && want_update_all)
         this.forEachRootFile(item => {
            this.forEachItem(fitem => { delete fitem._readobj; }, item);
            delete item._file;
         });

      return this.displayItems(allitems, options);
   }

   /** @summary Display all provided elements
     * @returns {Promise} when drawing finished
     * @private */
   HierarchyPainter.prototype.displayItems = function(items, options) {

      if (!items || (items.length == 0))
         return Promise.resolve(true);

      let h = this;

      if (!options) options = [];
      while (options.length < items.length)
         options.push("");

      if ((options.length == 1) && (options[0] == "iotest")) {
         this.clearHierarchy();
         d3.select("#" + this.disp_frameid).html("<h2>Start I/O test</h2>")

         let tm0 = new Date();
         return this.getObject(items[0]).then(() => {
            let tm1 = new Date();
            d3.select("#" + this.disp_frameid).append("h2").html("Item " + items[0] + " reading time = " + (tm1.getTime() - tm0.getTime()) + "ms");
            return Promise.resolve(true);
         });
      }

      let dropitems = new Array(items.length), dropopts = new Array(items.length), images = new Array(items.length);

      // First of all check that items are exists, look for cycle extension and plus sign
      for (let i = 0; i < items.length; ++i) {
         dropitems[i] = dropopts[i] = null;

         let item = items[i], can_split = true;

         if (item && item.indexOf("img:")==0) { images[i] = true; continue; }

         if (item && (item.length>1) && (item[0]=='\'') && (item[item.length-1]=='\'')) {
            items[i] = item.substr(1, item.length-2);
            can_split = false;
         }

         let elem = h.findItem({ name: items[i], check_keys: true });
         if (elem) { items[i] = h.itemFullName(elem); continue; }

         if (can_split && (items[i][0]=='[') && (items[i][items[i].length-1]==']')) {
            dropitems[i] = parseAsArray(items[i]);
            items[i] = dropitems[i].shift();
         } else
         if (can_split && (items[i].indexOf("+") > 0)) {
            dropitems[i] = items[i].split("+");
            items[i] = dropitems[i].shift();
         }

         if (dropitems[i] && dropitems[i].length > 0) {
            // allow to specify _same_ item in different file
            for (let j = 0; j < dropitems[i].length; ++j) {
               let pos = dropitems[i][j].indexOf("_same_");
               if ((pos>0) && (h.findItem(dropitems[i][j])===null))
                  dropitems[i][j] = dropitems[i][j].substr(0,pos) + items[i].substr(pos);

               elem = h.findItem({ name: dropitems[i][j], check_keys: true });
               if (elem) dropitems[i][j] = h.itemFullName(elem);
            }

            if ((options[i][0] == "[") && (options[i][options[i].length-1] == "]")) {
               dropopts[i] = parseAsArray(options[i]);
               options[i] = dropopts[i].shift();
            } else
            if (options[i].indexOf("+") > 0) {
               dropopts[i] = options[i].split("+");
               options[i] = dropopts[i].shift();
            } else {
               dropopts[i] = [];
            }

            while (dropopts[i].length < dropitems[i].length) dropopts[i].push("");
         }

         // also check if subsequent items has _same_, than use name from first item
         let pos = items[i].indexOf("_same_");
         if ((pos>0) && !h.findItem(items[i]) && (i>0))
            items[i] = items[i].substr(0,pos) + items[0].substr(pos);

         elem = h.findItem({ name: items[i], check_keys: true });
         if (elem) items[i] = h.itemFullName(elem);
      }

      // now check that items can be displayed
      for (let n = items.length-1; n>=0; --n) {
         if (images[n]) continue;
         let hitem = h.findItem(items[n]);
         if (!hitem || h.canDisplay(hitem, options[n])) continue;
         // try to expand specified item
         h.expandItem(items[n], null, true);
         items.splice(n, 1);
         options.splice(n, 1);
         dropitems.splice(n, 1);
      }

      if (items.length == 0)
         return Promise.resolve(true);

      let frame_names = new Array(items.length), items_wait = new Array(items.length);
      for (let n=0; n < items.length;++n) {
         items_wait[n] = 0;
         let fname = items[n], k = 0;
         if (items.indexOf(fname) < n) items_wait[n] = true; // if same item specified, one should wait first drawing before start next
         let p = options[n].indexOf("frameid:");
         if (p>=0) {
            fname = options[n].substr(p+8);
            options[n] = options[n].substr(0,p);
         } else {
            while (frame_names.indexOf(fname)>=0)
               fname = items[n] + "_" + k++;
         }
         frame_names[n] = fname;
      }

      // now check if several same items present - select only one for the drawing
      // if draw option includes 'main', such item will be drawn first
      for (let n=0; n<items.length;++n) {
         if (items_wait[n] !== 0) continue;
         let found_main = n;
         for (let k=0; k<items.length;++k)
            if ((items[n]===items[k]) && (options[k].indexOf('main')>=0)) found_main = k;
         for (let k=0; k<items.length;++k)
            if (items[n]===items[k]) items_wait[k] = (found_main != k);
      }

      return this.createDisplay().then(mdi => {
         if (!mdi) return false;

         // Than create empty frames for each item
         for (let i = 0; i < items.length; ++i)
            if (options[i].indexOf('update:')!==0) {
               mdi.createFrame(frame_names[i]);
               options[i] += "::_display_on_frame_::"+frame_names[i];
            }

         function DropNextItem(indx, painter) {
            if (painter && dropitems[indx] && (dropitems[indx].length > 0))
               return h.dropItem(dropitems[indx].shift(), painter.getDom(), dropopts[indx].shift()).then(() => DropNextItem(indx, painter));

            dropitems[indx] = null; // mark that all drop items are processed
            items[indx] = null; // mark item as ready

            for (let cnt = 0; cnt < items.length; ++cnt) {
               if (dropitems[cnt]) isany = true;
               if (items[cnt]===null) continue; // ignore completed item
               if (items_wait[cnt] && items.indexOf(items[cnt])===cnt) {
                  items_wait[cnt] = false;
                  return h.display(items[cnt], options[cnt]).then(painter => DropNextItem(cnt, painter));
               }
            }
         }

         let promises = [];

         // We start display of all items parallel, but only if they are not the same
         for (let i = 0; i < items.length; ++i)
            if (!items_wait[i])
               promises.push(h.display(items[i], options[i]).then(painter => DropNextItem(i, painter)));

         return Promise.all(promises);
      });
   }

   /** @summary Reload hierarchy and refresh html code
     * @returns {Promise} when completed */
   HierarchyPainter.prototype.reload = function() {
      if ('_online' in this.h)
         return this.openOnline(this.h._online).then(() => this.refreshHtml());
      return Promise.resolve(false);
   }

   HierarchyPainter.prototype.updateTreeNode = function() {
      // dummy function, will be redefined when jquery part loaded
   }

   /** @summary activate (select) specified item
     * @param {Array} items - array of items names
     * @param {boolean} [force] - if specified, all required sub-levels will be opened
     * @private */
   HierarchyPainter.prototype.activateItems = function(items, force) {

      if (typeof items == 'string') items = [ items ];

      let active = [], // array of elements to activate
          update = []; // array of elements to update
      this.forEachItem(item => { if (item._background) { active.push(item); delete item._background; } });

      let mark_active = () => {
         if (typeof this.updateBackground !== 'function') return;

         for (let n=update.length-1;n>=0;--n)
            this.updateTreeNode(update[n]);

         for (let n=0;n<active.length;++n)
            this.updateBackground(active[n], force);
      }

      let find_next = (itemname, prev_found) => {
         if (itemname === undefined) {
            // extract next element
            if (items.length == 0) return mark_active();
            itemname = items.shift();
         }

         let hitem = this.findItem(itemname);

         if (!hitem) {
            let d = this.findItem({ name: itemname, last_exists: true, check_keys: true, allow_index: true });
            if (!d || !d.last) return find_next();
            d.now_found = this.itemFullName(d.last);

            if (force) {

               // if after last expand no better solution found - skip it
               if ((prev_found!==undefined) && (d.now_found === prev_found)) return find_next();

               return this.expandItem(d.now_found).then(res => {
                  if (!res) return find_next();
                  let newname = this.itemFullName(d.last);
                  if (newname.length > 0) newname+="/";
                  find_next(newname + d.rest, d.now_found);
               });
            }
            hitem = d.last;
         }

         if (hitem) {
            // check that item is visible (opened), otherwise should enable parent

            let prnt = hitem._parent;
            while (prnt) {
               if (!prnt._isopen) {
                  if (force) {
                     prnt._isopen = true;
                     if (update.indexOf(prnt)<0) update.push(prnt);
                  } else {
                     hitem = prnt; break;
                  }
               }
               prnt = prnt._parent;
            }

            hitem._background = 'LightSteelBlue';
            if (active.indexOf(hitem)<0) active.push(hitem);
         }

         find_next();
      }

      if (force && this.brlayout) {
         if (!this.brlayout.browser_kind)
           return this.createBrowser('float', true).then(() => find_next());
         if (!this.brlayout.browser_visible)
            this.brlayout.toggleBrowserVisisbility();
      }

      // use recursion
      find_next();
   }

   /** @summary Check if item can be (potentially) expand
     * @private */
   HierarchyPainter.prototype.canExpandItem = function(item) {
      if (!item) return false;
      if (item._expand) return true;
      let handle = jsrp.getDrawHandle(item._kind, "::expand");
      return handle && (handle.expand_item || handle.expand);
   }

   /** @summary expand specified item
     * @param {String} itemname - item name
     * @returns {Promise} when ready */
   HierarchyPainter.prototype.expandItem = function(itemname, d3cont, silent) {
      let hitem = this.findItem(itemname);

      if (!hitem && d3cont)
         return Promise.resolve();

      let DoExpandItem = (_item, _obj) => {

         if (typeof _item._expand == 'string')
            _item._expand = JSROOT.findFunction(item._expand);

         if (typeof _item._expand !== 'function') {
            let handle = jsrp.getDrawHandle(_item._kind, "::expand");

            if (handle && handle.expand_item) {
               _obj = _obj[handle.expand_item];
              if (_obj && _obj._typename)
                 handle = jsrp.getDrawHandle("ROOT."+_obj._typename, "::expand");
            }

            if (handle && handle.expand) {
               if (typeof handle.expand == 'string')
                  return JSROOT.require(handle.prereq).then(() => {
                     _item._expand = handle.expand = JSROOT.findFunction(handle.expand);
                     return _item._expand ? DoExpandItem(_item, _obj) : true;
                  });
               _item._expand = handle.expand;
            }
         }

         // try to use expand function
         if (_obj && _item && (typeof _item._expand === 'function')) {
            if (_item._expand(_item, _obj)) {
               _item._isopen = true;
               if (_item._parent && !_item._parent._isopen) {
                  _item._parent._isopen = true; // also show parent
                  if (!silent) this.updateTreeNode(_item._parent);
               } else {
                  if (!silent) this.updateTreeNode(_item, d3cont);
               }
               return Promise.resolve(_item);
            }
         }

         if (_obj && objectHierarchy(_item, _obj)) {
            _item._isopen = true;
            if (_item._parent && !_item._parent._isopen) {
               _item._parent._isopen = true; // also show parent
               if (!silent) this.updateTreeNode(_item._parent);
            } else {
               if (!silent) this.updateTreeNode(_item, d3cont);
            }
            return Promise.resolve(_item);
         }

         return Promise.resolve(-1);
      }

      let promise = Promise.resolve(-1);

      if (hitem) {
         // item marked as it cannot be expanded, also top item cannot be changed
         if ((hitem._more === false) || (!hitem._parent && hitem._childs))
            return Promise.resolve();

         if (hitem._childs && hitem._isopen) {
            hitem._isopen = false;
            if (!silent) this.updateTreeNode(hitem, d3cont);
            return Promise.resolve();
         }

         if (hitem._obj) promise = DoExpandItem(hitem, hitem._obj, itemname);
      }

      return promise.then(res => {
         if (res !== -1) return res; // done

         jsrp.showProgress("Loading " + itemname);

         return this.getObject(itemname, "hierarchy_expand").then(res => {

            jsrp.showProgress();

            if (res.obj) return DoExpandItem(res.item, res.obj).then(res => { return (res !== -1) ? res : undefined; });
         });
      });

   }

   /** @summary Return main online item
     * @private */
   HierarchyPainter.prototype.getTopOnlineItem = function(item) {
      if (item) {
         while (item && (!('_online' in item))) item = item._parent;
         return item;
      }

      if (!this.h) return null;
      if ('_online' in this.h) return this.h;
      if (this.h._childs && ('_online' in this.h._childs[0])) return this.h._childs[0];
      return null;
   }

   /** @summary Call function for each item which corresponds to JSON file
     * @private */
   HierarchyPainter.prototype.forEachJsonFile = function(func) {
      if (!this.h) return;
      if ('_jsonfile' in this.h)
         return func(this.h);

      if (this.h._childs)
         for (let n = 0; n < this.h._childs.length; ++n) {
            let item = this.h._childs[n];
            if ('_jsonfile' in item) func(item);
         }
   }

   /** @summary Open JSON file
     * @param {string} filepath - URL to JSON file
     * @returns {Promise} when object ready */
   HierarchyPainter.prototype.openJsonFile = function(filepath) {
      let isfileopened = false;
      this.forEachJsonFile(item => { if (item._jsonfile==filepath) isfileopened = true; });
      if (isfileopened) return Promise.resolve();

      return JSROOT.httpRequest(filepath, 'object').then(res => {
         let h1 = { _jsonfile: filepath, _kind: "ROOT." + res._typename, _jsontmp: res, _name: filepath.split("/").pop() };
         if (res.fTitle) h1._title = res.fTitle;
         h1._get = function(item /* ,itemname */) {
            if (item._jsontmp)
               return Promise.resolve(item._jsontmp);
            return JSROOT.httpRequest(item._jsonfile, 'object')
                         .then(res => {
                             item._jsontmp = res;
                             return res;
                          });
         }
         if (!this.h) this.h = h1; else
         if (this.h._kind == 'TopFolder') this.h._childs.push(h1); else {
            let h0 = this.h, topname = ('_jsonfile' in h0) ? "Files" : "Items";
            this.h = { _name: topname, _kind: 'TopFolder', _childs : [h0, h1] };
         }

         return this.refreshHtml();
      });
   }

   /** @summary Call function for each item which corresponds to ROOT file
     * @private */
   HierarchyPainter.prototype.forEachRootFile = function(func) {
      if (!this.h) return;
      if ((this.h._kind == "ROOT.TFile") && this.h._file)
         return func(this.h);

      if (this.h._childs)
         for (let n = 0; n < this.h._childs.length; ++n) {
            let item = this.h._childs[n];
            if ((item._kind == 'ROOT.TFile') && ('_fullurl' in item))
               func(item);
         }
   }

   /** @summary Open ROOT file
     * @param {string} filepath - URL to ROOT file
     * @returns {Promise} when file is opened */
   HierarchyPainter.prototype.openRootFile = function(filepath) {
      // first check that file with such URL already opened

      let isfileopened = false;
      this.forEachRootFile(item => { if (item._fullurl===filepath) isfileopened = true; });
      if (isfileopened) return Promise.resolve();

      jsrp.showProgress("Opening " + filepath + " ...");

      return JSROOT.openFile(filepath).then(file => {

         let h1 = this.fileHierarchy(file);
         h1._isopen = true;
         if (!this.h) {
            this.h = h1;
            if (this._topname) h1._name = this._topname;
         } else if (this.h._kind == 'TopFolder') {
            this.h._childs.push(h1);
         }  else {
            let h0 = this.h, topname = (h0._kind == "ROOT.TFile") ? "Files" : "Items";
            this.h = { _name: topname, _kind: 'TopFolder', _childs : [h0, h1], _isopen: true };
         }

         return this.refreshHtml();
      }).catch(() => {
         // make CORS warning
         if (!d3.select("#gui_fileCORS").style("background","red").empty())
             setTimeout(function() { d3.select("#gui_fileCORS").style("background",''); }, 5000);
         return false;
      }).finally(() => jsrp.showProgress());
   }

   /** @summary Apply loaded TStyle object
     * @desc One also can specify item name of JSON file name where style is loaded
     * @param {object|string} style - either TStyle object of item name where object can be load */
   HierarchyPainter.prototype.applyStyle = function(style) {
      if (!style)
         return Promise.resolve();

      if (typeof style === 'object') {
         if (style._typename === "TStyle")
            JSROOT.extend(JSROOT.gStyle, style);
         return Promise.resolve();
      }

      if (typeof style === 'string') {

         let item = this.findItem({ name: style, allow_index: true, check_keys: true });

         if (item!==null)
            return this.getObject(item).then(res => this.applyStyle(res.obj));

         if (style.indexOf('.json') > 0)
            return JSROOT.httpRequest(style, 'object')
                         .then(res => this.applyStyle(res));
      }

      return Promise.resolve();
   }

   /** @summary Provides information abouf file item
     * @private */
   HierarchyPainter.prototype.getFileProp = function(itemname) {
      let item = this.findItem(itemname);
      if (!item) return null;

      let subname = item._name;
      while (item._parent) {
         item = item._parent;
         if ('_file' in item)
            return { kind: "file", fileurl: item._file.fURL, itemname: subname, localfile: !!item._file.fLocalFile };

         if ('_jsonfile' in item)
            return { kind: "json", fileurl: item._jsonfile, itemname: subname };

         subname = item._name + "/" + subname;
      }

      return null;
   }

   /** @summary Provides URL for online item
     * @desc Such URL can be used  to request data from the server
     * @returns string or null if item is not online
     * @private */
   HierarchyPainter.prototype.getOnlineItemUrl = function(item) {
      if (typeof item == "string") item = this.findItem(item);
      let prnt = item;
      while (prnt && (prnt._online===undefined)) prnt = prnt._parent;
      return prnt ? (prnt._online + this.itemFullName(item, prnt)) : null;
   }

   /** @summary Returns true if item is online
     * @private */
   HierarchyPainter.prototype.isOnlineItem = function(item) {
      return this.getOnlineItemUrl(item) !== null;
   }

   /** @summary method used to request object from the http server
     * @returns {Promise} with requested object
     * @private */
   HierarchyPainter.prototype.getOnlineItem = function(item, itemname, option) {

      let url = itemname, h_get = false, req = "", req_kind = "object", draw_handle = null;

      if (option === 'hierarchy_expand') { h_get = true; option = undefined; }

      if (item) {
         url = this.getOnlineItemUrl(item);
         let func = null;
         if ('_kind' in item) draw_handle = jsrp.getDrawHandle(item._kind);

         if (h_get) {
            req = 'h.json?compact=3';
            item._expand = onlineHierarchy; // use proper expand function
         } else if ('_make_request' in item) {
            func = JSROOT.findFunction(item._make_request);
         } else if ((draw_handle!=null) && ('make_request' in draw_handle)) {
            func = draw_handle.make_request;
         }

         if (typeof func == 'function') {
            // ask to make request
            let dreq = func(this, item, url, option);
            // result can be simple string or object with req and kind fields
            if (dreq!=null)
               if (typeof dreq == 'string') req = dreq; else {
                  if ('req' in dreq) req = dreq.req;
                  if ('kind' in dreq) req_kind = dreq.kind;
               }
         }

         if ((req.length==0) && (item._kind.indexOf("ROOT.")!=0))
           req = 'item.json.gz?compact=3';
      }

      if (!itemname && item && ('_cached_draw_object' in this) && (req.length == 0)) {
         // special handling for online draw when cashed
         let obj = this._cached_draw_object;
         delete this._cached_draw_object;
         return Promise.resolve(obj);
      }

      if (req.length == 0) req = 'root.json.gz?compact=23';

      if (url.length > 0) url += "/";
      url += req;

      return new Promise(resolveFunc => {

         let itemreq = JSROOT.NewHttpRequest(url, req_kind, obj => {

            let func = null;

            if (!h_get && item && ('_after_request' in item)) {
               func = JSROOT.findFunction(item._after_request);
            } else if (draw_handle && ('after_request' in draw_handle))
               func = draw_handle.after_request;

            if (typeof func == 'function') {
               let res = func(this, item, obj, option, itemreq);
               if (res && (typeof res == "object")) obj = res;
            }

            resolveFunc(obj);
         });

         itemreq.send(null);
      });
   }

   /** @summary Access THttpServer with provided address
     * @param {string} server_address - URL to server like "http://localhost:8090/"
     * @returns {Promise} when ready */
   HierarchyPainter.prototype.openOnline = function(server_address) {
      let AdoptHierarchy = result => {
         this.h = result;
         if (!result) return Promise.resolve(null);

         if (('_title' in this.h) && (this.h._title!='')) document.title = this.h._title;

         result._isopen = true;

         // mark top hierarchy as online data and
         this.h._online = server_address;

         this.h._get = (item, itemname, option) => {
            return this.getOnlineItem(item, itemname, option);
         }

         this.h._expand = onlineHierarchy;

         let scripts = [], modules = [];
         this.forEachItem(item => {
            if ('_childs' in item) item._expand = onlineHierarchy;

            if ('_autoload' in item) {
               let arr = item._autoload.split(";");
               for (let n = 0; n < arr.length; ++n)
                  if ((arr[n].length>3) &&
                      ((arr[n].lastIndexOf(".js")==arr[n].length-3) ||
                      (arr[n].lastIndexOf(".css")==arr[n].length-4))) {
                     if (!scripts.find(elem => elem == arr[n])) scripts.push(arr[n]);
                  } else {
                     if (arr[n] && !modules.find(elem => elem ==arr[n])) modules.push(arr[n]);
                  }
            }
         });

         return JSROOT.require(modules)
               .then(() => JSROOT.loadScript(scripts))
               .then(() => {
                  this.forEachItem(item => {
                     if (!('_drawfunc' in item) || !('_kind' in item)) return;
                     let typename = "kind:" + item._kind;
                     if (item._kind.indexOf('ROOT.')==0) typename = item._kind.slice(5);
                     let drawopt = item._drawopt;
                     if (!jsrp.canDraw(typename) || drawopt)
                        jsrp.addDrawFunc({ name: typename, func: item._drawfunc, script: item._drawscript, opt: drawopt });
                  });

                  return this;
               });
      }

      if (!server_address) server_address = "";

      if (typeof server_address == 'object') {
         let h = server_address;
         server_address = "";
         return AdoptHierarchy(h);
      }

      return JSROOT.httpRequest(server_address + "h.json?compact=3", 'object')
                   .then(hh => AdoptHierarchy(hh));
   }

   /** @summary Get properties for online item  - server name and relative name
     * @private */
   HierarchyPainter.prototype.getOnlineProp = function(itemname) {
      let item = this.findItem(itemname);
      if (!item) return null;

      let subname = item._name;
      while (item._parent) {
         item = item._parent;

         if ('_online' in item) {
            return {
               server: item._online,
               itemname: subname
            };
         }
         subname = item._name + "/" + subname;
      }

      return null;
   }

   /** @summary Fill context menu for online item
     * @private */
   HierarchyPainter.prototype.fillOnlineMenu = function(menu, onlineprop, itemname) {

      let node = this.findItem(itemname),
          sett = jsrp.getDrawSettings(node._kind, 'nosame;noinspect'),
          handle = jsrp.getDrawHandle(node._kind),
          root_type = (typeof node._kind == 'string') ? node._kind.indexOf("ROOT.") == 0 : false;

      if (sett.opts && (node._can_draw !== false)) {
         sett.opts.push('inspect');
         menu.addDrawMenu("Draw", sett.opts, arg => this.display(itemname, arg));
      }

      if (!node._childs && (node._more !== false) && (node._more || root_type || sett.expand))
         menu.add("Expand", () => this.expandItem(itemname));

      if (handle && ('execute' in handle))
         menu.add("Execute", () => this.executeCommand(itemname, menu.tree_node));

      if (sett.opts && (node._can_draw !== false))
         menu.addDrawMenu("Draw in new window", sett.opts,
                           arg => window.open(onlineprop.server + "?nobrowser&item=" + onlineprop.itemname +
                                              (this.isMonitoring() ? "&monitoring=" + this.getMonitoringInterval() : "") +
                                              (arg ? "&opt=" + arg : "")));

      if (sett.opts && (sett.opts.length > 0) && root_type && (node._can_draw !== false))
         menu.addDrawMenu("Draw as png", sett.opts,
                           arg => window.open(onlineprop.server + onlineprop.itemname + "/root.png?w=600&h=400" + (arg ? "&opt=" + arg : "")));

      if ('_player' in node)
         menu.add("Player", () => this.player(itemname));
   }

   /** @summary Assign existing hierarchy to the painter and refresh HTML code
     * @private */
   HierarchyPainter.prototype.setHierarchy = function(h) {
      this.h = h;
      this.refreshHtml();
   }

   /** @summary Configures monitoring interval
     * @param {number} interval - repetition interval in ms
     * @param {boolean} flag - initial monitoring state */
   HierarchyPainter.prototype.setMonitoring = function(interval, monitor_on) {

      this._runMonitoring("cleanup");

      if (interval) {
         interval = parseInt(interval);
         if (!isNaN(interval) && (interval > 0)) {
            this._monitoring_interval = Math.max(100,interval);
            monitor_on = true;
         } else {
            this._monitoring_interval = 3000;
         }
      }

      this._monitoring_on = monitor_on;

      if (this.isMonitoring())
         this._runMonitoring();
   }

   /** @summary Runs monitoring event loop
     * @private */
   HierarchyPainter.prototype._runMonitoring = function(arg) {
      if ((arg == "cleanup") || !this.isMonitoring()) {
         if (this._monitoring_handle) {
            clearTimeout(this._monitoring_handle);
            delete this._monitoring_handle;
         }

         if (this._monitoring_frame) {
            cancelAnimationFrame(this._monitoring_frame);
            delete this._monitoring_frame;
         }
         return;
      }

      if (arg == "frame") {
         // process of timeout, request animation frame
         delete this._monitoring_handle;
         this._monitoring_frame = requestAnimationFrame(this._runMonitoring.bind(this,"draw"));
         return;
      }

      if (arg == "draw") {
         delete this._monitoring_frame;
         this.updateItems();
      }

      this._monitoring_handle = setTimeout(this._runMonitoring.bind(this,"frame"), this.getMonitoringInterval());
   }

   /** @summary Returns configured monitoring interval in ms */
   HierarchyPainter.prototype.getMonitoringInterval = function() {
      return this._monitoring_interval || 3000;
   }

   /** @summary Enable/disable monitoring
     * @param {boolean} on - if monitoring enabled */
   HierarchyPainter.prototype.enableMonitoring = function(on) {
      this.setMonitoring(undefined, on);
   }

   /** @summary Returns true when monitoring is enabled */
   HierarchyPainter.prototype.isMonitoring = function() {
      return this._monitoring_on;
   }

   /** @summary Assign default layout and place where drawing will be performed
     * @param {string} layout - layout like "simple" or "grid2x2"
     * @param {string} frameid - DOM element id where object drawing will be performed */
   HierarchyPainter.prototype.setDisplay = function(layout, frameid) {
      if (!frameid && (typeof layout == 'object')) {
         this.disp = layout;
         this.disp_kind = 'custom';
         this.disp_frameid = null;
      } else {
         this.disp_kind = layout;
         this.disp_frameid = frameid;
      }

      if (!this.register_resize) {
         this.register_resize = true;
         jsrp.registerForResize(this);
      }
   }

   /** @summary Returns configured layout */
   HierarchyPainter.prototype.getLayout = function() {
      return this.disp_kind;
   }

   /** @summary Remove painter reference from hierarhcy
     * @private */
   HierarchyPainter.prototype.removePainter = function(obj_painter) {
      this.forEachItem(item => {
         if (item._painter === obj_painter) {
            // delete painter reference
            delete item._painter;
            // also clear data which could be associated with item
            if (typeof item.clear == 'function') item.clear();
         }
      });
   }

   /** @summary Cleanup all items in hierarchy
     * @private */
   HierarchyPainter.prototype.clearHierarchy = function(withbrowser) {
      if (this.disp) {
         this.disp.cleanup();
         delete this.disp;
      }

      let plainarr = [];

      this.forEachItem(item => {
         delete item._painter; // remove reference on the painter
         // when only display cleared, try to clear all browser items
         if (!withbrowser && (typeof item.clear == 'function')) item.clear();
         if (withbrowser) plainarr.push(item);
      });

      if (withbrowser) {
         // cleanup all monitoring loops
         this.enableMonitoring(false);
         // simplify work for javascript and delete all (ok, most of) cross-references
         this.selectDom().html("");
         plainarr.forEach(d => { delete d._parent; delete d._childs; delete d._obj; delete d._d3cont; });
         delete this.h;
      }
   }

   /** @summary Returns actual MDI display object
     * @desc It should an instance of {@link JSROOT.MDIDsiplay} class */
   HierarchyPainter.prototype.getDisplay = function() {
      return this.disp;
   }

   /** @summary method called when MDI element is cleaned up
     * @desc hook to perform extra actions when frame is cleaned
     * @private */
   HierarchyPainter.prototype.cleanupFrame = function(divid) {

      let lst = JSROOT.cleanup(divid);

      // we remove all painters references from items
      if (lst && (lst.length > 0))
         this.forEachItem(item => {
            if (item._painter && lst.indexOf(item._painter) >= 0)
               delete item._painter;
         });
   }

   /** @summary Creates configured JSROOT.MDIDisplay object
     * @returns {Promise} when ready
     * @private */
   HierarchyPainter.prototype.createDisplay = function() {

      if ('disp' in this) {
         if ((this.disp.numDraw() > 0) || (this.disp_kind == "custom"))
            return Promise.resolve(this.disp);
         this.disp.cleanup();
         delete this.disp;
      }

      // check that we can found frame where drawing should be done
      if (!document.getElementById(this.disp_frameid))
         return Promise.resolve(null);

      if ((this.disp_kind == "simple") ||
          ((this.disp_kind.indexOf("grid") == 0) && (this.disp_kind.indexOf("gridi") < 0)))
           this.disp = new GridDisplay(this.disp_frameid, this.disp_kind);
      else
         return JSROOT.require('jq2d').then(() => this.createDisplay());

      if (this.disp)
         this.disp.cleanupFrame = this.cleanupFrame.bind(this);

      return Promise.resolve(this.disp);
   }

   /** @summary If possible, creates custom JSROOT.MDIDisplay for given item
     * @param itemname - name of item, for which drawing is created
     * @param custom_kind - display kind
     * @returns {Promise} with mdi object created
     * @private */
   HierarchyPainter.prototype.createCustomDisplay = function(itemname, custom_kind) {

      if (this.disp_kind != "simple")
         return this.createDisplay();

      this.disp_kind = custom_kind;

      // check if display can be erased
      if (this.disp) {
         let num = this.disp.numDraw();
         if ((num > 1) || ((num == 1) && !this.disp.findFrame(itemname)))
            return this.createDisplay();
         this.disp.cleanup();
         delete this.disp;
      }

      return this.createDisplay();
   }

   /** @summary function updates object drawings for other painters
     * @private */
   HierarchyPainter.prototype.updateOnOtherFrames = function(painter, obj) {
      let mdi = this.disp, handle = null, isany = false;
      if (!mdi) return false;

      if (obj._typename) handle = jsrp.getDrawHandle("ROOT." + obj._typename);
      if (handle && handle.draw_field && obj[handle.draw_field])
         obj = obj[handle.draw_field];

      mdi.forEachPainter((p, frame) => {
         if ((p === painter) || (p.getItemName() != painter.getItemName())) return;
         mdi.activateFrame(frame);
         if ((typeof p.redrawObject == 'function') && p.redrawObject(obj)) isany = true;
      });
      return isany;
   }

   /** @summary Process resize event
     * @private */
   HierarchyPainter.prototype.checkResize = function(size) {
      if (this.disp) this.disp.checkMDIResize(null, size);
   }

   /** @summary Start GUI
     * @returns {Promise} when ready
     * @private */
   HierarchyPainter.prototype.startGUI = function(gui_div, url) {

      let d = JSROOT.decodeUrl(url);

      function GetOption(opt) {
         let res = d.get(opt, null);
         if ((res===null) && gui_div && !gui_div.empty() && gui_div.node().hasAttribute(opt))
            res = gui_div.attr(opt);
         return res;
      }

      function GetUrlOptionAsArray(opt) {

         let res = [];

         while (opt.length > 0) {
            let separ = opt.indexOf(";");
            let part = (separ>0) ? opt.substr(0, separ) : opt;

            if (separ > 0) opt = opt.substr(separ+1); else opt = "";

            let canarray = true;
            if (part[0]=='#') { part = part.substr(1); canarray = false; }

            let val = d.get(part,null);

            if (canarray)
               res = res.concat(parseAsArray(val));
            else if (val!==null)
               res.push(val);
         }
         return res;
      }

      let GetOptionAsArray = (opt) => {
         let res = GetUrlOptionAsArray(opt);
         if (res.length>0 || !gui_div || gui_div.empty()) return res;
         while (opt.length>0) {
            let separ = opt.indexOf(";");
            let part = separ>0 ? opt.substr(0, separ) : opt;
            if (separ>0) opt = opt.substr(separ+1); else opt = "";

            let canarray = true;
            if (part[0]=='#') { part = part.substr(1); canarray = false; }
            if (part==='files') continue; // special case for normal UI

            if (!gui_div.node().hasAttribute(part)) continue;

            let val = gui_div.attr(part);

            if (canarray) res = res.concat(parseAsArray(val));
            else if (val!==null) res.push(val);
         }
         return res;
      }

      let prereq = GetOption('prereq') || "",
          filesdir = d.get("path") || "", // path used in normal gui
          filesarr = GetOptionAsArray("#file;files"),
          localfile = GetOption("localfile"),
          jsonarr = GetOptionAsArray("#json;jsons"),
          expanditems = GetOptionAsArray("expand"),
          itemsarr = GetOptionAsArray("#item;items"),
          optionsarr = GetOptionAsArray("#opt;opts"),
          monitor = GetOption("monitoring"),
          layout = GetOption("layout"),
          style = GetOptionAsArray("#style"),
          statush = 0, status = GetOption("status"),
          browser_kind = GetOption("browser"),
          browser_configured = !!browser_kind,
          title = GetOption("title");

      if (monitor === null)
         monitor = 0;
      else if (monitor === "")
         monitor = 3000;
      else
         monitor = parseInt(monitor);

      if (GetOption("float")!==null) { browser_kind = 'float'; browser_configured = true; } else
      if (GetOption("fix")!==null) { browser_kind = 'fix'; browser_configured = true; }

      this.no_select = GetOption("noselect");

      if (GetOption('files_monitoring')!==null) this.files_monitoring = true;

      if (title) document.title = title;

      let load = GetOption("load");

      if (expanditems.length==0 && (GetOption("expand")==="")) expanditems.push("");

      if (filesdir) {
         for (let i=0;i<filesarr.length;++i) filesarr[i] = filesdir + filesarr[i];
         for (let i=0;i<jsonarr.length;++i) jsonarr[i] = filesdir + jsonarr[i];
      }

      if ((itemsarr.length == 0) && GetOption("item") === "") itemsarr.push("");

      if ((jsonarr.length == 1) && (itemsarr.length == 0) && (expanditems.length==0)) itemsarr.push("");

      if (!this.disp_kind) {
         if ((typeof layout == "string") && (layout.length > 0))
            this.disp_kind = layout;
         else
         switch (itemsarr.length) {
           case 0:
           case 1: this.disp_kind = 'simple'; break;
           case 2: this.disp_kind = 'vert2'; break;
           case 3: this.disp_kind = 'vert21'; break;
           case 4: this.disp_kind = 'vert22'; break;
           case 5: this.disp_kind = 'vert32'; break;
           case 6: this.disp_kind = 'vert222'; break;
           case 7: this.disp_kind = 'vert322'; break;
           case 8: this.disp_kind = 'vert332'; break;
           case 9: this.disp_kind = 'vert333'; break;
           default: this.disp_kind = 'flex';
         }
      }

      if (status==="no")
         status = null;
      else if (status==="off") {
         this.status_disabled = true;
         status = null;
      } else if (status==="on")
         status = true;
      else if (status!==null) {
         statush = parseInt(status);
         if (isNaN(statush) || (statush < 5)) statush = 0;
         status = true;
      }
      if (this.no_select === "") this.no_select = true;

      if (!browser_kind)
         browser_kind = "fix";
      else if (browser_kind === "no")
         browser_kind = "";
      else if (browser_kind==="off") {
         browser_kind = "";
         status = null;
         this.exclude_browser = true;
      }
      if (GetOption("nofloat")!==null) this.float_browser_disabled = true;

      if (this.start_without_browser) browser_kind = "";

      if ((status || browser_kind) && !JSROOT.batch_mode) prereq += "jq2d;";

      this._topname = GetOption("topname");

      let openAllFiles = () => {
         let promise;

         if (prereq) {
            promise = JSROOT.require(prereq); prereq = "";
         } else if (load) {
            promise = JSROOT.loadScript(load.split(";")); load = "";
         } else if (browser_kind) {
            promise = this.createBrowser(browser_kind); browser_kind = "";
         } else if (status!==null) {
            promise = this.createStatusLine(statush, status); status = null;
         } else if (jsonarr.length > 0)
            promise = this.openJsonFile(jsonarr.shift());
         else if (filesarr.length > 0)
            promise = this.openRootFile(filesarr.shift());
         else if ((localfile!==null) && (typeof this.selectLocalFile == 'function')) {
            localfile = null; promise = this.selectLocalFile();
         } else if (expanditems.length > 0)
            promise = this.expandItem(expanditems.shift());
         else if (style.length > 0)
            promise = this.applyStyle(style.shift());
         else
            return this.refreshHtml()
                   .then(() => this.displayItems(itemsarr, optionsarr))
                   .then(() => {
                      this.setMonitoring(monitor);
                      return itemsarr ? this.refreshHtml() : this; // this is final return
                   });

         return promise.then(openAllFiles);
      }

      let h0 = null;
      if (this.is_online) {
         if (typeof GetCachedHierarchy == 'function')
            h0 = GetCachedHierarchy();
         if (typeof h0 !== 'object') h0 = "";

         if ((this.is_online == "draw") && !itemsarr.length)
            itemsarr.push("");
      }

      if (h0 !== null)
         return this.openOnline(h0).then(() => {
            // check if server enables monitoring
            if (!this.exclude_browser && !browser_configured && ('_browser' in this.h)) {
               browser_kind = this.h._browser;
               if (browser_kind==="no") browser_kind = ""; else
               if (browser_kind==="off") { browser_kind = ""; status = null; this.exclude_browser = true; }
            }

            if (('_monitoring' in this.h) && !monitor)
               monitor = this.h._monitoring;

            if (('_loadfile' in this.h) && (filesarr.length==0))
               filesarr = parseAsArray(this.h._loadfile);

            if (('_drawitem' in this.h) && (itemsarr.length==0)) {
               itemsarr = parseAsArray(this.h._drawitem);
               optionsarr = parseAsArray(this.h._drawopt);
            }

            if (('_layout' in this.h) && !layout && ((this.is_online != "draw") || (itemsarr.length > 1)))
               this.disp_kind = this.h._layout;

            if (('_toptitle' in this.h) && this.exclude_browser && document)
               document.title = this.h._toptitle;

            if (gui_div)
               this.prepareGuiDiv(gui_div, this.disp_kind);

            return openAllFiles();
         });

      if (gui_div)
         this.prepareGuiDiv(gui_div, this.disp_kind);

      return openAllFiles();
   }

   /** @summary Prepare div element - create layout and buttons
     * @private */
   HierarchyPainter.prototype.prepareGuiDiv = function(myDiv, layout) {

      this.gui_div = myDiv.attr('id');

      this.brlayout = new BrowserLayout(this.gui_div, this);

      this.brlayout.create(!this.exclude_browser);

      if (!this.exclude_browser) {
         let btns = this.brlayout.createBrowserBtns();

         JSROOT.require(['interactive']).then(inter => {
            inter.ToolbarIcons.createSVG(btns, inter.ToolbarIcons.diamand, 15, "toggle fix-pos browser")
                               .style("margin","3px").on("click", () => this.createBrowser("fix", true));

            if (!this.float_browser_disabled)
               inter.ToolbarIcons.createSVG(btns, inter.ToolbarIcons.circle, 15, "toggle float browser")
                                  .style("margin","3px").on("click", () => this.createBrowser("float", true));

            if (!this.status_disabled)
               inter.ToolbarIcons.createSVG(btns, inter.ToolbarIcons.three_circles, 15, "toggle status line")
                                  .style("margin","3px").on("click", () => this.createStatusLine(0, "toggle"));
          });
      }

      this.setDisplay(layout, this.brlayout.drawing_divid());
   }

   /** @summary Create status line
     * @param {number} [height] - size of the status line
     * @param [mode] - false / true / "toggle"
     * @returns {Promise} when ready */
   HierarchyPainter.prototype.createStatusLine = function(height, mode) {
      if (this.status_disabled || !this.gui_div || !this.brlayout)
         return Promise.resolve("");
      return this.brlayout.createStatusLine(height, mode);
   }

   /** @summary Create browser layout
     * @private */
   HierarchyPainter.prototype.createBrowser = function(browser_kind, update_html) {
      if (!this.gui_div)
         return Promise.resolve(false);

      return JSROOT.require('jq2d').then(() => this.createBrowser(browser_kind, update_html));
   }

   /** @summary Redraw hierarchy
     * @desc works only when inspector or streamer info is displayed
     * @private */
   HierarchyPainter.prototype.redrawObject = function(obj) {
      if (!this._inspector && !this._streamer_info) return false;
      if (this._streamer_info)
         this.h = createStreamerInfoContent(obj)
      else
         this.h = createInspectorContent(obj);
      return this.refreshHtml().then(() => { this.setTopPainter(); });
   }

   // ======================================================================================

   /** @summary tag item in hierarchy painter as streamer info
     * @desc this function used on THttpServer to mark streamer infos list
     * as fictional TStreamerInfoList class, which has special draw function
     * @private */
   JSROOT.markAsStreamerInfo = function(h,item,obj) {
      if (obj && (obj._typename=='TList'))
         obj._typename = 'TStreamerInfoList';
   }

   /** @summary Build gui without visisble hierarchy browser
     * @desc avoid loading of jquery part
     * @private */
   JSROOT.buildNobrowserGUI = function(gui_element, gui_kind) {

      let myDiv = (typeof gui_element == 'string') ? d3.select('#' + gui_element) : d3.select(gui_element);
      if (myDiv.empty()) return alert('no div for simple nobrowser gui found');

      let online = false, drawing = false;
      if (gui_kind == 'online')
         online = true;
      else if (gui_kind == 'draw')
         online = drawing = true;

      if (myDiv.attr("ignoreurl") === "true")
         JSROOT.settings.IgnoreUrlOptions = true;

      jsrp.readStyleFromURL();

      let d = JSROOT.decodeUrl(), guisize = d.get("divsize");
      if (guisize) {
         guisize = guisize.split("x");
         if (guisize.length != 2) guisize = null;
      }

      if (guisize) {
         myDiv.style('position',"relative").style('width', guisize[0] + "px").style('height', guisize[1] + "px");
      } else {
         d3.select('html').style('height','100%');
         d3.select('body').style('min-height','100%').style('margin',0).style('overflow',"hidden");
         myDiv.style('position',"absolute").style('left',0).style('top',0).style('bottom',0).style('right',0).style('padding',1);
      }

      let hpainter = new HierarchyPainter('root', null);

      if (online) hpainter.is_online = drawing ? "draw" : "online";
      if (drawing) hpainter.exclude_browser = true;

      hpainter.start_without_browser = true; // indicate that browser not required at the beginning

      return hpainter.startGUI(myDiv, () => {
         if (!drawing) return hpainter;
         let func = JSROOT.findFunction('GetCachedObject');
         let obj = (typeof func == 'function') ? JSROOT.parse(func()) : null;
         if (obj) hpainter._cached_draw_object = obj;
         let opt = d.get("opt", "");

         if (d.has("websocket")) opt+=";websocket";
         console.log('try to draw first object');

         return hpainter.display("", opt).then(() => { return hpainter; });
      });
   }

   /** @summary Display streamer info
     * @private */
   jsrp.drawStreamerInfo = function(divid, lst) {
      let painter = new HierarchyPainter('sinfo', divid, 'white');

      painter._streamer_info = true;
      painter.h = createStreamerInfoContent(lst);

      // painter.selectDom().style('overflow','auto');

      return painter.refreshHtml().then(() => {
         painter.setTopPainter();
         return painter;
      });
   }

   // ======================================================================================

   /** @summary Display inspector
     * @private */
   jsrp.drawInspector = function(divid, obj) {

      JSROOT.cleanup(divid);
      let painter = new HierarchyPainter('inspector', divid, 'white');

      // in batch mode HTML drawing is not possible, just keep object reference for a minute
      if (JSROOT.batch_mode) {
         painter.selectDom().property("_json_object_", obj);
         return Promise.resolve(painter);
      }

      painter.default_by_click = "expand"; // default action
      painter.with_icons = false;
      painter._inspector = true; // keep

      if (painter.selectDom().classed("jsroot_inspector"))
         painter.removeInspector = function() {
            this.selectDom().remove();
         }

      painter.fill_context = function(menu, hitem) {
         let sett = jsrp.getDrawSettings(hitem._kind, 'nosame');
         if (sett.opts)
            menu.addDrawMenu("nosub:Draw", sett.opts, function(arg) {
               if (!hitem || !hitem._obj) return;
               let obj = hitem._obj, dom = this.selectDom();
               if (this.removeInspector) {
                  dom = dom.node().parentNode;
                  this.removeInspector();
                  if (arg == "inspect")
                     return this.showInspector(obj);
               }
               JSROOT.cleanup(dom);
               JSROOT.draw(dom, obj, arg);
            });
      }

      painter.h = createInspectorContent(obj);

      return painter.refreshHtml().then(() => {
         painter.setTopPainter();
         return painter;
      });
   }

   // ================================================================

   /**
    * @summary Base class to manage multiple document interface for drawings
    *
    * @class
    * @memberof JSROOT
    * @extends JSROOT.BasePainter
    * @private
    */

   class MDIDisplay extends JSROOT.BasePainter {
      constructor(frameid) {
         super();
         this.frameid = frameid;
         if (frameid != "$batch$") {
            this.setDom(frameid);
            this.selectDom().property('mdi', this);
         }
         this.cleanupFrame = JSROOT.cleanup; // use standard cleanup function by default
         this.active_frame_title = ""; // keep title of active frame
      }

      /** @summary method called before new frame is created */
      beforeCreateFrame(title) { this.active_frame_title = title; }

      /** @summary method dedicated to iterate over existing panels
        * @param {function} userfunc is called with arguments (frame)
        * @param {boolean} only_visible let select only visible frames */
      forEachFrame(/* userfunc, only_visible */) { console.warn("forEachFrame not implemented in MDIDisplay"); }

      /** @summary method dedicated to iterate over existing panles
        * @param {function} userfunc is called with arguments (painter, frame)
        * @param {boolean} only_visible let select only visible frames */
      forEachPainter(userfunc, only_visible) {
         this.forEachFrame(frame => {
            let dummy = new JSROOT.ObjectPainter(frame);
            dummy.forEachPainter(painter => userfunc(painter, frame));
         }, only_visible);
      }

      /** @summary Returns total number of drawings */
      numDraw() {
         let cnt = 0;
         this.forEachFrame(() => ++cnt);
         return cnt;
      }

      /** @summary Serach for the frame using item name */
      findFrame(searchtitle, force) {
         let found_frame = null;

         this.forEachFrame(frame => {
            if (d3.select(frame).attr('frame_title') == searchtitle)
               found_frame = frame;
         });

         if (!found_frame && force)
            found_frame = this.createFrame(searchtitle);

         return found_frame;
      }

      /** @summary Activate frame */
      activateFrame(frame) { this.active_frame_title = d3.select(frame).attr('frame_title'); }

      /** @summary Return active frame */
      getActiveFrame() { return this.findFrame(this.active_frame_title); }

      /** @summary perform resize for each frame
        * @protected */
      checkMDIResize(only_frame_id, size) {

         let resized_frame = null;

         this.forEachPainter((painter, frame) => {

            if (only_frame_id && (d3.select(frame).attr('id') != only_frame_id)) return;

            if ((painter.getItemName()!==null) && (typeof painter.checkResize == 'function')) {
               // do not call resize for many painters on the same frame
               if (resized_frame === frame) return;
               painter.checkResize(size);
               resized_frame = frame;
            }
         });
      }

      /** @summary Cleanup all drawings */
      cleanup() {
         this.active_frame_title = "";

         this.forEachFrame(this.cleanupFrame);

         this.selectDom().html("").property('mdi', null);
      }

   } // class MDIDisplay


   // ==================================================

   /**
    * @summary Custom MDI display
    *
    * @class
    * @memberof JSROOT
    * @desc All HTML frames should be created before and add via {@link CustomDisplay.addFrame} calls
    * @private
    */

   class CustomDisplay extends MDIDisplay {
      constructor() {
         super("dummy");
         this.frames = {}; // array of configured frames
      }

      addFrame(divid, itemname) {
         if (!(divid in this.frames)) this.frames[divid] = "";

         this.frames[divid] += (itemname + ";");
      }

      forEachFrame(userfunc /* ,  only_visible */) {
         let ks = Object.keys(this.frames);
         for (let k = 0; k < ks.length; ++k) {
            let node = d3.select("#"+ks[k]);
            if (!node.empty())
               userfunc(node.node());
         }
      }

      createFrame(title) {
         this.beforeCreateFrame(title);

         let ks = Object.keys(this.frames);
         for (let k = 0; k < ks.length; ++k) {
            let items = this.frames[ks[k]];
            if (items.indexOf(title+";") >= 0)
               return d3.select("#"+ks[k]).node();
         }
         return null;
      }

      cleanup() {
         super.cleanup();
         this.forEachFrame(frame => d3.select(frame).html(""));
      }

   } // class CustomDisplay

   // ================================================

   /**
    * @summary Generic grid MDI display
    *
    * @class
    * @memberof JSROOT
    * @private
    */


   class GridDisplay extends MDIDisplay {

    /** @summary Create GridDisplay instance
      * @param {string} frameid - where grid display is created
      * @param {string} kind - kind of grid
      * @desc  following kinds are supported
      *    - vertical or horizontal - only first letter matters, defines basic orientation
      *    - 'x' in the name disable interactive separators
      *    - v4 or h4 - 4 equal elements in specified direction
      *    - v231 -  created 3 vertical elements, first divided on 2, second on 3 and third on 1 part
      *    - v23_52 - create two vertical elements with 2 and 3 subitems, size ratio 5:2
      *    - gridNxM - normal grid layout without interactive separators
      *    - gridiNxM - grid layout with interactive separators
      *    -  simple - no layout, full frame used for object drawings */

      constructor(frameid, kind, kind2) {

         super(frameid);

         this.framecnt = 0;
         this.getcnt = 0;
         this.groups = [];
         this.vertical = kind && (kind[0] == 'v');
         this.use_separarators = !kind || (kind.indexOf("x")<0);
         this.simple_layout = false;

         this.selectDom().style('overflow','hidden');

         if (kind === "simple") {
            this.simple_layout = true;
            this.use_separarators = false;
            this.framecnt = 1;
            return;
         }

         let num = 2, arr = undefined, sizes = undefined;

         if ((kind.indexOf("grid") == 0) || kind2) {
            if (kind2) kind = kind + "x" + kind2;
                  else kind = kind.substr(4).trim();
            this.use_separarators = false;
            if (kind[0] === "i") {
               this.use_separarators = true;
               kind = kind.substr(1);
            }

            let separ = kind.indexOf("x"), sizex = 3, sizey = 3;

            if (separ > 0) {
               sizey = parseInt(kind.substr(separ + 1));
               sizex = parseInt(kind.substr(0, separ));
            } else {
               sizex = sizey = parseInt(kind);
            }

            if (isNaN(sizex)) sizex = 3;
            if (isNaN(sizey)) sizey = 3;

            if (sizey > 1) {
               this.vertical = true;
               num = sizey;
               if (sizex>1) {
                  arr = new Array(num);
                  for (let k = 0; k < num; ++k) arr[k] = sizex;
               }
            } else if (sizex > 1) {
               this.vertical = false;
               num = sizex;
            } else {
               this.simple_layout = true;
               this.use_separarators = false;
               this.framecnt = 1;
               return;
            }
            kind = "";
         }

         if (kind && kind.indexOf("_")>0) {
            let arg = parseInt(kind.substr(kind.indexOf("_")+1), 10);
            if (!isNaN(arg) && (arg>10)) {
               kind = kind.substr(0, kind.indexOf("_"));
               sizes = [];
               while (arg>0) {
                  sizes.unshift(Math.max(arg % 10, 1));
                  arg = Math.round((arg-sizes[0])/10);
                  if (sizes[0]===0) sizes[0]=1;
               }
            }
         }

         kind = kind ? parseInt(kind.replace( /^\D+/g, ''), 10) : 0;
         if (kind && (kind>1)) {
            if (kind<10) {
               num = kind;
            } else {
               arr = [];
               while (kind>0) {
                  arr.unshift(kind % 10);
                  kind = Math.round((kind-arr[0])/10);
                  if (arr[0]==0) arr[0]=1;
               }
               num = arr.length;
            }
         }

         if (sizes && (sizes.length!==num)) sizes = undefined;

         if (!this.simple_layout)
            this.createGroup(this, this.selectDom(), num, arr, sizes);
      }

      createGroup(handle, main, num, childs, sizes) {

         if (!sizes) sizes = new Array(num);
         let sum1 = 0, sum2 = 0;
         for (let n=0;n<num;++n) sum1 += (sizes[n] || 1);
         for (let n=0;n<num;++n) {
            sizes[n] = Math.round(100 * (sizes[n] || 1) / sum1);
            sum2 += sizes[n];
            if (n==num-1) sizes[n] += (100-sum2); // make 100%
         }

         for (let cnt = 0; cnt<num; ++cnt) {
            let group = { id: cnt, drawid: -1, position: 0, size: sizes[cnt] };
            if (cnt>0) group.position = handle.groups[cnt-1].position + handle.groups[cnt-1].size;
            group.position0 = group.position;

            if (!childs || !childs[cnt] || childs[cnt]<2) group.drawid = this.framecnt++;

            handle.groups.push(group);

            let elem = main.append("div").attr('groupid', group.id);

            if (handle.vertical)
               elem.style('float', 'bottom').style('height',group.size+'%').style('width','100%');
            else
               elem.style('float', 'left').style('width',group.size+'%').style('height','100%');

            if (group.drawid>=0) {
               elem.classed('jsroot_newgrid', true);
               if (typeof this.frameid === 'string')
                  elem.attr('id', this.frameid + "_" + group.drawid);
            } else {
               elem.style('display','flex').style('flex-direction', handle.vertical ? "row" : "column");
            }

            if (childs && (childs[cnt]>1)) {
               group.vertical = !handle.vertical;
               group.groups = [];
               elem.style('overflow','hidden');
               this.createGroup(group, elem, childs[cnt]);
            }
         }

         if (this.use_separarators && this.createSeparator)
            for (let cnt=1;cnt<num;++cnt)
               this.createSeparator(handle, main, handle.groups[cnt]);
      }

      forEachFrame(userfunc /*, only_visible */) {
         if (this.simple_layout)
            userfunc(this.getGridFrame());
         else
            this.selectDom().selectAll('.jsroot_newgrid').each(function() {
               userfunc(d3.select(this).node());
            });
      }

      /** @summary Returns active frame */
      getActiveFrame() {
         if (this.simple_layout) return this.getGridFrame();

         let found = super.getActiveFrame();
         if (found) return found;

         this.forEachFrame(frame => { if (!found) found = frame; }, true);

         return found;
      }

      /** @summary Returns number of frames in grid layout */
      numGridFrames() { return this.framecnt; }

      /** @summary Return grid frame by its id */
      getGridFrame(id) {
         if (this.simple_layout)
            return this.selectDom('origin').node();
         let res = null;
         this.selectDom().selectAll('.jsroot_newgrid').each(function() {
            if (id-- === 0) res = this;
         });
         return res;
      }

      /** @summary Create new frame */
      createFrame(title) {
         this.beforeCreateFrame(title);

         let frame = null, maxloop = this.framecnt || 2;

         while (!frame && maxloop--) {
            frame = this.getGridFrame(this.getcnt);
            if (!this.simple_layout && this.framecnt)
               this.getcnt = (this.getcnt+1) % this.framecnt;

            if (d3.select(frame).classed("jsroot_fixed_frame")) frame = null;
         }

         if (frame) {
            this.cleanupFrame(frame);
            d3.select(frame).attr('frame_title', title);
         }

         return frame;
      }

   } // class GridDisplay

   // ==================================================

   /**
    * @summary Batch MDI display
    *
    * @class
    * @memberof JSROOT
    * @desc Can be used together with hierarchy painter in node.js
    * @private
    */

   class BatchDisplay extends MDIDisplay {
      constructor(width, height) {
         super("$batch$");
         this.frames = []; // array of configured frames
         this.width = width || 1200;
         this.height = height || 800;
      }

      forEachFrame(userfunc) {
         this.frames.forEach(userfunc)
      }

      createFrame(title) {
         this.beforeCreateFrame(title);

         let frame;
         if (!JSROOT.nodejs) {
            frame = d3.select('body').append("div").style("visible", "hidden");
         } else {
            if (!JSROOT._.nodejs_document) {
             // use eval while old minifier is not able to parse newest Node.js syntax
               const { JSDOM } = require("jsdom");
              JSROOT._.nodejs_window = (new JSDOM("<!DOCTYPE html>hello")).window;
              JSROOT._.nodejs_document = JSROOT._.nodejs_window.document; // used with three.js
              JSROOT._.nodejs_window.d3 = d3.select(JSROOT._.nodejs_document); //get d3 into the dom
            }
            frame = JSROOT._.nodejs_window.d3.select('body').append('div');
         }

         if (this.frames.length == 0)
            JSROOT._.svg_3ds = undefined;

         frame.attr("width", this.width).attr("height", this.height);
         frame.style("width", this.width + "px").style("height", this.height + "px");
         frame.attr("id","jsroot_batch_" + this.frames.length);
         frame.attr('frame_title', title);
         this.frames.push(frame.node());

         return frame.node();
      }

      /** @summary Returns number of created frames */
      numFrames() { return this.frames.length; }

      /** @summary returns JSON representation if any
        * @desc Now works only for inspector, can be called once */
      makeJSON(id, spacing) {
         let frame = this.frames[id];
         if (!frame) return;
         let obj = d3.select(frame).property('_json_object_');
         if (obj) {
            d3.select(frame).property('_json_object_', null);
            return JSROOT.toJSON(obj, spacing);
         }
      }

      /** @summary Create SVG for specified frame id */
      makeSVG(id) {
         let frame = this.frames[id];
         if (!frame) return;
         let main = d3.select(frame);
         let has_workarounds = JSROOT._.svg_3ds && jsrp.processSvgWorkarounds;
         main.select('svg')
             .attr("xmlns", "http://www.w3.org/2000/svg")
             .attr("xmlns:xlink", "http://www.w3.org/1999/xlink")
             .attr("width", this.width)
             .attr("height", this.height)
             .attr("title", null).attr("style", null).attr("class", null).attr("x", null).attr("y", null);

         let svg = main.html();
         if (has_workarounds)
            svg = jsrp.processSvgWorkarounds(svg, id != this.frames.length-1);

         svg = jsrp.compressSVG(svg);

         main.remove();
         return svg;
      }
   } // class BatchDisplay


   // export all functions and classes

   jsrp.drawList = drawList;

   jsrp.folderHierarchy = folderHierarchy;
   jsrp.taskHierarchy = taskHierarchy;
   jsrp.listHierarchy = listHierarchy;
   jsrp.objectHierarchy = objectHierarchy;
   jsrp.keysHierarchy = keysHierarchy;

   JSROOT.BrowserLayout = BrowserLayout;
   JSROOT.HierarchyPainter = HierarchyPainter;

   JSROOT.MDIDisplay = MDIDisplay;
   JSROOT.CustomDisplay = CustomDisplay;
   JSROOT.GridDisplay = GridDisplay;
   JSROOT.BatchDisplay = BatchDisplay;

   return JSROOT;

});
