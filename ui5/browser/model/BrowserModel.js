sap.ui.define([
    "sap/ui/model/json/JSONModel",
    "rootui5/browser/model/BrowserListBinding"
], function(JSONModel, BrowserListBinding) {
   "use strict";

    let hRootModel = JSONModel.extend("rootui5.browser.model.BrowserModel", {

        constructor: function() {
            JSONModel.apply(this);
            this.setProperty("/", {
                nodes: {}, // nodes shown in the TreeTable, should be flat list
                length: 0  // total number of elements
            });

            // this is true hierarchy, created on the client side and used for creation of flat list
            this.h = {
               name: "__Holder__",
               expanded: true
            };

            this.useIndexSuffix = true; // use index suffix in path for handling name duplication

            this.loadDataCounter = 0; // counter of number of nodes

            this.sortMethod = "name"; // "name", "size"
            this.reverseOrder = false;
            this.itemsFilter = "";
            this.showHidden = false;

            this.threshold = 100; // default threshold to prefetch items
        },

        assignTreeTable: function(t) {
           this.treeTable = t;
        },

        /** @summary Set sort method */
        setSortMethod: function(arg) { this.sortMethod = arg; },

        /** @summary Get sort method */
        getSortMethod: function() { return this.sortMethod; },

        /** @summary Set show hidden flag */
        setShowHidden: function(flag) { this.showHidden = flag; },

        /** @summary Is show hidden flag set */
        isShowHidden: function() { return this.showHidden; },

        /** @summary Set reverse order */
        setReverseOrder: function(on) { this.reverseOrder = on; },

        /** @summary Is reverse order */
        isReverseOrder: function() { return this.reverseOrder; },

        /** @summary Set full model
          * @desc Method can be used when complete hierarchy is ready and can be used directly */
        setFullModel: function(topnode) {
           this.fullModel = true;
           if (topnode.length) {
              this.h.nchilds = topnode.length;
              this.h.childs = topnode;
           } else {
              this.h.nchilds = 1;
              this.h.childs = [ topnode ];
           }
           delete this.h._requested; // reply on top element can be full description

           if (!this.mainModel) {
              this.mainModel = this.h.childs;
              this.mainFullModel = true;
           }

           topnode.expanded = true;
           this.reset_nodes = true;
           delete this.noData;
           this.scanShifts();
           if (this.oBinding)
              this.oBinding.checkUpdate(true);
        },

        /** @summary Clear full model */
        clearFullModel: function() {
           if (!this.fullModel) return;

           delete this.h.childs;
           delete this.h.nchilds;
           delete this.fullModel;
           delete this.noData;
           this.reset_nodes = true;
           if (this.oBinding)
              this.oBinding.checkUpdate(true);
        },

        setNoData: function(on) {
           this.noData = on;
           if (this.oBinding)
              this.oBinding.checkUpdate(true);
        },

        bindTree: function(sPath, oContext, aFilters, mParameters, aSorters) {
           console.log('BINDING TREE!!!!!!!!!!!!! ' + sPath);

           this.oBinding = new BrowserListBinding(this, sPath, oContext, aFilters, mParameters, aSorters);
           return this.oBinding;
        },

        getLength: function() {
           if (this.noData) return 0;
           return this.getProperty("/length");
        },

        /** @summary Code index as string */
        codeIndex: function(indx) {
           return this.useIndexSuffix ? "###" + indx + "$$$" : "";
        },

        /** @summary Extract index from string */
        extractIndex: function(name) {
           if (!name || typeof name != "string") return "";
           let p1 = name.lastIndexOf("###"), p2 = name.lastIndexOf("$$$");
           if ((p1 < 0) || (p2 < 0) || (p1 >= p2) || (p2 != name.length - 3)) return name;

           let indx = parseInt(name.substr(p1+3, p2-p1-3));
           if (isNaN(indx) || (indx < 0)) return name;
           return { n: name.substr(0, p1), i: indx };
        },

        /** @summary Get node by path which is array of strings, optionally includes indicies */
        getNodeByPath: function(path) {
           let curr = this.h;
           if (!path || (path.length == 0)) return curr;

           let names = path.slice(); // make copy to avoid changes in argument

           while (names.length > 0) {
              let name = this.extractIndex(names.shift()), find = false, indx = -1;

              if (typeof name != "string") {
                 indx = name.i;
                 name = name.n;
              }

              if (indx >= 0) {
                 indx -= (curr.first || 0);
                 if (curr.childs[indx] && (curr.childs[indx].name == name)) {
                    curr = curr.childs[indx];
                    find = true;
                    console.log(`Match index ${indx} with name ${name}`);
                 }
              }

              for (let k = 0; !find && (k < curr.childs.length); ++k) {
                 if (curr.childs[k].name == name) {
                    curr = curr.childs[k];
                    find = true;
                 }
              }

              if (!find) return null;
           }
           return curr;
        },

        /** @summary expand node by given path
          * @desc When path not exists - try to send request */
        expandNodeByPath: function(path) {
           if (!path || (typeof path !== "string") || (path == "/")) return -1;

           let names = path.split("/"), curr = this.h, currpath = [];

           while (names.length > 0) {
              let name = names.shift(), find = false;
              if (!name) continue; // ignore start or stop slash

              if (!curr.childs) {
                 // request childs for current element
                 // TODO: we do not know child index, but simply can suply search child as argument
                 if (!this.fullModel && curr.nchilds && (curr.nchilds > 0)) {
                    curr.expanded = true;
                    this.reset_nodes = true;
                    this._expanding_path = path;
                    this.submitRequest(false, curr, currpath, "expanding");
                    break;
                 }
                 return -1;
              }

              for (let k=0;k<curr.childs.length;++k) {
                 if (curr.childs[k].name == name) {
                    this.reset_nodes = true;
                    curr.expanded = true;
                    curr = curr.childs[k];
                    find = true;
                    break;
                 }
              }

              if (!find) return -1;

              currpath.push(curr.name);
           }

           return this.scanShifts(curr);
        },

        sendFirstRequest: function(websocket) {
           this._websocket = websocket;
           // submit top-level request already when construct model
           this.submitRequest(false, this.h);
        },

        /** @summary Reload main model
          * One can force to submit new request while settings were changed
          * One also can force to reload items on the server side if they can be potentially changed */
        reloadMainModel: function(force_request, force_reload, path) {
           if (this.mainModel && !force_request && !force_reload) {
              this.h.nchilds = this.mainModel.length;
              this.h.childs = this.mainModel;
              this.h.expanded = true;
              this.reset_nodes = true;
              this.fullModel = this.mainFullModel;
              delete this.noData;
              this.scanShifts();
              if (this.oBinding)
                 this.oBinding.checkUpdate(true);
           } else if (!this.fullModel) {
              // send request, content will be reassigned
              this.submitRequest(force_reload, this.h, path);
           }

        },

        /** @summary submit next request to the server
          * @param {Array} path - path as array of strings
          * @desc directly use web socket, later can be dedicated channel */
        submitRequest: function(force_reload, elem, path, first, number) {
           if (first === "expanding") {
              first = 0;
           } else {
              delete this._expanding_path;
           }

           if (!this._websocket || elem._requested || this.fullModel) return;
           elem._requested = true;

           this.loadDataCounter++;

           let request = {
              path: path || [],
              first: first || 0,
              number: number || this.threshold || 100,
              sort: this.sortMethod || "",
              reverse: this.reverseOrder || false,
              hidden: this.showHidden ? true : false,
              reload: force_reload ? true : false,  // rescan items by server even when path was not changed
              regex: this.itemsFilter ? "^(" + this.itemsFilter + ".*)$" : ""
           };
           this._websocket.send("BRREQ:" + JSON.stringify(request));
        },

        // process reply from server
        // In the future reply will be received via websocket channel
        processResponse: function(reply) {

           this.loadDataCounter--;

           let elem = this.getNodeByPath(reply.path);

           if (!elem) { console.error('DID NOT FOUND ' + reply.path); return; }

           if (!elem._requested) console.error('ELEMENT WAS NOT REQUESTED!!!!', reply.path);
           delete elem._requested;

           let smart_merge = false;

           if ((elem.nchilds === reply.nchilds) && elem.childs && reply.nodes) {
              if (elem.first + elem.childs.length == reply.first) {
                 elem.childs = elem.childs.concat(reply.nodes);
                 smart_merge = true;
              } else if (reply.first + reply.nodes.length == elem.first) {
                 elem.first = reply.first;
                 elem.childs = reply.nodes.concat(elem.childs);
                 smart_merge = true;
              }
           }

           if (!smart_merge) {
              elem.nchilds = reply.nchilds;
              elem.childs = reply.nodes;
              elem.first = reply.first || 0;
           }

           // remember main model
           if ((reply.path.length === 0) && !this.mainModel) {
              this.mainModel = elem.childs;
              this.mainFullModel = false;
           }

           this.scanShifts();

           // reset existing nodes if reply does not match with expectation
           if (!smart_merge)
              this.reset_nodes = true;

           if (this.loadDataCounter == 0)
              if (this.oBinding)
                 this.oBinding.checkUpdate(true);

           if (this._expanding_path) {
              let d = this._expanding_path;
              delete this._expanding_path;
              let index = this.expandNodeByPath(d);
              if ((index > 0) && this.treeTable)
                 this.treeTable.setFirstVisibleRow(Math.max(0, index - Math.round(this.treeTable.getVisibleRowCount()/2)));
           }

        },

        getNodeByIndex: function(indx) {
           let nodes = this.getProperty("/nodes");
           return nodes ? nodes[indx] : null;
        },

        // return element of hierarchical structure by TreeTable index
        getElementByIndex: function(indx) {
           let node = this.getNodeByIndex(indx);
           return node ? node._elem : null;
        },

        // function used to calculate all ids shifts and total number of elements
        // if element specified - returns index of that element
        scanShifts: function(for_elem) {

           let id = 0, full = this.fullModel, res = -1;

           function scan(lvl, elem) {

              if (elem === for_elem) res = id;

              if (lvl >= 0) id++;

              let before_id = id;

              if (elem.expanded) {
                 if (elem.childs === undefined) {
                    // do nothing, childs are not visible as long as we do not have any list

                    // id += 0;
                 } else {

                    // gap at the begin
                    if (!full && elem.first)
                       id += elem.first;

                    // jump over all childs
                    for (let k = 0; k < elem.childs.length; ++k)
                       scan(lvl+1, elem.childs[k]);

                    // gap at the end
                    if (!full) {
                       let _last = (elem.first || 0) + elem.childs.length,
                           _remains = elem.nchilds  - _last;
                       if (_remains > 0) id += _remains;
                    }
                 }
              }

              // this shift can be later applied to jump over all elements
              elem._shift = id - before_id;
           }

           scan(-1, this.h);

           this.setProperty("/length", id);

           return for_elem ? res : id;
        },

        // main  method to create flat list of nodes - only whose which are specified in selection
        // following arguments:
        //    args.begin     - first visisble element from flat list
        //    args.end       - first not-visisble element
        //    args.threshold - extra elements (before/after) which probably should be prefetched
        // returns holder object with all existing nodes
        buildFlatNodes: function(args) {

           if (this.noData) return null;

           let id = 0,            // current id, at the same time number of items
               threshold = args.threshold || this.threshold || 100,
               threshold2 = Math.round(threshold/2),
               nodes = this.reset_nodes ? {} : this.getProperty("/nodes");

           // main method to scan through all existing sub-folders
           let scan = (lvl, elem, path) => {

              // create elements with safety margin
              if ((lvl >= 0) && (nodes !== null) && !nodes[id] && (id >= args.begin - threshold2) && (id < args.end + threshold2)) {
                 nodes[id] = {
                    name: elem.name,
                    path: path.slice(), // make array copy
                    index: id,
                    _elem: elem,
                    isLeaf: !elem.nchilds,
                    // these are required by list binding, should be eliminated in the future
                    type: elem.nchilds ? "folder" : "file",
                    level: lvl,
                    context: this.getContext("/nodes/" + id),
                    nodeState: {
                       expanded: !!elem.expanded,
                       selected: !!elem.selected,
                       sum: false // ????
                    }
                 };
                 if (typeof this.addNodeAttributes == 'function')
                    this.addNodeAttributes(nodes[id], elem);
              }

              if (lvl >= 0) id++;

              if (!elem.expanded) return;

              if (elem.childs === undefined) {
                 // add new request - can we check if only special part of childs is required?

                 // TODO: probably one could guess more precise request
                 if ((elem.nchilds === undefined) || (elem.nchilds !== 0))
                   this.submitRequest(false, elem, path);

                 return;
              }

              // check if scan is required
              if (((id + elem._shift) < args.begin - threshold2) || (id >= args.end + threshold2)) {
                 id += elem._shift;
                 return;
              }

              // when not all childs from very beginning is loaded, but may be required
              if (elem.first && !this.fullModel) {

                 // check if requests are needed to load part in the begin of the list
                 if (args.begin - id - threshold2 < elem.first) {

                    let first = Math.max(args.begin - id - threshold2, 0),
                        number = Math.min(elem.first - first, threshold);

                    this.submitRequest(false, elem, path, first, number);
                 }

                 id += elem.first;
              }


              let subpath = path.slice(), // make copy to avoid changes in argument
                  subindx = subpath.push("") - 1;
              for (let k = 0; k < elem.childs.length; ++k) {
                 subpath[subindx] = elem.childs[k].name + this.codeIndex((elem.first || 0) + k);
                 scan(lvl+1, elem.childs[k], subpath);
              }

              // check if more elements are required

              if (!this.fullModel) {
                 let _last = (elem.first || 0) + elem.childs.length;
                 let _remains = elem.nchilds  - _last;

                 if (_remains > 0) {
                    if (args.end + threshold2 > id) {

                       let first = _last, number = args.end + threshold2 - id;
                       if (number < threshold) number = threshold; // always request much
                       if (number > _remains) number = _remains; // but not too much
                       if (number > threshold) {
                          first += (number - threshold);
                          number = threshold;
                       }

                       this.submitRequest(false, elem, path, first, number);
                    }

                    id += _remains;
                 }
              }
           }

           // start scan from very top
           scan(-1, this.h, []);

           if (this.getProperty("/length") != id) {
              // console.error('LENGTH MISMATCH', this.getProperty("/length"), id);
              this.setProperty("/length", id); // update length property
           }

           if (this.reset_nodes) {
              this.setProperty("/nodes", nodes);
              delete this.reset_nodes;
           }

           return nodes;
        },

        // toggle expand state of specified node
        toggleNode: function(index) {

           let node = this.getNodeByIndex(index),
               elem = node ? node._elem : null;

           if (!node || !elem) return;

           if (elem.expanded) {
              delete elem.expanded;
              if (!this.fullModel)
                 delete elem.childs; // TODO: for the future keep childs but make request if expand once again

              // close folder - reassign shifts
              this.reset_nodes = true;
              this.scanShifts();

              return true;

           } else if (elem.nchilds || !elem.index) {

              elem.expanded = true;
              // structure is changing but not immediately

              if (this.fullModel) {
                 this.reset_nodes = true;
                 this.scanShifts();
              }

              return true;
           }
        },

        changeItemsFilter: function(newValue) {
           if (newValue === undefined)
              newValue = this.getProperty("/itemsFilter") || "";

           // ignore same value
           if (newValue === this.itemsFilter)
              return;

           this.itemsFilter = newValue;

           // now we should request values once again
           this.submitRequest(false, this.h);
        }

    });

    return hRootModel;

});
