sap.ui.define([
    "sap/ui/model/json/JSONModel",
    "rootui5/geom/model/BrowserListBinding",
    "sap/base/Log"
], function(JSONModel, BrowserListBinding, Log) {
   "use strict";

    var hRootModel = JSONModel.extend("rootui5.eve7.model.BrowserModel", {

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

            this.loadDataCounter = 0; // counter of number of nodes

            this.sortOrder = "";

            this.threshold = 100; // default threshold to prefetch items
        },

        /* Method can be used when complete hierarchy is ready and can be used directly */
        setFullModel: function(topnode) {
           this.fullModel = true;
           this.h.nchilds = 1;
           this.h.childs = [ topnode ];
           topnode.expanded = true;
           this.scanShifts();
           if (this.oBinding)
              this.oBinding.checkUpdate(true);
        },

        clearFullModel: function() {
           delete this.h.childs;
           delete this.h.nchilds;
           delete this.fullModel;
        },


        bindTree: function(sPath, oContext, aFilters, mParameters, aSorters) {
           Log.warning("root.model.hModel#bindTree() " + sPath);

           console.log('BINDING TREE!!!!!!!!!!!!! ' + sPath);

           this.oBinding = new BrowserListBinding(this, sPath, oContext, aFilters, mParameters, aSorters);
           return this.oBinding;
        },

        getLength: function() {
           return this.getProperty("/length");
        },

        getNodeByPath: function(path) {
           var curr = this.h;
           if (!path || (typeof path !== "string") || (path == "/")) return curr;

           var names = path.split("/");

           while (names.length > 0) {
              var name = names.shift(), find = false;
              if (!name) continue;

              for (var k=0;k<curr.childs.length;++k) {
                 if (curr.childs[k].name == name) {
                    curr = curr.childs[k];
                    find = true;
                    break;
                 }
              }

              if (!find) return null;
           }
           return curr;
        },


        sendFirstRequest: function(websocket) {
           this._websocket = websocket;
           // submit top-level request already when construct model
           this.submitRequest(this.h, "/");
        },

        // submit next request to the server
        // directly use web socket, later can be dedicated channel
        submitRequest: function(elem, path, first, number) {

           if (!this._websocket || elem._requested || this.fullModel) return;
           elem._requested = true;

           this.loadDataCounter++;

           var request = {
              _typename: "ROOT::Experimental::RBrowserRequest",
              path: path,
              first: first || 0,
              number: number || this.threshold || 100,
              sort: this.sortOrder || ""
           };

           this._websocket.Send("BRREQ:" + JSON.stringify(request));
        },

        // process reply from server
        // In the future reply will be received via websocket channel
        processResponse: function(reply) {

           this.loadDataCounter--;

           // console.log('PROCESS BR RESPONSE', reply.path, reply);

           var elem = this.getNodeByPath(reply.path);

           if (!elem) { console.error('DID NOT FOUND ' + reply.path); return; }

           if (!elem._requested) console.error('ELEMENT WAS NOT REQUESTED!!!!', reply.path);
           delete elem._requested;

           var smart_merge = false;

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

           this.scanShifts();

           // reset existing nodes if reply does not match with expectation
           if (!smart_merge)
              this.reset_nodes = true;

           if (this.loadDataCounter == 0)
              if (this.oBinding)
                 this.oBinding.checkUpdate(true);
        },

        // return element of hierarchical structure by TreeTable index
        getElementByIndex: function(indx) {
           var nodes = this.getProperty("/nodes"),
               node = nodes ? nodes[indx] : null;

           return node ? node._elem : null;
        },

        // function used to calculate all ids shifts and total number of elements
        scanShifts: function() {

           var id = 0, full = this.fullModel;

           function scan(lvl, elem) {

              if (lvl >= 0) id++;

              var before_id = id;

              if (elem.expanded) {
                 if (elem.childs === undefined) {
                    // do nothing, childs are not visible as long as we do not have any list

                    // id += 0;
                 } else {

                    // gap at the begin
                    if (!full && elem.first)
                       id += elem.first;

                    // jump over all childs
                    for (var k=0;k<elem.childs.length;++k)
                       scan(lvl+1, elem.childs[k]);

                    // gap at the end
                    if (!full) {
                       var _last = (elem.first || 0) + elem.childs.length;
                       var _remains = elem.nchilds  - _last;
                       if (_remains > 0) id += _remains;
                    }
                 }
              }

              // this shift can be later applied to jump over all elements
              elem._shift = id - before_id;
           }

           scan(-1, this.h);

           this.setProperty("/length", id);

           return id;
        },

        // main  method to create flat list of nodes - only whose which are specified in selection
        // following arguments:
        //    args.begin     - first visisble element from flat list
        //    args.end       - first not-visisble element
        //    args.threshold - extra elements (before/after) which probably should be prefetched
        // returns total number of nodes
        buildFlatNodes: function(args) {

           var pthis = this,
               id = 0,            // current id, at the same time number of items
               threshold = args.threshold || this.threshold || 100,
               threshold2 = Math.round(threshold/2),
               nodes = this.reset_nodes ? {} : this.getProperty("/nodes");

           // main method to scan through all existing sub-folders
           function scan(lvl, elem, path) {

              // create elements with safety margin
              if ((lvl >= 0) && (nodes !== null) && !nodes[id] && (id >= args.begin - threshold2) && (id < args.end + threshold2) )
                 nodes[id] = {
                    name: elem.name,
                    level: lvl,
                    index: id,
                    _elem: elem,
                    // these are optional, should be eliminated in the future
                    type: elem.nchilds || (id == 0) ? "folder" : "file",
                    isLeaf: !elem.nchilds,
                    expanded: !!elem.expanded
                 };

              if (lvl >= 0) id++;

              if (!elem.expanded) return;

              if (elem.childs === undefined) {
                 // add new request - can we check if only special part of childs is required?

                 // TODO: probably one could guess more precise request
                 pthis.submitRequest(elem, path);

                 return;
              }

              // check if scan is required
              if (((id + elem._shift) < args.begin - threshold2) || (id >= args.end + threshold2)) {
                 id += elem._shift;
                 return;
              }

              // when not all childs from very beginning is loaded, but may be required
              if (elem.first && !pthis.fullModel) {

                 // check if requests are needed to load part in the begin of the list
                 if (args.begin - id - threshold2 < elem.first) {

                    var first = Math.max(args.begin - id - threshold2, 0),
                        number = Math.min(elem.first - first, threshold);

                    pthis.submitRequest(elem, path, first, number);
                 }

                 id += elem.first;
              }


              for (var k=0;k<elem.childs.length;++k)
                 scan(lvl+1, elem.childs[k], path + elem.childs[k].name + "/");

              // check if more elements are required

              if (!pthis.fullModel) {
                 var _last = (elem.first || 0) + elem.childs.length;
                 var _remains = elem.nchilds  - _last;

                 if (_remains > 0) {
                    if (args.end + threshold2 > id) {

                       var first = _last, number = args.end + threshold2 - id;
                       if (number < threshold) number = threshold; // always request much
                       if (number > _remains) number = _remains; // but not too much
                       if (number > threshold) {
                          first += (number - threshold);
                          number = threshold;
                       }

                       console.log('submit request for last', path, first,number)

                       pthis.submitRequest(elem, path, first, number);
                    }

                    id += _remains;
                 }
              }
           }

           // start scan from very top
           scan(-1, this.h, "/");

           if (this.getProperty("/length") != id) {
              // console.error('LENGTH MISMATCH', this.getProperty("/length"), id);
              this.setProperty("/length", id); // update length property
           }

           if (this.reset_nodes) {
              this.setProperty("/nodes", nodes);
              delete this.reset_nodes;
           }

           return id;
        },

        // toggle expand state of specified node
        toggleNode: function(index) {

           var elem = this.getElementByIndex(index);
           if (!elem) return;

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

        // change sorting method, for now server supports default, "direct" and "reverse"
        changeSortOrder: function(newValue) {
           if (newValue === undefined)
               newValue = this.getProperty("/sortOrder") || "";

           if ((newValue !== "") && (newValue !=="direct") && (newValue !== "reverse")) {
              console.error('WRONG sorting order ', newValue, 'use default');
              newValue = "";
           }

           // ignore same value
           if (newValue === this.sortOrder)
              return;


           this.sortOrder = newValue;

           // now we should request values once again

           this.submitRequest(this.h, "/");

        }

    });

    return hRootModel;

});
