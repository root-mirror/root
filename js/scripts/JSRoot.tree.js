// Collect all TTree-relevant methods for reading, drawing and hierarchy browsing

JSROOT.define(['io', 'math'], (jsrio, jsrmath) => {

   "use strict";

   JSROOT.BranchType = {
      kLeafNode: 0, kBaseClassNode: 1, kObjectNode: 2, kClonesNode: 3,
      kSTLNode: 4, kClonesMemberNode: 31, kSTLMemberNode: 41
   };

   jsrio.BranchBits = {
      kDoNotProcess: JSROOT.BIT(10), // Active bit for branches
      kIsClone: JSROOT.BIT(11), // to indicate a TBranchClones
      kBranchObject: JSROOT.BIT(12), // branch is a TObject*
      kBranchAny: JSROOT.BIT(17), // branch is an object*
      kAutoDelete: JSROOT.BIT(15),
      kDoNotUseBufferMap: JSROOT.BIT(22) // If set, at least one of the entry in the branch will use the buffer's map of classname and objects.
   }

   /**
    * @summary Class to read data from TTree
    *
    * @class
    * @memberof JSROOT
    * @desc Instance of TSelector can be used to access TTree data
    */

   function TSelector() {
      this.branches = []; // list of branches to read
      this.names = []; // list of member names for each branch in tgtobj
      this.directs = []; // indication if only branch without any children should be read
      this.break_execution = 0;
      this.tgtobj = {};
   }

   /** @summary Add branch to the selector
    * @desc Either branch name or branch itself should be specified
    * Second parameter defines member name in the tgtobj
    * If selector.addBranch("px", "read_px") is called,
    * branch will be read into selector.tgtobj.read_px member
    * If second parameter not specified, branch name (here "px") will be used
    * If branch object specified as first parameter and second parameter missing,
    * then member like "br0", "br1" and so on will be assigned
    * @param {string|Object} branch - name of branch (or branch object itself}
    * @param {string} [name] - member name in tgtobj where data will be read */
   TSelector.prototype.addBranch = function(branch, name, direct) {
      if (!name)
         name = (typeof branch === 'string') ? branch : ("br" + this.branches.length);
      this.branches.push(branch);
      this.names.push(name);
      this.directs.push(direct);
      return this.branches.length - 1;
   }

   /** @summary returns index of branch
     * @private */
   TSelector.prototype.indexOfBranch = function(branch) {
      return this.branches.indexOf(branch);
   }

   /** @summary returns name of branch
     * @private */
   TSelector.prototype.nameOfBranch = function(indx) {
      return this.names[indx];
   }

   /** @summary function called during TTree processing
    * @abstract
    * @param {number} progress - current value between 0 and 1 */
   TSelector.prototype.ShowProgress = function(/* progress */) {
   }

   /** @summary call this function to abort processing */
   TSelector.prototype.Abort = function() {
      this.break_execution = -1111;
   }

   /** @summary function called before start processing
    * @abstract
    * @param {object} tree - tree object */
   TSelector.prototype.Begin = function(/* tree */) {
   }

   /** @summary function called when next entry extracted from the tree
    * @abstract
    * @param {number} entry - read entry number */
   TSelector.prototype.Process = function(/* entry */) {
   }

   /** @summary function called at the very end of processing
    * @abstract
    * @param {boolean} res - true if all data were correctly processed */
   TSelector.prototype.Terminate = function(/* res */) {
   }

   // =================================================================

   /** @summary Checks array kind
     * @desc return 0 when not array
     * 1 - when arbitrary array
     * 2 - when plain (1-dim) array with same-type content */
   function checkArrayPrototype(arr, check_content) {
      if (typeof arr !== 'object') return 0;
      let proto = Object.prototype.toString.apply(arr);
      if (proto.indexOf('[object') !== 0) return 0;
      const pos = proto.indexOf('Array]');
      if (pos < 0) return 0;
      if (pos > 8) return 2; // this is typed array like Int32Array

      if (!check_content) return 1; //
      let typ, plain = true;
      for (let k = 0; k < arr.length; ++k) {
         let sub = typeof arr[k];
         if (!typ) typ = sub;
         if (sub !== typ) { plain = false; break; }
         if ((sub == "object") && checkArrayPrototype(arr[k])) { plain = false; break; }
      }

      return plain ? 2 : 1;
   }

   /**
    * @summary Class to iterate over array elements
    *
    * @class
    * @private
    */

   function ArrayIterator(arr, select, tgtobj) {
      this.object = arr;
      this.value = 0; // value always used in iterator
      this.arr = []; // all arrays
      this.indx = []; // all indexes
      this.cnt = -1; // current index counter
      this.tgtobj = tgtobj;

      if (typeof select === 'object')
         this.select = select; // remember indexes for selection
      else
         this.select = []; // empty array, undefined for each dimension means iterate over all indexes
   }

   /** @summary next element */
   ArrayIterator.prototype.next = function() {
      let obj, typ, cnt = this.cnt;

      if (cnt >= 0) {

         if (++this.fastindx < this.fastlimit) {
            this.value = this.fastarr[this.fastindx];
            return true;
         }

         while (--cnt >= 0) {
            if ((this.select[cnt] === undefined) && (++this.indx[cnt] < this.arr[cnt].length)) break;
         }
         if (cnt < 0) return false;
      }

      while (true) {

         obj = (cnt < 0) ? this.object : (this.arr[cnt])[this.indx[cnt]];

         typ = obj ? typeof obj : "any";

         if (typ === "object") {
            if (obj._typename !== undefined) {
               if (JSROOT.isRootCollection(obj)) { obj = obj.arr; typ = "array"; }
               else typ = "any";
            } else if (!isNaN(obj.length) && (checkArrayPrototype(obj) > 0)) {
               typ = "array";
            } else {
               typ = "any";
            }
         }

         if (this.select[cnt + 1] == "$self$") {
            this.value = obj;
            this.fastindx = this.fastlimit = 0;
            this.cnt = cnt + 1;
            return true;
         }

         if ((typ == "any") && (typeof this.select[cnt + 1] === "string")) {
            // this is extraction of the member from arbitrary class
            this.arr[++cnt] = obj;
            this.indx[cnt] = this.select[cnt]; // use member name as index
            continue;
         }

         if ((typ === "array") && ((obj.length > 0) || (this.select[cnt + 1] === "$size$"))) {
            this.arr[++cnt] = obj;
            switch (this.select[cnt]) {
               case undefined: this.indx[cnt] = 0; break;
               case "$last$": this.indx[cnt] = obj.length - 1; break;
               case "$size$":
                  this.value = obj.length;
                  this.fastindx = this.fastlimit = 0;
                  this.cnt = cnt;
                  return true;
               default:
                  if (!isNaN(this.select[cnt])) {
                     this.indx[cnt] = this.select[cnt];
                     if (this.indx[cnt] < 0) this.indx[cnt] = obj.length - 1;
                  } else {
                     // this is compile variable as array index - can be any expression
                     this.select[cnt].produce(this.tgtobj);
                     this.indx[cnt] = Math.round(this.select[cnt].get(0));
                  }
            }
         } else {

            if (cnt < 0) return false;

            this.value = obj;
            if (this.select[cnt] === undefined) {
               this.fastarr = this.arr[cnt];
               this.fastindx = this.indx[cnt];
               this.fastlimit = this.fastarr.length;
            } else {
               this.fastindx = this.fastlimit = 0; // no any iteration on that level
            }

            this.cnt = cnt;
            return true;
         }
      }

      // unreachable code
      // return false;
   }

   /** @summary reset iterator */
   ArrayIterator.prototype.reset = function() {
      this.arr = [];
      this.indx = [];
      delete this.fastarr;
      this.cnt = -1;
      this.value = 0;
   }

   // ============================================================================

   /**
    * @summary object with single variable in TTree::Draw expression
    *
    * @class
    * @memberof JSROOT
    * @private
    */

   function TDrawVariable(globals) {
      this.globals = globals;

      this.code = "";
      this.brindex = []; // index of used branches from selector
      this.branches = []; // names of branches in target object
      this.brarray = []; // array specifier for each branch
      this.func = null; // generic function for variable calculation

      this.kind = undefined;
      this.buf = []; // buffer accumulates temporary values
   }

   /** @summary Parse variable
     * @desc when only_branch specified, its placed in the front of the expression */
   TDrawVariable.prototype.parse = function(tree, selector, code, only_branch, branch_mode) {

      function is_start_symbol(symb) {
         if ((symb >= "A") && (symb <= "Z")) return true;
         if ((symb >= "a") && (symb <= "z")) return true;
         return (symb === "_");
      }

      function is_next_symbol(symb) {
         if (is_start_symbol(symb)) return true;
         if ((symb >= "0") && (symb <= "9")) return true;
         return false;
      }

      if (!code) code = ""; // should be empty string at least

      this.code = (only_branch ? only_branch.fName : "") + code;

      let pos = 0, pos2 = 0, br = null;
      while ((pos < code.length) || only_branch) {

         let arriter = [];

         if (only_branch) {
            br = only_branch;
            only_branch = undefined;
         } else {
            // first try to find branch
            pos2 = pos;
            while ((pos2 < code.length) && (is_next_symbol(code[pos2]) || code[pos2] === ".")) pos2++;
            if (code[pos2] == "$") {
               let repl = "";
               switch (code.substr(pos, pos2 - pos)) {
                  case "LocalEntry":
                  case "Entry": repl = "arg.$globals.entry"; break;
                  case "Entries": repl = "arg.$globals.entries"; break;
               }
               if (repl) {
                  code = code.substr(0, pos) + repl + code.substr(pos2 + 1);
                  pos = pos + repl.length;
                  continue;
               }
            }

            br = tree.FindBranchComplex(code.substr(pos, pos2 - pos));
            if (!br) { pos = pos2 + 1; continue; }

            // when full id includes branch name, replace only part of extracted expression
            if (br.branch && (br.rest !== undefined)) {
               pos2 -= br.rest.length;
               branch_mode = undefined; // maybe selection of the sub-object done
               br = br.branch;
            }

            // when code ends with the point - means object itself will be accessed
            // sometime branch name itself ends with the point
            if ((pos2 >= code.length - 1) && (code[code.length - 1] === ".")) {
               arriter.push("$self$");
               pos2 = code.length;
            }
         }

         // now extract all levels of iterators
         while (pos2 < code.length) {

            if ((code[pos2] === "@") && (code.substr(pos2, 5) == "@size") && (arriter.length == 0)) {
               pos2 += 5;
               branch_mode = true;
               break;
            }

            if (code[pos2] === ".") {
               // this is object member
               let prev = ++pos2;

               if ((code[prev] === "@") && (code.substr(prev, 5) === "@size")) {
                  arriter.push("$size$");
                  pos2 += 5;
                  break;
               }

               if (!is_start_symbol(code[prev])) {
                  arriter.push("$self$"); // last point means extraction of object itself
                  break;
               }

               while ((pos2 < code.length) && is_next_symbol(code[pos2])) pos2++;

               // this is looks like function call - do not need to extract member with
               if (code[pos2] == "(") { pos2 = prev - 1; break; }

               // this is selection of member, but probably we need to activate iterator for ROOT collection
               if ((arriter.length === 0) && br) {
                  // TODO: if selected member is simple data type - no need to make other checks - just break here
                  if ((br.fType === JSROOT.BranchType.kClonesNode) || (br.fType === JSROOT.BranchType.kSTLNode)) {
                     arriter.push(undefined);
                  } else {
                     let objclass = getBranchObjectClass(br, tree, false, true);
                     if (objclass && JSROOT.isRootCollection(null, objclass)) arriter.push(undefined);
                  }
               }
               arriter.push(code.substr(prev, pos2 - prev));
               continue;
            }

            if (code[pos2] !== "[") break;

            // simple []
            if (code[pos2 + 1] == "]") { arriter.push(undefined); pos2 += 2; continue; }

            let prev = pos2++, cnt = 0;
            while ((pos2 < code.length) && ((code[pos2] != "]") || (cnt > 0))) {
               if (code[pos2] == '[') cnt++; else if (code[pos2] == ']') cnt--;
               pos2++;
            }
            let sub = code.substr(prev + 1, pos2 - prev - 1);
            switch (sub) {
               case "":
               case "$all$": arriter.push(undefined); break;
               case "$last$": arriter.push("$last$"); break;
               case "$size$": arriter.push("$size$"); break;
               case "$first$": arriter.push(0); break;
               default:
                  if (!isNaN(parseInt(sub))) {
                     arriter.push(parseInt(sub));
                  } else {
                     // try to compile code as draw variable
                     let subvar = new TDrawVariable(this.globals);
                     if (!subvar.parse(tree, selector, sub)) return false;
                     arriter.push(subvar);
                  }
            }
            pos2++;
         }

         if (arriter.length === 0) arriter = undefined; else
            if ((arriter.length === 1) && (arriter[0] === undefined)) arriter = true;

         let indx = selector.indexOfBranch(br);
         if (indx < 0) indx = selector.addBranch(br, undefined, branch_mode);

         branch_mode = undefined;

         this.brindex.push(indx);
         this.branches.push(selector.nameOfBranch(indx));
         this.brarray.push(arriter);

         // this is simple case of direct usage of the branch
         if ((pos === 0) && (pos2 === code.length) && (this.branches.length === 1)) {
            this.direct_branch = true; // remember that branch read as is
            return true;
         }

         let replace = "arg.var" + (this.branches.length - 1);

         code = code.substr(0, pos) + replace + code.substr(pos2);

         pos = pos + replace.length;
      }

      // support usage of some standard TMath functions
      code = code.replace(/TMath::Exp\(/g, 'Math.exp(')
                 .replace(/TMath::Abs\(/g, 'Math.abs(')
                 .replace(/TMath::Prob\(/g, 'arg.$math.Prob(')
                 .replace(/TMath::Gaus\(/g, 'arg.$math.Gaus(');

      this.func = new Function("arg", "return (" + code + ")");

      return true;
   }

   /** @summary Check if it is dummy variable */
   TDrawVariable.prototype.is_dummy = function() {
      return (this.branches.length === 0) && !this.func;
   }

   /** @summary Produce variable
     * @desc after reading tree braches into the object, calculate variable value */
   TDrawVariable.prototype.produce = function(obj) {

      this.length = 1;
      this.isarray = false;

      if (this.is_dummy()) {
         this.value = 1.; // used as dummy weight variable
         this.kind = "number";
         return;
      }

      let arg = { $globals: this.globals, $math: jsrmath }, usearrlen = -1, arrs = [];
      for (let n = 0; n < this.branches.length; ++n) {
         let name = "var" + n;
         arg[name] = obj[this.branches[n]];

         // try to check if branch is array and need to be iterated
         if (this.brarray[n] === undefined)
            this.brarray[n] = (checkArrayPrototype(arg[name]) > 0) || JSROOT.isRootCollection(arg[name]);

         // no array - no pain
         if (this.brarray[n] === false) continue;

         // check if array can be used as is - one dimension and normal values
         if ((this.brarray[n] === true) && (checkArrayPrototype(arg[name], true) === 2)) {
            // plain array, can be used as is
            arrs[n] = arg[name];
         } else {
            let iter = new ArrayIterator(arg[name], this.brarray[n], obj);
            arrs[n] = [];
            while (iter.next()) arrs[n].push(iter.value);
         }
         if ((usearrlen < 0) || (usearrlen < arrs[n].length)) usearrlen = arrs[n].length;
      }

      if (usearrlen < 0) {
         this.value = this.direct_branch ? arg.var0 : this.func(arg);
         if (!this.kind) this.kind = typeof this.value;
         return;
      }

      if (usearrlen == 0) {
         // empty array - no any histogram should be filled
         this.length = 0;
         this.value = 0;
         return;
      }

      this.length = usearrlen;
      this.isarray = true;

      if (this.direct_branch) {
         this.value = arrs[0]; // just use array
      } else {
         this.value = new Array(usearrlen);

         for (let k = 0; k < usearrlen; ++k) {
            for (let n = 0; n < this.branches.length; ++n) {
               if (arrs[n]) arg["var" + n] = arrs[n][k];
            }
            this.value[k] = this.func(arg);
         }
      }

      if (!this.kind) this.kind = typeof this.value[0];
   }

   /** @summary Get variable */
   TDrawVariable.prototype.get = function(indx) {
      return this.isarray ? this.value[indx] : this.value;
   }

   /** @summary Append array to the buffer */
   TDrawVariable.prototype.appendArray = function(tgtarr) {
      this.buf = this.buf.concat(tgtarr[this.branches[0]]);
   }

   // =============================================================================

   /**
    * @summary Selector class for TTree::Draw function
    *
    * @class
    * @memberof JSROOT
    * @extends JSROOT.TSelector
    * @private
    */

   function TDrawSelector() {
      TSelector.call(this);

      this.ndim = 0;
      this.vars = []; // array of expression variables
      this.cut = null; // cut variable
      this.hist = null;
      this.histo_drawopt = "";
      this.hist_name = "$htemp";
      this.hist_title = "Result of TTree::Draw";
      this.graph = false;
      this.hist_args = []; // arguments for histogram creation
      this.arr_limit = 1000;  // number of accumulated items before create histogram
      this.htype = "F";
      this.monitoring = 0;
      this.globals = {}; // object with global parameters, which could be used in any draw expression
      this.last_progress = 0;
      this.aver_diff = 0;
   }

   TDrawSelector.prototype = Object.create(TSelector.prototype);

   /** @summary Set draw selector callbacks */
   TDrawSelector.prototype.setCallback = function(result_callback, progress_callback) {
      this.result_callback = result_callback;
      this.progress_callback = progress_callback;
   }

   /** @summary Parse parameters */
   TDrawSelector.prototype.parseParameters = function(tree, args, expr) {

      if (!expr || (typeof expr !== "string")) return "";

      // parse parameters which defined at the end as expression;par1name:par1value;par2name:par2value
      let pos = expr.lastIndexOf(";");
      while (pos >= 0) {
         let parname = expr.substr(pos + 1), parvalue = undefined;
         expr = expr.substr(0, pos);
         pos = expr.lastIndexOf(";");

         let separ = parname.indexOf(":");
         if (separ > 0) { parvalue = parname.substr(separ + 1); parname = parname.substr(0, separ); }

         let intvalue = parseInt(parvalue);
         if (!parvalue || isNaN(intvalue)) intvalue = undefined;

         switch (parname) {
            case "num":
            case "entries":
            case "numentries":
               if (parvalue === "all") args.numentries = tree.fEntries; else
                  if (parvalue === "half") args.numentries = Math.round(tree.fEntries / 2); else
                     if (intvalue !== undefined) args.numentries = intvalue;
               break;
            case "first":
               if (intvalue !== undefined) args.firstentry = intvalue;
               break;
            case "mon":
            case "monitor":
               args.monitoring = (intvalue !== undefined) ? intvalue : 5000;
               break;
            case "player":
               args.player = true;
               break;
            case "dump":
               args.dump = true;
               break;
            case "maxseg":
            case "maxrange":
               if (intvalue) tree.$file.fMaxRanges = intvalue;
               break;
            case "accum":
               if (intvalue) this.arr_limit = intvalue;
               break;
            case "htype":
               if (parvalue && (parvalue.length === 1)) {
                  this.htype = parvalue.toUpperCase();
                  if ((this.htype !== "C") && (this.htype !== "S") && (this.htype !== "I")
                     && (this.htype !== "F") && (this.htype !== "L") && (this.htype !== "D")) this.htype = "F";
               }
               break;
            case "hbins":
               this.hist_nbins = parseInt(parvalue);
               if (isNaN(this.hist_nbins) || (this.hist_nbins <= 3)) delete this.hist_nbins;
               break;
            case "drawopt":
               args.drawopt = parvalue;
               break;
            case "graph":
               args.graph = intvalue || true;
               break;
         }
      }

      pos = expr.lastIndexOf(">>");
      if (pos >= 0) {
         let harg = expr.substr(pos + 2).trim();
         expr = expr.substr(0, pos).trim();
         pos = harg.indexOf("(");
         if (pos > 0) {
            this.hist_name = harg.substr(0, pos);
            harg = harg.substr(pos);
         }
         if (harg === "dump") {
            args.dump = true;
         } else if (harg.indexOf("Graph") == 0) {
            args.graph = true;
         } else if (pos < 0) {
            this.hist_name = harg;
         } else if ((harg[0] == "(") && (harg[harg.length - 1] == ")")) {
            harg = harg.substr(1, harg.length - 2).split(",");
            let isok = true;
            for (let n = 0; n < harg.length; ++n) {
               harg[n] = (n % 3 === 0) ? parseInt(harg[n]) : parseFloat(harg[n]);
               if (isNaN(harg[n])) isok = false;
            }
            if (isok) this.hist_args = harg;
         }
      }

      if (args.dump) {
         this.dump_values = true;
         args.reallocate_objects = true;
         if (args.numentries === undefined) args.numentries = 10;
      }

      return expr;
   }

   /** @summary Parse draw expression */
   TDrawSelector.prototype.parseDrawExpression = function(tree, args) {

      // parse complete expression
      let expr = this.parseParameters(tree, args, args.expr), cut = "";

      // parse option for histogram creation
      this.hist_title = "drawing '" + expr + "' from " + tree.fName;

      let pos = 0;
      if (args.cut) {
         cut = args.cut;
      } else {
         pos = expr.replace(/TMath::/g, 'TMath__').lastIndexOf("::"); // avoid confusion due-to :: in the namespace
         if (pos > 0) {
            cut = expr.substr(pos + 2).trim();
            expr = expr.substr(0, pos).trim();
         }
      }

      args.parse_expr = expr;
      args.parse_cut = cut;

      // let names = expr.split(":"); // to allow usage of ? operator, we need to handle : as well
      let names = [], nbr1 = 0, nbr2 = 0, prev = 0;
      for (pos = 0; pos < expr.length; ++pos) {
         switch (expr[pos]) {
            case "(": nbr1++; break;
            case ")": nbr1--; break;
            case "[": nbr2++; break;
            case "]": nbr2--; break;
            case ":":
               if (expr[pos + 1] == ":") { pos++; continue; }
               if (!nbr1 && !nbr2 && (pos > prev)) names.push(expr.substr(prev, pos - prev));
               prev = pos + 1;
               break;
         }
      }
      if (!nbr1 && !nbr2 && (pos > prev)) names.push(expr.substr(prev, pos - prev));

      if ((names.length < 1) || (names.length > 3)) return false;

      this.ndim = names.length;

      let is_direct = !cut;

      for (let n = 0; n < this.ndim; ++n) {
         this.vars[n] = new TDrawVariable(this.globals);
         if (!this.vars[n].parse(tree, this, names[n])) return false;
         if (!this.vars[n].direct_branch) is_direct = false;
      }

      this.cut = new TDrawVariable(this.globals);
      if (cut)
         if (!this.cut.parse(tree, this, cut)) return false;

      if (!this.branches.length) {
         console.warn('no any branch is selected');
         return false;
      }

      if (is_direct) this.ProcessArrays = this.ProcessArraysFunc;

      this.monitoring = args.monitoring;

      this.graph = args.graph;

      if (args.drawopt !== undefined)
         this.histo_drawopt = args.drawopt;
      else
         this.histo_drawopt = (this.ndim === 2) ? "col" : "";

      return true;
   }

   /** @summary Draw only specified branch */
   TDrawSelector.prototype.drawOnlyBranch = function(tree, branch, expr, args) {
      this.ndim = 1;

      if (expr.indexOf("dump") == 0) expr = ";" + expr;

      expr = this.parseParameters(tree, args, expr);

      this.monitoring = args.monitoring;

      if (args.dump) {
         this.dump_values = true;
         args.reallocate_objects = true;
      }

      if (this.dump_values) {

         this.hist = []; // array of dump objects

         this.leaf = args.leaf;

         // branch object remains, therefore we need to copy fields to see them all
         this.copy_fields = ((args.branch.fLeaves && (args.branch.fLeaves.arr.length > 1)) ||
            (args.branch.fBranches && (args.branch.fBranches.arr.length > 0))) && !args.leaf;

         this.addBranch(branch, "br0", args.direct_branch); // add branch

         this.Process = this.ProcessDump;

         return true;
      }

      this.vars[0] = new TDrawVariable(this.globals);
      if (!this.vars[0].parse(tree, this, expr, branch, args.direct_branch)) return false;
      this.hist_title = "drawing branch '" + branch.fName + (expr ? "' expr:'" + expr : "") + "'  from " + tree.fName;

      this.cut = new TDrawVariable(this.globals);

      if (this.vars[0].direct_branch) this.ProcessArrays = this.ProcessArraysFunc;

      return true;
   }

   /** @summary Begin processing */
   TDrawSelector.prototype.Begin = function(tree) {
      this.globals.entries = tree.fEntries;

      if (this.monitoring)
         this.lasttm = new Date().getTime();
   }

   /** @summary Show progress */
   TDrawSelector.prototype.ShowProgress = function(value) {
      if (typeof document == 'undefined' || !JSROOT.Painter) return;

      if ((value === undefined) || isNaN(value)) return JSROOT.Painter.showProgress();

      if (this.last_progress !== value) {
         let diff = value - this.last_progress;
         if (!this.aver_diff) this.aver_diff = diff;
         this.aver_diff = diff * 0.3 + this.aver_diff * 0.7;
      }

      let ndig = 0;
      if (this.aver_diff <= 0) ndig = 0; else
         if (this.aver_diff < 0.0001) ndig = 3; else
            if (this.aver_diff < 0.001) ndig = 2; else
               if (this.aver_diff < 0.01) ndig = 1;

      let main_box = document.createElement("p"),
         text_node = document.createTextNode("TTree draw " + (value * 100).toFixed(ndig) + " %  "),
         selector = this;

      main_box.appendChild(text_node);
      main_box.title = "Click on element to break drawing";

      main_box.onclick = function() {
         if (++selector.break_execution < 3) {
            main_box.title = "Tree draw will break after next I/O operation";
            return text_node.nodeValue = "Breaking ... ";
         }
         selector.Abort();
         JSROOT.Painter.showProgress();
      }

      JSROOT.Painter.showProgress(main_box);
      this.last_progress = value;
   }

   /** @summary Get bins for bits histogram */
   TDrawSelector.prototype.getBitsBins = function(nbits, res) {
      res.nbins = res.max = nbits;
      res.fLabels = JSROOT.create("THashList");
      for (let k = 0; k < nbits; ++k) {
         let s = JSROOT.create("TObjString");
         s.fString = k.toString();
         s.fUniqueID = k + 1;
         res.fLabels.Add(s);
      }
      return res;
   }

   /** @summary Get min.max bins */
   TDrawSelector.prototype.getMinMaxBins = function(axisid, nbins) {

      let res = { min: 0, max: 0, nbins: nbins, k: 1., fLabels: null, title: "" };

      if (axisid >= this.ndim) return res;

      let arr = this.vars[axisid].buf;

      res.title = this.vars[axisid].code || "";

      if (this.vars[axisid].kind === "object") {
         // this is any object type
         let typename, similar = true, maxbits = 8;
         for (let k = 0; k < arr.length; ++k) {
            if (!arr[k]) continue;
            if (!typename) typename = arr[k]._typename;
            if (typename !== arr[k]._typename) similar = false; // check all object types
            if (arr[k].fNbits) maxbits = Math.max(maxbits, arr[k].fNbits + 1);
         }

         if (typename && similar) {
            if ((typename === "TBits") && (axisid === 0)) {
               this.fill1DHistogram = this.fillTBitsHistogram;
               if (maxbits % 8) maxbits = (maxbits & 0xfff0) + 8;

               if ((this.hist_name === "bits") && (this.hist_args.length == 1) && this.hist_args[0])
                  maxbits = this.hist_args[0];

               return this.getBitsBins(maxbits, res);
            }
         }
      }

      if (this.vars[axisid].kind === "string") {
         res.lbls = []; // all labels

         for (let k = 0; k < arr.length; ++k)
            if (res.lbls.indexOf(arr[k]) < 0)
               res.lbls.push(arr[k]);

         res.lbls.sort();
         res.max = res.nbins = res.lbls.length;

         res.fLabels = JSROOT.create("THashList");
         for (let k = 0; k < res.lbls.length; ++k) {
            let s = JSROOT.create("TObjString");
            s.fString = res.lbls[k];
            s.fUniqueID = k + 1;
            if (s.fString === "") s.fString = "<empty>";
            res.fLabels.Add(s);
         }
      } else if ((axisid === 0) && (this.hist_name === "bits") && (this.hist_args.length <= 1)) {
         this.fill1DHistogram = this.FillBitsHistogram;
         return this.getBitsBins(this.hist_args[0] || 32, res);
      } else if (axisid * 3 + 2 < this.hist_args.length) {
         res.nbins = this.hist_args[axisid * 3];
         res.min = this.hist_args[axisid * 3 + 1];
         res.max = this.hist_args[axisid * 3 + 2];
      } else {

         res.min = Math.min.apply(null, arr);
         res.max = Math.max.apply(null, arr);

         if (this.hist_nbins)
            nbins = res.nbins = this.hist_nbins;

         res.isinteger = (Math.round(res.min) === res.min) && (Math.round(res.max) === res.max);
         if (res.isinteger)
            for (let k = 0; k < arr.length; ++k)
               if (arr[k] !== Math.round(arr[k])) { res.isinteger = false; break; }

         if (res.isinteger) {
            res.min = Math.round(res.min);
            res.max = Math.round(res.max);
            if (res.max - res.min < nbins * 5) {
               res.min -= 1;
               res.max += 2;
               res.nbins = Math.round(res.max - res.min);
            } else {
               let range = (res.max - res.min + 2), step = Math.floor(range / nbins);
               while (step * nbins < range) step++;
               res.max = res.min + nbins * step;
            }
         } else if (res.min >= res.max) {
            res.max = res.min;
            if (Math.abs(res.min) < 100) { res.min -= 1; res.max += 1; } else
               if (res.min > 0) { res.min *= 0.9; res.max *= 1.1; } else { res.min *= 1.1; res.max *= 0.9; }
         } else {
            res.max += (res.max - res.min) / res.nbins;
         }
      }

      res.k = res.nbins / (res.max - res.min);

      res.GetBin = function(value) {
         const bin = this.lbls ? this.lbls.indexOf(value) : Math.floor((value - this.min) * this.k);
         return (bin < 0) ? 0 : ((bin > this.nbins) ? this.nbins + 1 : bin + 1);
      }

      return res;
   }

   /** @summary Create histogram */
   TDrawSelector.prototype.createHistogram = function() {
      if (this.hist || !this.vars[0].buf) return;

      if (this.dump_values) {
         // just create array where dumped valus will be collected
         this.hist = [];

         // reassign fill method
         this.fill1DHistogram = this.fill2DHistogram = this.fill3DHistogram = this.dumpValues;
      } else if (this.graph) {
         let N = this.vars[0].buf.length;

         if (this.ndim == 1) {
            // A 1-dimensional graph will just have the x axis as an index
            this.hist = JSROOT.createTGraph(N, Array.from(Array(N).keys()), this.vars[0].buf);
         } else if (this.ndim == 2) {
            this.hist = JSROOT.createTGraph(N, this.vars[0].buf, this.vars[1].buf);
            delete this.vars[1].buf;
         }

         this.hist.fTitle = this.hist_title;
         this.hist.fName = "Graph";

      } else {

         this.x = this.getMinMaxBins(0, (this.ndim > 1) ? 50 : 200);

         this.y = this.getMinMaxBins(1, 50);

         this.z = this.getMinMaxBins(2, 50);

         switch (this.ndim) {
            case 1: this.hist = JSROOT.createHistogram("TH1" + this.htype, this.x.nbins); break;
            case 2: this.hist = JSROOT.createHistogram("TH2" + this.htype, this.x.nbins, this.y.nbins); break;
            case 3: this.hist = JSROOT.createHistogram("TH3" + this.htype, this.x.nbins, this.y.nbins, this.z.nbins); break;
         }

         this.hist.fXaxis.fTitle = this.x.title;
         this.hist.fXaxis.fXmin = this.x.min;
         this.hist.fXaxis.fXmax = this.x.max;
         this.hist.fXaxis.fLabels = this.x.fLabels;

         if (this.ndim > 1) this.hist.fYaxis.fTitle = this.y.title;
         this.hist.fYaxis.fXmin = this.y.min;
         this.hist.fYaxis.fXmax = this.y.max;
         this.hist.fYaxis.fLabels = this.y.fLabels;

         if (this.ndim > 2) this.hist.fZaxis.fTitle = this.z.title;
         this.hist.fZaxis.fXmin = this.z.min;
         this.hist.fZaxis.fXmax = this.z.max;
         this.hist.fZaxis.fLabels = this.z.fLabels;

         this.hist.fName = this.hist_name;
         this.hist.fTitle = this.hist_title;
         this.hist.fOption = this.histo_drawopt;
         this.hist.$custom_stat = (this.hist_name == "$htemp") ? 111110 : 111111;
      }

      let var0 = this.vars[0].buf, cut = this.cut.buf, len = var0.length;

      if (!this.graph) {
         switch (this.ndim) {
            case 1: {
               for (let n = 0; n < len; ++n)
                  this.fill1DHistogram(var0[n], cut ? cut[n] : 1.);
               break;
            }
            case 2: {
               let var1 = this.vars[1].buf;
               for (let n = 0; n < len; ++n)
                  this.fill2DHistogram(var0[n], var1[n], cut ? cut[n] : 1.);
               delete this.vars[1].buf;
               break;
            }
            case 3: {
               let var1 = this.vars[1].buf, var2 = this.vars[2].buf;
               for (let n = 0; n < len; ++n)
                  this.fill3DHistogram(var0[n], var1[n], var2[n], cut ? cut[n] : 1.);
               delete this.vars[1].buf;
               delete this.vars[2].buf;
               break;
            }
         }
      }

      delete this.vars[0].buf;
      delete this.cut.buf;
   }

   /** @summary Fill TBits histogram */
   TDrawSelector.prototype.fillTBitsHistogram = function(xvalue, weight) {
      if (!weight || !xvalue || !xvalue.fNbits || !xvalue.fAllBits) return;

      const sz = Math.min(xvalue.fNbits + 1, xvalue.fNbytes * 8);

      for (let bit = 0, mask = 1, b = 0; bit < sz; ++bit) {
         if (xvalue.fAllBits[b] && mask) {
            if (bit <= this.x.nbins)
               this.hist.fArray[bit + 1] += weight;
            else
               this.hist.fArray[this.x.nbins + 1] += weight;
         }

         mask *= 2;
         if (mask >= 0x100) { mask = 1; ++b; }
      }
   }

   /** @summary Fill bits histogram */
   TDrawSelector.prototype.FillBitsHistogram = function(xvalue, weight) {
      if (!weight) return;

      for (let bit = 0, mask = 1; bit < this.x.nbins; ++bit) {
         if (xvalue & mask) this.hist.fArray[bit + 1] += weight;
         mask *= 2;
      }
   }

   /** @summary Fill 1D histogram */
   TDrawSelector.prototype.fill1DHistogram = function(xvalue, weight) {
      let bin = this.x.GetBin(xvalue);
      this.hist.fArray[bin] += weight;

      if (!this.x.lbls) {
         this.hist.fTsumw += weight;
         this.hist.fTsumwx += weight * xvalue;
         this.hist.fTsumwx2 += weight * xvalue * xvalue;
      }
   }

   /** @summary Fill 2D histogram */
   TDrawSelector.prototype.fill2DHistogram = function(xvalue, yvalue, weight) {
      let xbin = this.x.GetBin(xvalue),
         ybin = this.y.GetBin(yvalue);

      this.hist.fArray[xbin + (this.x.nbins + 2) * ybin] += weight;
      if (!this.x.lbls && !this.y.lbls) {
         this.hist.fTsumw += weight;
         this.hist.fTsumwx += weight * xvalue;
         this.hist.fTsumwy += weight * yvalue;
         this.hist.fTsumwx2 += weight * xvalue * xvalue;
         this.hist.fTsumwxy += weight * xvalue * yvalue;
         this.hist.fTsumwy2 += weight * yvalue * yvalue;
      }
   }

   /** @summary Fill 3D histogram */
   TDrawSelector.prototype.fill3DHistogram = function(xvalue, yvalue, zvalue, weight) {
      let xbin = this.x.GetBin(xvalue),
         ybin = this.y.GetBin(yvalue),
         zbin = this.z.GetBin(zvalue);

      this.hist.fArray[xbin + (this.x.nbins + 2) * (ybin + (this.y.nbins + 2) * zbin)] += weight;
      if (!this.x.lbls && !this.y.lbls && !this.z.lbls) {
         this.hist.fTsumw += weight;
         this.hist.fTsumwx += weight * xvalue;
         this.hist.fTsumwy += weight * yvalue;
         this.hist.fTsumwz += weight * zvalue;
         this.hist.fTsumwx2 += weight * xvalue * xvalue;
         this.hist.fTsumwy2 += weight * yvalue * yvalue;
         this.hist.fTsumwz2 += weight * zvalue * zvalue;
         this.hist.fTsumwxy += weight * xvalue * yvalue;
         this.hist.fTsumwxz += weight * xvalue * zvalue;
         this.hist.fTsumwyz += weight * yvalue * zvalue;
      }
   }

   /** @summary Dump values */
   TDrawSelector.prototype.dumpValues = function(v1, v2, v3, v4) {
      let obj;
      switch (this.ndim) {
         case 1: obj = { x: v1, weight: v2 }; break;
         case 2: obj = { x: v1, y: v2, weight: v3 }; break;
         case 3: obj = { x: v1, y: v2, z: v3, weight: v4 }; break;
      }

      if (this.cut.is_dummy()) {
         if (this.ndim === 1) obj = v1; else delete obj.weight;
      }

      this.hist.push(obj);
   }

    /** @summary function used when all branches can be read as array
      * @desc most typical usage - histogramming of single branch */
   TDrawSelector.prototype.ProcessArraysFunc = function(/*entry*/) {

      if (this.arr_limit || this.graph) {
         let var0 = this.vars[0], len = this.tgtarr.br0.length,
            var1 = this.vars[1], var2 = this.vars[2];
         if ((var0.buf.length === 0) && (len >= this.arr_limit) && !this.graph) {
            // special use case - first array large enough to create histogram directly base on it
            var0.buf = this.tgtarr.br0;
            if (var1) var1.buf = this.tgtarr.br1;
            if (var2) var2.buf = this.tgtarr.br2;
         } else {
            for (let k = 0; k < len; ++k) {
               var0.buf.push(this.tgtarr.br0[k]);
               if (var1) var1.buf.push(this.tgtarr.br1[k]);
               if (var2) var2.buf.push(this.tgtarr.br2[k]);
            }
         }
         var0.kind = "number";
         if (var1) var1.kind = "number";
         if (var2) var2.kind = "number";
         this.cut.buf = null; // do not create buffer for cuts
         if (!this.graph && (var0.buf.length >= this.arr_limit)) {
            this.createHistogram();
            this.arr_limit = 0;
         }
      } else {
         let br0 = this.tgtarr.br0, len = br0.length;
         switch (this.ndim) {
            case 1: {
               for (let k = 0; k < len; ++k)
                  this.fill1DHistogram(br0[k], 1.);
               break;
            }
            case 2: {
               let br1 = this.tgtarr.br1;
               for (let k = 0; k < len; ++k)
                  this.fill2DHistogram(br0[k], br1[k], 1.);
               break;
            }
            case 3: {
               let br1 = this.tgtarr.br1, br2 = this.tgtarr.br2;
               for (let k = 0; k < len; ++k)
                  this.fill3DHistogram(br0[k], br1[k], br2[k], 1.);
               break;
            }
         }
      }
   }

   /** @summary simple dump of the branch - no need to analyze something */
   TDrawSelector.prototype.ProcessDump = function(/* entry */) {

      let res = this.leaf ? this.tgtobj.br0[this.leaf] : this.tgtobj.br0;

      if (res && this.copy_fields) {
         if (checkArrayPrototype(res) === 0) {
            this.hist.push(JSROOT.extend({}, res));
         } else {
            this.hist.push(res);
         }
      } else {
         this.hist.push(res);
      }
   }

   /** @summary Normal TSelector Process handler */
   TDrawSelector.prototype.Process = function(entry) {

      this.globals.entry = entry; // can be used in any expression

      this.cut.produce(this.tgtobj);
      if (!this.dump_values && !this.cut.value) return;

      for (let n = 0; n < this.ndim; ++n)
         this.vars[n].produce(this.tgtobj);

      let var0 = this.vars[0], var1 = this.vars[1], var2 = this.vars[2], cut = this.cut;

      if (this.graph || this.arr_limit) {
         switch (this.ndim) {
            case 1:
               for (let n0 = 0; n0 < var0.length; ++n0) {
                  var0.buf.push(var0.get(n0));
                  cut.buf.push(cut.value);
               }
               break;
            case 2:
               for (let n0 = 0; n0 < var0.length; ++n0)
                  for (let n1 = 0; n1 < var1.length; ++n1) {
                     var0.buf.push(var0.get(n0));
                     var1.buf.push(var1.get(n1));
                     cut.buf.push(cut.value);
                  }
               break;
            case 3:
               for (let n0 = 0; n0 < var0.length; ++n0)
                  for (let n1 = 0; n1 < var1.length; ++n1)
                     for (let n2 = 0; n2 < var2.length; ++n2) {
                        var0.buf.push(var0.get(n0));
                        var1.buf.push(var1.get(n1));
                        var2.buf.push(var2.get(n2));
                        cut.buf.push(cut.value);
                     }
               break;
         }
         if (!this.graph && var0.buf.length >= this.arr_limit) {
            this.createHistogram();
            this.arr_limit = 0;
         }
      } else if (this.hist) {
         switch (this.ndim) {
            case 1:
               for (let n0 = 0; n0 < var0.length; ++n0)
                  this.fill1DHistogram(var0.get(n0), cut.value);
               break;
            case 2:
               for (let n0 = 0; n0 < var0.length; ++n0)
                  for (let n1 = 0; n1 < var1.length; ++n1)
                     this.fill2DHistogram(var0.get(n0), var1.get(n1), cut.value);
               break;
            case 3:
               for (let n0 = 0; n0 < var0.length; ++n0)
                  for (let n1 = 0; n1 < var1.length; ++n1)
                     for (let n2 = 0; n2 < var2.length; ++n2)
                        this.fill3DHistogram(var0.get(n0), var1.get(n1), var2.get(n2), cut.value);
               break;
         }
      }

      if (this.monitoring && this.hist && !this.dump_values) {
         let now = new Date().getTime();
         if (now - this.lasttm > this.monitoring) {
            this.lasttm = now;
            if (this.progress_callback)
               this.progress_callback(this.hist);
         }
      }
   }

   /** @summary Normal TSelector Terminate handler */
   TDrawSelector.prototype.Terminate = function(res) {
      if (res && !this.hist) this.createHistogram();

      this.ShowProgress();

      if (this.result_callback)
         this.result_callback(this.hist);
   }

   // ======================================================================

   /** @summary return TStreamerElement associated with the branch - if any
     * @desc unfortunately, branch.fID is not number of element in streamer info
     * @private */
   function findBrachStreamerElement(branch, file) {
      if (!branch || !file || (branch._typename !== "TBranchElement") || (branch.fID < 0) || (branch.fStreamerType < 0)) return null;

      let s_i = file.findStreamerInfo(branch.fClassName, branch.fClassVersion, branch.fCheckSum),
         arr = (s_i && s_i.fElements) ? s_i.fElements.arr : null;
      if (!arr) return null;

      let match_name = branch.fName,
         pos = match_name.indexOf("[");
      if (pos > 0) match_name = match_name.substr(0, pos);
      pos = match_name.lastIndexOf(".");
      if (pos > 0) match_name = match_name.substr(pos + 1);

      function match_elem(elem) {
         if (!elem) return false;
         if (elem.fName !== match_name) return false;
         if (elem.fType === branch.fStreamerType) return true;
         if ((elem.fType === jsrio.kBool) && (branch.fStreamerType === jsrio.kUChar)) return true;
         if (((branch.fStreamerType === jsrio.kSTL) || (branch.fStreamerType === jsrio.kSTL + jsrio.kOffsetL) ||
            (branch.fStreamerType === jsrio.kSTLp) || (branch.fStreamerType === jsrio.kSTLp + jsrio.kOffsetL))
            && (elem.fType === jsrio.kStreamer)) return true;
         console.warn('Should match element', elem.fType, 'with branch', branch.fStreamerType);
         return false;
      }

      // first check branch fID - in many cases gut guess
      if (match_elem(arr[branch.fID])) return arr[branch.fID];

      // console.warn('Missmatch with branch name and extracted element', branch.fName, match_name, (s_elem ? s_elem.fName : "---"));

      for (let k = 0; k < arr.length; ++k)
         if ((k !== branch.fID) && match_elem(arr[k])) return arr[k];

      console.error('Did not found/match element for branch', branch.fName, 'class', branch.fClassName);

      return null;
   }

   /** @summary return type name of given member in the class
     * @private */
   function defineMemberTypeName(file, parent_class, member_name) {
      let s_i = file.findStreamerInfo(parent_class),
         arr = (s_i && s_i.fElements) ? s_i.fElements.arr : null,
         elem = null;
      if (!arr) return "";

      for (let k = 0; k < arr.length; ++k) {
         if (arr[k].fTypeName === "BASE") {
            let res = defineMemberTypeName(file, arr[k].fName, member_name);
            if (res) return res;
         } else
            if (arr[k].fName === member_name) { elem = arr[k]; break; }
      }

      if (!elem) return "";

      let clname = elem.fTypeName;
      if (clname[clname.length - 1] === "*") clname = clname.substr(0, clname.length - 1);

      return clname;
   }

   /** @summary return class name of the object, stored in the branch
     * @private */
   function getBranchObjectClass(branch, tree, with_clones, with_leafs) {

      if (!branch || (branch._typename !== "TBranchElement")) return "";

      if ((branch.fType === JSROOT.BranchType.kLeafNode) && (branch.fID === -2) && (branch.fStreamerType === -1)) {
         // object where all sub-branches will be collected
         return branch.fClassName;
      }

      if (with_clones && branch.fClonesName && ((branch.fType === JSROOT.BranchType.kClonesNode) || (branch.fType === JSROOT.BranchType.kSTLNode)))
         return branch.fClonesName;

      let s_elem = findBrachStreamerElement(branch, tree.$file);

      if ((branch.fType === JSROOT.BranchType.kBaseClassNode) && s_elem && (s_elem.fTypeName === "BASE"))
         return s_elem.fName;

      if (branch.fType === JSROOT.BranchType.kObjectNode) {
         if (s_elem && ((s_elem.fType === jsrio.kObject) || (s_elem.fType === jsrio.kAny)))
            return s_elem.fTypeName;
         return "TObject";
      }

      if ((branch.fType === JSROOT.BranchType.kLeafNode) && s_elem && with_leafs) {
         if ((s_elem.fType === jsrio.kObject) || (s_elem.fType === jsrio.kAny)) return s_elem.fTypeName;
         if (s_elem.fType === jsrio.kObjectp) return s_elem.fTypeName.substr(0, s_elem.fTypeName.length - 1);
      }

      return "";
   }

   /** @summary create fast list to assign all methods to the object
     * @private */
   function makeMethodsList(typename) {
      let methods = JSROOT.getMethods(typename);

      let res = {
         names: [],
         values: [],
         Create: function() {
            let obj = {};
            for (let n = 0; n < this.names.length; ++n)
               obj[this.names[n]] = this.values[n];
            return obj;
         }
      }

      res.names.push("_typename"); res.values.push(typename);
      for (let key in methods) {
         res.names.push(key);
         res.values.push(methods[key]);
      }
      return res;
   }

   /** @summary try to define classname for the branch member, scanning list of branches
     * @private */
   function detectBranchMemberClass(brlst, prefix, start) {
      let clname = "";
      for (let kk = (start || 0); kk < brlst.arr.length; ++kk)
         if ((brlst.arr[kk].fName.indexOf(prefix) === 0) && brlst.arr[kk].fClassName)
            clname = brlst.arr[kk].fClassName;
      return clname;
   }

   /**
    * @summary JSROOT methods for ROOT TTree class
    *
    * @class
    * @alias JSROOT.TTree
    * @hideconstructor
    * @desc TTree only can be read from the existing ROOT file, there is no possibility to create and fill tree
    * @example
    * JSROOT.openFile("https://root.cern/js/files/hsimple.root")
    *       .then(file => file.readObject("ntuple;1"))
    *       .then(tree => JSROOT.draw("drawing", tree, "px:py::pz>5"));
    */
   let TTreeMethods = {};

   /** @summary Process selector
     * @param {object} selector - instance of {@link JSROOT.TSelector} class
     * @param {object} [args] - different arguments
     * @param {number} [args.firstentry] - first entry to process, 0 when not specified
     * @param {number} [args.numentries] - number of entries to process, all when not specified
     * @returns {Promise} with TSelector instance */
   TTreeMethods.Process = function(selector, args) {
      // function similar to the TTree::Process

      if (!args) args = {};

      if (!selector || !this.$file || !selector.branches) {
         if (selector) selector.Terminate(false);
         return Promise.reject(Error("required parameter missing for TTree::Process"));
      }

      // central handle with all information required for reading
      let handle = {
         tree: this, // keep tree reference
         file: this.$file, // keep file reference
         selector: selector, // reference on selector
         arr: [], // list of branches
         curr: -1,  // current entry ID
         current_entry: -1, // current processed entry
         simple_read: true, // all baskets in all used branches are in sync,
         process_arrays: true // one can process all branches as arrays
      };

      let namecnt = 0;

      function CreateLeafElem(leaf, name) {
         // function creates TStreamerElement which corresponds to the elementary leaf
         let datakind = 0;
         switch (leaf._typename) {
            case 'TLeafF': datakind = jsrio.kFloat; break;
            case 'TLeafD': datakind = jsrio.kDouble; break;
            case 'TLeafO': datakind = jsrio.kBool; break;
            case 'TLeafB': datakind = leaf.fIsUnsigned ? jsrio.kUChar : jsrio.kChar; break;
            case 'TLeafS': datakind = leaf.fIsUnsigned ? jsrio.kUShort : jsrio.kShort; break;
            case 'TLeafI': datakind = leaf.fIsUnsigned ? jsrio.kUInt : jsrio.kInt; break;
            case 'TLeafL': datakind = leaf.fIsUnsigned ? jsrio.kULong64 : jsrio.kLong64; break;
            case 'TLeafC': datakind = jsrio.kTString; break; // datakind = leaf.fIsUnsigned ? jsrio.kUChar : jsrio.kChar; break;
            default: return null;
         }
         return jsrio.createStreamerElement(name || leaf.fName, datakind);
      }

      function FindInHandle(branch) {
         for (let k = 0; k < handle.arr.length; ++k)
            if (handle.arr[k].branch === branch) return handle.arr[k];
         return null;
      }

      function AddBranchForReading(branch, target_object, target_name, read_mode) {
         // central method to add branch for reading
         // read_mode == true - read only this branch
         // read_mode == '$child$' is just member of object from for STL or clonesarray
         // read_mode == '<any class name>' is sub-object from STL or clonesarray, happens when such new object need to be created
         // read_mode == '.member_name' select only reading of member_name instead of complete object

         if (typeof branch === 'string')
            branch = handle.tree.FindBranch(branch);

         if (!branch) { console.error('Did not found branch'); return null; }

         let item = FindInHandle(branch);

         if (item) {
            console.error('Branch already configured for reading', branch.fName);
            if (item.tgt !== target_object) console.error('Target object differs');
            return null;
         }

         if (!branch.fEntries) {
            console.warn('Branch ', branch.fName, ' does not have entries');
            return null;
         }

         // console.log('Add branch', branch.fName);

         item = {
            branch: branch,
            tgt: target_object, // used target object - can be differ for object members
            name: target_name,
            index: -1, // index in the list of read branches
            member: null, // member to read branch
            type: 0, // keep type identifier
            curr_entry: -1, // last processed entry
            raw: null, // raw buffer for reading
            basket: null, // current basket object
            curr_basket: 0,  // number of basket used for processing
            read_entry: -1,  // last entry which is already read
            staged_entry: -1, // entry which is staged for reading
            first_readentry: -1, // first entry to read
            staged_basket: 0,  // last basket staged for reading
            numentries: branch.fEntries,
            numbaskets: branch.fWriteBasket, // number of baskets which can be read from the file
            counters: null, // branch indexes used as counters
            ascounter: [], // list of other branches using that branch as counter
            baskets: [], // array for read baskets,
            staged_prev: 0, // entry limit of previous I/O request
            staged_now: 0, // entry limit of current I/O request
            progress_showtm: 0, // last time when progress was showed
            GetBasketEntry: function(k) {
               if (!this.branch || (k > this.branch.fMaxBaskets)) return 0;
               let res = (k < this.branch.fMaxBaskets) ? this.branch.fBasketEntry[k] : 0;
               if (res) return res;
               let bskt = (k > 0) ? this.branch.fBaskets.arr[k - 1] : null;
               return bskt ? (this.branch.fBasketEntry[k - 1] + bskt.fNevBuf) : 0;
            },
            GetTarget: function(tgtobj) {
               // returns target object which should be used for the branch reading
               if (!this.tgt) return tgtobj;
               for (let k = 0; k < this.tgt.length; ++k) {
                  let sub = this.tgt[k];
                  if (!tgtobj[sub.name]) tgtobj[sub.name] = sub.lst.Create();
                  tgtobj = tgtobj[sub.name];
               }
               return tgtobj;
            },
            GetEntry: function(entry) {
               // This should be equivalent to TBranch::GetEntry() method
               let shift = entry - this.first_entry, off;
               if (!this.branch.TestBit(jsrio.BranchBits.kDoNotUseBufferMap))
                  this.raw.clearObjectMap();
               if (this.basket.fEntryOffset) {
                  off = this.basket.fEntryOffset[shift];
                  if (this.basket.fDisplacement)
                     this.raw.fDisplacement = this.basket.fDisplacement[shift];
               } else {
                  off = this.basket.fKeylen + this.basket.fNevBufSize * shift;
               }
               this.raw.locate(off - this.raw.raw_shift);

               // this.member.func(this.raw, this.GetTarget(tgtobj));
            }
         };

         // last basket can be stored directly with the branch
         while (item.GetBasketEntry(item.numbaskets + 1)) item.numbaskets++;

         // check all counters if we
         let nb_leaves = branch.fLeaves ? branch.fLeaves.arr.length : 0,
            leaf = (nb_leaves > 0) ? branch.fLeaves.arr[0] : null,
            elem = null, // TStreamerElement used to create reader
            member = null, // member for actual reading of the branch
            is_brelem = (branch._typename === "TBranchElement"),
            child_scan = 0, // scan child branches after main branch is appended
            item_cnt = null, item_cnt2 = null, object_class = "";

         if (branch.fBranchCount) {

            item_cnt = FindInHandle(branch.fBranchCount);

            if (!item_cnt)
               item_cnt = AddBranchForReading(branch.fBranchCount, target_object, "$counter" + namecnt++, true);

            if (!item_cnt) { console.error('Cannot add counter branch', branch.fBranchCount.fName); return null; }

            let BranchCount2 = branch.fBranchCount2;

            if (!BranchCount2 && (branch.fBranchCount.fStreamerType === jsrio.kSTL) &&
               ((branch.fStreamerType === jsrio.kStreamLoop) || (branch.fStreamerType === jsrio.kOffsetL + jsrio.kStreamLoop))) {
               // special case when count member from kStreamLoop not assigned as fBranchCount2
               let elemd = findBrachStreamerElement(branch, handle.file),
                  arrd = branch.fBranchCount.fBranches.arr;

               if (elemd && elemd.fCountName && arrd)
                  for (let k = 0; k < arrd.length; ++k)
                     if (arrd[k].fName === branch.fBranchCount.fName + "." + elemd.fCountName) {
                        BranchCount2 = arrd[k];
                        break;
                     }

               if (!BranchCount2) console.error('Did not found branch for second counter of kStreamLoop element');
            }

            if (BranchCount2) {
               item_cnt2 = FindInHandle(BranchCount2);

               if (!item_cnt2) item_cnt2 = AddBranchForReading(BranchCount2, target_object, "$counter" + namecnt++, true);

               if (!item_cnt2) { console.error('Cannot add counter branch2', BranchCount2.fName); return null; }
            }
         } else
            if (nb_leaves === 1 && leaf && leaf.fLeafCount) {
               let br_cnt = handle.tree.FindBranch(leaf.fLeafCount.fName);

               if (br_cnt) {
                  item_cnt = FindInHandle(br_cnt);

                  if (!item_cnt) item_cnt = AddBranchForReading(br_cnt, target_object, "$counter" + namecnt++, true);

                  if (!item_cnt) { console.error('Cannot add counter branch', br_cnt.fName); return null; }
               }
            }

         function ScanBranches(lst, master_target, chld_kind) {
            if (!lst || !lst.arr.length) return true;

            let match_prefix = branch.fName;
            if (match_prefix[match_prefix.length - 1] === ".") match_prefix = match_prefix.substr(0, match_prefix.length - 1);
            if ((typeof read_mode === "string") && (read_mode[0] == ".")) match_prefix += read_mode;
            match_prefix += ".";

            for (let k = 0; k < lst.arr.length; ++k) {
               let br = lst.arr[k];
               if ((chld_kind > 0) && (br.fType !== chld_kind)) continue;

               if (br.fType === JSROOT.BranchType.kBaseClassNode) {
                  if (!ScanBranches(br.fBranches, master_target, chld_kind)) return false;
                  continue;
               }

               let elem = findBrachStreamerElement(br, handle.file);
               if (elem && (elem.fTypeName === "BASE")) {
                  // if branch is data of base class, map it to original target
                  if (br.fTotBytes && !AddBranchForReading(br, target_object, target_name, read_mode)) return false;
                  if (!ScanBranches(br.fBranches, master_target, chld_kind)) return false;
                  continue;
               }

               let subname = br.fName, chld_direct = 1;

               if (br.fName.indexOf(match_prefix) === 0) {
                  subname = subname.substr(match_prefix.length);
               } else {
                  if (chld_kind > 0) continue; // for defined children names prefix must be present
               }

               let p = subname.indexOf('[');
               if (p > 0) subname = subname.substr(0, p);
               p = subname.indexOf('<');
               if (p > 0) subname = subname.substr(0, p);

               if (chld_kind > 0) {
                  chld_direct = "$child$";
                  let pp = subname.indexOf(".");
                  if (pp > 0) chld_direct = detectBranchMemberClass(lst, branch.fName + "." + subname.substr(0, pp + 1), k) || "TObject";
               }

               if (!AddBranchForReading(br, master_target, subname, chld_direct)) return false;
            }

            return true;
         }

         if (branch._typename === "TBranchObject") {
            member = {
               name: target_name,
               typename: branch.fClassName,
               virtual: leaf.fVirtual,
               func: function(buf, obj) {
                  let clname = this.typename;
                  if (this.virtual) clname = buf.readFastString(buf.ntou1() + 1);
                  obj[this.name] = buf.classStreamer({}, clname);
               }
            };
         } else

            if ((branch.fType === JSROOT.BranchType.kClonesNode) || (branch.fType === JSROOT.BranchType.kSTLNode)) {

               elem = jsrio.createStreamerElement(target_name, jsrio.kInt);

               if (!read_mode || ((typeof read_mode === "string") && (read_mode[0] === ".")) || (read_mode === 1)) {
                  handle.process_arrays = false;

                  member = {
                     name: target_name,
                     conttype: branch.fClonesName || "TObject",
                     reallocate: args.reallocate_objects,
                     func: function(buf, obj) {
                        let size = buf.ntoi4(), n = 0, arr = obj[this.name];
                        if (!arr || this.reallocate) {
                           arr = obj[this.name] = new Array(size);
                        } else {
                           n = arr.length;
                           arr.length = size; // reallocate array
                        }

                        while (n < size) arr[n++] = this.methods.Create(); // create new objects
                     }
                  }

                  if ((typeof read_mode === "string") && (read_mode[0] === ".")) {
                     member.conttype = detectBranchMemberClass(branch.fBranches, branch.fName + read_mode);
                     if (!member.conttype) {
                        console.error('Cannot select object', read_mode, "in the branch", branch.fName);
                        return null;
                     }
                  }

                  member.methods = makeMethodsList(member.conttype);

                  child_scan = (branch.fType === JSROOT.BranchType.kClonesNode) ? JSROOT.BranchType.kClonesMemberNode : JSROOT.BranchType.kSTLMemberNode;
               }
            } else

               if ((object_class = getBranchObjectClass(branch, handle.tree))) {

                  if (read_mode === true) {
                     console.warn('Object branch ' + object_class + ' can not have data to be read directly');
                     return null;
                  }

                  handle.process_arrays = false;

                  let newtgt = new Array(target_object ? (target_object.length + 1) : 1);
                  for (let l = 0; l < newtgt.length - 1; ++l) newtgt[l] = target_object[l];
                  newtgt[newtgt.length - 1] = { name: target_name, lst: makeMethodsList(object_class) };

                  if (!ScanBranches(branch.fBranches, newtgt, 0)) return null;

                  return item; // this kind of branch does not have baskets and not need to be read

               } else if (is_brelem && (nb_leaves === 1) && (leaf.fName === branch.fName) && (branch.fID == -1)) {

                  elem = jsrio.createStreamerElement(target_name, branch.fClassName);

                  if (elem.fType === jsrio.kAny) {

                     let streamer = handle.file.getStreamer(branch.fClassName, { val: branch.fClassVersion, checksum: branch.fCheckSum });
                     if (!streamer) { elem = null; console.warn('not found streamer!'); } else
                        member = {
                           name: target_name,
                           typename: branch.fClassName,
                           streamer: streamer,
                           func: function(buf, obj) {
                              let res = { _typename: this.typename };
                              for (let n = 0; n < this.streamer.length; ++n)
                                 this.streamer[n].func(buf, res);
                              obj[this.name] = res;
                           }
                        };
                  }

                  // elem.fType = jsrio.kAnyP;

                  // only STL containers here
                  // if (!elem.fSTLtype) elem = null;
               } else if (is_brelem && (nb_leaves <= 1)) {

                  elem = findBrachStreamerElement(branch, handle.file);

                  // this is basic type - can try to solve problem differently
                  if (!elem && branch.fStreamerType && (branch.fStreamerType < 20))
                     elem = jsrio.createStreamerElement(target_name, branch.fStreamerType);

               } else if (nb_leaves === 1) {
                  // no special constrains for the leaf names

                  elem = CreateLeafElem(leaf, target_name);

               } else if ((branch._typename === "TBranch") && (nb_leaves > 1)) {
                  // branch with many elementary leaves

                  let arr = new Array(nb_leaves), isok = true;
                  for (let l = 0; l < nb_leaves; ++l) {
                     arr[l] = CreateLeafElem(branch.fLeaves.arr[l]);
                     arr[l] = jsrio.createMember(arr[l], handle.file);
                     if (!arr[l]) isok = false;
                  }

                  if (isok)
                     member = {
                        name: target_name,
                        leaves: arr,
                        func: function(buf, obj) {
                           let tgt = obj[this.name], l = 0;
                           if (!tgt) obj[this.name] = tgt = {};
                           while (l < this.leaves.length)
                              this.leaves[l++].func(buf, tgt);
                        }
                     }
               }

         if (!elem && !member) {
            console.warn('Not supported branch kind', branch.fName, branch._typename);
            return null;
         }

         if (!member) {
            member = jsrio.createMember(elem, handle.file);

            if ((member.base !== undefined) && member.basename) {
               // when element represent base class, we need handling which differ from normal IO
               member.func = function(buf, obj) {
                  if (!obj[this.name]) obj[this.name] = { _typename: this.basename };
                  buf.classStreamer(obj[this.name], this.basename);
               };
            }
         }

         if (item_cnt && (typeof read_mode === "string")) {

            member.name0 = item_cnt.name;

            let snames = target_name.split(".");

            if (snames.length === 1) {
               // no point in the name - just plain array of objects
               member.get = function(arr, n) { return arr[n]; }
            } else if (read_mode === "$child$") {
               console.error('target name contains point, but suppose to be direct child', target_name);
               return null;
            } else if (snames.length === 2) {
               target_name = member.name = snames[1];
               member.name1 = snames[0];
               member.subtype1 = read_mode;
               member.methods1 = makeMethodsList(member.subtype1);
               member.get = function(arr, n) {
                  let obj1 = arr[n][this.name1];
                  if (!obj1) obj1 = arr[n][this.name1] = this.methods1.Create();
                  return obj1;
               }
            } else {
               // very complex task - we need to reconstruct several embedded members with their types
               // try our best - but not all data types can be reconstructed correctly
               // while classname is not enough - there can be different versions

               if (!branch.fParentName) {
                  console.error('Not possible to provide more than 2 parts in the target name', target_name);
                  return null;
               }

               target_name = member.name = snames.pop(); // use last element
               member.snames = snames; // remember all sub-names
               member.smethods = []; // and special handles to create missing objects

               let parent_class = branch.fParentName; // unfortunately, without version

               for (let k = 0; k < snames.length; ++k) {
                  let chld_class = defineMemberTypeName(handle.file, parent_class, snames[k]);
                  member.smethods[k] = makeMethodsList(chld_class || "AbstractClass");
                  parent_class = chld_class;
               }
               member.get = function(arr, n) {
                  let obj1 = arr[n][this.snames[0]];
                  if (!obj1) obj1 = arr[n][this.snames[0]] = this.smethods[0].Create();
                  for (let k = 1; k < this.snames.length; ++k) {
                     let obj2 = obj1[this.snames[k]];
                     if (!obj2) obj2 = obj1[this.snames[k]] = this.smethods[k].Create();
                     obj1 = obj2;
                  }
                  return obj1;
               }
            }

            // case when target is sub-object and need to be created before


            if (member.objs_branch_func) {
               // STL branch provides special function for the reading
               member.func = member.objs_branch_func;
            } else {
               member.func0 = member.func;

               member.func = function(buf, obj) {
                  let arr = obj[this.name0], n = 0; // objects array where reading is done
                  while (n < arr.length)
                     this.func0(buf, this.get(arr, n++)); // read all individual object with standard functions
               }
            }

         } else if (item_cnt) {

            handle.process_arrays = false;

            if ((elem.fType === jsrio.kDouble32) || (elem.fType === jsrio.kFloat16)) {
               // special handling for compressed floats

               member.stl_size = item_cnt.name;
               member.func = function(buf, obj) {
                  obj[this.name] = this.readarr(buf, obj[this.stl_size]);
               }

            } else
               if (((elem.fType === jsrio.kOffsetP + jsrio.kDouble32) || (elem.fType === jsrio.kOffsetP + jsrio.kFloat16)) && branch.fBranchCount2) {
                  // special handling for variable arrays of compressed floats in branch - not tested

                  member.stl_size = item_cnt.name;
                  member.arr_size = item_cnt2.name;
                  member.func = function(buf, obj) {
                     let sz0 = obj[this.stl_size], sz1 = obj[this.arr_size], arr = new Array(sz0);
                     for (let n = 0; n < sz0; ++n)
                        arr[n] = (buf.ntou1() === 1) ? this.readarr(buf, sz1[n]) : [];
                     obj[this.name] = arr;
                  }

               } else
                  // special handling of simple arrays
                  if (((elem.fType > 0) && (elem.fType < jsrio.kOffsetL)) || (elem.fType === jsrio.kTString) ||
                     (((elem.fType > jsrio.kOffsetP) && (elem.fType < jsrio.kOffsetP + jsrio.kOffsetL)) && branch.fBranchCount2)) {

                     member = {
                        name: target_name,
                        stl_size: item_cnt.name,
                        type: elem.fType,
                        func: function(buf, obj) {
                           obj[this.name] = buf.readFastArray(obj[this.stl_size], this.type);
                        }
                     };

                     if (branch.fBranchCount2) {
                        member.type -= jsrio.kOffsetP;
                        member.arr_size = item_cnt2.name;
                        member.func = function(buf, obj) {
                           let sz0 = obj[this.stl_size], sz1 = obj[this.arr_size], arr = new Array(sz0);
                           for (let n = 0; n < sz0; ++n)
                              arr[n] = (buf.ntou1() === 1) ? buf.readFastArray(sz1[n], this.type) : [];
                           obj[this.name] = arr;
                        }
                     }

                  } else
                     if ((elem.fType > jsrio.kOffsetP) && (elem.fType < jsrio.kOffsetP + jsrio.kOffsetL) && member.cntname) {

                        member.cntname = item_cnt.name;
                     } else
                        if (elem.fType == jsrio.kStreamer) {
                           // with streamers one need to extend existing array

                           if (item_cnt2)
                              throw new Error('Second branch counter not supported yet with jsrio.kStreamer');

                           // function provided by normal I/O
                           member.func = member.branch_func;
                           member.stl_size = item_cnt.name;
                        } else
                           if ((elem.fType === jsrio.kStreamLoop) || (elem.fType === jsrio.kOffsetL + jsrio.kStreamLoop)) {
                              if (item_cnt2) {
                                 // special solution for kStreamLoop
                                 member.stl_size = item_cnt.name;
                                 member.cntname = item_cnt2.name;
                                 member.func = member.branch_func; // this is special function, provided by base I/O
                              } else {
                                 member.cntname = item_cnt.name;
                              }
                           } else {

                              member.name = "$stl_member";

                              let loop_size_name;

                              if (item_cnt2) {
                                 if (member.cntname) {
                                    loop_size_name = item_cnt2.name;
                                    member.cntname = "$loop_size";
                                 } else {
                                    throw new Error('Second branch counter not used - very BAD');
                                 }
                              }

                              let stlmember = {
                                 name: target_name,
                                 stl_size: item_cnt.name,
                                 loop_size: loop_size_name,
                                 member0: member,
                                 func: function(buf, obj) {
                                    let cnt = obj[this.stl_size], arr = new Array(cnt);
                                    for (let n = 0; n < cnt; ++n) {
                                       if (this.loop_size) obj.$loop_size = obj[this.loop_size][n];
                                       this.member0.func(buf, obj);
                                       arr[n] = obj.$stl_member;
                                    }
                                    delete obj.$stl_member;
                                    delete obj.$loop_size;
                                    obj[this.name] = arr;
                                 }
                              };

                              member = stlmember;
                           }
         }

         // set name used to store result
         member.name = target_name;

         item.member = member; // member for reading
         if (elem) item.type = elem.fType;
         item.index = handle.arr.length; // index in the global list of branches

         if (item_cnt) {
            item.counters = [item_cnt.index];
            item_cnt.ascounter.push(item.index);

            if (item_cnt2) {
               item.counters.push(item_cnt2.index);
               item_cnt2.ascounter.push(item.index);
            }
         }

         handle.arr.push(item);

         // now one should add all other child branches
         if (child_scan)
            if (!ScanBranches(branch.fBranches, target_object, child_scan)) return null;

         return item;
      }

      // main loop to add all branches from selector for reading
      for (let nn = 0; nn < selector.branches.length; ++nn) {

         let item = AddBranchForReading(selector.branches[nn], undefined, selector.names[nn], selector.directs[nn]);

         if (!item) {
            selector.Terminate(false);
            return Promise.reject(Error("Fail to add branch ${selector.branches[nn]}"));
         }
      }

      // check if simple reading can be performed and there are direct data in branch

      for (let h = 1; (h < handle.arr.length) && handle.simple_read; ++h) {

         let item = handle.arr[h], item0 = handle.arr[0];

         if ((item.numentries !== item0.numentries) || (item.numbaskets !== item0.numbaskets)) handle.simple_read = false;
         for (let n = 0; n < item.numbaskets; ++n)
            if (item.GetBasketEntry(n) !== item0.GetBasketEntry(n)) handle.simple_read = false;
      }

      // now calculate entries range

      handle.firstentry = handle.lastentry = 0;
      for (let nn = 0; nn < handle.arr.length; ++nn) {
         let branch = handle.arr[nn].branch, e1 = branch.fFirstEntry;
         if (e1 === undefined) e1 = (branch.fBasketBytes[0] ? branch.fBasketEntry[0] : 0);
         handle.firstentry = Math.max(handle.firstentry, e1);
         handle.lastentry = (nn === 0) ? (e1 + branch.fEntries) : Math.min(handle.lastentry, e1 + branch.fEntries);
      }

      if (handle.firstentry >= handle.lastentry) {
         selector.Terminate(false);
         return Promise.reject(Error('No any common events for selected branches'));
      }

      handle.process_min = handle.firstentry;
      handle.process_max = handle.lastentry;

      let resolveFunc, rejectFunc; // Promise methods

      if (!isNaN(args.firstentry) && (args.firstentry > handle.firstentry) && (args.firstentry < handle.lastentry))
         handle.process_min = args.firstentry;

      handle.current_entry = handle.staged_now = handle.process_min;

      if (!isNaN(args.numentries) && (args.numentries > 0)) {
         let max = handle.process_min + args.numentries;
         if (max < handle.process_max) handle.process_max = max;
      }

      if ((typeof selector.ProcessArrays === 'function') && handle.simple_read) {
         // this is indication that selector can process arrays of values
         // only strictly-matched tree structure can be used for that

         for (let k = 0; k < handle.arr.length; ++k) {
            let elem = handle.arr[k];
            if ((elem.type <= 0) || (elem.type >= jsrio.kOffsetL) || (elem.type === jsrio.kCharStar)) handle.process_arrays = false;
         }

         if (handle.process_arrays) {
            // create other members for fast processing

            selector.tgtarr = {}; // object with arrays

            for (let nn = 0; nn < handle.arr.length; ++nn) {
               let item = handle.arr[nn],
                  elem = jsrio.createStreamerElement(item.name, item.type);

               elem.fType = item.type + jsrio.kOffsetL;
               elem.fArrayLength = 10;
               elem.fArrayDim = 1;
               elem.fMaxIndex[0] = 10; // 10 if artificial number, will be replaced during reading

               item.arrmember = jsrio.createMember(elem, handle.file);
            }
         }
      } else {
         handle.process_arrays = false;
      }

      /** read basket with tree data, selecting different files */
      function ReadBaskets(bitems) {

         function ExtractPlaces() {
            // extract places to read and define file name

            let places = [], filename = "";

            for (let n = 0; n < bitems.length; ++n) {
               if (bitems[n].done) continue;

               let branch = bitems[n].branch;

               if (places.length === 0)
                  filename = branch.fFileName;
               else if (filename !== branch.fFileName)
                  continue;

               bitems[n].selected = true; // mark which item was selected for reading

               places.push(branch.fBasketSeek[bitems[n].basket], branch.fBasketBytes[bitems[n].basket]);
            }

            return places.length > 0 ? { places: places, filename: filename } : null;
         }

         function ReadProgress(value) {

            if ((handle.staged_prev === handle.staged_now) ||
               (handle.process_max <= handle.process_min)) return;

            let tm = new Date().getTime();

            if (tm - handle.progress_showtm < 500) return; // no need to show very often

            handle.progress_showtm = tm;

            let portion = (handle.staged_prev + value * (handle.staged_now - handle.staged_prev)) /
               (handle.process_max - handle.process_min);

            handle.selector.ShowProgress(portion);
         }

         function ProcessBlobs(blobs, places) {
            if (!blobs || ((places.length > 2) && (blobs.length * 2 !== places.length)))
               return Promise.resolve(null);

            if (places.length == 2) blobs = [ blobs ];

            function DoProcessing(k) {

               for (; k < bitems.length; ++k) {
                  if (!bitems[k].selected) continue;

                  bitems[k].selected = false;
                  bitems[k].done = true;

                  let blob = blobs.shift(),
                     buf = new JSROOT.TBuffer(blob, 0, handle.file),
                     basket = buf.classStreamer({}, "TBasket");

                  if (basket.fNbytes !== bitems[k].branch.fBasketBytes[bitems[k].basket])
                     console.error('mismatch in read basket sizes', bitems[k].branch.fBasketBytes[bitems[k].basket]);

                  // items[k].obj = basket; // keep basket object itself if necessary

                  bitems[k].bskt_obj = basket; // only number of entries in the basket are relevant for the moment

                  if (basket.fKeylen + basket.fObjlen === basket.fNbytes) {
                     // use data from original blob
                     buf.raw_shift = 0;

                     bitems[k].raw = buf; // here already unpacked buffer

                    if (bitems[k].branch.fEntryOffsetLen > 0)
                        buf.readBasketEntryOffset(basket, buf.raw_shift);

                    continue;
                  }

                  // unpack data and create new blob
                  return jsrio.R__unzip(blob, basket.fObjlen, false, buf.o).then(objblob => {

                     if (objblob) {
                        buf = new JSROOT.TBuffer(objblob, 0, handle.file);
                        buf.raw_shift = basket.fKeylen;
                        buf.fTagOffset = basket.fKeylen;
                     } else {
                        throw new Error('FAIL TO UNPACK');
                     }

                     bitems[k].raw = buf; // here already unpacked buffer

                     if (bitems[k].branch.fEntryOffsetLen > 0)
                        buf.readBasketEntryOffset(basket, buf.raw_shift);

                     return DoProcessing(k+1);  // continue processing
                  });
               }

               let req = ExtractPlaces();

               if (req)
                  return handle.file.readBuffer(req.places, req.filename, ReadProgress).then(blobs => ProcessBlobs(blobs)).catch(() => { return null; });

               return Promise.resolve(bitems);
             }

             return DoProcessing(0);
         }

         let req = ExtractPlaces();

         // extract places where to read
         if (req)
            return handle.file.readBuffer(req.places, req.filename, ReadProgress).then(blobs => ProcessBlobs(blobs, req.places)).catch(() => { return null; });

         return Promise.resolve(null);
      }

      function ReadNextBaskets() {

         let totalsz = 0, bitems = [], isany = true, is_direct = false, min_staged = handle.process_max;

         while ((totalsz < 1e6) && isany) {
            isany = false;
            // very important, loop over branches in reverse order
            // let check counter branch after reading of normal branch is prepared
            for (let n = handle.arr.length - 1; n >= 0; --n) {
               let elem = handle.arr[n];

               while (elem.staged_basket < elem.numbaskets) {

                  let k = elem.staged_basket++;

                  // no need to read more baskets, process_max is not included
                  if (elem.GetBasketEntry(k) >= handle.process_max) break;

                  // check which baskets need to be read
                  if (elem.first_readentry < 0) {
                     let lmt = elem.GetBasketEntry(k + 1),
                        not_needed = (lmt <= handle.process_min);

                     //for (let d=0;d<elem.ascounter.length;++d) {
                     //   let dep = handle.arr[elem.ascounter[d]]; // dependent element
                     //   if (dep.first_readentry < lmt) not_needed = false; // check that counter provide required data
                     //}

                     if (not_needed) continue; // if that basket not required, check next

                     elem.curr_basket = k; // basket where reading will start

                     elem.first_readentry = elem.GetBasketEntry(k); // remember which entry will be read first
                  }

                  // check if basket already loaded in the branch

                  let bitem = {
                     id: n, // to find which element we are reading
                     branch: elem.branch,
                     basket: k,
                     raw: null // here should be result
                  };

                  let bskt = elem.branch.fBaskets.arr[k];
                  if (bskt) {
                     bitem.raw = bskt.fBufferRef;
                     if (bitem.raw)
                        bitem.raw.locate(0); // reset pointer - same branch may be read several times
                     else
                        bitem.raw = new JSROOT.TBuffer(null, 0, handle.file); // create dummy buffer - basket has no data
                     bitem.raw.raw_shift = bskt.fKeylen;

                     if (bskt.fBufferRef && (elem.branch.fEntryOffsetLen > 0))
                        bitem.raw.readBasketEntryOffset(bskt, bitem.raw.raw_shift);

                     bitem.bskt_obj = bskt;
                     is_direct = true;
                     elem.baskets[k] = bitem;
                  } else {
                     bitems.push(bitem);
                     totalsz += elem.branch.fBasketBytes[k];
                     isany = true;
                  }

                  elem.staged_entry = elem.GetBasketEntry(k + 1);

                  min_staged = Math.min(min_staged, elem.staged_entry);

                  break;
               }
            }
         }

         if ((totalsz === 0) && !is_direct) {
            handle.selector.Terminate(true);
            return resolveFunc(handle.selector);
         }

         handle.staged_prev = handle.staged_now;
         handle.staged_now = min_staged;

         let portion = 0;
         if (handle.process_max > handle.process_min)
            portion = (handle.staged_prev - handle.process_min) / (handle.process_max - handle.process_min);

         handle.selector.ShowProgress(portion);

         handle.progress_showtm = new Date().getTime();

         if (totalsz > 0) return ReadBaskets(bitems).then(ProcessBaskets);

         if (is_direct) return ProcessBaskets([]); // directly process baskets

         throw new Error("No any data is requested - never come here");
      }

      function ProcessBaskets(bitems) {
         // this is call-back when next baskets are read

         if ((handle.selector.break_execution !== 0) || (bitems === null)) {
            handle.selector.Terminate(false);
            return resolveFunc(handle.selector);
         }

         // redistribute read baskets over branches
         for (let n = 0; n < bitems.length; ++n)
            handle.arr[bitems[n].id].baskets[bitems[n].basket] = bitems[n];

         // now process baskets

         let isanyprocessed = false;

         while (true) {

            let loopentries = 100000000, n, elem;

            // first loop used to check if all required data exists
            for (n = 0; n < handle.arr.length; ++n) {

               elem = handle.arr[n];

               if (!elem.raw || !elem.basket || (elem.first_entry + elem.basket.fNevBuf <= handle.current_entry)) {
                  delete elem.raw;
                  delete elem.basket;

                  if ((elem.curr_basket >= elem.numbaskets)) {
                     if (n == 0) {
                        handle.selector.Terminate(true);
                        return resolveFunc(handle.selector);
                     }
                     continue; // ignore non-master branch
                  }

                  // this is single response from the tree, includes branch, bakset number, raw data
                  let bitem = elem.baskets[elem.curr_basket];

                  // basket not read
                  if (!bitem) {
                     // no data, but no any event processed - problem
                     if (!isanyprocessed) {
                        handle.selector.Terminate(false);
                        return rejectFunc(Error(`no data for ${elem.branch.fName} basket ${elem.curr_basket}`));
                     }

                     // try to read next portion of tree data
                     return ReadNextBaskets();
                  }

                  elem.raw = bitem.raw;
                  elem.basket = bitem.bskt_obj;
                  // elem.nev = bitem.fNevBuf; // number of entries in raw buffer
                  elem.first_entry = elem.GetBasketEntry(bitem.basket);

                  bitem.raw = null; // remove reference on raw buffer
                  bitem.branch = null; // remove reference on the branch
                  bitem.bskt_obj = null; // remove reference on the branch
                  elem.baskets[elem.curr_basket++] = undefined; // remove from array
               }

               // define how much entries can be processed before next raw buffer will be finished
               loopentries = Math.min(loopentries, elem.first_entry + elem.basket.fNevBuf - handle.current_entry);
            }

            // second loop extracts all required data

            // do not read too much
            if (handle.current_entry + loopentries > handle.process_max)
               loopentries = handle.process_max - handle.current_entry;

            if (handle.process_arrays && (loopentries > 1)) {
               // special case - read all data from baskets as arrays

               for (n = 0; n < handle.arr.length; ++n) {
                  elem = handle.arr[n];

                  elem.GetEntry(handle.current_entry);

                  elem.arrmember.arrlength = loopentries;
                  elem.arrmember.func(elem.raw, handle.selector.tgtarr);

                  elem.raw = null;
               }

               handle.selector.ProcessArrays(handle.current_entry);

               handle.current_entry += loopentries;

               isanyprocessed = true;
            } else

               // main processing loop
               while (loopentries--) {

                  for (n = 0; n < handle.arr.length; ++n) {
                     elem = handle.arr[n];

                     // locate buffer offset at proper place
                     elem.GetEntry(handle.current_entry);

                     elem.member.func(elem.raw, elem.GetTarget(handle.selector.tgtobj));
                  }

                  handle.selector.Process(handle.current_entry);

                  handle.current_entry++;

                  isanyprocessed = true;
               }

            if (handle.current_entry >= handle.process_max) {
               handle.selector.Terminate(true);
               return resolveFunc(handle.selector);
            }
         }
      }

      return new Promise((resolve, reject) => {

         resolveFunc = resolve;
         rejectFunc = reject;

         // call begin before first entry is read
         handle.selector.Begin(this);

         ReadNextBaskets();

      });

   }

   /** @summary Special branch search
    * @desc Name can include extra part, which will be returned in the result
    * @param {string} name - name of the branch
    * @returns {Object} with "branch" and "rest" members
    * @private */
   TTreeMethods.FindBranchComplex = function(name, lst, only_search) {

      let top_search = false, search = name, res = null;

      if (lst === undefined) {
         top_search = true;
         lst = this.fBranches;
         let pos = search.indexOf("[");
         if (pos > 0) search = search.substr(0, pos);
      }

      if (!lst || (lst.arr.length === 0)) return null;

      for (let n = 0; n < lst.arr.length; ++n) {
         let brname = lst.arr[n].fName;
         if (brname[brname.length - 1] == "]")
            brname = brname.substr(0, brname.indexOf("["));

         // special case when branch name includes STL map name
         if ((search.indexOf(brname) !== 0) && (brname.indexOf("<") > 0)) {
            let p1 = brname.indexOf("<"), p2 = brname.lastIndexOf(">");
            brname = brname.substr(0, p1) + brname.substr(p2 + 1);
         }

         if (brname === search) { res = { branch: lst.arr[n], rest: "" }; break; }

         if (search.indexOf(brname) !== 0) continue;

         // this is a case when branch name is in the begin of the search string

         // check where point is
         let pnt = brname.length;
         if (brname[pnt - 1] === '.') pnt--;
         if (search[pnt] !== '.') continue;

         res = this.FindBranchComplex(search, lst.arr[n].fBranches);
         if (!res) res = this.FindBranchComplex(search.substr(pnt + 1), lst.arr[n].fBranches);

         if (!res) res = { branch: lst.arr[n], rest: search.substr(pnt) };

         break;
      }

      if (top_search && !only_search && !res && (search.indexOf("br_") == 0)) {
         let p = 3;
         while ((p < search.length) && (search[p] >= '0') && (search[p] <= '9')) ++p;
         let br = (p > 3) ? this.GetBranch(parseInt(search.slice(3,p))) : null;
         if (br) res = { branch: br, rest: search.substr(p) };
      }

      if (!top_search || !res) return res;

      if (name.length > search.length) res.rest += name.substr(search.length);

      return res;
   }

   /** @summary Search branch with specified name
    * @param {string} name - name of the branch
    * @returns {Object} found branch */
   TTreeMethods.FindBranch = function(name) {
      let res = this.FindBranchComplex(name, undefined, true);
      return (!res || (res.rest.length > 0)) ? null : res.branch;
   }

   /** @summary Returns number of branches in the TTree
    * @desc Checks also sub-branches in the branches
    * @returns {number} number of branches */
   TTreeMethods.GetNumBranches = function() {
      function Count(obj) {
         if (!obj || !obj.fBranches) return 0;
         let nchld = 0;
         obj.fBranches.arr.forEach(sub => nchld += Count(sub));
         return obj.fBranches.arr.length + nchld;
      }

      return Count(this);
   }

   /** @summary Get branch with specified id
     * @desc All sub-branches checked as well
     * @returns {Object} branch */
   TTreeMethods.GetBranch = function(id) {
      if ((id === undefined) || isNaN(id)) return;
      let res, seq = 0;
      function Scan(obj) {
         if (obj && obj.fBranches)
            obj.fBranches.arr.forEach(br => {
               if (seq++ === id) res = br;
               if (!res) Scan(br);
            });
      }

      Scan(this);
      return res;
   }

   /** @summary  JSROOT implementation of TTree::Draw
     * @disc in callback returns histogram and draw options
     * @param {object|string} args - different setting or simply draw expression
     * @param {string} args.expr - draw expression
     * @param {string} [args.cut=undefined]   - cut expression (also can be part of 'expr' after '::')
     * @param {string} [args.drawopt=undefined] - draw options for result histogram
     * @param {number} [args.firstentry=0] - first entry to process
     * @param {number} [args.numentries=undefined] - number of entries to process, all by default
     * @param {object} [args.branch=undefined] - TBranch object from TTree itself for the direct drawing
     * @param {function} [args.progress=undefined] - function called during histogram accumulation with argument { obj: draw_object, opt: draw_options }
     * @returns {Promise} with object like { obj: draw_object, opt: draw_options }
     */
   TTreeMethods.Draw = function(args) {

      if (typeof args === 'string') args = { expr: args };

      if (!args.expr) args.expr = "";

      // special debugging code
      if (args.expr === "testio")
         return this.IOTest(args);

      let selector = new TDrawSelector();

      if (args.branch) {
         if (!selector.drawOnlyBranch(this, args.branch, args.expr, args)) selector = null;
      } else {
         if (!selector.parseDrawExpression(this, args)) selector = null;
      }

      if (!selector)
         return Promise.reject(Error("Fail to create selector for specified expression"));

      return new Promise(resolve => {
         selector.setCallback(resolve, args.progress);
         this.Process(selector, args);
      });
   }

   /** @summary Performs generic I/O test for all branches in the TTree
     * @desc Used when "testio" draw option for TTree is specified
     * @private */
   TTreeMethods.IOTest = function(args) {
      args.branches = [];
      args.names = [];
      args.nchilds = [];
      args.nbr = 0;

      function CollectBranches(obj, prntname) {
         if (!obj || !obj.fBranches) return 0;

         let cnt = 0;

         for (let n = 0; n < obj.fBranches.arr.length; ++n) {
            let br = obj.fBranches.arr[n],
               name = (prntname ? prntname + "/" : "") + br.fName;
            args.branches.push(br);
            args.names.push(name);
            args.nchilds.push(0);
            let pos = args.nchilds.length - 1;
            cnt += br.fLeaves ? br.fLeaves.arr.length : 0;
            let nchld = CollectBranches(br, name);

            cnt += nchld;
            args.nchilds[pos] = nchld;

         }
         return cnt;
      }

      let numleaves = CollectBranches(this);

      args.names.push("Total are " + args.branches.length + " branches with " + numleaves + " leaves");

      args.lasttm = new Date().getTime();
      args.lastnbr = 0;

      let TestNextBranch = () => {

         let selector = new TSelector;

         selector.addBranch(args.branches[args.nbr], "br0");

         selector.Process = function() {
            if (this.tgtobj.br0 === undefined)
               this.fail = true;
         }

         selector.Terminate = function(res) {
            if (typeof res !== 'string')
               res = (!res || this.fails) ? "FAIL" : "ok";

            args.names[args.nbr] = res + " " + args.names[args.nbr];
            args.nbr++;

            if (args.nbr >= args.branches.length) {
               if (JSROOT.Painter) JSROOT.Painter.showProgress();
               return args.resolveFunc(args.names);
            }

            let now = new Date().getTime();

            if ((now - args.lasttm > 5000) || (args.nbr - args.lastnbr > 50)) {
               args.lasttm = now;
               args.lastnbr = args.nbr;
               setTimeout(TestNextBranch, 1); // use timeout to avoid deep recursion
            } else
               TestNextBranch();
         }

         if (JSROOT.Painter)
            JSROOT.Painter.showProgress("br " + args.nbr + "/" + args.branches.length + " " + args.names[args.nbr]);

         let br = args.branches[args.nbr],
            object_class = getBranchObjectClass(br, this),
            num = br.fEntries,
            skip_branch = (!br.fLeaves || (br.fLeaves.arr.length === 0));

         if (object_class) skip_branch = (args.nchilds[args.nbr] > 100);

         // skip_branch = args.nchilds[args.nbr]>1;

         if (skip_branch || (num <= 0)) {
            // ignore empty branches or objects with too-many subbranch
            // if (object_class) console.log('Ignore branch', br.fName, 'class', object_class, 'with', args.nchilds[args.nbr],'subbranches');
            selector.Terminate("ignore");
         } else {

            let drawargs = { numentries: 10 },
               first = br.fFirstEntry || 0,
               last = br.fEntryNumber || (first + num);

            if (num < drawargs.numentries) {
               drawargs.numentries = num;
            } else {
               // select randomly first entry to test I/O
               drawargs.firstentry = first + Math.round((last - first - drawargs.numentries) * Math.random());
            }

            // keep console output for debug purposes
            console.log('test branch', br.fName, 'first', (drawargs.firstentry || 0), "num", drawargs.numentries);

            this.Process(selector, drawargs);
         }
      }

      return new Promise(resolve => {
         args.resolveFunc = resolve;
         TestNextBranch();
      });
   }

   /** @summary Create hierarchy of TTree object
     * @private */
   JSROOT.treeHierarchy = function(node, obj) {
      if (obj._typename != 'TTree' && obj._typename != 'TNtuple' && obj._typename != 'TNtupleD' ) return false;

      function CreateBranchItem(node, branch, tree, parent_branch) {
         if (!node || !branch) return false;

         let nb_branches = branch.fBranches ? branch.fBranches.arr.length : 0,
             nb_leaves = branch.fLeaves ? branch.fLeaves.arr.length : 0;

         function ClearName(arg) {
            let pos = arg.indexOf("[");
            if (pos>0) arg = arg.substr(0, pos);
            if (parent_branch && arg.indexOf(parent_branch.fName)==0) {
               arg = arg.substr(parent_branch.fName.length);
               if (arg[0]===".") arg = arg.substr(1);
            }
            return arg;
         }

         branch.$tree = tree; // keep tree pointer, later do it more smart

         let subitem = {
               _name : ClearName(branch.fName),
               _kind : "ROOT." + branch._typename,
               _title : branch.fTitle,
               _obj : branch
         };

         if (!node._childs) node._childs = [];

         node._childs.push(subitem);

         if (branch._typename==='TBranchElement')
            subitem._title += " from " + branch.fClassName + ";" + branch.fClassVersion;

         if (nb_branches > 0) {
            subitem._more = true;
            subitem._expand = function(bnode,bobj) {
               // really create all sub-branch items
               if (!bobj) return false;

               if (!bnode._childs) bnode._childs = [];

               if (bobj.fLeaves && (bobj.fLeaves.arr.length === 1) &&
                   ((bobj.fType === JSROOT.BranchType.kClonesNode) || (bobj.fType === JSROOT.BranchType.kSTLNode))) {
                    bobj.fLeaves.arr[0].$branch = bobj;
                    bnode._childs.push({
                       _name: "@size",
                       _title: "container size",
                       _kind: "ROOT.TLeafElement",
                       _icon: "img_leaf",
                       _obj: bobj.fLeaves.arr[0],
                       _more : false
                    });
                 }

               for (let i=0; i<bobj.fBranches.arr.length; ++i)
                  CreateBranchItem(bnode, bobj.fBranches.arr[i], bobj.$tree, bobj);

               let object_class = getBranchObjectClass(bobj, bobj.$tree, true),
                   methods = object_class ? JSROOT.getMethods(object_class) : null;

               if (methods && (bobj.fBranches.arr.length>0))
                  for (let key in methods) {
                     if (typeof methods[key] !== 'function') continue;
                     let s = methods[key].toString();
                     if ((s.indexOf("return")>0) && (s.indexOf("function ()")==0))
                        bnode._childs.push({
                           _name: key+"()",
                           _title: "function " + key + " of class " + object_class,
                           _kind: "ROOT.TBranchFunc", // fictional class, only for drawing
                           _obj: { _typename: "TBranchFunc", branch: bobj, func: key },
                           _more : false
                        });

                  }

               return true;
            }
            return true;
         } else if (nb_leaves === 1) {
            subitem._icon = "img_leaf";
            subitem._more = false;
         } else if (nb_leaves > 1) {
            subitem._childs = [];
            for (let j = 0; j < nb_leaves; ++j) {
               branch.fLeaves.arr[j].$branch = branch; // keep branch pointer for drawing
               let leafitem = {
                  _name : ClearName(branch.fLeaves.arr[j].fName),
                  _kind : "ROOT." + branch.fLeaves.arr[j]._typename,
                  _obj: branch.fLeaves.arr[j]
               }
               subitem._childs.push(leafitem);
            }
         }

         return true;
      }

      node._childs = [];
      node._tree = obj;  // set reference, will be used later by TTree::Draw

      for (let i=0; i < obj.fBranches.arr.length; ++i)
         CreateBranchItem(node, obj.fBranches.arr[i], obj);

      return true;
   }

   /** @summary function called from JSROOT.draw()
    * @desc just envelope for real TTree::Draw method which do the main job
    * Can be also used for the branch and leaf object
    * @private */
   JSROOT.drawTree = function(divid, obj, opt) {

      if (!JSROOT.ObjectPainter) {
         console.error('JSROOT.drawTree called without loading of JSRoot.painter.js - how it can happen?');
         return null;
      }

      let painter = new JSROOT.ObjectPainter(divid, obj),
          tree = obj, args = opt;

      if (obj._typename == "TBranchFunc") {
         // fictional object, created only in browser
         args = { expr: "." + obj.func + "()", branch: obj.branch };
         if (opt && opt.indexOf("dump")==0)
            args.expr += ">>" + opt;
         else if (opt)
            args.expr += opt;
         tree = obj.branch.$tree;
      } else if (obj.$branch) {
         // this is drawing of the single leaf from the branch
         args = { expr: "." + obj.fName + (opt || ""), branch: obj.$branch };
         if ((args.branch.fType === JSROOT.BranchType.kClonesNode) || (args.branch.fType === JSROOT.BranchType.kSTLNode)) {
            // special case of size
            args.expr = opt;
            args.direct_branch = true;
         }

         tree = obj.$branch.$tree;
      } else if (obj.$tree) {
         // this is drawing of the branch

         // if generic object tried to be drawn without specifying any options, it will be just dump
         if (!opt && obj.fStreamerType && (obj.fStreamerType !== JSROOT.IO.kTString) &&
             (obj.fStreamerType >= JSROOT.IO.kObject) && (obj.fStreamerType <= JSROOT.IO.kAnyP)) opt = "dump";

         args = { expr: opt, branch: obj };
         tree = obj.$tree;
      } else {

         if ((args==='player') || !args) {
            return JSROOT.require("jq2d").then(() => {
               JSROOT.createTreePlayer(painter);
               painter.configureTree(tree);
               painter.Show();
               return painter;
            });
         }

         if (typeof args === 'string') args = { expr: args };
      }

      if (!tree) {
         console.error('No TTree object available for TTree::Draw');
         return null;
      }

      JSROOT.cleanup(divid);

      let create_player = 0, finalResolve;

      let process_result = function(obj, intermediate) {

         let drawid;

         if (!args.player)
            drawid = divid;
         else if (create_player === 2)
            drawid = painter.drawid;

         if (drawid)
            return JSROOT.redraw(drawid, obj); // return painter for histogram

         if (create_player === 1)
            return intermediate ? null : new Promise(resolve => { finalResolve = resolve; });

         // redirect drawing to the player
         create_player = 1;

         return JSROOT.require("jq2d").then(() => {
            JSROOT.createTreePlayer(painter);
            painter.configureTree(tree);
            painter.Show(args);
            create_player = 2;
            return JSROOT.redraw(painter.drawid, obj).then(objpainter => {
               painter.setItemName("TreePlayer"); // item name used by MDI when process resize
               if (finalResolve) finalResolve(objpainter);
               return objpainter; // return painter for histogram
            });
         });
      }

      args.progress = obj => process_result(obj, true);

      // use in result handling same function as for progress handling
      return tree.Draw(args).then(obj => process_result(obj));
   }

   JSROOT.TSelector = TSelector;
   JSROOT.TDrawVariable = TDrawVariable;
   JSROOT.TDrawSelector = TDrawSelector;
   jsrio.TTreeMethods = TTreeMethods;

   return JSROOT;

});
