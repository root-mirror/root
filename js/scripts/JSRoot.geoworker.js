JSROOT = {}; // just place holder for JSROOT.GEO functions

JSROOT.BIT = function(n) { return 1 << (n); }

let jscsg;

// emulate behavour, only required for csg and geobase
JSROOT.define = function(args, factory) {
   if (args.length == 1)
      jscsg = factory(THREE);
   else
      factory(THREE, jscsg);
}

importScripts("three.min.js", "JSRoot.csg.js", "JSRoot.geobase.js");

if (console) console.log('geoworker started ' + THREE.REVISION);

let clones = null;

onmessage = function(e) {

   if (typeof e.data == 'string') {
      console.log('Worker get message ' + e.data);
      return;
   }

   if (typeof e.data != 'object') return;

   e.data.tm1 = new Date().getTime();

   if (e.data.init) {
      // console.log('start worker ' +  (e.data.tm1 -  e.data.tm0));

      let nodes = e.data.clones;
      if (nodes) {
         // console.log('get clones ' + nodes.length);
         clones = new JSROOT.GEO.ClonedNodes(null, nodes);
         clones.setVisLevel(e.data.vislevel);
         clones.setMaxVisNodes(e.data.maxvisnodes);
         delete e.data.clones;
         clones.sortmap = e.data.sortmap;
      }

      // used in composite shape
      JSROOT.browser = e.data.browser;

      e.data.tm2 = new Date().getTime();

      return postMessage(e.data);
   }

   if (e.data.shapes) {
      // this is task to create geometries in the worker

      let shapes = e.data.shapes, transferables = [];

      // build all shapes up to specified limit, also limit execution time
      for (let n=0;n<100;++n) {
         let res = clones.buildShapes(shapes, e.data.limit, 1000);
         if (res.done) break;
         postMessage({ progress: "Worker creating: " + res.shapes + " / " + shapes.length + " shapes,  "  + res.faces + " faces" });
      }

      for (let n=0;n<shapes.length;++n) {
         let item = shapes[n];

         if (item.geom) {
            let bufgeom;
            if (item.geom instanceof THREE.BufferGeometry) {
               bufgeom = item.geom;
            } else {
               let bufgeom = new THREE.BufferGeometry();
               bufgeom.fromGeometry(item.geom);
            }

            item.buf_pos = bufgeom.attributes.position.array;
            item.buf_norm = bufgeom.attributes.normal.array;

            // use nice feature of HTML workers with transferable
            // we allow to take ownership of buffer from local array
            // therefore buffer content not need to be copied
            transferables.push(item.buf_pos.buffer, item.buf_norm.buffer);

            delete item.geom;
         }

         delete item.shape; // no need to send back shape
      }

      e.data.tm2 = new Date().getTime();

      return postMessage(e.data, transferables);
   }

   if (e.data.collect !== undefined) {
      // this is task to collect visible nodes using camera position

      // first mark all visible flags
      clones.setVisibleFlags(e.data.flags);
      delete e.data.flags;

      clones.produceIdShifts();

      let matrix = null;
      if (e.data.matrix)
         matrix = new THREE.Matrix4().fromArray(e.data.matrix);
      delete e.data.matrix;

      let res = clones.collectVisibles(e.data.collect, JSROOT.GEO.createFrustum(matrix));

      e.data.new_nodes = res.lst;
      e.data.complete = res.complete; // inform if all nodes are selected

      e.data.tm2 = new Date().getTime();

      // console.log('Collect visibles in worker ' + e.data.new_nodes.length + ' takes ' + (e.data.tm2-e.data.tm1));

      return postMessage(e.data);
   }

}
