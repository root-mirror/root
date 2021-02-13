/// @file JSRoot.geobase.js
/// Basic functions for work with TGeo classes

JSROOT.define(['three', 'csg'], (THREE, ThreeBSP) => {

   "use strict";

   /** @summary Collection of TGeo-related methods, loaded with ```JSROOT.require('geo').then(geo => ...)```
     * @namespace
     * @alias JSROOT.GEO
     */
   let geo = {
      GradPerSegm: 6,      // grad per segment in cylinder/spherical symmetry shapes
      CompressComp: true,  // use faces compression in composite shapes
      CompLimit: 20        // maximal number of components in composite shape
   };

   const kindGeo = 0,    // TGeoNode / TGeoShape
         kindEve = 1,    // TEveShape / TEveGeoShapeExtract
         kindShape = 2;  // special kind for single shape handling

   /** @summary TGeo-related bits
     * @private */
   geo.BITS = {
         kVisOverride     : JSROOT.BIT(0),           // volume's vis. attributes are overwritten
         kVisNone         : JSROOT.BIT(1),           // the volume/node is invisible, as well as daughters
         kVisThis         : JSROOT.BIT(2),           // this volume/node is visible
         kVisDaughters    : JSROOT.BIT(3),           // all leaves are visible
         kVisOneLevel     : JSROOT.BIT(4),           // first level daughters are visible (not used)
         kVisStreamed     : JSROOT.BIT(5),           // true if attributes have been streamed
         kVisTouched      : JSROOT.BIT(6),           // true if attributes are changed after closing geom
         kVisOnScreen     : JSROOT.BIT(7),           // true if volume is visible on screen
         kVisContainers   : JSROOT.BIT(12),          // all containers visible
         kVisOnly         : JSROOT.BIT(13),          // just this visible
         kVisBranch       : JSROOT.BIT(14),          // only a given branch visible
         kVisRaytrace     : JSROOT.BIT(15)           // raytracing flag
      };

   /** @summary Test fGeoAtt bits
     * @private */
   geo.TestBit = function(volume, f) {
      let att = volume.fGeoAtt;
      return att === undefined ? false : ((att & f) !== 0);
   }

   /** @summary Set fGeoAtt bit
     * @private */
   geo.SetBit = function(volume, f, value) {
      if (volume.fGeoAtt === undefined) return;
      volume.fGeoAtt = value ? (volume.fGeoAtt | f) : (volume.fGeoAtt & ~f);
   }

   /** @summary Toggle fGeoAttBit
     * @private */
   geo.ToggleBit = function(volume, f) {
      if (volume.fGeoAtt !== undefined)
         volume.fGeoAtt = volume.fGeoAtt ^ (f & 0xffffff);
   }

   /** @summary Implementation of TGeoVolume::InvisibleAll
     * @private */
   geo.setInvisibleAll = function(volume, flag) {
      if (flag===undefined) flag = true;

      geo.SetBit(volume, geo.BITS.kVisThis, !flag);
      // geo.SetBit(this, geo.BITS.kVisDaughters, !flag);

      if (volume.fNodes)
         for (let n = 0; n < volume.fNodes.arr.length; ++n) {
            let sub = volume.fNodes.arr[n].fVolume;
            geo.SetBit(sub, geo.BITS.kVisThis, !flag);
            // geo.SetBit(sub, geo.BITS.kVisDaughters, !flag);
         }
   }

   /** @summary method used to avoid duplication of warnings
    * @private */
   geo.warn = function(msg) {
      if (geo._warn_msgs === undefined) geo._warn_msgs = {};
      if (geo._warn_msgs[msg] !== undefined) return;
      geo._warn_msgs[msg] = true;
      console.warn(msg);
   }

   /** @summary Analyze TGeo node kind
    *  @desc  0 - TGeoNode
    *         1 - TEveGeoNode
    *        -1 - unsupported
    * @returns detected node kind
    * @private */
   geo.getNodeKind = function(obj) {
      if ((obj === undefined) || (obj === null) || (typeof obj !== 'object')) return -1;
      return ('fShape' in obj) && ('fTrans' in obj) ? 1 : 0;
   }

   /** @summary Returns number of shapes
     * @desc Used to count total shapes number in composites
     * @private */
   geo.countNumShapes = function(shape) {
      if (!shape) return 0;
      // if (shape._typename=="TGeoHalfSpace") geo.HalfSpace = true;
      if (shape._typename !== 'TGeoCompositeShape') return 1;
      return geo.countNumShapes(shape.fNode.fLeft) + geo.countNumShapes(shape.fNode.fRight);
   }

   /** @summary Create normal to plane, defined with three points
     * @private */
   geo.produceNormal = function(x1,y1,z1,x2,y2,z2,x3,y3,z3) {

      let pA = new THREE.Vector3(x1,y1,z1),
          pB = new THREE.Vector3(x2,y2,z2),
          pC = new THREE.Vector3(x3,y3,z3),
          cb = new THREE.Vector3(),
          ab = new THREE.Vector3();

      cb.subVectors( pC, pB );
      ab.subVectors( pA, pB );
      cb.cross(ab );

      return cb;
   }

   // ==========================================================================

   /**
     * @summary Helper class for geometry creation
     *
     * @class
     * @memberof JSROOT.GEO
     * @private
     */

   function GeometryCreator(numfaces) {
      this.nfaces = numfaces;
      this.indx = 0;
      this.pos = new Float32Array(numfaces*9);
      this.norm = new Float32Array(numfaces*9);

      return this;
   }

   /** @summary Add face with 3 vertices
     * @private */
   GeometryCreator.prototype.addFace3 = function(x1,y1,z1,
                                                 x2,y2,z2,
                                                 x3,y3,z3) {
      let indx = this.indx, pos = this.pos;
      pos[indx] = x1;
      pos[indx+1] = y1;
      pos[indx+2] = z1;
      pos[indx+3] = x2;
      pos[indx+4] = y2;
      pos[indx+5] = z2;
      pos[indx+6] = x3;
      pos[indx+7] = y3;
      pos[indx+8] = z3;
      this.last4 = false;
      this.indx = indx + 9;
   }

   /** @summary Start polygon
     * @private */
   GeometryCreator.prototype.startPolygon = function() {}

   /** @summary Stop polygon
     * @private */
   GeometryCreator.prototype.stopPolygon = function() {}

   /** @summary Add face with 4 vertices
     * @desc From four vertices one normally creates two faces (1,2,3) and (1,3,4)
     * if (reduce==1), first face is reduced
     * if (reduce==2), second face is reduced
     * @private */
   GeometryCreator.prototype.addFace4 = function(x1,y1,z1,
                                                 x2,y2,z2,
                                                 x3,y3,z3,
                                                 x4,y4,z4,
                                                 reduce) {
      let indx = this.indx, pos = this.pos;

      if (reduce !== 1) {
         pos[indx] = x1;
         pos[indx+1] = y1;
         pos[indx+2] = z1;
         pos[indx+3] = x2;
         pos[indx+4] = y2;
         pos[indx+5] = z2;
         pos[indx+6] = x3;
         pos[indx+7] = y3;
         pos[indx+8] = z3;
         indx+=9;
      }

      if (reduce !== 2) {
         pos[indx] = x1;
         pos[indx+1] = y1;
         pos[indx+2] = z1;
         pos[indx+3] = x3;
         pos[indx+4] = y3;
         pos[indx+5] = z3;
         pos[indx+6] = x4;
         pos[indx+7] = y4;
         pos[indx+8] = z4;
         indx+=9;
      }

      this.last4 = (indx !== this.indx + 9);
      this.indx = indx;
   }

   /** @summary Specify normal for face with 4 vertices
     * @desc same as addFace4, assign normals for each individual vertex
     * reduce has same meaning and should be the same
     * @private */
   GeometryCreator.prototype.setNormal4 = function(nx1,ny1,nz1,
                                                   nx2,ny2,nz2,
                                                   nx3,ny3,nz3,
                                                   nx4,ny4,nz4,
                                                   reduce) {
      if (this.last4 && reduce)
         return console.error('missmatch between addFace4 and setNormal4 calls');

      let indx = this.indx - (this.last4 ? 18 : 9), norm = this.norm;

      if (reduce!==1) {
         norm[indx] = nx1;
         norm[indx+1] = ny1;
         norm[indx+2] = nz1;
         norm[indx+3] = nx2;
         norm[indx+4] = ny2;
         norm[indx+5] = nz2;
         norm[indx+6] = nx3;
         norm[indx+7] = ny3;
         norm[indx+8] = nz3;
         indx+=9;
      }

      if (reduce!==2) {
         norm[indx] = nx1;
         norm[indx+1] = ny1;
         norm[indx+2] = nz1;
         norm[indx+3] = nx3;
         norm[indx+4] = ny3;
         norm[indx+5] = nz3;
         norm[indx+6] = nx4;
         norm[indx+7] = ny4;
         norm[indx+8] = nz4;
      }
   }

   /** @summary Recalculate Z with provided func
     * @private */
   GeometryCreator.prototype.recalcZ = function(func) {
      let pos = this.pos,
          last = this.indx,
          indx = last - (this.last4 ? 18 : 9);

      while (indx < last) {
         pos[indx+2] = func(pos[indx], pos[indx+1], pos[indx+2]);
         indx+=3;
      }
   }

   /** @summary Caclualte normal
     * @private */
   GeometryCreator.prototype.calcNormal = function() {
      if (!this.cb) {
         this.pA = new THREE.Vector3();
         this.pB = new THREE.Vector3();
         this.pC = new THREE.Vector3();
         this.cb = new THREE.Vector3();
         this.ab = new THREE.Vector3();
      }

      this.pA.fromArray( this.pos, this.indx - 9 );
      this.pB.fromArray( this.pos, this.indx - 6 );
      this.pC.fromArray( this.pos, this.indx - 3 );

      this.cb.subVectors( this.pC, this.pB );
      this.ab.subVectors( this.pA, this.pB );
      this.cb.cross( this.ab );

      this.setNormal(this.cb.x, this.cb.y, this.cb.z);
   }

   /** @summary Set normal
     * @private */
   GeometryCreator.prototype.setNormal = function(nx,ny,nz) {
      let indx = this.indx - 9, norm = this.norm;

      norm[indx]   = norm[indx+3] = norm[indx+6] = nx;
      norm[indx+1] = norm[indx+4] = norm[indx+7] = ny;
      norm[indx+2] = norm[indx+5] = norm[indx+8] = nz;

      if (this.last4) {
         indx -= 9;
         norm[indx]   = norm[indx+3] = norm[indx+6] = nx;
         norm[indx+1] = norm[indx+4] = norm[indx+7] = ny;
         norm[indx+2] = norm[indx+5] = norm[indx+8] = nz;
      }
   }

   /** @summary Set normal
     * @desc special shortcut, when same normals can be applied for 1-2 point and 3-4 point
     * @private */
   GeometryCreator.prototype.setNormal_12_34 = function(nx12,ny12,nz12,nx34,ny34,nz34,reduce) {
      if (reduce===undefined) reduce = 0;

      let indx = this.indx - ((reduce>0) ? 9 : 18), norm = this.norm;

      if (reduce!==1) {
         norm[indx]   = nx12;
         norm[indx+1] = ny12;
         norm[indx+2] = nz12;
         norm[indx+3] = nx12;
         norm[indx+4] = ny12;
         norm[indx+5] = nz12;
         norm[indx+6] = nx34;
         norm[indx+7] = ny34;
         norm[indx+8] = nz34;
         indx+=9;
      }

      if (reduce!==2) {
         norm[indx]   = nx12;
         norm[indx+1] = ny12;
         norm[indx+2] = nz12;
         norm[indx+3] = nx34;
         norm[indx+4] = ny34;
         norm[indx+5] = nz34;
         norm[indx+6] = nx34;
         norm[indx+7] = ny34;
         norm[indx+8] = nz34;
         indx+=9;
      }
   }

   /** @summary Create geometry
     * @private */
   GeometryCreator.prototype.create = function() {
      if (this.nfaces !== this.indx/9)
         console.error('Mismatch with created ' + this.nfaces + ' and filled ' + this.indx/9 + ' number of faces');

      let geometry = new THREE.BufferGeometry();
      geometry.setAttribute( 'position', new THREE.BufferAttribute( this.pos, 3 ) );
      geometry.setAttribute( 'normal', new THREE.BufferAttribute( this.norm, 3 ) );
      return geometry;
   }

   // ================================================================================

   /** @summary Helper class for ThreeBSP geometry creation
     *
     * @class
     * @memberof JSROOT.GEO
     * @private
     */

   function PolygonsCreator() {
      this.polygons = [];
   }

   /** @summary Start polygon
     * @private */
   PolygonsCreator.prototype.startPolygon = function(normal) {
      this.multi = 1;
      this.mnormal = normal;
   }

   /** @summary Stop polygon
     * @private */
   PolygonsCreator.prototype.stopPolygon = function() {
      if (!this.multi) return;
      this.multi = 0;
      console.error('Polygon should be already closed at this moment');
   }

   /** @summary Add face with 3 vertices
     * @private */
   PolygonsCreator.prototype.addFace3 = function(x1,y1,z1,
                                                 x2,y2,z2,
                                                 x3,y3,z3) {
      this.addFace4(x1,y1,z1,x2,y2,z2,x3,y3,z3,x3,y3,z3,2);
   }

   /** @summary Add face with 4 vertices
     * @desc From four vertices one normally creates two faces (1,2,3) and (1,3,4)
     * if (reduce==1), first face is reduced
     * if (reduce==2), second face is reduced
     * @private */
   PolygonsCreator.prototype.addFace4 = function(x1,y1,z1,
                                                 x2,y2,z2,
                                                 x3,y3,z3,
                                                 x4,y4,z4,
                                                 reduce) {
      if (reduce === undefined) reduce = 0;

      this.v1 = new ThreeBSP.Vertex( x1, y1, z1, 0, 0, 0 );
      this.v2 = (reduce===1) ? null : new ThreeBSP.Vertex( x2, y2, z2, 0, 0, 0 );
      this.v3 = new ThreeBSP.Vertex( x3, y3, z3, 0, 0, 0 );
      this.v4 = (reduce===2) ? null : new ThreeBSP.Vertex( x4, y4, z4, 0, 0, 0 );

      this.reduce = reduce;

      if (this.multi) {

         if (reduce!==2) console.error('polygon not supported for not-reduced faces');

         let polygon;

         if (this.multi++ === 1) {
            polygon = new ThreeBSP.Polygon;

            polygon.vertices.push(this.mnormal ? this.v2 : this.v3);
            this.polygons.push(polygon);
         } else {
            polygon = this.polygons[this.polygons.length-1];
            // check that last vertice equals to v2
            let last = this.mnormal ? polygon.vertices[polygon.vertices.length-1] : polygon.vertices[0],
                comp = this.mnormal ? this.v2 : this.v3;

            if (comp.diff(last) > 1e-12)
               console.error('vertex missmatch when building polygon');
         }

         let first = this.mnormal ? polygon.vertices[0] : polygon.vertices[polygon.vertices.length-1],
             next = this.mnormal ? this.v3 : this.v2;

         if (next.diff(first) < 1e-12) {
            //console.log('polygon closed!!!', polygon.vertices.length);
            this.multi = 0;
         } else
         if (this.mnormal) {
            polygon.vertices.push(this.v3);
         } else {
            polygon.vertices.unshift(this.v2);
         }

         return;

      }

      let polygon = new ThreeBSP.Polygon;

      switch (reduce) {
         case 0: polygon.vertices.push(this.v1, this.v2, this.v3, this.v4); break;
         case 1: polygon.vertices.push(this.v1, this.v3, this.v4); break;
         case 2: polygon.vertices.push(this.v1, this.v2, this.v3); break;
      }

      this.polygons.push(polygon);
   }

   /** @summary Specify normal for face with 4 vertices
     * @desc same as addFace4, assign normals for each individual vertex
     * reduce has same meaning and should be the same
     * @private */
   PolygonsCreator.prototype.setNormal4 = function(nx1,ny1,nz1,
                                                   nx2,ny2,nz2,
                                                   nx3,ny3,nz3,
                                                   nx4,ny4,nz4) {
      this.v1.setnormal(nx1,ny1,nz1);
      if (this.v2) this.v2.setnormal(nx2,ny2,nz2);
      this.v3.setnormal(nx3,ny3,nz3);
      if (this.v4) this.v4.setnormal(nx4,ny4,nz4);
   }

   /** @summary Set normal
     * @desc special shortcut, when same normals can be applied for 1-2 point and 3-4 point
     * @private */
   PolygonsCreator.prototype.setNormal_12_34 = function(nx12,ny12,nz12,nx34,ny34,nz34) {
      this.v1.setnormal(nx12,ny12,nz12);
      if (this.v2) this.v2.setnormal(nx12,ny12,nz12);
      this.v3.setnormal(nx34,ny34,nz34);
      if (this.v4) this.v4.setnormal(nx34,ny34,nz34);
   }

   /** @summary Calculate normal
     * @private */
   PolygonsCreator.prototype.calcNormal = function() {

      if (!this.cb) {
         this.pA = new THREE.Vector3();
         this.pB = new THREE.Vector3();
         this.pC = new THREE.Vector3();
         this.cb = new THREE.Vector3();
         this.ab = new THREE.Vector3();
      }

      this.pA.set( this.v1.x, this.v1.y, this.v1.z);

      if (this.reduce!==1) {
         this.pB.set( this.v2.x, this.v2.y, this.v2.z);
         this.pC.set( this.v3.x, this.v3.y, this.v3.z);
      } else {
         this.pB.set( this.v3.x, this.v3.y, this.v3.z);
         this.pC.set( this.v4.x, this.v4.y, this.v4.z);
      }

      this.cb.subVectors( this.pC, this.pB );
      this.ab.subVectors( this.pA, this.pB );
      this.cb.cross( this.ab );

      this.setNormal(this.cb.x, this.cb.y, this.cb.z);
   }

   /** @summary Set normal
     * @private */
   PolygonsCreator.prototype.setNormal = function(nx,ny,nz) {
      this.v1.setnormal(nx,ny,nz);
      if (this.v2) this.v2.setnormal(nx,ny,nz);
      this.v3.setnormal(nx,ny,nz);
      if (this.v4) this.v4.setnormal(nx,ny,nz);
   }

   /** @summary Recalculate Z with provided func
     * @private */
   PolygonsCreator.prototype.recalcZ = function(func) {
      this.v1.z = func(this.v1.x, this.v1.y, this.v1.z);
      if (this.v2) this.v2.z = func(this.v2.x, this.v2.y, this.v2.z);
      this.v3.z = func(this.v3.x, this.v3.y, this.v3.z);
      if (this.v4) this.v4.z = func(this.v4.x, this.v4.y, this.v4.z);
   }

   /** @summary Create geometry
     * @private */
   PolygonsCreator.prototype.create = function() {
      return { polygons: this.polygons };
   }

   // ================= all functions to create geometry ===================================

   /** @summary Creates cube geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createCubeBuffer(shape, faces_limit) {

      if (faces_limit < 0) return 12;

      let dx = shape.fDX, dy = shape.fDY, dz = shape.fDZ;

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(12);

      creator.addFace4(dx,dy,dz, dx,-dy,dz, dx,-dy,-dz, dx,dy,-dz); creator.setNormal(1,0,0);

      creator.addFace4(-dx,dy,-dz, -dx,-dy,-dz, -dx,-dy,dz, -dx,dy,dz); creator.setNormal(-1,0,0);

      creator.addFace4(-dx,dy,-dz, -dx,dy,dz, dx,dy,dz, dx,dy,-dz); creator.setNormal(0,1,0);

      creator.addFace4(-dx,-dy,dz, -dx,-dy,-dz, dx,-dy,-dz, dx,-dy,dz); creator.setNormal(0,-1,0);

      creator.addFace4(-dx,dy,dz, -dx,-dy,dz, dx,-dy,dz, dx,dy,dz); creator.setNormal(0,0,1);

      creator.addFace4(dx,dy,-dz, dx,-dy,-dz, -dx,-dy,-dz, -dx,dy,-dz); creator.setNormal(0,0,-1);

      return creator.create();
   }

   /** @summary Creates 8 edges geometry
     * @memberof JSROOT.GEO
     * @private */
   function create8edgesBuffer( v, faces_limit ) {

      let indicies = [ 4,7,6,5,  0,3,7,4,  4,5,1,0,  6,2,1,5,  7,3,2,6,  1,2,3,0 ];

      let creator = (faces_limit > 0) ? new PolygonsCreator : new GeometryCreator(12);

      for (let n=0;n<indicies.length;n+=4) {
         let i1 = indicies[n]*3,
             i2 = indicies[n+1]*3,
             i3 = indicies[n+2]*3,
             i4 = indicies[n+3]*3;
         creator.addFace4(v[i1], v[i1+1], v[i1+2], v[i2], v[i2+1], v[i2+2],
                          v[i3], v[i3+1], v[i3+2], v[i4], v[i4+1], v[i4+2]);
         if (n===0)
            creator.setNormal(0,0,1);
         else if (n===20)
            creator.setNormal(0,0,-1);
         else
            creator.calcNormal();
      }

      return creator.create();
   }

   /** @summary Creates PARA geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createParaBuffer( shape, faces_limit ) {

      if (faces_limit < 0) return 12;

      let txy = shape.fTxy, txz = shape.fTxz, tyz = shape.fTyz;

      let v = [
          -shape.fZ*txz-txy*shape.fY-shape.fX, -shape.fY-shape.fZ*tyz,  -shape.fZ,
          -shape.fZ*txz+txy*shape.fY-shape.fX,  shape.fY-shape.fZ*tyz,  -shape.fZ,
          -shape.fZ*txz+txy*shape.fY+shape.fX,  shape.fY-shape.fZ*tyz,  -shape.fZ,
          -shape.fZ*txz-txy*shape.fY+shape.fX, -shape.fY-shape.fZ*tyz,  -shape.fZ,
           shape.fZ*txz-txy*shape.fY-shape.fX, -shape.fY+shape.fZ*tyz,   shape.fZ,
           shape.fZ*txz+txy*shape.fY-shape.fX,  shape.fY+shape.fZ*tyz,   shape.fZ,
           shape.fZ*txz+txy*shape.fY+shape.fX,  shape.fY+shape.fZ*tyz,   shape.fZ,
           shape.fZ*txz-txy*shape.fY+shape.fX, -shape.fY+shape.fZ*tyz,   shape.fZ ];

      return create8edgesBuffer(v, faces_limit );
   }

   /** @summary Creates Ttrapezoid geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createTrapezoidBuffer( shape, faces_limit ) {

      if (faces_limit < 0) return 12;

      let y1, y2;
      if (shape._typename == "TGeoTrd1") {
         y1 = y2 = shape.fDY;
      } else {
         y1 = shape.fDy1; y2 = shape.fDy2;
      }

      let v = [
            -shape.fDx1,  y1, -shape.fDZ,
             shape.fDx1,  y1, -shape.fDZ,
             shape.fDx1, -y1, -shape.fDZ,
            -shape.fDx1, -y1, -shape.fDZ,
            -shape.fDx2,  y2,  shape.fDZ,
             shape.fDx2,  y2,  shape.fDZ,
             shape.fDx2, -y2,  shape.fDZ,
            -shape.fDx2, -y2,  shape.fDZ
         ];

      return create8edgesBuffer(v, faces_limit );
   }


   /** @summary Creates arb8 geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createArb8Buffer( shape, faces_limit ) {

      if (faces_limit < 0) return 12;

      let vertices = [
            shape.fXY[0][0], shape.fXY[0][1], -shape.fDZ,
            shape.fXY[1][0], shape.fXY[1][1], -shape.fDZ,
            shape.fXY[2][0], shape.fXY[2][1], -shape.fDZ,
            shape.fXY[3][0], shape.fXY[3][1], -shape.fDZ,
            shape.fXY[4][0], shape.fXY[4][1],  shape.fDZ,
            shape.fXY[5][0], shape.fXY[5][1],  shape.fDZ,
            shape.fXY[6][0], shape.fXY[6][1],  shape.fDZ,
            shape.fXY[7][0], shape.fXY[7][1],  shape.fDZ
         ];
      const indicies = [
            4,7,6,   6,5,4,   3,7,4,   4,0,3,
            5,1,0,   0,4,5,   6,2,1,   1,5,6,
            7,3,2,   2,6,7,   1,2,3,   3,0,1 ];

      // detect same vertices on both Z-layers
      for (let side=0;side<vertices.length;side += vertices.length/2)
         for (let n1 = side; n1 < side + vertices.length/2 - 3 ; n1+=3)
            for (let n2 = n1+3; n2 < side + vertices.length/2 ; n2+=3)
               if ((vertices[n1] === vertices[n2]) &&
                   (vertices[n1+1] === vertices[n2+1]) &&
                   (vertices[n1+2] === vertices[n2+2])) {
                      for (let k=0;k<indicies.length;++k)
                        if (indicies[k] === n2/3) indicies[k] = n1/3;
                  }


      let map = [], // list of existing faces (with all rotations)
          numfaces = 0;

      for (let k=0;k<indicies.length;k+=3) {
         let id1 = indicies[k]*100   + indicies[k+1]*10 + indicies[k+2],
             id2 = indicies[k+1]*100 + indicies[k+2]*10 + indicies[k],
             id3 = indicies[k+2]*100 + indicies[k]*10   + indicies[k+1];

         if ((indicies[k] == indicies[k+1]) || (indicies[k] == indicies[k+2]) || (indicies[k+1] == indicies[k+2]) ||
             (map.indexOf(id1)>=0) || (map.indexOf(id2)>=0) || (map.indexOf(id3)>=0)) {
            indicies[k] = indicies[k+1] = indicies[k+2] = -1;
         } else {
            map.push(id1,id2,id3);
            numfaces++;
         }
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      // let creator = new GeometryCreator(numfaces);

      for (let n=0; n < indicies.length; n+=6) {
         let i1 = indicies[n]   * 3,
             i2 = indicies[n+1] * 3,
             i3 = indicies[n+2] * 3,
             i4 = indicies[n+3] * 3,
             i5 = indicies[n+4] * 3,
             i6 = indicies[n+5] * 3,
             norm = null;

         if ((i1>=0) && (i4>=0) && faces_limit) {
            // try to identify two faces with same normal - very useful if one can create face4
            if (n===0) norm = new THREE.Vector3(0,0,1); else
            if (n===30) norm = new THREE.Vector3(0,0,-1); else {
               let norm1 = geo.produceNormal(vertices[i1], vertices[i1+1], vertices[i1+2],
                                                vertices[i2], vertices[i2+1], vertices[i2+2],
                                                vertices[i3], vertices[i3+1], vertices[i3+2]);

               norm1.normalize();

               let norm2 = geo.produceNormal(vertices[i4], vertices[i4+1], vertices[i4+2],
                                                vertices[i5], vertices[i5+1], vertices[i5+2],
                                                vertices[i6], vertices[i6+1], vertices[i6+2]);

               norm2.normalize();

               if (norm1.distanceToSquared(norm2) < 1e-12) norm = norm1;
            }
         }

         if (norm !== null) {
            creator.addFace4(vertices[i1], vertices[i1+1], vertices[i1+2],
                             vertices[i2], vertices[i2+1], vertices[i2+2],
                             vertices[i3], vertices[i3+1], vertices[i3+2],
                             vertices[i5], vertices[i5+1], vertices[i5+2]);
            creator.setNormal(norm.x, norm.y, norm.z);
         }  else {
            if (i1>=0) {
               creator.addFace3(vertices[i1], vertices[i1+1], vertices[i1+2],
                                vertices[i2], vertices[i2+1], vertices[i2+2],
                                vertices[i3], vertices[i3+1], vertices[i3+2]);
               creator.calcNormal();
            }
            if (i4>=0) {
               creator.addFace3(vertices[i4], vertices[i4+1], vertices[i4+2],
                                vertices[i5], vertices[i5+1], vertices[i5+2],
                                vertices[i6], vertices[i6+1], vertices[i6+2]);
               creator.calcNormal();
            }
         }
      }

      return creator.create();
   }

   /** @summary Creates sphere geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createSphereBuffer( shape, faces_limit ) {
      let radius = [shape.fRmax, shape.fRmin],
          phiStart = shape.fPhi1,
          phiLength = shape.fPhi2 - shape.fPhi1,
          thetaStart = shape.fTheta1,
          thetaLength = shape.fTheta2 - shape.fTheta1,
          widthSegments = shape.fNseg,
          heightSegments = shape.fNz,
          noInside = (radius[1] <= 0);

      if (faces_limit > 0) {
         let fact = (noInside ? 2 : 4) * widthSegments * heightSegments / faces_limit;

         if (fact > 1.) {
            widthSegments = Math.max(4, Math.floor(widthSegments/Math.sqrt(fact)));
            heightSegments = Math.max(4, Math.floor(heightSegments/Math.sqrt(fact)));
         }
      }

      let numoutside = widthSegments * heightSegments * 2,
          numtop = widthSegments * 2,
          numbottom = widthSegments * 2,
          numcut = phiLength === 360 ? 0 : heightSegments * (noInside ? 2 : 4),
          epsilon = 1e-10;

      if (noInside) numbottom = numtop = widthSegments;

      if (faces_limit < 0) return numoutside * (noInside ? 1 : 2) + numtop + numbottom + numcut;

      let _sinp = new Float32Array(widthSegments+1),
          _cosp = new Float32Array(widthSegments+1),
          _sint = new Float32Array(heightSegments+1),
          _cost = new Float32Array(heightSegments+1);

      for (let n=0;n<=heightSegments;++n) {
         let theta = (thetaStart + thetaLength/heightSegments*n)*Math.PI/180;
         _sint[n] = Math.sin(theta);
         _cost[n] = Math.cos(theta);
      }

      for (let n=0;n<=widthSegments;++n) {
         let phi = (phiStart + phiLength/widthSegments*n)*Math.PI/180;
         _sinp[n] = Math.sin(phi);
         _cosp[n] = Math.cos(phi);
      }

      if (Math.abs(_sint[0]) <= epsilon) { numoutside -= widthSegments; numtop = 0; }
      if (Math.abs(_sint[heightSegments]) <= epsilon) { numoutside -= widthSegments; numbottom = 0; }

      let numfaces = numoutside * (noInside ? 1 : 2) + numtop + numbottom + numcut;

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      // let creator = new GeometryCreator(numfaces);

      for (let side=0;side<2;++side) {
         if ((side===1) && noInside) break;

         let r = radius[side],
             s = (side===0) ? 1 : -1,
             d1 = 1 - side, d2 = 1 - d1;

         // use direct algorithm for the sphere - here normals and position can be calculated directly
         for (let k=0;k<heightSegments;++k) {

            let k1 = k + d1, k2 = k + d2;

            let skip = 0;
            if (Math.abs(_sint[k1]) <= epsilon) skip = 1; else
            if (Math.abs(_sint[k2]) <= epsilon) skip = 2;

            for (let n=0;n<widthSegments;++n) {
               creator.addFace4(
                     r*_sint[k1]*_cosp[n],   r*_sint[k1] *_sinp[n],   r*_cost[k1],
                     r*_sint[k1]*_cosp[n+1], r*_sint[k1] *_sinp[n+1], r*_cost[k1],
                     r*_sint[k2]*_cosp[n+1], r*_sint[k2] *_sinp[n+1], r*_cost[k2],
                     r*_sint[k2]*_cosp[n],   r*_sint[k2] *_sinp[n],   r*_cost[k2],
                     skip);
               creator.setNormal4(
                     s*_sint[k1]*_cosp[n],   s*_sint[k1] *_sinp[n],   s*_cost[k1],
                     s*_sint[k1]*_cosp[n+1], s*_sint[k1] *_sinp[n+1], s*_cost[k1],
                     s*_sint[k2]*_cosp[n+1], s*_sint[k2] *_sinp[n+1], s*_cost[k2],
                     s*_sint[k2]*_cosp[n],   s*_sint[k2] *_sinp[n],   s*_cost[k2],
                     skip);
            }
         }
      }

      // top/bottom
      for (let side=0; side<=heightSegments; side+=heightSegments)
         if (Math.abs(_sint[side]) >= epsilon) {
            let ss = _sint[side], cc = _cost[side],
                d1 = (side===0) ? 0 : 1, d2 = 1 - d1;
            for (let n=0;n<widthSegments;++n) {
               creator.addFace4(
                     radius[1] * ss * _cosp[n+d1], radius[1] * ss * _sinp[n+d1], radius[1] * cc,
                     radius[0] * ss * _cosp[n+d1], radius[0] * ss * _sinp[n+d1], radius[0] * cc,
                     radius[0] * ss * _cosp[n+d2], radius[0] * ss * _sinp[n+d2], radius[0] * cc,
                     radius[1] * ss * _cosp[n+d2], radius[1] * ss * _sinp[n+d2], radius[1] * cc,
                     noInside ? 2 : 0);
               creator.calcNormal();
            }
         }

      // cut left/right sides
      if (phiLength < 360) {
         for (let side=0;side<=widthSegments;side+=widthSegments) {
            let ss = _sinp[side], cc = _cosp[side],
                d1 = (side === 0) ? 1 : 0, d2 = 1 - d1;

            for (let k=0;k<heightSegments;++k) {
               creator.addFace4(
                     radius[1] * _sint[k+d1] * cc, radius[1] * _sint[k+d1] * ss, radius[1] * _cost[k+d1],
                     radius[0] * _sint[k+d1] * cc, radius[0] * _sint[k+d1] * ss, radius[0] * _cost[k+d1],
                     radius[0] * _sint[k+d2] * cc, radius[0] * _sint[k+d2] * ss, radius[0] * _cost[k+d2],
                     radius[1] * _sint[k+d2] * cc, radius[1] * _sint[k+d2] * ss, radius[1] * _cost[k+d2],
                     noInside ? 2 : 0);
               creator.calcNormal();
            }
         }
      }

      return creator.create();
   }

   /** @summary Creates tube geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createTubeBuffer( shape, faces_limit) {
      let outerR, innerR; // inner/outer tube radius
      if ((shape._typename == "TGeoCone") || (shape._typename == "TGeoConeSeg")) {
         outerR = [ shape.fRmax2, shape.fRmax1 ];
         innerR = [ shape.fRmin2, shape.fRmin1 ];
      } else {
         outerR = [ shape.fRmax, shape.fRmax ];
         innerR = [ shape.fRmin, shape.fRmin ];
      }

      let hasrmin = (innerR[0] > 0) || (innerR[1] > 0),
          thetaStart = 0, thetaLength = 360;

      if ((shape._typename == "TGeoConeSeg") || (shape._typename == "TGeoTubeSeg") || (shape._typename == "TGeoCtub")) {
         thetaStart = shape.fPhi1;
         thetaLength = shape.fPhi2 - shape.fPhi1;
      }

      let radiusSegments = Math.max(4, Math.round(thetaLength/geo.GradPerSegm));

      // external surface
      let numfaces = radiusSegments * (((outerR[0] <= 0) || (outerR[1] <= 0)) ? 1 : 2);

      // internal surface
      if (hasrmin)
         numfaces += radiusSegments * (((innerR[0] <= 0) || (innerR[1] <= 0)) ? 1 : 2);

      // upper cap
      if (outerR[0] > 0) numfaces += radiusSegments * ((innerR[0]>0) ? 2 : 1);
      // bottom cup
      if (outerR[1] > 0) numfaces += radiusSegments * ((innerR[1]>0) ? 2 : 1);

      if (thetaLength < 360)
         numfaces += ((outerR[0] > innerR[0]) ? 2 : 0) + ((outerR[1] > innerR[1]) ? 2 : 0);

      if (faces_limit < 0) return numfaces;

      let phi0 = thetaStart*Math.PI/180,
          dphi = thetaLength/radiusSegments*Math.PI/180,
          _sin = new Float32Array(radiusSegments+1),
          _cos = new Float32Array(radiusSegments+1);

      for (let seg=0; seg<=radiusSegments; ++seg) {
         _cos[seg] = Math.cos(phi0+seg*dphi);
         _sin[seg] = Math.sin(phi0+seg*dphi);
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      let calcZ;
      if (shape._typename == "TGeoCtub")
         calcZ = function(x,y,z) {
            let arr = (z<0) ? shape.fNlow : shape.fNhigh;
            return ((z<0) ? -shape.fDz : shape.fDz) - (x*arr[0] + y*arr[1]) / arr[2];
         }

      // create outer/inner tube
      for (let side = 0; side<2; ++side) {
         if ((side === 1) && !hasrmin) break;

         let R = (side === 0) ? outerR : innerR,
             d1 = side, d2 = 1 - side, nxy = 1., nz = 0;

         if (R[0] !== R[1]) {
            let angle = Math.atan2((R[1]-R[0]), 2*shape.fDZ);
            nxy = Math.cos(angle);
            nz = Math.sin(angle);
         }

         if (side === 1) { nxy *= -1; nz *= -1; };

         let reduce = 0;
         if (R[0] <= 0) reduce = 2; else
         if (R[1] <= 0) reduce = 1;

         for (let seg=0;seg<radiusSegments;++seg) {
            creator.addFace4(
                  R[0] * _cos[seg+d1], R[0] * _sin[seg+d1],  shape.fDZ,
                  R[1] * _cos[seg+d1], R[1] * _sin[seg+d1], -shape.fDZ,
                  R[1] * _cos[seg+d2], R[1] * _sin[seg+d2], -shape.fDZ,
                  R[0] * _cos[seg+d2], R[0] * _sin[seg+d2],  shape.fDZ,
                  reduce );

            if (calcZ) creator.recalcZ(calcZ);

            creator.setNormal_12_34(nxy*_cos[seg+d1], nxy*_sin[seg+d1], nz,
                                    nxy*_cos[seg+d2], nxy*_sin[seg+d2], nz,
                                    reduce);
         }
      }

      // create upper/bottom part
      for (let side = 0; side<2; ++side) {
         if (outerR[side] <= 0) continue;

         let d1 = side, d2 = 1- side,
             sign = (side == 0) ? 1 : -1,
             reduce = (innerR[side] <= 0) ? 2 : 0;
         if ((reduce==2) && (thetaLength === 360) && !calcZ) creator.startPolygon(side===0);
         for (let seg=0;seg<radiusSegments;++seg) {
            creator.addFace4(
                  innerR[side] * _cos[seg+d1], innerR[side] * _sin[seg+d1], sign*shape.fDZ,
                  outerR[side] * _cos[seg+d1], outerR[side] * _sin[seg+d1], sign*shape.fDZ,
                  outerR[side] * _cos[seg+d2], outerR[side] * _sin[seg+d2], sign*shape.fDZ,
                  innerR[side] * _cos[seg+d2], innerR[side] * _sin[seg+d2], sign*shape.fDZ,
                  reduce);
            if (calcZ) {
               creator.recalcZ(calcZ);
               creator.calcNormal();
            } else {
               creator.setNormal(0,0,sign);
            }
         }

         creator.stopPolygon();
      }

      // create cut surfaces
      if (thetaLength < 360) {
         creator.addFace4(innerR[1] * _cos[0], innerR[1] * _sin[0], -shape.fDZ,
                          outerR[1] * _cos[0], outerR[1] * _sin[0], -shape.fDZ,
                          outerR[0] * _cos[0], outerR[0] * _sin[0],  shape.fDZ,
                          innerR[0] * _cos[0], innerR[0] * _sin[0],  shape.fDZ,
                          (outerR[0] === innerR[0]) ? 2 : ((innerR[1]===outerR[1]) ? 1 : 0) );
         if (calcZ) creator.recalcZ(calcZ);
         creator.calcNormal();

         creator.addFace4(innerR[0] * _cos[radiusSegments], innerR[0] * _sin[radiusSegments],  shape.fDZ,
                          outerR[0] * _cos[radiusSegments], outerR[0] * _sin[radiusSegments],  shape.fDZ,
                          outerR[1] * _cos[radiusSegments], outerR[1] * _sin[radiusSegments], -shape.fDZ,
                          innerR[1] * _cos[radiusSegments], innerR[1] * _sin[radiusSegments], -shape.fDZ,
                          (outerR[0] === innerR[0]) ? 1 : ((innerR[1]===outerR[1]) ? 2 : 0));

         if (calcZ) creator.recalcZ(calcZ);
         creator.calcNormal();
      }

      return creator.create();
   }

   /** @summary Creates eltu geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createEltuBuffer( shape , faces_limit ) {
      let radiusSegments = Math.max(4, Math.round(360/geo.GradPerSegm));

      if (faces_limit < 0) return radiusSegments*4;

      // calculate all sin/cos tables in advance
      let x = new Float32Array(radiusSegments+1),
          y = new Float32Array(radiusSegments+1);
      for (let seg=0; seg<=radiusSegments; ++seg) {
          let phi = seg/radiusSegments*2*Math.PI;
          x[seg] = shape.fRmin*Math.cos(phi);
          y[seg] = shape.fRmax*Math.sin(phi);
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(radiusSegments*4),
          nx1 = 1, ny1 = 0, nx2 = 1, ny2 = 0;

      // create tube faces
      for (let seg=0; seg<radiusSegments; ++seg) {
         creator.addFace4(x[seg],   y[seg],   +shape.fDZ,
                          x[seg],   y[seg],   -shape.fDZ,
                          x[seg+1], y[seg+1], -shape.fDZ,
                          x[seg+1], y[seg+1],  shape.fDZ);

         // calculate normals ourself
         nx1 = nx2; ny1 = ny2;
         nx2 = x[seg+1] * shape.fRmax / shape.fRmin;
         ny2 = y[seg+1] * shape.fRmin / shape.fRmax;
         let dist = Math.sqrt(nx2*nx2 + ny2*ny2);
         nx2 = nx2 / dist; ny2 = ny2/dist;

         creator.setNormal_12_34(nx1,ny1,0,nx2,ny2,0);
      }

      // create top/bottom sides
      for (let side=0;side<2;++side) {
         let sign = (side===0) ? 1 : -1, d1 = side, d2 = 1 - side;
         for (let seg=0; seg<radiusSegments; ++seg) {
            creator.addFace3(0,          0,          sign*shape.fDZ,
                             x[seg+d1],  y[seg+d1],  sign*shape.fDZ,
                             x[seg+d2],  y[seg+d2],  sign*shape.fDZ);
            creator.setNormal(0, 0, sign);
         }
      }

      return creator.create();
   }

   /** @summary Creates torus geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createTorusBuffer( shape, faces_limit ) {
      let radius = shape.fR,
          radialSegments = Math.max(6, Math.round(360/geo.GradPerSegm)),
          tubularSegments = Math.max(8, Math.round(shape.fDphi/geo.GradPerSegm));

      let numfaces = (shape.fRmin > 0 ? 4 : 2) * radialSegments * (tubularSegments + (shape.fDphi !== 360 ? 1 : 0));

      if (faces_limit < 0) return numfaces;

      if ((faces_limit > 0) && (numfaces > faces_limit)) {
         radialSegments = Math.floor(radialSegments/Math.sqrt(numfaces / faces_limit));
         tubularSegments = Math.floor(tubularSegments/Math.sqrt(numfaces / faces_limit));
         numfaces = (shape.fRmin > 0 ? 4 : 2) * radialSegments * (tubularSegments + (shape.fDphi !== 360 ? 1 : 0));
      }

      let _sinr = new Float32Array(radialSegments+1),
          _cosr = new Float32Array(radialSegments+1),
          _sint = new Float32Array(tubularSegments+1),
          _cost = new Float32Array(tubularSegments+1);

      for (let n=0;n<=radialSegments;++n) {
         _sinr[n] = Math.sin(n/radialSegments*2*Math.PI);
         _cosr[n] = Math.cos(n/radialSegments*2*Math.PI);
      }

      for (let t=0;t<=tubularSegments;++t) {
         let angle = (shape.fPhi1 + shape.fDphi*t/tubularSegments)/180*Math.PI;
         _sint[t] = Math.sin(angle);
         _cost[t] = Math.cos(angle);
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      // use vectors for normals calculation
      let p1 = new THREE.Vector3(), p2 = new THREE.Vector3(), p3 = new THREE.Vector3(), p4 = new THREE.Vector3(),
          n1 = new THREE.Vector3(), n2 = new THREE.Vector3(), n3 = new THREE.Vector3(), n4 = new THREE.Vector3(),
          center1 = new THREE.Vector3(), center2 = new THREE.Vector3();

      for (let side=0;side<2;++side) {
         if ((side > 0) && (shape.fRmin <= 0)) break;
         let tube = (side > 0) ? shape.fRmin : shape.fRmax,
             d1 = 1 - side, d2 = 1 - d1, ns = side>0 ? -1 : 1;

         for (let t=0;t<tubularSegments;++t) {
            let t1 = t + d1, t2 = t + d2;
            center1.x = radius * _cost[t1]; center1.y = radius * _sint[t1];
            center2.x = radius * _cost[t2]; center2.y = radius * _sint[t2];

            for (let n=0;n<radialSegments;++n) {
               p1.x = (radius + tube * _cosr[n])   * _cost[t1]; p1.y = (radius + tube * _cosr[n])   * _sint[t1]; p1.z = tube*_sinr[n];
               p2.x = (radius + tube * _cosr[n+1]) * _cost[t1]; p2.y = (radius + tube * _cosr[n+1]) * _sint[t1]; p2.z = tube*_sinr[n+1];
               p3.x = (radius + tube * _cosr[n+1]) * _cost[t2]; p3.y = (radius + tube * _cosr[n+1]) * _sint[t2]; p3.z = tube*_sinr[n+1];
               p4.x = (radius + tube * _cosr[n])   * _cost[t2]; p4.y = (radius + tube * _cosr[n])   * _sint[t2]; p4.z = tube*_sinr[n];

               creator.addFace4(p1.x, p1.y, p1.z,
                                p2.x, p2.y, p2.z,
                                p3.x, p3.y, p3.z,
                                p4.x, p4.y, p4.z);

               n1.subVectors( p1, center1 ).normalize();
               n2.subVectors( p2, center1 ).normalize();
               n3.subVectors( p3, center2 ).normalize();
               n4.subVectors( p4, center2 ).normalize();

               creator.setNormal4(ns*n1.x, ns*n1.y, ns*n1.z,
                                  ns*n2.x, ns*n2.y, ns*n2.z,
                                  ns*n3.x, ns*n3.y, ns*n3.z,
                                  ns*n4.x, ns*n4.y, ns*n4.z);
            }
         }
      }

      if (shape.fDphi !== 360)
         for (let t=0;t<=tubularSegments;t+=tubularSegments) {
            let tube1 = shape.fRmax, tube2 = shape.fRmin,
                d1 = (t>0) ? 0 : 1, d2 = 1 - d1,
                skip = (shape.fRmin) > 0 ?  0 : 1,
                nsign = t>0 ? 1 : -1;
            for (let n=0;n<radialSegments;++n) {
               creator.addFace4((radius + tube1 * _cosr[n+d1]) * _cost[t], (radius + tube1 * _cosr[n+d1]) * _sint[t], tube1*_sinr[n+d1],
                                (radius + tube2 * _cosr[n+d1]) * _cost[t], (radius + tube2 * _cosr[n+d1]) * _sint[t], tube2*_sinr[n+d1],
                                (radius + tube2 * _cosr[n+d2]) * _cost[t], (radius + tube2 * _cosr[n+d2]) * _sint[t], tube2*_sinr[n+d2],
                                (radius + tube1 * _cosr[n+d2]) * _cost[t], (radius + tube1 * _cosr[n+d2]) * _sint[t], tube1*_sinr[n+d2], skip);
               creator.setNormal(-nsign* _sint[t], nsign * _cost[t], 0);
            }
         }

      return creator.create();
   }


   /** @summary Creates polygon geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createPolygonBuffer( shape, faces_limit ) {
      let thetaStart = shape.fPhi1,
          thetaLength = shape.fDphi,
          radiusSegments = 60, factor = 1;

      if (shape._typename == "TGeoPgon") {
         radiusSegments = shape.fNedges;
         factor = 1. / Math.cos(Math.PI/180 * thetaLength / radiusSegments / 2);
      } else {
         radiusSegments = Math.max(5, Math.round(thetaLength/geo.GradPerSegm));
      }

      let usage = new Int16Array(2*shape.fNz), numusedlayers = 0, hasrmin = false;

      for (let layer=0; layer < shape.fNz; ++layer)
         if (shape.fRmin[layer] > 0) hasrmin = true;

      // return very rough estimation, number of faces may be much less
      if (faces_limit < 0) return (hasrmin ? 4 : 2) * radiusSegments * (shape.fNz-1);

      // coordinate of point on cut edge (x,z)
      let pnts = (thetaLength === 360) ? null : [];

      // first analyse levels - if we need to create all of them
      for (let side = 0; side < 2; ++side) {
         let rside = (side === 0) ? 'fRmax' : 'fRmin';

         for (let layer=0; layer < shape.fNz; ++layer) {

            // first create points for the layer
            let layerz = shape.fZ[layer], rad = shape[rside][layer];

            usage[layer*2+side] = 0;

            if ((layer > 0) && (layer < shape.fNz-1))
               if (((shape.fZ[layer-1] === layerz) && (shape[rside][layer-1] === rad)) ||
                   ((shape[rside][layer+1] === rad) && (shape[rside][layer-1] === rad))) {

                  // same Z and R as before - ignore
                  // or same R before and after

                  continue;
               }

            if ((layer>0) && ((side === 0) || hasrmin)) {
               usage[layer*2+side] = 1;
               numusedlayers++;
            }

            if (pnts !== null) {
               if (side === 0) {
                  pnts.push(new THREE.Vector2(factor*rad, layerz));
               } else if (rad < shape.fRmax[layer]) {
                  pnts.unshift(new THREE.Vector2(factor*rad, layerz));
               }
            }
         }
      }

      let numfaces = numusedlayers*radiusSegments*2;
      if (shape.fRmin[0] !== shape.fRmax[0]) numfaces += radiusSegments * (hasrmin ? 2 : 1);
      if (shape.fRmin[shape.fNz-1] !== shape.fRmax[shape.fNz-1]) numfaces += radiusSegments * (hasrmin ? 2 : 1);

      let cut_faces = null;

      if (pnts!==null) {
         if (pnts.length === shape.fNz * 2) {
            // special case - all layers are there, create faces ourself
            cut_faces = [];
            for (let layer = shape.fNz-1; layer>0; --layer) {
               if (shape.fZ[layer] === shape.fZ[layer-1]) continue;
               let right = 2*shape.fNz - 1 - layer;
               cut_faces.push([right, layer - 1, layer]);
               cut_faces.push([right, right + 1, layer-1]);
            }

         } else {
            // let three.js calculate our faces
            // console.log('triangulate polygon ' + shape.fShapeId);
            cut_faces = THREE.ShapeUtils.triangulateShape(pnts, []);
         }
         numfaces += cut_faces.length*2;
      }

      let phi0 = thetaStart*Math.PI/180, dphi = thetaLength/radiusSegments*Math.PI/180;

      // calculate all sin/cos tables in advance
      let _sin = new Float32Array(radiusSegments+1),
          _cos = new Float32Array(radiusSegments+1);
      for (let seg=0;seg<=radiusSegments;++seg) {
         _cos[seg] = Math.cos(phi0+seg*dphi);
         _sin[seg] = Math.sin(phi0+seg*dphi);
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      // add sides
      for (let side = 0; side < 2; ++side) {
         let rside = (side === 0) ? 'fRmax' : 'fRmin',
             z1 = shape.fZ[0], r1 = factor*shape[rside][0],
             d1 = 1 - side, d2 = side;

         for (let layer=0; layer < shape.fNz; ++layer) {

            if (usage[layer*2+side] === 0) continue;

            let z2 = shape.fZ[layer], r2 = factor*shape[rside][layer],
                nxy = 1, nz = 0;

            if ((r2 !== r1)) {
               let angle = Math.atan2((r2-r1), (z2-z1));
               nxy = Math.cos(angle);
               nz = Math.sin(angle);
            }

            if (side>0) { nxy*=-1; nz*=-1; }

            for (let seg=0;seg < radiusSegments;++seg) {
               creator.addFace4(r1 * _cos[seg+d1], r1 * _sin[seg+d1], z1,
                                r2 * _cos[seg+d1], r2 * _sin[seg+d1], z2,
                                r2 * _cos[seg+d2], r2 * _sin[seg+d2], z2,
                                r1 * _cos[seg+d2], r1 * _sin[seg+d2], z1);
               creator.setNormal_12_34(nxy*_cos[seg+d1], nxy*_sin[seg+d1], nz, nxy*_cos[seg+d2], nxy*_sin[seg+d2], nz);
            }

            z1 = z2; r1 = r2;
         }
      }

      // add top/bottom
      for (let layer=0; layer < shape.fNz; layer += (shape.fNz-1)) {

         let rmin = factor*shape.fRmin[layer], rmax = factor*shape.fRmax[layer];

         if (rmin === rmax) continue;

         let layerz = shape.fZ[layer],
             d1 = (layer===0) ? 1 : 0, d2 = 1 - d1,
             normalz = (layer===0) ? -1: 1;

         if (!hasrmin && !cut_faces) creator.startPolygon(layer>0);

         for (let seg=0;seg < radiusSegments;++seg) {
            creator.addFace4(rmin * _cos[seg+d1], rmin * _sin[seg+d1], layerz,
                             rmax * _cos[seg+d1], rmax * _sin[seg+d1], layerz,
                             rmax * _cos[seg+d2], rmax * _sin[seg+d2], layerz,
                             rmin * _cos[seg+d2], rmin * _sin[seg+d2], layerz,
                             hasrmin ? 0 : 2);
            creator.setNormal(0, 0, normalz);
         }

         creator.stopPolygon();
      }

      if (cut_faces)
         for (let seg = 0; seg <= radiusSegments; seg += radiusSegments) {
            let d1 = (seg === 0) ? 1 : 2, d2 = 3 - d1;
            for (let n=0;n<cut_faces.length;++n) {
               let a = pnts[cut_faces[n][0]],
                   b = pnts[cut_faces[n][d1]],
                   c = pnts[cut_faces[n][d2]];

               creator.addFace3(a.x * _cos[seg], a.x * _sin[seg], a.y,
                                b.x * _cos[seg], b.x * _sin[seg], b.y,
                                c.x * _cos[seg], c.x * _sin[seg], c.y);

               creator.calcNormal();
            }
         }

      return creator.create();
   }

   /** @summary Creates xtru geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createXtruBuffer( shape, faces_limit ) {
      let nfaces = (shape.fNz-1) * shape.fNvert * 2;

      if (faces_limit < 0) return nfaces + shape.fNvert*3;

      // create points
      let pnts = [];
      for (let vert = 0; vert < shape.fNvert; ++vert)
         pnts.push(new THREE.Vector2(shape.fX[vert], shape.fY[vert]));

      // console.log('triangulate Xtru ' + shape.fShapeId);
      let faces = THREE.ShapeUtils.triangulateShape(pnts , []);
      if (faces.length < pnts.length-2) {
         geo.warn('Problem with XTRU shape ' +shape.fName + ' with ' + pnts.length + ' vertices');
         faces = [];
      } else {
         nfaces += faces.length * 2;
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(nfaces);

      for (let layer = 0; layer < shape.fNz-1; ++layer) {
         let z1 = shape.fZ[layer], scale1 = shape.fScale[layer],
             z2 = shape.fZ[layer+1], scale2 = shape.fScale[layer+1],
             x01 = shape.fX0[layer], x02 = shape.fX0[layer+1],
             y01 = shape.fY0[layer], y02 = shape.fY0[layer+1];

         for (let vert1 = 0; vert1 < shape.fNvert; ++vert1) {
            let vert2 = (vert1+1) % shape.fNvert;
            creator.addFace4(scale1 * shape.fX[vert1] + x01, scale1 * shape.fY[vert1] + y01, z1,
                             scale2 * shape.fX[vert1] + x02, scale2 * shape.fY[vert1] + y02, z2,
                             scale2 * shape.fX[vert2] + x02, scale2 * shape.fY[vert2] + y02, z2,
                             scale1 * shape.fX[vert2] + x01, scale1 * shape.fY[vert2] + y01, z1);
            creator.calcNormal();
         }
      }

      for (let layer = 0; layer <= shape.fNz-1; layer+=(shape.fNz-1)) {
         let z = shape.fZ[layer], scale = shape.fScale[layer],
             x0 = shape.fX0[layer], y0 = shape.fY0[layer];

         for (let n=0;n<faces.length;++n) {
            let face = faces[n],
                pnt1 = pnts[face[0]],
                pnt2 = pnts[face[(layer===0) ? 2 : 1]],
                pnt3 = pnts[face[(layer===0) ? 1 : 2]];

            creator.addFace3(scale * pnt1.x + x0, scale * pnt1.y + y0, z,
                             scale * pnt2.x + x0, scale * pnt2.y + y0, z,
                             scale * pnt3.x + x0, scale * pnt3.y + y0, z);
            creator.setNormal(0,0,layer===0 ? -1 : 1);
         }
      }

      return creator.create();
   }

   /** @summary Creates para geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createParaboloidBuffer( shape, faces_limit ) {

      let radiusSegments = Math.max(4, Math.round(360/geo.GradPerSegm)),
          heightSegments = 30;

      if (faces_limit > 0) {
         let fact = 2*radiusSegments*(heightSegments+1) / faces_limit;
         if (fact > 1.) {
            radiusSegments = Math.max(5, Math.floor(radiusSegments/Math.sqrt(fact)));
            heightSegments = Math.max(5, Math.floor(heightSegments/Math.sqrt(fact)));
         }
      }

      let zmin = -shape.fDZ, zmax = shape.fDZ, rmin = shape.fRlo, rmax = shape.fRhi;

      // if no radius at -z, find intersection
      if (shape.fA >= 0) {
         if (shape.fB > zmin) zmin = shape.fB;
      } else {
         if (shape.fB < zmax) zmax = shape.fB;
      }

      let ttmin = Math.atan2(zmin, rmin), ttmax = Math.atan2(zmax, rmax);

      let numfaces = (heightSegments+1)*radiusSegments*2;
      if (rmin===0) numfaces -= radiusSegments*2; // complete layer
      if (rmax===0) numfaces -= radiusSegments*2; // complete layer

      if (faces_limit < 0) return numfaces;

      // calculate all sin/cos tables in advance
      let _sin = new Float32Array(radiusSegments+1),
          _cos = new Float32Array(radiusSegments+1);
      for (let seg=0;seg<=radiusSegments;++seg) {
         _cos[seg] = Math.cos(seg/radiusSegments*2*Math.PI);
         _sin[seg] = Math.sin(seg/radiusSegments*2*Math.PI);
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      let lastz = zmin, lastr = 0, lastnxy = 0, lastnz = -1;

      for (let layer = 0; layer <= heightSegments + 1; ++layer) {

         let layerz = 0, radius = 0, nxy = 0, nz = -1;

         if ((layer === 0) && (rmin===0)) continue;

         if ((layer === heightSegments + 1) && (lastr === 0)) break;

         switch (layer) {
            case 0: layerz = zmin; radius = rmin; break;
            case heightSegments: layerz = zmax; radius = rmax; break;
            case heightSegments + 1: layerz = zmax; radius = 0; break;
            default: {
               let tt = Math.tan(ttmin + (ttmax-ttmin) * layer / heightSegments);
               let delta = tt*tt - 4*shape.fA*shape.fB; // should be always positive (a*b<0)
               radius = 0.5*(tt+Math.sqrt(delta))/shape.fA;
               if (radius < 1e-6) radius = 0;
               layerz = radius*tt;
            }
         }

         nxy = shape.fA * radius;
         nz = (shape.fA > 0) ? -1 : 1;

         let skip = 0;
         if (lastr === 0) skip = 1; else
         if (radius === 0) skip = 2;

         for (let seg=0; seg<radiusSegments; ++seg) {
            creator.addFace4(radius*_cos[seg],   radius*_sin[seg], layerz,
                             lastr*_cos[seg],    lastr*_sin[seg], lastz,
                             lastr*_cos[seg+1],  lastr*_sin[seg+1], lastz,
                             radius*_cos[seg+1], radius*_sin[seg+1], layerz, skip);

            // use analytic normal values when open/closing paraboloid around 0
            // cut faces (top or bottom) set with simple normal
            if ((skip===0) || ((layer===1) && (rmin===0)) || ((layer===heightSegments+1) && (rmax===0)))
               creator.setNormal4(nxy*_cos[seg],       nxy*_sin[seg],       nz,
                                  lastnxy*_cos[seg],   lastnxy*_sin[seg],   lastnz,
                                  lastnxy*_cos[seg+1], lastnxy*_sin[seg+1], lastnz,
                                  nxy*_cos[seg+1],     nxy*_sin[seg+1],     nz, skip);
            else
               creator.setNormal(0, 0, (layer < heightSegments) ? -1 : 1);
         }

         lastz = layerz; lastr = radius;
         lastnxy = nxy; lastnz = nz;
      }

      return creator.create();
   }

   /** @summary Creates hype geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createHypeBuffer( shape, faces_limit ) {

      if ((shape.fTin===0) && (shape.fTout===0))
         return createTubeBuffer(shape, faces_limit);

      let radiusSegments = Math.max(4, Math.round(360/geo.GradPerSegm)),
          heightSegments = 30;

      let numfaces = radiusSegments * (heightSegments + 1) * ((shape.fRmin > 0) ? 4 : 2);

      if (faces_limit < 0) return numfaces;

      if ((faces_limit > 0) && (faces_limit > numfaces)) {
         radiusSegments = Math.max(4, Math.floor(radiusSegments/Math.sqrt(numfaces/faces_limit)));
         heightSegments = Math.max(4, Math.floor(heightSegments/Math.sqrt(numfaces/faces_limit)));
         numfaces = radiusSegments * (heightSegments + 1) * ((shape.fRmin > 0) ? 4 : 2);
      }

      // calculate all sin/cos tables in advance
      let _sin = new Float32Array(radiusSegments+1), _cos = new Float32Array(radiusSegments+1);
      for (let seg=0;seg<=radiusSegments;++seg) {
         _cos[seg] = Math.cos(seg/radiusSegments*2*Math.PI);
         _sin[seg] = Math.sin(seg/radiusSegments*2*Math.PI);
      }

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      // in-out side
      for (let side=0;side<2;++side) {
         if ((side > 0) && (shape.fRmin <= 0)) break;

         let r0 = (side > 0) ? shape.fRmin : shape.fRmax,
             tsq = (side > 0) ? shape.fTinsq : shape.fToutsq,
             d1 = 1- side, d2 = 1 - d1;

         // vertical layers
         for (let layer=0;layer<heightSegments;++layer) {
            let z1 = -shape.fDz + layer/heightSegments*2*shape.fDz,
                z2 = -shape.fDz + (layer+1)/heightSegments*2*shape.fDz,
                r1 = Math.sqrt(r0*r0+tsq*z1*z1),
                r2 = Math.sqrt(r0*r0+tsq*z2*z2);

            for (let seg=0; seg<radiusSegments; ++seg) {
               creator.addFace4(r1 * _cos[seg+d1], r1 * _sin[seg+d1], z1,
                                r2 * _cos[seg+d1], r2 * _sin[seg+d1], z2,
                                r2 * _cos[seg+d2], r2 * _sin[seg+d2], z2,
                                r1 * _cos[seg+d2], r1 * _sin[seg+d2], z1);
               creator.calcNormal();
            }
         }
      }

      // add caps
      for(let layer=0; layer<2; ++layer) {
         let z = (layer === 0) ? shape.fDz : -shape.fDz,
             r1 = Math.sqrt(shape.fRmax*shape.fRmax + shape.fToutsq*z*z),
             r2 = (shape.fRmin > 0) ? Math.sqrt(shape.fRmin*shape.fRmin + shape.fTinsq*z*z) : 0,
             skip = (shape.fRmin > 0) ? 0 : 1,
             d1 = 1 - layer, d2 = 1 - d1;
          for (let seg=0; seg<radiusSegments; ++seg) {
             creator.addFace4(r1 * _cos[seg+d1], r1 * _sin[seg+d1], z,
                              r2 * _cos[seg+d1], r2 * _sin[seg+d1], z,
                              r2 * _cos[seg+d2], r2 * _sin[seg+d2], z,
                              r1 * _cos[seg+d2], r1 * _sin[seg+d2], z, skip);
             creator.setNormal(0,0, (layer===0) ? 1 : -1)
          }

      }

      return creator.create();
   }

   /** @summary Creates tessalated geometrey
     * @memberof JSROOT.GEO
     * @private */
   function createTessellatedBuffer( shape, faces_limit) {
      let numfaces = 0;

      for (let i = 0; i < shape.fFacets.length; ++i) {
         let f = shape.fFacets[i];
         if (f.fNvert == 4) numfaces += 2;
                       else numfaces += 1;
      }

      if (faces_limit < 0) return numfaces;

      let creator = faces_limit ? new PolygonsCreator : new GeometryCreator(numfaces);

      for (let i = 0; i < shape.fFacets.length; ++i) {
         let f = shape.fFacets[i],
             v0 = shape.fVertices[f.fIvert[0]].fVec,
             v1 = shape.fVertices[f.fIvert[1]].fVec,
             v2 = shape.fVertices[f.fIvert[2]].fVec;

         if (f.fNvert == 4) {
            let v3 = shape.fVertices[f.fIvert[3]].fVec;
            creator.addFace4(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]);
            creator.calcNormal();
         } else {
            creator.addFace3(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
            creator.calcNormal();
         }
      }

      return creator.create();
   }

   /** @summary Creates THREE.Matrix4 from TGeoMatrix
     * @memberof JSROOT.GEO
     * @private */
   geo.createMatrix = function(matrix) {

      if (!matrix) return null;

      let translation, rotation, scale;

      switch (matrix._typename) {
         case 'TGeoTranslation': translation = matrix.fTranslation; break;
         case 'TGeoRotation': rotation = matrix.fRotationMatrix; break;
         case 'TGeoScale': scale = matrix.fScale; break;
         case 'TGeoGenTrans':
            scale = matrix.fScale; // no break, translation and rotation follows
         case 'TGeoCombiTrans':
            translation = matrix.fTranslation;
            if (matrix.fRotation) rotation = matrix.fRotation.fRotationMatrix;
            break;
         case 'TGeoHMatrix':
            translation = matrix.fTranslation;
            rotation = matrix.fRotationMatrix;
            scale = matrix.fScale;
            break;
         case 'TGeoIdentity':
            break;
         default:
            console.warn('unsupported matrix ' + matrix._typename);
      }

      if (!translation && !rotation && !scale) return null;

      let res = new THREE.Matrix4();

      if (rotation)
         res.set(rotation[0], rotation[1], rotation[2],  0,
                 rotation[3], rotation[4], rotation[5],  0,
                 rotation[6], rotation[7], rotation[8],  0,
                           0,           0,           0,  1);

      if (translation)
         res.setPosition(translation[0], translation[1], translation[2]);

      if (scale)
         res.scale(new THREE.Vector3(scale[0], scale[1], scale[2]));

      return res;
   }

   /** @summary Creates transformation matrix for TGeoNode
     * @desc created after node visibility flag is checked and volume cut is performed
     * @memberof JSROOT.GEO
     * @private */
   geo.getNodeMatrix = function(kind, node) {

      let matrix = null;

      if (kind === kindEve) {
         // special handling for EVE nodes

         matrix = new THREE.Matrix4();

         if (node.fTrans) {
            matrix.set(node.fTrans[0],  node.fTrans[4],  node.fTrans[8],  0,
                       node.fTrans[1],  node.fTrans[5],  node.fTrans[9],  0,
                       node.fTrans[2],  node.fTrans[6],  node.fTrans[10], 0,
                                    0,               0,               0,  1);
            // second - set position with proper sign
            matrix.setPosition(node.fTrans[12], node.fTrans[13], node.fTrans[14]);
         }
      } else if (node.fMatrix) {
         matrix = geo.createMatrix(node.fMatrix);
      } else if ((node._typename == "TGeoNodeOffset") && node.fFinder) {
         let kPatternReflected = JSROOT.BIT(14);
         if ((node.fFinder.fBits & kPatternReflected) !== 0)
            geo.warn('Unsupported reflected pattern ' + node.fFinder._typename);

         // if (node.fFinder._typename === 'TGeoPatternCylR') { }
         // if (node.fFinder._typename === 'TGeoPatternSphR') { }
         // if (node.fFinder._typename === 'TGeoPatternSphTheta') { }
         // if (node.fFinder._typename === 'TGeoPatternSphPhi') { }
         // if (node.fFinder._typename === 'TGeoPatternHoneycomb') { }
         switch(node.fFinder._typename) {
           case 'TGeoPatternX':
           case 'TGeoPatternY':
           case 'TGeoPatternZ':
           case 'TGeoPatternParaX':
           case 'TGeoPatternParaY':
           case 'TGeoPatternParaZ':
              let _shift = node.fFinder.fStart + (node.fIndex + 0.5) * node.fFinder.fStep;

              matrix = new THREE.Matrix4();

              switch (node.fFinder._typename[node.fFinder._typename.length-1]) {
                 case 'X': matrix.setPosition(_shift, 0, 0); break;
                 case 'Y': matrix.setPosition(0, _shift, 0); break;
                 case 'Z': matrix.setPosition(0, 0, _shift); break;
              }
              break;

           case 'TGeoPatternCylPhi':
              let phi = (Math.PI/180)*(node.fFinder.fStart+(node.fIndex+0.5)*node.fFinder.fStep),
                  _cos = Math.cos(phi), _sin = Math.sin(phi);

              matrix = new THREE.Matrix4();

              matrix.set(_cos, -_sin,  0,  0,
                         _sin,  _cos,  0,  0,
                            0,     0,  1,  0,
                            0,     0,  0,  1);
              break;

           case 'TGeoPatternCylR':
               // seems to be, require no transformation
               matrix = new THREE.Matrix4();
               break;

           case 'TGeoPatternTrapZ':
              let dz = node.fFinder.fStart + (node.fIndex+0.5)*node.fFinder.fStep;
              matrix = new THREE.Matrix4();
              matrix.setPosition(node.fFinder.fTxz*dz, node.fFinder.fTyz*dz, dz);
              break;

           default:
              geo.warn('Unsupported pattern type ' + node.fFinder._typename);
              break;
         }
      }

      return matrix;
   }

   let createGeometry; // will be function to create geometry

   /** @summary Returns number of faces for provided geometry
    * @param {Object} geom  - can be THREE.Geometry, THREE.BufferGeometry, ThreeBSP.Geometry or interim array of polygons
    * @private */
   geo.numGeometryFaces = function(geom) {
      if (!geom) return 0;

      if (geom instanceof ThreeBSP.Geometry)
         return geom.tree.numPolygons();

      if (geom.type == 'BufferGeometry') {
         var attr = geom.getAttribute('position');
         return attr && attr.count ? Math.round(attr.count / 3) : 0;
      }

      // special array of polygons
      if (geom.polygons)
         return geom.polygons.length;

      return geom.faces.length;
   }

   /** @summary Returns number of faces for provided geometry
    * @param {Object} geom  - can be THREE.Geometry, THREE.BufferGeometry, ThreeBSP.Geometry or interim array of polygons
    * @private */
   geo.numGeometryVertices = function(geom) {
      if (!geom) return 0;

      if (geom instanceof ThreeBSP.Geometry)
         return geom.tree.numPolygons() * 3;

      if (geom.type == 'BufferGeometry') {
         var attr = geom.getAttribute('position');
         return attr ? attr.count : 0;
      }

      if (geom.polygons)
         return geom.polygons.length * 4;

      return geom.vertices.length;
   }

   /** @summary Returns geometry bounding box
    * @memberof JSROOT.GEO
    * @private */
   function geomBoundingBox(geom) {
      if (!geom) return null;

      let polygons = null;

      if (geom instanceof ThreeBSP.Geometry)
         polygons = geom.tree.collectPolygons([]);
      else if (geom.polygons)
         polygons = geom.polygons;

      if (polygons!==null) {
         let box = new THREE.Box3();
         for (let n=0;n<polygons.length;++n) {
            let polygon = polygons[n], nvert = polygon.vertices.length;
            for (let k=0;k<nvert;++k)
               box.expandByPoint(polygon.vertices[k]);
         }
         return box;
      }

      if (!geom.boundingBox) geom.computeBoundingBox();

      return geom.boundingBox.clone();
   }

   /** @summary Creates half-space geometry for given shape
    * @desc Just big-enough triangle to make BSP calculations
    * @memberof JSROOT.GEO
    * @private */
   function createHalfSpace(shape, geom) {
      if (!shape || !shape.fN || !shape.fP) return null;

      // shape.fP = [0,0,15]; shape.fN = [0,1,1];

      let vertex = new THREE.Vector3(shape.fP[0], shape.fP[1], shape.fP[2]);

      let normal = new THREE.Vector3(shape.fN[0], shape.fN[1], shape.fN[2]);
      normal.normalize();

      let sz = 1e10;
      if (geom) {
         // using real size of other geometry, we probably improve precision
         let box = geomBoundingBox(geom);
         if (box) sz = box.getSize(new THREE.Vector3()).length() * 1000;
      }

      // console.log('normal', normal, 'vertex', vertex, 'size', sz);

      let v1 = new THREE.Vector3(-sz, -sz/2, 0),
          v2 = new THREE.Vector3(0, sz, 0),
          v3 = new THREE.Vector3(sz, -sz/2, 0),
          v4 = new THREE.Vector3(0, 0, -sz);

      let geometry = new THREE.Geometry();

      geometry.vertices.push(v1, v2, v3, v4);

      geometry.faces.push( new THREE.Face3( 0, 2, 1 ) );
      geometry.faces.push( new THREE.Face3( 0, 1, 3 ) );
      geometry.faces.push( new THREE.Face3( 1, 2, 3 ) );
      geometry.faces.push( new THREE.Face3( 2, 0, 3 ) );

      geometry.lookAt(normal);
      geometry.computeFaceNormals();

      v1.add(vertex);
      v2.add(vertex);
      v3.add(vertex);
      v4.add(vertex);
      return geometry;
   }

   /** @summary Returns number of faces for provided geometry
     * @param geom  - can be THREE.Geometry, THREE.BufferGeometry, ThreeBSP.Geometry or interim array of polygons
     * @memberof JSROOT.GEO
     * @private */
   function countGeometryFaces(geom) {
      if (!geom) return 0;

      if (geom instanceof ThreeBSP.Geometry)
         return geom.tree.numPolygons();

      if (geom.type == 'BufferGeometry') {
         let attr = geom.getAttribute('position');
         return attr && attr.count ? Math.round(attr.count / 3) : 0;
      }

      // special array of polygons
      if (geom.polygons)
         return geom.polygons.length;

      return geom.faces.length;
   }

   /** @summary Creates geometrey for composite shape
     * @memberof JSROOT.GEO
     * @private */
   function createComposite( shape, faces_limit ) {

      if (faces_limit < 0)
         return createGeometry(shape.fNode.fLeft, -10) +
                createGeometry(shape.fNode.fRight, -10);

      let geom1, geom2, bsp1, bsp2, return_bsp = false,
          matrix1 = geo.createMatrix(shape.fNode.fLeftMat),
          matrix2 = geo.createMatrix(shape.fNode.fRightMat);

      if (faces_limit === 0) faces_limit = 4000;
                        else return_bsp = true;

      if (matrix1 && (matrix1.determinant() < -0.9))
         geo.warn('Axis reflection in left composite shape - not supported');

      if (matrix2 && (matrix2.determinant() < -0.9))
         geo.warn('Axis reflections in right composite shape - not supported');

      if (shape.fNode.fLeft._typename == "TGeoHalfSpace") {
         geom1 = createHalfSpace(shape.fNode.fLeft);
      } else {
         geom1 = createGeometry(shape.fNode.fLeft, faces_limit);
      }

      if (!geom1) return null;

      let n1 = countGeometryFaces(geom1), n2 = 0;
      if (geom1._exceed_limit) n1 += faces_limit;

      if (n1 < faces_limit) {

         if (shape.fNode.fRight._typename == "TGeoHalfSpace") {
            geom2 = createHalfSpace(shape.fNode.fRight, geom1);
         } else {
            geom2 = createGeometry(shape.fNode.fRight, faces_limit);
         }

         n2 = countGeometryFaces(geom2);
      }

      if ((n1 + n2 >= faces_limit) || !geom2) {
         if (geom1.polygons) {
            geom1 = ThreeBSP.CreateBufferGeometry(geom1.polygons);
            n1 = countGeometryFaces(geom1);
         }
         if (matrix1) geom1.applyMatrix4(matrix1);
         // if (!geom1._exceed_limit) console.log('reach faces limit', faces_limit, 'got', n1, n2);
         geom1._exceed_limit = true;
         return geom1;
      }

      bsp1 = new ThreeBSP.Geometry(geom1, matrix1, geo.CompressComp ? 0 : undefined);

      bsp2 = new ThreeBSP.Geometry(geom2, matrix2, bsp1.maxid);

      // take over maxid from both geometries
      bsp1.maxid = bsp2.maxid;

      switch(shape.fNode._typename) {
         case 'TGeoIntersection': bsp1.direct_intersect(bsp2);  break; // "*"
         case 'TGeoUnion': bsp1.direct_union(bsp2); break;   // "+"
         case 'TGeoSubtraction': bsp1.direct_subtract(bsp2); break; // "/"
         default:
            geo.warn('unsupported bool operation ' + shape.fNode._typename + ', use first geom');
      }

      if (countGeometryFaces(bsp1) === 0) {
         geo.warn('Zero faces in comp shape'
               + ' left: ' + shape.fNode.fLeft._typename +  ' ' + countGeometryFaces(geom1) + ' faces'
               + ' right: ' + shape.fNode.fRight._typename + ' ' + countGeometryFaces(geom2) + ' faces'
               + '  use first');
         bsp1 = new ThreeBSP.Geometry(geom1, matrix1);
      }

      return return_bsp ? { polygons: bsp1.toPolygons() } : bsp1.toBufferGeometry();
   }

   /** @summary Try to create projected geometry
     * @memberof JSROOT.GEO
     * @private */
   function projectGeometry(geom, matrix, projection, position, flippedMesh) {

      if (!geom.boundingBox) geom.computeBoundingBox();

      let box = geom.boundingBox.clone();

      box.applyMatrix4(matrix);

      if (!position) position = 0;

      if (((box.min[projection]>=position) && (box.max[projection]>=position)) ||
          ((box.min[projection]<=position) && (box.max[projection]<=position))) {
         return null; // not interesting
      }

      let bsp1 = new ThreeBSP.Geometry(geom, matrix, 0, flippedMesh),
          sizex = 2*Math.max(Math.abs(box.min.x), Math.abs(box.max.x)),
          sizey = 2*Math.max(Math.abs(box.min.y), Math.abs(box.max.y)),
          sizez = 2*Math.max(Math.abs(box.min.z), Math.abs(box.max.z)),
          size = 10000;

      switch (projection) {
         case "x": size = Math.max(sizey,sizez); break;
         case "y": size = Math.max(sizex,sizez); break;
         case "z": size = Math.max(sizex,sizey); break;
      }

      let bsp2 = ThreeBSP.CreateNormal(projection, position, size);

      bsp1.cut_from_plane(bsp2);

      return bsp2.toBufferGeometry();
   }

   /** @summary Creates geometry model for the provided shape
    * @desc
    *  - if limit === 0 (or undefined) returns THREE.BufferGeometry
    *  - if limit < 0 just returns estimated number of faces
    *  - if limit > 0 return list of ThreeBSP polygons (used only for composite shapes)
    * @param {Object} shape - instance of TGeoShape object
    * @param {Number} limit - defines return value, see details
    * @memberof JSROOT.GEO
    * @private */
   createGeometry = function( shape, limit ) {
      if (limit === undefined) limit = 0;

      try {
         switch (shape._typename) {
            case "TGeoBBox": return createCubeBuffer( shape, limit );
            case "TGeoPara": return createParaBuffer( shape, limit );
            case "TGeoTrd1":
            case "TGeoTrd2": return createTrapezoidBuffer( shape, limit );
            case "TGeoArb8":
            case "TGeoTrap":
            case "TGeoGtra": return createArb8Buffer( shape, limit );
            case "TGeoSphere": return createSphereBuffer( shape , limit );
            case "TGeoCone":
            case "TGeoConeSeg":
            case "TGeoTube":
            case "TGeoTubeSeg":
            case "TGeoCtub": return createTubeBuffer( shape, limit );
            case "TGeoEltu": return createEltuBuffer( shape, limit );
            case "TGeoTorus": return createTorusBuffer( shape, limit );
            case "TGeoPcon":
            case "TGeoPgon": return createPolygonBuffer( shape, limit );
            case "TGeoXtru": return createXtruBuffer( shape, limit );
            case "TGeoParaboloid": return createParaboloidBuffer( shape, limit );
            case "TGeoHype": return createHypeBuffer( shape, limit );
            case "TGeoTessellated": return createTessellatedBuffer( shape, limit );
            case "TGeoCompositeShape": return createComposite( shape, limit );
            case "TGeoShapeAssembly": break;
            case "TGeoScaledShape": {
               let res = createGeometry(shape.fShape, limit);
               if (shape.fScale && (limit>=0) && (typeof res === 'object') && (typeof res.scale === 'function'))
                  res.scale(shape.fScale.fScale[0],shape.fScale.fScale[1],shape.fScale.fScale[2]);
               return res;
            }
            case "TGeoHalfSpace":
               if (limit < 0) return 1; // half space if just plane used in composite
               // no break here - warning should appear
            default: geo.warn('unsupported shape type ' + shape._typename);
         }
      } catch(e) {
         let place = "";
         if (e.stack !== undefined) {
            place = e.stack.split("\n")[0];
            if (place.indexOf(e.message) >= 0) place = e.stack.split("\n")[1];
                                          else place = " at: " + place;
         }
         geo.warn(shape._typename + " err: " + e.message + place);
      }

      return limit < 0 ? 0 : null;
   }

   /** @summary Provides info about geo object, used for tooltip info
     * @param {Object} obj - any kind of TGeo-related object like shape or node or volume
     * @memberof JSROOT.GEO
     * @private */
   function provideObjectInfo(obj) {
      let info = [], shape = null;

      if (obj.fVolume !== undefined) shape = obj.fVolume.fShape; else
      if (obj.fShape !== undefined) shape = obj.fShape; else
      if ((obj.fShapeBits !== undefined) && (obj.fShapeId !== undefined)) shape = obj;

      if (!shape) {
         info.push(obj._typename);
         return info;
      }

      let sz = Math.max(shape.fDX, shape.fDY, shape.fDZ);
      let useexp = (sz>1e7) || (sz<1e-7);

      function conv(v) {
         if (v===undefined) return "???";
         if ((v==Math.round(v) && v<1e7)) return Math.round(v);
         return useexp ? v.toExponential(4) : v.toPrecision(7);
      }

      info.push(shape._typename);

      info.push("DX="+conv(shape.fDX) + " DY="+conv(shape.fDY) + " DZ="+conv(shape.fDZ));

      switch (shape._typename) {
         case "TGeoBBox": break;
         case "TGeoPara": info.push("Alpha=" + shape.fAlpha + " Phi=" + shape.fPhi + " Theta=" + shape.fTheta); break;
         case "TGeoTrd2": info.push("Dy1=" + conv(shape.fDy1) + " Dy2=" + conv(shape.fDy1));
         case "TGeoTrd1": info.push("Dx1=" + conv(shape.fDx1) + " Dx2=" + conv(shape.fDx1)); break;
         case "TGeoArb8": break;
         case "TGeoTrap": break;
         case "TGeoGtra": break;
         case "TGeoSphere":
            info.push("Rmin=" + conv(shape.fRmin) + " Rmax=" + conv(shape.fRmax));
            info.push("Phi1=" + shape.fPhi1 + " Phi2=" + shape.fPhi2);
            info.push("Theta1=" + shape.fTheta1 + " Theta2=" + shape.fTheta2);
            break;
         case "TGeoConeSeg":
            info.push("Phi1=" + shape.fPhi1 + " Phi2=" + shape.fPhi2);
         case "TGeoCone":
            info.push("Rmin1=" + conv(shape.fRmin1) + " Rmax1=" + conv(shape.fRmax1));
            info.push("Rmin2=" + conv(shape.fRmin2) + " Rmax2=" + conv(shape.fRmax2));
            break;
         case "TGeoCtub":
         case "TGeoTubeSeg":
            info.push("Phi1=" + shape.fPhi1 + " Phi2=" + shape.fPhi2);
         case "TGeoEltu":
         case "TGeoTube":
            info.push("Rmin=" + conv(shape.fRmin) + " Rmax=" + conv(shape.fRmax));
            break;
         case "TGeoTorus":
            info.push("Rmin=" + conv(shape.fRmin) + " Rmax=" + conv(shape.fRmax));
            info.push("Phi1=" + shape.fPhi1 + " Dphi=" + shape.fDphi);
            break;
         case "TGeoPcon":
         case "TGeoPgon": break;
         case "TGeoXtru": break;
         case "TGeoParaboloid":
            info.push("Rlo=" + conv(shape.fRlo) + " Rhi=" + conv(shape.fRhi));
            info.push("A=" + conv(shape.fA) + " B=" + conv(shape.fB));
            break;
         case "TGeoHype":
            info.push("Rmin=" + conv(shape.fRmin) + " Rmax=" + conv(shape.fRmax));
            info.push("StIn=" + conv(shape.fStIn) + " StOut=" + conv(shape.fStOut));
            break;
         case "TGeoCompositeShape": break;
         case "TGeoShapeAssembly": break;
         case "TGeoScaledShape":
            info = provideObjectInfo(shape.fShape);
            if (shape.fScale)
               info.unshift('Scale X=' + shape.fScale.fScale[0] + " Y=" + shape.fScale.fScale[1] + " Z=" + shape.fScale.fScale[2]);
            break;
      }

      return info;
   }

   /** @summary Creates projection matrix for the camera
     * @memberof JSROOT.GEO
     * @private */
   function createProjectionMatrix(camera) {
      let cameraProjectionMatrix = new THREE.Matrix4();

      camera.updateMatrixWorld();
      camera.matrixWorldInverse.getInverse( camera.matrixWorld );
      cameraProjectionMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse);

      return cameraProjectionMatrix;
   }

   /** @summary Creates frustum
     * @memberof JSROOT.GEO
     * @private */
   function createFrustum(source) {
      if (!source) return null;

      if (source instanceof THREE.PerspectiveCamera)
         source = createProjectionMatrix(source);

      let frustum = new THREE.Frustum();
      frustum.setFromProjectionMatrix(source);

      frustum.corners = new Float32Array([
          1,  1,  1,
          1,  1, -1,
          1, -1,  1,
          1, -1, -1,
         -1,  1,  1,
         -1,  1, -1,
         -1, -1,  1,
         -1, -1, -1,
          0,  0,  0 // also check center of the shape
      ]);

      frustum.test = new THREE.Vector3(0,0,0);

      frustum.CheckShape = function(matrix, shape) {
         let pnt = this.test, len = this.corners.length, corners = this.corners, i;

         for (i = 0; i < len; i+=3) {
            pnt.x = corners[i] * shape.fDX;
            pnt.y = corners[i+1] * shape.fDY;
            pnt.z = corners[i+2] * shape.fDZ;
            if (this.containsPoint(pnt.applyMatrix4(matrix))) return true;
        }

        return false;
      }

      frustum.CheckBox = function(box) {
         let pnt = this.test, cnt = 0;
         pnt.set(box.min.x, box.min.y, box.min.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.min.x, box.min.y, box.max.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.min.x, box.max.y, box.min.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.min.x, box.max.y, box.max.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.max.x, box.max.y, box.max.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.max.x, box.min.y, box.max.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.max.x, box.max.y, box.min.z);
         if (this.containsPoint(pnt)) cnt++;
         pnt.set(box.max.x, box.max.y, box.max.z);
         if (this.containsPoint(pnt)) cnt++;
         return cnt>5; // only if 6 edges and more are seen, we think that box is fully visible
      }

      return frustum;
   }

   /** @summary Compares two stacks.
     * @returns {number} length where stacks are the same
     * @private */
   geo.compareStacks = function(stack1, stack2) {
      if (!stack1 || !stack2) return 0;
      if (stack1 === stack2) return stack1.length;
      let len = Math.min(stack1.length, stack2.length);
      for (let k=0;k<len;++k)
         if (stack1[k] !== stack2[k]) return k;
      return len;
   }

   /** @summary Checks if two stack arrays are identical
     * @private */
   geo.isSameStack = function(stack1, stack2) {
      if (!stack1 || !stack2) return false;
      if (stack1 === stack2) return true;
      if (stack1.length !== stack2.length) return false;
      for (let k=0;k<stack1.length;++k)
         if (stack1[k] !== stack2[k]) return false;
      return true;
   }


   // ====================================================================

   /**
     * @summary class for working with cloned nodes
     *
     * @class
     * @memberof JSROOT.GEO
     * @private
     */

   function ClonedNodes(obj, clones) {
      this.toplevel = true; // indicate if object creates top-level structure with Nodes and Volumes folder
      this.name_prefix = ""; // name prefix used for nodes names
      this.maxdepth = 1;  // maximal hierarchy depth, required for transparency
      this.vislevel = 4;  // maximal depth of nodes visibility aka gGeoManager->SetVisLevel, same default
      this.maxnodes = 10000; // maximal number of visisble nodes aka gGeoManager->fMaxVisNodes

      if (obj) {
         if (obj.$geoh) this.toplevel = false;
         this.createClones(obj);
      } else if (clones) {
         this.nodes = clones;
      }
   }

   /** @summary Set maximal depth for nodes visibility */
   ClonedNodes.prototype.setVisLevel = function(lvl) {
      this.vislevel = lvl && !isNaN(lvl) ? lvl : 4;
   }

   /** @summary Returns maximal depth for nodes visibility */
   ClonedNodes.prototype.getVisLevel = function() {
      return this.vislevel;
   }

   /** @summary Set maximal number of visible nodes */
   ClonedNodes.prototype.setMaxVisNodes = function(v) {
      this.maxnodes = !isNaN(v) ? v : 10000;
   }

   /** @summary Returns configured maximal number of visible nodes */
   ClonedNodes.prototype.getMaxVisNodes = function() {
      return this.maxnodes;
   }

   /** @summary Insert node into existing array */
   ClonedNodes.prototype.updateNode = function(node) {
      if (node && !isNaN(node.id) && (node.id < this.nodes.length))
         this.nodes[node.id] = node;
   }

   /** @summary Returns TGeoShape for element with given indx */
   ClonedNodes.prototype.getNodeShape = function(indx) {
      if (!this.origin || !this.nodes) return null;
      let obj = this.origin[indx], clone = this.nodes[indx];
      if (!obj || !clone) return null;
      if (clone.kind === kindGeo) {
         if (obj.fVolume) return obj.fVolume.fShape;
      } else {
         return obj.fShape;
      }
      return null;
   }

   /** @summary function to cleanup as much as possible structures
     * @desc Provided parameters drawnodes and drawshapes are arrays created during building of geometry */
   ClonedNodes.prototype.cleanup = function(drawnodes, drawshapes) {

      if (drawnodes) {
         for (let n=0;n<drawnodes.length;++n) {
            delete drawnodes[n].stack;
            drawnodes[n] = undefined;
         }
      }

      if (drawshapes) {
         for (let n=0;n<drawshapes.length;++n) {
            delete drawshapes[n].geom;
            drawshapes[n] = undefined;
         }
      }

      if (this.nodes) {
         for (let n=0;n<this.nodes.length;++n) {
            if (this.nodes[n])
               delete this.nodes[n].chlds;
         }
      }

      delete this.nodes;
      delete this.origin;

      delete this.sortmap;
   }

   /** @summary Create complete description for provided Geo object
     * @private */
   ClonedNodes.prototype.createClones = function(obj, sublevel, kind) {
       if (!sublevel) {

          if (obj && obj._typename == "$$Shape$$")
             return this.createClonesForShape(obj);

          this.origin = [];
          sublevel = 1;
          kind = geo.getNodeKind(obj);
       }

       if ((kind < 0) || !obj || ('_refid' in obj)) return;

       obj._refid = this.origin.length;
       this.origin.push(obj);
       if (sublevel>this.maxdepth) this.maxdepth = sublevel;

       let chlds = null;
       if (kind === kindGeo)
          chlds = (obj.fVolume && obj.fVolume.fNodes) ? obj.fVolume.fNodes.arr : null;
       else
          chlds = obj.fElements ? obj.fElements.arr : null;

       if (chlds !== null) {
          geo.checkDuplicates(obj, chlds);
          for (let i = 0; i < chlds.length; ++i)
             this.createClones(chlds[i], sublevel+1, kind);
       }

       if (sublevel > 1) return;

       this.nodes = [];

       let sortarr = [];

       // first create nodes objects
       for (let n=0; n<this.origin.length; ++n) {
          // let obj = this.origin[n];
          let node = { id: n, kind: kind, vol: 0, nfaces: 0 };
          this.nodes.push(node);
          sortarr.push(node); // array use to produce sortmap
       }

       // than fill children lists
       for (let n=0;n<this.origin.length;++n) {
          let obj = this.origin[n], clone = this.nodes[n];

          let chlds = null, shape = null;

          if (kind === kindEve) {
             shape = obj.fShape;
             if (obj.fElements) chlds = obj.fElements.arr;
          } else if (obj.fVolume) {
             shape = obj.fVolume.fShape;
             if (obj.fVolume.fNodes) chlds = obj.fVolume.fNodes.arr;
          }

          let matrix = geo.getNodeMatrix(kind, obj);
          if (matrix) {
             clone.matrix = matrix.elements; // take only matrix elements, matrix will be constructed in worker
             if (clone.matrix[0] === 1) {
                let issimple = true;
                for (let k = 1; (k < clone.matrix.length) && issimple; ++k)
                   issimple = (clone.matrix[k] === ((k===5) || (k===10) || (k===15) ? 1 : 0));
                if (issimple) delete clone.matrix;
             }
             if (clone.matrix && (kind == kindEve)) clone.abs_matrix = true;
          }
          if (shape) {
             clone.fDX = shape.fDX;
             clone.fDY = shape.fDY;
             clone.fDZ = shape.fDZ;
             clone.vol = shape.fDX*shape.fDY*shape.fDZ;
             if (shape.$nfaces === undefined)
                shape.$nfaces = createGeometry(shape, -1);
             clone.nfaces = shape.$nfaces;
             if (clone.nfaces <= 0) clone.vol = 0;
          }

          if (!chlds) continue;

          // in cloned object children is only list of ids
          clone.chlds = new Array(chlds.length);
          for (let k=0;k<chlds.length;++k)
             clone.chlds[k] = chlds[k]._refid;
       }

       // remove _refid identifiers from original objects
       for (let n=0;n<this.origin.length;++n)
          delete this.origin[n]._refid;

       // do sorting once
       sortarr.sort(function(a,b) { return b.vol - a.vol; });

       // remember sort map and also sortid
       this.sortmap = new Array(this.nodes.length);
       for (let n=0;n<this.nodes.length;++n) {
          this.sortmap[n] = sortarr[n].id;
          sortarr[n].sortid = n;
       }
   }

   /** @summary Create elementary item with single already existing shape
     * @desc used by details view of geometry shape
     * @private */
   ClonedNodes.prototype.createClonesForShape = function(obj) {
      this.origin = [];

      // indicate that just plain shape is used
      this.plain_shape = obj;

      let node = {
            id: 0, sortid: 0, kind: kindShape,
            name: "Shape",
            nfaces: obj.nfaces,
            fDX: 1, fDY: 1, fDZ: 1, vol: 1,
            vis: true
         };

      this.nodes = [ node ];
   }

   /** @summary Count all visisble nodes */
   ClonedNodes.prototype.countVisibles = function() {
      let cnt = 0;
      if (this.nodes)
         for (let k=0;k<this.nodes.length;++k)
            if (this.nodes[k].vis)
               cnt++;
      return cnt;
   }

   /** @summary Mark visisble nodes.
     * @desc Set only basic flags, actual visibility depends from hierarchy */
   ClonedNodes.prototype.markVisibles = function(on_screen, copy_bits, hide_top_volume) {
      if (this.plain_shape) return 1;
      if (!this.origin || !this.nodes) return 0;

      let res = 0;

      for (let n=0;n<this.nodes.length;++n) {
         let clone = this.nodes[n],
             obj = this.origin[n];

         clone.vis = 0; // 1 - only with last level
         delete clone.nochlds;

         if (clone.kind === kindGeo) {
            if (obj.fVolume) {
               if (on_screen) {
                  // on screen bits used always, childs always checked
                  clone.vis = geo.TestBit(obj.fVolume, geo.BITS.kVisOnScreen) ? 99 : 0;

                  if ((n==0) && clone.vis && hide_top_volume) clone.vis = 0;

                  if (copy_bits) {
                     geo.SetBit(obj.fVolume, geo.BITS.kVisNone, false);
                     geo.SetBit(obj.fVolume, geo.BITS.kVisThis, (clone.vis > 0));
                     geo.SetBit(obj.fVolume, geo.BITS.kVisDaughters, true);
                  }
               } else {
                  clone.vis = !geo.TestBit(obj.fVolume, geo.BITS.kVisNone) &&
                               geo.TestBit(obj.fVolume, geo.BITS.kVisThis) ? 99 : 0;

                  if (!geo.TestBit(obj, geo.BITS.kVisDaughters) ||
                      !geo.TestBit(obj.fVolume, geo.BITS.kVisDaughters)) clone.nochlds = true;

                  // node with childs only shown in case if it is last level in hierarchy
                  if ((clone.vis > 0) && clone.chlds && !clone.nochlds) clone.vis = 1;

                  // special handling for top node
                  if (n==0) {
                     if (hide_top_volume) clone.vis = 0;
                     delete clone.nochlds;
                  }
               }
            }
         } else {
            clone.vis = obj.fRnrSelf ? 99 : 0;

            // when the only node is selected, draw it
            if ((n===0) && (this.nodes.length===1)) clone.vis = 99;

            this.vislevel = 9999; // automatically take all volumes
         }

         // shape with zero volume or without faces will not be observed
         if ((clone.vol <= 0) || (clone.nfaces <= 0)) clone.vis = 0;

         if (clone.vis) res++;
      }

      return res;
   }

   /** @summary After visibility flags is set, produce idshift for all nodes as it would be maximum level */
   ClonedNodes.prototype.produceIdShifts = function() {
      for (let k=0;k<this.nodes.length;++k)
         this.nodes[k].idshift = -1;

      function scan_func(nodes, node) {
         if (node.idshift < 0) {
            node.idshift = 0;
            if (node.chlds)
               for(let k = 0; k<node.chlds.length; ++k)
                  node.idshift += scan_func(nodes, nodes[node.chlds[k]]);
         }

         return node.idshift + 1;
      }

      scan_func(this.nodes, this.nodes[0]);
   }

   /** @summary Extract only visibility flags
     * @desc Used to transfer them to the worker */
   ClonedNodes.prototype.getVisibleFlags = function() {
      let res = new Array(this.nodes.length);
      for (let n=0;n<this.nodes.length;++n)
         res[n] = { vis: this.nodes[n].vis, nochlds: this.nodes[n].nochlds };
      return res;
   }

   /** @summary Assign only visibility flags, extracted with getVisibleFlags */
   ClonedNodes.prototype.setVisibleFlags = function(flags) {
      if (!this.nodes || !flags || !flags.length != this.nodes.length)
         return 0;

      let res = 0;
      for (let n=0;n<this.nodes.length;++n) {
         let clone = this.nodes[n];

         clone.vis = flags[n].vis;
         clone.nochlds = flags[n].nochlds;
         if (clone.vis) res++;
      }

      return res;
   }

   /** @summary Scan visible nodes in hierarchy, starting from nodeid
     * @desc Each entry in hierarchy get its unique id, which is not changed with visibility flags */
   ClonedNodes.prototype.scanVisible = function(arg, vislvl) {

      if (!this.nodes) return 0;

      if (vislvl === undefined) {
         if (!arg) arg = {};

         vislvl = arg.vislvl || this.vislevel || 4; // default 3 in ROOT
         if (vislvl > 88) vislvl = 88;

         arg.stack = new Array(100); // current stack
         arg.nodeid = 0;
         arg.counter = 0; // sequence ID of the node, used to identify it later
         arg.last = 0;
         arg.CopyStack = function(factor) {
            let entry = { nodeid: this.nodeid, seqid: this.counter, stack: new Array(this.last) };
            if (factor) entry.factor = factor; // factor used to indicate importance of entry, will be build as first
            for (let n=0;n<this.last;++n) entry.stack[n] = this.stack[n+1]; // copy stack
            return entry;
         }

         if (arg.domatrix) {
            arg.matrices = [];
            arg.mpool = [ new THREE.Matrix4() ]; // pool of Matrix objects to avoid permanent creation
            arg.getmatrix = function() { return this.matrices[this.last]; }
         }
      }

      let res = 0, node = this.nodes[arg.nodeid];

      if (arg.domatrix) {
         if (!arg.mpool[arg.last+1])
            arg.mpool[arg.last+1] = new THREE.Matrix4();

         let prnt = (arg.last > 0) ? arg.matrices[arg.last-1] : new THREE.Matrix4();
         if (node.matrix) {
            arg.matrices[arg.last] = arg.mpool[arg.last].fromArray(prnt.elements);
            arg.matrices[arg.last].multiply(arg.mpool[arg.last+1].fromArray(node.matrix));
         } else {
            arg.matrices[arg.last] = prnt;
         }
      }

      if (node.nochlds) vislvl = 0;

      if (node.vis > vislvl) {
         if (!arg.func || arg.func(node)) res++;
      }

      arg.counter++;

      if ((vislvl > 0) && node.chlds) {
         arg.last++;
         for (let i = 0; i < node.chlds.length; ++i) {
            arg.nodeid = node.chlds[i];
            arg.stack[arg.last] = i; // in the stack one store index of child, it is path in the hierarchy
            res += this.scanVisible(arg, vislvl-1);
         }
         arg.last--;
      } else {
         arg.counter += (node.idshift || 0);
      }

      if (arg.last === 0) {
         delete arg.last;
         delete arg.stack;
         delete arg.CopyStack;
         delete arg.counter;
         delete arg.matrices;
         delete arg.mpool;
         delete arg.getmatrix;
      }

      return res;
   }

   /** @summary Return node name with given id.
    * @desc Either original object or description is used */
   ClonedNodes.prototype.getNodeName = function(nodeid) {
      if (this.origin) {
         let obj = this.origin[nodeid];
         return obj ? geo.getObjectName(obj) : "";
      }
      let node = this.nodes[nodeid];
      return node ? node.name : "";
   }

   /** @summary Returns description for provide stack */
   ClonedNodes.prototype.resolveStack = function(stack, withmatrix) {

      let res = { id: 0, obj: null, node: this.nodes[0], name: this.name_prefix };

      // if (!this.toplevel || (this.nodes.length === 1) || (res.node.kind === 1)) res.name = "";

      if (withmatrix) {
         res.matrix = new THREE.Matrix4();
         if (res.node.matrix) res.matrix.fromArray(res.node.matrix);
      }

      if (this.origin)
         res.obj = this.origin[0];

      //if (!res.name)
      //   res.name = this.getNodeName(0);

      if (stack)
         for(let lvl=0;lvl<stack.length;++lvl) {
            res.id = res.node.chlds[stack[lvl]];
            res.node = this.nodes[res.id];

            if (this.origin)
               res.obj = this.origin[res.id];

            let subname = this.getNodeName(res.id);
            if (subname) {
               if (res.name) res.name+="/";
               res.name += subname;
            }

            if (withmatrix && res.node.matrix)
               res.matrix.multiply(new THREE.Matrix4().fromArray(res.node.matrix));
         }

      return res;
   }

   /** @summary Create stack array based on nodes ids array.
    * @desc Ids list should correspond to existing nodes hierarchy */
   ClonedNodes.prototype.buildStackByIds = function(ids) {
      if (!ids) return null;

      if (ids[0] !== 0) {
         console.error('wrong ids - first should be 0');
         return null;
      }

      let node = this.nodes[0], stack = [];

      for (let k=1;k<ids.length;++k) {
         let nodeid = ids[k];
         if (!node) return null;
         let chindx = node.chlds.indexOf(nodeid);
         if (chindx < 0) {
            console.error('wrong nodes ids ' + ids[k] + ' is not child of ' + ids[k-1]);
            return null;
         }

         stack.push(chindx);
         node = this.nodes[nodeid];
      }

      return stack;
   }

   /** @summary Retuns ids array which correspond to the stack */
   ClonedNodes.prototype.buildIdsByStack = function(stack) {
      if (!stack) return null;
      let node = this.nodes[0], ids = [0];
      for (let k=0;k<stack.length;++k) {
         let id = node.chlds[stack[k]];
         ids.push(id);
         node = this.nodes[id];
      }
      return ids;
   }

   /** @summary Returns true if stack includes at any place provided nodeid */
   ClonedNodes.prototype.isIdInStack = function(nodeid, stack) {

      if (!nodeid) return true;

      let node = this.nodes[0], id = 0;

      for(let lvl = 0; lvl < stack.length; ++lvl) {
         id = node.chlds[stack[lvl]];
         if (id == nodeid) return true;
         node = this.nodes[id];
      }

      return false;
   }

   /** @summary Find stack by name which include names of all parents */
   ClonedNodes.prototype.findStackByName = function(fullname) {

      let names = fullname.split('/'), currid = 0, stack = [];

      if (this.getNodeName(currid) !== names[0]) return null;

      for (let n=1;n<names.length;++n) {
         let node = this.nodes[currid];
         if (!node.chlds) return null;

         for (let k=0;k<node.chlds.length;++k) {
            let chldid = node.chlds[k];
            if (this.getNodeName(chldid) === names[n]) { stack.push(k); currid = chldid; break; }
         }

         // no new entry - not found stack
         if (stack.length === n - 1) return null;
      }

      return stack;
   }

   /** @summary Set usage of default ROOT colors */
   ClonedNodes.prototype.setDefaultColors = function(on) {
      this.use_dflt_colors = on;
      if (this.use_dflt_colors && !this.dflt_table) {

         let dflt = { kWhite:0,  kBlack:1, kGray:920,
               kRed:632, kGreen:416, kBlue:600, kYellow:400, kMagenta:616, kCyan:432,
               kOrange:800, kSpring:820, kTeal:840, kAzure:860, kViolet:880, kPink:900 };

         let nmax = 110, col = [];
         for (let i=0;i<nmax;i++) col.push(dflt.kGray);

         //  here we should create a new TColor with the same rgb as in the default
         //  ROOT colors used below
         col[ 3] = dflt.kYellow-10;
         col[ 4] = col[ 5] = dflt.kGreen-10;
         col[ 6] = col[ 7] = dflt.kBlue-7;
         col[ 8] = col[ 9] = dflt.kMagenta-3;
         col[10] = col[11] = dflt.kRed-10;
         col[12] = dflt.kGray+1;
         col[13] = dflt.kBlue-10;
         col[14] = dflt.kOrange+7;
         col[16] = dflt.kYellow+1;
         col[20] = dflt.kYellow-10;
         col[24] = col[25] = col[26] = dflt.kBlue-8;
         col[29] = dflt.kOrange+9;
         col[79] = dflt.kOrange-2;

         this.dflt_table = col;
      }
   }

   /** @summary Provide different properties of draw entry nodeid
    * @desc Only if node visible, material will be created*/
   ClonedNodes.prototype.getDrawEntryProperties = function(entry) {

      let clone = this.nodes[entry.nodeid];
      let visible = true;

      if (clone.kind === kindShape) {
         let prop = { name: clone.name, nname: clone.name, shape: null, material: null, chlds: null };
         let _opacity = entry.opacity || 1;
         prop.fillcolor = new THREE.Color( entry.color ? "rgb(" + entry.color + ")" : "blue" );
         prop.material = new THREE.MeshLambertMaterial( { transparent: _opacity < 1,
                          opacity: _opacity, wireframe: false, color: prop.fillcolor,
                          side: THREE.FrontSide /* THREE.DoubleSide*/, vertexColors: THREE.NoColors /*THREE.VertexColors */,
                          depthWrite: _opacity == 1 } );
         prop.material.inherentOpacity = _opacity;

         return prop;
      }

      if (!this.origin) {
         console.error('origin not there - kind', clone.kind, entry.nodeid, clone);
         return null;
      }

      let node = this.origin[entry.nodeid];

      if (clone.kind === kindEve) {
         // special handling for EVE nodes

         let prop = { name: geo.getObjectName(node), nname: geo.getObjectName(node), shape: node.fShape, material: null, chlds: null };

         if (node.fElements !== null) prop.chlds = node.fElements.arr;

         if (visible) {
            let _opacity = Math.min(1, node.fRGBA[3]);
            prop.fillcolor = new THREE.Color( node.fRGBA[0], node.fRGBA[1], node.fRGBA[2] );
            prop.material = new THREE.MeshLambertMaterial( { transparent: _opacity < 1,
                             opacity: _opacity, wireframe: false, color: prop.fillcolor,
                             side: THREE.FrontSide /* THREE.DoubleSide*/, vertexColors: THREE.NoColors /*THREE.VertexColors */,
                             depthWrite: _opacity == 1 } );
            prop.material.inherentOpacity = _opacity;
         }

         return prop;
      }

      let volume = node.fVolume;

      let prop = { name: geo.getObjectName(volume), nname: geo.getObjectName(node), volume: node.fVolume, shape: volume.fShape, material: null, chlds: null };

      if (node.fVolume.fNodes !== null) prop.chlds = node.fVolume.fNodes.arr;

      if (volume) prop.linewidth = volume.fLineWidth;

      if (visible) {

         // TODO: maybe correctly extract ROOT colors here?
         let _opacity = 1.0, jsrp = JSROOT.Painter,
             root_colors = jsrp ? jsrp.root_colors : ['white', 'black', 'red', 'green', 'blue', 'yellow', 'magenta', 'cyan'];

         if (entry.custom_color)
            prop.fillcolor = entry.custom_color;
         else if ((volume.fFillColor > 1) && (volume.fLineColor == 1))
            prop.fillcolor = root_colors[volume.fFillColor];
         else if (volume.fLineColor >= 0)
            prop.fillcolor = root_colors[volume.fLineColor];

         if (volume.fMedium && volume.fMedium.fMaterial) {
            let mat = volume.fMedium.fMaterial,
                fillstyle = mat.fFillStyle,
                transparency = (fillstyle < 3000 || fillstyle > 3100) ? 0 : fillstyle - 3000;

            if (this.use_dflt_colors) {
               let matZ = Math.round(mat.fZ),
                   icol = this.dflt_table[matZ];
               prop.fillcolor = root_colors[icol];
               if (mat.fDensity < 0.1) transparency = 60;
            }

            if (transparency > 0)
               _opacity = (100.0 - transparency) / 100.0;
            if (prop.fillcolor === undefined)
               prop.fillcolor = root_colors[mat.fFillColor];
         }
         if (prop.fillcolor === undefined)
            prop.fillcolor = "lightgrey";

         prop.material = new THREE.MeshLambertMaterial( { transparent: _opacity < 1,
                              opacity: _opacity, wireframe: false, color: prop.fillcolor,
                              side: THREE.FrontSide /* THREE.DoubleSide */, vertexColors: THREE.NoColors /*THREE.VertexColors*/,
                              depthWrite: _opacity == 1 } );
         prop.material.inherentOpacity = _opacity;

      }

      return prop;
   }

   /** @summary Creates hierarchy of Object3D for given stack entry
     * @desc Such hierarchy repeats hierarchy of TGeoNodes and set matrix for the objects drawing
     * also set renderOrder, required to handle transparency */
   ClonedNodes.prototype.createObject3D = function(stack, toplevel, options) {

      let node = this.nodes[0], three_prnt = toplevel, draw_depth = 0,
          force = (typeof options == 'object') || (options==='force');

      for(let lvl=0; lvl<=stack.length; ++lvl) {
         let nchld = (lvl > 0) ? stack[lvl-1] : 0;
         // extract current node
         if (lvl>0)  node = this.nodes[node.chlds[nchld]];

         let obj3d = undefined;

         if (three_prnt.children)
            for (let i = 0; i < three_prnt.children.length; ++i) {
               if (three_prnt.children[i].nchld === nchld) {
                  obj3d = three_prnt.children[i];
                  break;
               }
            }

         if (obj3d) {
            three_prnt = obj3d;
            if (obj3d.$jsroot_drawable) draw_depth++;
            continue;
         }

         if (!force) return null;

         obj3d = new THREE.Object3D();

         if (node.abs_matrix) {
            obj3d.absMatrix = new THREE.Matrix4();
            obj3d.absMatrix.fromArray(node.matrix);
         } else if (node.matrix) {
            obj3d.matrix.fromArray(node.matrix);
            obj3d.matrix.decompose( obj3d.position, obj3d.quaternion, obj3d.scale );
         }

         // this.accountNodes(obj3d);
         obj3d.nchld = nchld; // mark index to find it again later

         // add the mesh to the scene
         three_prnt.add(obj3d);

         // this is only for debugging - test inversion of whole geometry
         if ((lvl==0) && (typeof options == 'object') && options.scale) {
            if ((options.scale.x < 0) || (options.scale.y < 0) || (options.scale.z < 0)) {
               obj3d.scale.copy(options.scale);
               obj3d.updateMatrix();
            }
         }

         obj3d.updateMatrixWorld();

         three_prnt = obj3d;
      }

      if ((options === 'mesh') || (options === 'delete_mesh')) {
         let mesh = null;
         if (three_prnt)
            for (let n=0; (n<three_prnt.children.length) && !mesh;++n) {
               let chld = three_prnt.children[n];
               if ((chld.type === 'Mesh') && (chld.nchld === undefined)) mesh = chld;
            }

         if ((options === 'mesh') || !mesh) return mesh;

         let res = three_prnt;
         while (mesh && (mesh !== toplevel)) {
            three_prnt = mesh.parent;
            three_prnt.remove(mesh);
            mesh = (three_prnt.children.length == 0) ? three_prnt : null;
         }

         return res;
      }

      if (three_prnt) {
         three_prnt.$jsroot_drawable = true;
         three_prnt.$jsroot_depth = draw_depth;
      }

      return three_prnt;
   }

   /** @summary Get volume boundary
     * @private */
   ClonedNodes.prototype.getVolumeBoundary = function(viscnt, facelimit, nodeslimit) {

      let result = { min: 0, max: 1, sortidcut: 0 };

      if (!this.sortmap) {
         console.error('sorting map do not exist');
         return result;
      }

      let maxNode, currNode, cnt=0, facecnt=0;

      for (let n = 0; (n < this.sortmap.length) && (cnt < nodeslimit) && (facecnt < facelimit); ++n) {
         let id = this.sortmap[n];
         if (viscnt[id] === 0) continue;
         currNode = this.nodes[id];
         if (!maxNode) maxNode = currNode;
         cnt += viscnt[id];
         facecnt += viscnt[id] * currNode.nfaces;
      }

      if (!currNode) {
         console.error('no volumes selected');
         return result;
      }

      // console.log('Volume boundary ' + currNode.vol + '  cnt ' + cnt + '  faces ' + facecnt);
      result.max = maxNode.vol;
      result.min = currNode.vol;
      result.sortidcut = currNode.sortid; // latest node is not included
      return result;
   }

   /** @summary Collects visible nodes, using maxlimit
     * @desc One can use map to define cut based on the volume or serious of cuts */
   ClonedNodes.prototype.collectVisibles = function(maxnumfaces, frustum) {

      // in simple case shape as it is
      if (this.plain_shape)
         return { lst: [ { nodeid: 0, seqid: 0, stack: [], factor: 1, shapeid: 0, server_shape: this.plain_shape } ], complete: true };

      let arg = {
         facecnt: 0,
         viscnt: new Array(this.nodes.length), // counter for each node
         vislvl: this.getVisLevel(),
         reset: function() {
            this.total = 0;
            this.facecnt = 0;
            for (let n=0;n<this.viscnt.length;++n) this.viscnt[n] = 0;
         },
         // nodes: this.nodes,
         func: function(node) {
            this.total++;
            this.facecnt += node.nfaces;
            this.viscnt[node.id]++;
            return true;
         }
      };

      arg.reset();

      let total = this.scanVisible(arg),
          maxnumnodes = this.getMaxVisNodes();

      if (maxnumnodes > 0) {
         while ((total > maxnumnodes) && (arg.vislvl > 1)) {
            arg.vislvl--;
            arg.reset();
            total = this.scanVisible(arg);
         }
      }

      this.actual_level = arg.vislvl; // not used, can be shown somewhere in the gui

      let minVol = 0, maxVol = 0, camVol = -1, camFact = 10, sortidcut = this.nodes.length + 1;

      console.log('Total visible nodes ' + total + ' numfaces ' + arg.facecnt);

      if (arg.facecnt > maxnumfaces) {

         let bignumfaces = maxnumfaces * (frustum ? 0.8 : 1.0),
             bignumnodes = maxnumnodes * (frustum ? 0.8 : 1.0);

         // define minimal volume, which always to shown
         let boundary = this.getVolumeBoundary(arg.viscnt, bignumfaces, bignumnodes);

         minVol = boundary.min;
         maxVol = boundary.max;
         sortidcut = boundary.sortidcut;

         if (frustum) {
             arg.domatrix = true;
             arg.frustum = frustum;
             arg.totalcam = 0;
             arg.func = function(node) {
                if (node.vol <= minVol) // only small volumes are interesting
                   if (this.frustum.CheckShape(this.getmatrix(), node)) {
                      this.viscnt[node.id]++;
                      this.totalcam += node.nfaces;
                   }

                return true;
             }

             for (let n=0;n<arg.viscnt.length;++n) arg.viscnt[n] = 0;

             this.scanVisible(arg);

             if (arg.totalcam > maxnumfaces*0.2)
                camVol = this.getVolumeBoundary(arg.viscnt, maxnumfaces*0.2, maxnumnodes*0.2).min;
             else
                camVol = 0;

             camFact = maxVol / ((camVol>0) ? (camVol>0) : minVol);

             // console.log('Limit for camera ' + camVol + '  faces in camera view ' + arg.totalcam);
         }
      }

      arg.items = [];

      arg.func = function(node) {
         if (node.sortid < sortidcut) {
            this.items.push(this.CopyStack());
         } else if ((camVol >= 0) && (node.vol > camVol)) {
            if (this.frustum.CheckShape(this.getmatrix(), node))
               this.items.push(this.CopyStack(camFact));
         }
         return true;
      }

      this.scanVisible(arg);

      return { lst: arg.items, complete: minVol === 0 };
   }

   /** @summary Merge list of drawn objects
     * @desc In current list we should mark if object already exists
     * from previous list we should collect objects which are not there */
   ClonedNodes.prototype.mergeVisibles = function(current, prev) {

      let indx2 = 0, del = [];
      for (let indx1=0; (indx1<current.length) && (indx2<prev.length); ++indx1) {

         while ((indx2 < prev.length) && (prev[indx2].seqid < current[indx1].seqid)) {
            del.push(prev[indx2++]); // this entry should be removed
         }

         if ((indx2 < prev.length) && (prev[indx2].seqid === current[indx1].seqid)) {
            if (prev[indx2].done) current[indx1].done = true; // copy ready flag
            indx2++;
         }
      }

      // remove rest
      while (indx2<prev.length)
         del.push(prev[indx2++]);

      return del; //
   }

   /** @summary Collect all uniques shapes which should be build
    *  @desc Check if same shape used many time for drawing */
   ClonedNodes.prototype.collectShapes = function(lst) {

      // nothing else - just that single shape
      if (this.plain_shape)
         return [ this.plain_shape ];

      let shapes = [];

      for (let i=0;i<lst.length;++i) {
         let entry = lst[i];
         let shape = this.getNodeShape(entry.nodeid);

         if (!shape) continue; // strange, but avoid misleading

         if (shape._id === undefined) {
            shape._id = shapes.length;

            shapes.push({ id: shape._id, shape: shape, vol: this.nodes[entry.nodeid].vol, refcnt: 1, factor: 1, ready: false });

            // shapes.push( { obj: shape, vol: this.nodes[entry.nodeid].vol });
         } else {
            shapes[shape._id].refcnt++;
         }

         entry.shape = shapes[shape._id]; // remember shape used

         // use maximal importance factor to push element to the front
         if (entry.factor && (entry.factor>entry.shape.factor))
            entry.shape.factor = entry.factor;
      }

      // now sort shapes in volume decrease order
      shapes.sort(function(a,b) { return b.vol*b.factor - a.vol*a.factor; })

      // now set new shape ids according to the sorted order and delete temporary field
      for (let n=0;n<shapes.length;++n) {
         let item = shapes[n];
         item.id = n; // set new ID
         delete item.shape._id; // remove temporary field
      }

      // as last action set current shape id to each entry
      for (let i=0;i<lst.length;++i) {
         let entry = lst[i];
         if (entry.shape) {
            entry.shapeid = entry.shape.id; // keep only id for the entry
            delete entry.shape; // remove direct references
         }
      }

      return shapes;
   }

   /** @summary Merge shape lists */
   ClonedNodes.prototype.mergeShapesLists = function(oldlst, newlst) {

      if (!oldlst) return newlst;

      // set geometry to shape object itself
      for (let n=0;n<oldlst.length;++n) {
         let item = oldlst[n];

         item.shape._geom = item.geom;
         delete item.geom;

         if (item.geomZ!==undefined) {
            item.shape._geomZ = item.geomZ;
            delete item.geomZ;
         }
      }

      // take from shape (if match)
      for (let n=0;n<newlst.length;++n) {
         let item = newlst[n];

         if (item.shape._geom !== undefined) {
            item.geom = item.shape._geom;
            delete item.shape._geom;
         }

         if (item.shape._geomZ !== undefined) {
            item.geomZ = item.shape._geomZ;
            delete item.shape._geomZ;
         }
      }

      // now delete all unused geometries
      for (let n=0;n<oldlst.length;++n) {
         let item = oldlst[n];
         delete item.shape._geom;
         delete item.shape._geomZ;
      }

      return newlst;
   }

   /** @summary Build shapes */
   ClonedNodes.prototype.buildShapes = function(lst, limit, timelimit) {

      let created = 0,
          tm1 = new Date().getTime(),
          res = { done: false, shapes: 0, faces: 0, notusedshapes: 0 };

      for (let n=0;n<lst.length;++n) {
         let item = lst[n];

         // if enough faces are produced, nothing else is required
         if (res.done) { item.ready = true; continue; }

         if (!item.ready) {
            item._typename = "$$Shape$$"; // let reuse item for direct drawing
            item.ready = true;
            if (item.geom === undefined) {
               item.geom = createGeometry(item.shape);
               if (item.geom) created++; // indicate that at least one shape was created
            }
            item.nfaces = countGeometryFaces(item.geom);
         }

         res.shapes++;
         if (!item.used) res.notusedshapes++;
         res.faces += item.nfaces*item.refcnt;

         if (res.faces >= limit) {
            res.done = true;
         } else if ((created > 0.01*lst.length) && (timelimit!==undefined)) {
            let tm2 = new Date().getTime();
            if (tm2-tm1 > timelimit) return res;
         }
      }

      res.done = true;

      return res;
   }

   /// =====================================================================

   /** @summary Returns geo object name
     * @desc Can appens some special suffixes
     * @private */
   geo.getObjectName = function(obj) {
      if (!obj || !obj.fName) return "";
      return obj.fName + (obj.$geo_suffix ? obj.$geo_suffix : "");
   }

   /** @summary Check duplicates
     * @private */
   geo.checkDuplicates = function(parent, chlds) {
      if (parent) {
         if (parent.$geo_checked) return;
         parent.$geo_checked = true;
      }

      let names = [], cnts = [];
      for (let k=0;k<chlds.length;++k) {
         let chld = chlds[k];
         if (!chld || !chld.fName) continue;
         if (!chld.$geo_suffix) {
            let indx = names.indexOf(chld.fName);
            if (indx>=0) {
               let cnt = cnts[indx] || 1;
               while(names.indexOf(chld.fName+"#"+cnt)>=0) ++cnt;
               chld.$geo_suffix = "#" + cnt;
               cnts[indx] = cnt+1;
            }
         }
         names.push(geo.getObjectName(chld));
      }
   }

  /** @summary Create flipped mesh for the shape
    * @desc When transformation matrix includes one or several inversion of axis,
    * one should inverse geometry object, otherwise THREE.js cannot correctly draw it
    * @param {Object} shape - TGeoShape object
    * @param {Object} material - material
    * @memberof SJROOT.GEO
    * @private */

   function createFlippedMesh(shape, material) {

      let flip =  new THREE.Vector3(1,1,-1);

      if (shape.geomZ === undefined) {

         if (shape.geom.type == 'BufferGeometry') {

            let pos = shape.geom.getAttribute('position').array,
                norm = shape.geom.getAttribute('normal').array,
                index = shape.geom.getIndex();

            if (index) {
               // we need to unfold all points to
               let arr = index.array,
                   i0 = shape.geom.drawRange.start,
                   ilen = shape.geom.drawRange.count;
               if (i0 + ilen > arr.length) ilen = arr.length - i0;

               let dpos = new Float32Array(ilen*3), dnorm = new Float32Array(ilen*3);
               for (let ii = 0; ii < ilen; ++ii) {
                  let k = arr[i0 + ii];
                  if ((k<0) || (k*3>=pos.length)) console.log('strange index', k*3, pos.length);
                  dpos[ii*3] = pos[k*3];
                  dpos[ii*3+1] = pos[k*3+1];
                  dpos[ii*3+2] = pos[k*3+2];
                  dnorm[ii*3] = norm[k*3];
                  dnorm[ii*3+1] = norm[k*3+1];
                  dnorm[ii*3+2] = norm[k*3+2];
               }

               pos = dpos; norm = dnorm;
            }

            let len = pos.length, n, shift = 0,
                newpos = new Float32Array(len),
                newnorm = new Float32Array(len);

            // we should swap second and third point in each face
            for (n=0; n<len; n+=3) {
               newpos[n]   = pos[n+shift];
               newpos[n+1] = pos[n+1+shift];
               newpos[n+2] = -pos[n+2+shift];

               newnorm[n]   = norm[n+shift];
               newnorm[n+1] = norm[n+1+shift];
               newnorm[n+2] = -norm[n+2+shift];

               shift+=3; if (shift===6) shift=-3; // values 0,3,-3
            }

            shape.geomZ = new THREE.BufferGeometry();
            shape.geomZ.setAttribute( 'position', new THREE.BufferAttribute( newpos, 3 ) );
            shape.geomZ.setAttribute( 'normal', new THREE.BufferAttribute( newnorm, 3 ) );
            // normals are calculated with normal geometry and correctly scaled
            // geom.computeVertexNormals();

         } else {

            shape.geomZ = shape.geom.clone();

            shape.geomZ.scale(flip.x, flip.y, flip.z);

            let face, d, n = 0;
            while(n < shape.geomZ.faces.length) {
               face = geom.faces[n++];
               d = face.b; face.b = face.c; face.c = d;
            }

            // normals are calculated with normal geometry and correctly scaled
            // geom.computeFaceNormals();
         }
      }

      let mesh = new THREE.Mesh( shape.geomZ, material );
      mesh.scale.copy(flip);
      mesh.updateMatrix();

      mesh._flippedMesh = true;

      return mesh;
   }

   /** @summary extract code of Box3.expandByObject
     * @desc Major difference - do not traverse hierarchy
     * @memberof JSROOT.GEO
     * @private */
   function getBoundingBox(node, box3, local_coordinates) {
      if (!node || !node.geometry) return box3;

      if (!box3) { box3 = new THREE.Box3(); box3.makeEmpty(); }

      if (!local_coordinates) node.updateMatrixWorld();

      let v1 = new THREE.Vector3(),
          geometry = node.geometry;

      if ( geometry.isGeometry ) {
         let vertices = geometry.vertices;
         for (let i = 0, l = vertices.length; i < l; i ++ ) {
            v1.copy( vertices[ i ] );
            if (!local_coordinates) v1.applyMatrix4( node.matrixWorld );
            box3.expandByPoint( v1 );
         }
      } else if ( geometry.isBufferGeometry ) {
         let attribute = geometry.attributes.position;
         if ( attribute !== undefined ) {
            for (let i = 0, l = attribute.count; i < l; i ++ ) {
               // v1.fromAttribute( attribute, i ).applyMatrix4( node.matrixWorld );
               v1.fromBufferAttribute( attribute, i );
               if (!local_coordinates) v1.applyMatrix4( node.matrixWorld );
               box3.expandByPoint( v1 );
            }
         }
      }

      return box3;
   }

   /** @summary Cleanup shape entity
    * @private */
   geo.cleanupShape = function(shape) {
      if (!shape) return;

      if (shape.geom && (typeof shape.geom.dispose == 'funciton'))
         shape.geom.dispose();

      if (shape.geomZ && (typeof shape.geomZ.dispose == 'funciton'))
         shape.geomZ.dispose();

      delete shape.geom;
      delete shape.geomZ;
   }

   /** @summary Set rendering order for created hierarchy
    * @desc depending from provided method sort differently objects
    * @param toplevel - top element
    * @param origin - camera position used to provide sorting
    * @param method - name of sorting method like "pnt", "ray", "size", "dflt"  */
   geo.produceRenderOrder = function(toplevel, origin, method, clones) {

      let raycast = new THREE.Raycaster();

      function setdefaults(top) {
         if (!top) return;
         top.traverse(function(obj) {
            obj.renderOrder = 0;
            if (obj.material) obj.material.depthWrite = true; // by default depthWriting enabled
         });
      }

      function traverse(obj, lvl, arr) {
         // traverse hierarchy and extract all children of given level
         // if (obj.$jsroot_depth===undefined) return;

         if (!obj.children) return;

         for (let k=0;k<obj.children.length;++k) {
            let chld = obj.children[k];
            if (chld.$jsroot_order === lvl) {
               if (chld.material) {
                  if (chld.material.transparent) {
                     chld.material.depthWrite = false; // disable depth writing for transparent
                     arr.push(chld);
                  } else {
                     setdefaults(chld);
                  }
               }
            } else if ((obj.$jsroot_depth===undefined) || (obj.$jsroot_depth < lvl)) {
               traverse(chld, lvl, arr);
            }
         }
      }

      function sort(arr, minorder, maxorder) {
         // resort meshes using ray caster and camera position
         // idea to identify meshes which are in front or behind

         if (arr.length > 300) {
            // too many of them, just set basic level and exit
            for (let i=0;i<arr.length;++i) arr[i].renderOrder = (minorder + maxorder)/2;
            return false;
         }

         let tmp_vect = new THREE.Vector3();

         // first calculate distance to the camera
         // it gives preliminary order of volumes

         for (let i=0;i<arr.length;++i) {
            let mesh = arr[i],
                box3 = mesh.$jsroot_box3;

            if (!box3)
               mesh.$jsroot_box3 = box3 = getBoundingBox(mesh);

            if (method === 'size') {
               mesh.$jsroot_distance = box3.getSize(new THREE.Vector3());
               continue;
            }

            if (method === "pnt") {
               mesh.$jsroot_distance = origin.distanceTo(box3.getCenter(tmp_vect));
               continue;
            }

            let dist = Math.min(origin.distanceTo(box3.min), origin.distanceTo(box3.max));

            let pnt = new THREE.Vector3(box3.min.x, box3.min.y, box3.max.z);
            dist = Math.min(dist, origin.distanceTo(pnt));
            pnt.set(box3.min.x, box3.max.y, box3.min.z)
            dist = Math.min(dist, origin.distanceTo(pnt));
            pnt.set(box3.max.x, box3.min.y, box3.min.z)
            dist = Math.min(dist, origin.distanceTo(pnt));

            pnt.set(box3.max.x, box3.max.y, box3.min.z)
            dist = Math.min(dist, origin.distanceTo(pnt));

            pnt.set(box3.max.x, box3.min.y, box3.max.z)
            dist = Math.min(dist, origin.distanceTo(pnt));

            pnt.set(box3.min.x, box3.max.y, box3.max.z)
            dist = Math.min(dist, origin.distanceTo(pnt));

            mesh.$jsroot_distance = dist;
         }

         arr.sort(function(a,b) { return a.$jsroot_distance - b.$jsroot_distance; });

         let resort = new Array(arr.length);

         for (let i=0;i<arr.length;++i) {
            arr[i].$jsroot_index = i;
            resort[i] = arr[i];
         }

         if (method==="ray")
         for (let i=arr.length-1;i>=0;--i) {
            let mesh = arr[i],
                box3 = mesh.$jsroot_box3,
                direction = box3.getCenter(tmp_vect);

            for(let ntry=0; ntry<2;++ntry) {

               direction.sub(origin).normalize();

               raycast.set( origin, direction );

               let intersects = raycast.intersectObjects(arr, false); // only plain array

               let unique = [];

               for (let k1=0;k1<intersects.length;++k1) {
                  if (unique.indexOf(intersects[k1].object)<0) unique.push(intersects[k1].object);
                  // if (intersects[k1].object === mesh) break; // trace until object itself
               }

               intersects = unique;

               if ((intersects.indexOf(mesh)<0) && (ntry>0))
                  console.log('MISS', clones ? clones.resolveStack(mesh.stack).name : "???");

               if ((intersects.indexOf(mesh)>=0) || (ntry>0)) break;

               let pos = mesh.geometry.attributes.position.array;

               direction = new THREE.Vector3((pos[0]+pos[3]+pos[6])/3, (pos[1]+pos[4]+pos[7])/3, (pos[2]+pos[5]+pos[8])/3);

               direction.applyMatrix4(mesh.matrixWorld);
            }

            // now push first object in intersects to the front
            for (let k1=0;k1<intersects.length-1;++k1) {
               let mesh1 = intersects[k1], mesh2 = intersects[k1+1],
                   i1 = mesh1.$jsroot_index, i2 = mesh2.$jsroot_index;
               if (i1<i2) continue;
               for (let ii=i2;ii<i1;++ii) {
                  resort[ii] = resort[ii+1];
                  resort[ii].$jsroot_index = ii;
               }
               resort[i1] = mesh2;
               mesh2.$jsroot_index = i1;
            }

         }

         for (let i=0;i<resort.length;++i) {
            resort[i].renderOrder = maxorder - (i+1) / (resort.length+1) * (maxorder-minorder);
            delete resort[i].$jsroot_index;
            delete resort[i].$jsroot_distance;
         }

         return true;
      }

      function process(obj, lvl, minorder, maxorder) {
         let arr = [], did_sort = false;

         traverse(obj, lvl, arr);

         if (!arr.length) return;

         if (minorder === maxorder) {
            for (let k=0;k<arr.length;++k)
               arr[k].renderOrder = minorder;
         } else {
           did_sort = sort(arr, minorder, maxorder);
           if (!did_sort) minorder = maxorder = (minorder + maxorder) / 2;
         }

         for (let k=0;k<arr.length;++k) {
            let next = arr[k].parent, min = minorder, max = maxorder;

            if (did_sort) {
               max = arr[k].renderOrder;
               min = max - (maxorder - minorder) / (arr.length + 2);
            }

            process(next, lvl+1, min, max);
         }
      }

      if (!method || (method==="dflt"))
         setdefaults(toplevel);
      else
         process(toplevel, 0, 1, 1000000);
   }

   /** @summary Build three.js model for given geometry object
     * @param {Object} obj - TGeo-related object
     * @param {Object} [opt] - options
     * @param {Number} [opt.vislevel] - visibility level like TGeoManager, when not specified - show all
     * @param {Number} [opt.numnodes=1000] - maximal number of visible nodes
     * @param {Number} [opt.numfaces=100000] - approx maximal number of created triangles
     * @param {boolean} [opt.doubleside=false] - use double-side material
     * @param {boolean} [opt.wireframe=false] - show wireframe for created shapes
     * @param {boolean} [opt.dflt_colors=false] - use default ROOT colors
     * @returns {object} THREE.Object3D with created model
     * @example
     * JSROOT.require('geom')
     *       .then(geo => {
     *           let obj3d = geo.build(obj);
     *           // this is three.js object and can be now inserted in the scene
     *        });
     */
   geo.build = function(obj, opt) {

      if (!obj) return null;

      if (!opt) opt = {};
      if (!opt.numfaces) opt.numfaces = 100000;
      if (!opt.numnodes) opt.numnodes = 1000;
      if (!opt.frustum) opt.frustum = null;

      opt.res_mesh = opt.res_faces = 0;

      let shape = null, hide_top = false;

      if (('fShapeBits' in obj) && ('fShapeId' in obj)) {
         shape = obj; obj = null;
      } else if ((obj._typename === 'TGeoVolumeAssembly') || (obj._typename === 'TGeoVolume')) {
         shape = obj.fShape;
      } else if ((obj._typename === "TEveGeoShapeExtract") || (obj._typename === "ROOT::Experimental::REveGeoShapeExtract")) {
         shape = obj.fShape;
      } else if (obj._typename === 'TGeoManager') {
         obj = obj.fMasterVolume;
         hide_top = !opt.showtop;
         shape = obj.fShape;
      } else if (obj.fVolume) {
         shape = obj.fVolume.fShape;
      } else {
         obj = null;
      }

      if (opt.composite && shape && (shape._typename == 'TGeoCompositeShape') && shape.fNode)
         obj = geo.buildCompositeVolume(shape);

      if (!obj && shape)
         obj = JSROOT.extend(JSROOT.create("TEveGeoShapeExtract"),
                   { fTrans: null, fShape: shape, fRGBA: [0, 1, 0, 1], fElements: null, fRnrSelf: true });

      if (!obj) return null;

      if (obj._typename.indexOf('TGeoVolume') === 0)
         obj = { _typename:"TGeoNode", fVolume: obj, fName: obj.fName, $geoh: obj.$geoh, _proxy: true };

      let clones = new ClonedNodes(obj);
      clones.setVisLevel(opt.vislevel);
      clones.setMaxVisNodes(opt.numnodes);

      if (opt.dflt_colors)
         clones.setDefaultColors(true);

      let uniquevis = opt.no_screen ? 0 : clones.markVisibles(true);
      if (uniquevis <= 0)
         uniquevis = clones.markVisibles(false, false, hide_top);
      else
         uniquevis = clones.markVisibles(true, true, hide_top); // copy bits once and use normal visibility bits

      clones.produceIdShifts();

      // collect visible nodes
      let res = clones.collectVisibles(opt.numfaces, opt.frustum);

      let draw_nodes = res.lst;

      // collect shapes
      let shapes = clones.collectShapes(draw_nodes);

      clones.buildShapes(shapes, opt.numfaces);

      let toplevel = new THREE.Object3D();

      for (let n=0; n < draw_nodes.length;++n) {
         let entry = draw_nodes[n];
         if (entry.done) continue;

         let shape = shapes[entry.shapeid];
         if (!shape.ready) {
            console.warn('shape marked as not ready when should');
            break;
         }
         entry.done = true;
         shape.used = true; // indicate that shape was used in building

         if (!shape.geom || (shape.nfaces === 0)) {
            // node is visible, but shape does not created
            clones.createObject3D(entry.stack, toplevel, 'delete_mesh');
            continue;
         }

         let prop = clones.getDrawEntryProperties(entry);

         opt.res_mesh++;
         opt.res_faces += shape.nfaces;

         let obj3d = clones.createObject3D(entry.stack, toplevel, opt);

         prop.material.wireframe = opt.wireframe;

         prop.material.side = opt.doubleside ? THREE.DoubleSide : THREE.FrontSide;

         let mesh = null, matrix = obj3d.absMatrix || obj3d.matrixWorld;

         if (matrix.determinant() > -0.9) {
            mesh = new THREE.Mesh(shape.geom, prop.material);
         } else {
            mesh = createFlippedMesh(shape, prop.material);
         }

         obj3d.add(mesh);

         if (obj3d.absMatrix) {
            mesh.matrix.copy(obj3d.absMatrix);
            mesh.matrix.decompose( obj3d.position, obj3d.quaternion, obj3d.scale );
            mesh.updateMatrixWorld();
         }

         // specify rendering order, required for transparency handling
         //if (obj3d.$jsroot_depth !== undefined)
         //   mesh.renderOrder = clones.maxdepth - obj3d.$jsroot_depth;
         //else
         //   mesh.renderOrder = clones.maxdepth - entry.stack.length;
      }

      return toplevel;
   }

   geo.projectGeometry = projectGeometry;
   geo.countGeometryFaces = countGeometryFaces;
   geo.createGeometry = createGeometry;
   geo.createProjectionMatrix = createProjectionMatrix;
   geo.createFrustum = createFrustum;
   geo.createFlippedMesh = createFlippedMesh;
   geo.getBoundingBox = getBoundingBox;
   geo.provideObjectInfo = provideObjectInfo;

   geo.ClonedNodes = ClonedNodes;
   JSROOT.GEO = geo;

   if (JSROOT.nodejs) module.exports = geo;

   return geo;
});
