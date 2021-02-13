/// @file JSRoot.csg.js
/// CSG library for THREE.js

JSROOT.define(['three'], function(THREE) {

   "use strict";

   const EPSILON = 1e-5,
         COPLANAR = 0,
         FRONT = 1,
         BACK = 2,
         SPANNING = 3;

   function Vertex(x, y, z, nx, ny, nz) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.nx = nx;
      this.ny = ny;
      this.nz = nz;
   }

   Vertex.prototype.setnormal = function ( nx, ny, nz ) {
      this.nx = nx;
      this.ny = ny;
      this.nz = nz;
   }

   Vertex.prototype.clone = function() {
      return new Vertex( this.x, this.y, this.z, this.nx, this.ny, this.nz);
   }

   Vertex.prototype.add = function( vertex ) {
      this.x += vertex.x;
      this.y += vertex.y;
      this.z += vertex.z;
      return this;
   }

   Vertex.prototype.subtract = function( vertex ) {
      this.x -= vertex.x;
      this.y -= vertex.y;
      this.z -= vertex.z;
      return this;
   }

   Vertex.prototype.multiplyScalar = function( scalar ) {
      this.x *= scalar;
      this.y *= scalar;
      this.z *= scalar;
      return this;
   }

   Vertex.prototype.cross = function( vertex ) {
      let x = this.x,
          y = this.y,
          z = this.z;

      this.x = y * vertex.z - z * vertex.y;
      this.y = z * vertex.x - x * vertex.z;
      this.z = x * vertex.y - y * vertex.x;

      return this;
   }

   Vertex.prototype.normalize = function() {
      let length = Math.sqrt( this.x * this.x + this.y * this.y + this.z * this.z );

      this.x /= length;
      this.y /= length;
      this.z /= length;

      return this;
   }

   Vertex.prototype.dot = function( vertex ) {
      return this.x*vertex.x + this.y*vertex.y + this.z*vertex.z;
   }

   Vertex.prototype.diff = function( vertex ) {
      let dx = (this.x - vertex.x),
          dy = (this.y - vertex.y),
          dz = (this.z - vertex.z),
          len2 = this.x*this.x + this.y*this.y + this.z*this.z;

      return (dx*dx + dy*dy + dz*dz) / (len2>0 ? len2 : 1e-10);
   }

/*
   Vertex.prototype.lerp = function( a, t ) {
      this.add(
         a.clone().subtract( this ).multiplyScalar( t )
      );

      this.normal.add(
         a.normal.clone().sub( this.normal ).multiplyScalar( t )
      );

      //this.uv.add(
      //   a.uv.clone().sub( this.uv ).multiplyScalar( t )
      //);

      return this;
   };
   Vertex.prototype.interpolate = function( other, t ) {
      return this.clone().lerp( other, t );
   };
*/

   Vertex.prototype.interpolate = function( a, t ) {
      let t1 = 1-t;
      return new Vertex(this.x*t1 + a.x*t, this.y*t1 + a.y*t, this.z*t1 + a.z*t,
                                 this.nx*t1 + a.nx*t, this.ny*t1 + a.ny*t, this.nz*t1 + a.nz*t);
   }

   Vertex.prototype.applyMatrix4 = function ( m ) {

      // input: THREE.Matrix4 affine matrix

      let x = this.x, y = this.y, z = this.z, e = m.elements;

      this.x = e[0] * x + e[4] * y + e[8]  * z + e[12];
      this.y = e[1] * x + e[5] * y + e[9]  * z + e[13];
      this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

      x = this.nx; y = this.ny; z = this.nz;

      this.nx = e[0] * x + e[4] * y + e[8]  * z;
      this.ny = e[1] * x + e[5] * y + e[9]  * z;
      this.nz = e[2] * x + e[6] * y + e[10] * z;

      return this;
   }

   // =======================================================================================

   function Polygon(vertices) {
      this.vertices = vertices || [];
      this.nsign = 1;
      if ( this.vertices.length > 0 )
         this.calculateProperties();
   }

   Polygon.prototype.copyProperties = function(parent, more) {
      this.normal = parent.normal; // .clone();
      this.w = parent.w;
      this.nsign = parent.nsign;
      if (more && (parent.id !== undefined)) {
         this.id = parent.id;
         this.parent = parent;
      }
      return this;
   }

   Polygon.prototype.calculateProperties = function() {
      if (this.normal) return;

      let a = this.vertices[0],
          b = this.vertices[1],
          c = this.vertices[2];

      this.nsign = 1;

      this.normal = b.clone().subtract( a ).cross(
         c.clone().subtract( a )
      ).normalize();

      this.w = this.normal.clone().dot( a );
      return this;
   }

   Polygon.prototype.clone = function() {
      let vertice_count = this.vertices.length,
          polygon = new Polygon;

      for (let i = 0; i < vertice_count; ++i )
         polygon.vertices.push( this.vertices[i].clone() );

      return polygon.copyProperties(this);
   }

   Polygon.prototype.flip = function() {

      /// normal is not changed, only sign variable
      //this.normal.multiplyScalar( -1 );
      //this.w *= -1;

      this.nsign *= -1;

      this.vertices.reverse();

      return this;
   }

   Polygon.prototype.classifyVertex = function( vertex ) {
      let side_value = this.nsign * (this.normal.dot( vertex ) - this.w);

      if ( side_value < -EPSILON ) return BACK;
      if ( side_value > EPSILON ) return FRONT;
      return COPLANAR;
   }

   Polygon.prototype.classifySide = function( polygon ) {
      let i, classification,
          num_positive = 0, num_negative = 0,
          vertice_count = polygon.vertices.length;

      for ( i = 0; i < vertice_count; ++i ) {
         classification = this.classifyVertex( polygon.vertices[i] );
         if ( classification === FRONT ) {
            ++num_positive;
         } else if ( classification === BACK ) {
            ++num_negative;
         }
      }

      if ( num_positive > 0 && num_negative === 0 ) return FRONT;
      if ( num_positive === 0 && num_negative > 0 ) return BACK;
      if ( num_positive === 0 && num_negative === 0 ) return COPLANAR;
      return SPANNING;
   }

   Polygon.prototype.splitPolygon = function( polygon, coplanar_front, coplanar_back, front, back ) {
      let classification = this.classifySide( polygon );

      if ( classification === COPLANAR ) {

         ( (this.nsign * polygon.nsign * this.normal.dot( polygon.normal ) > 0) ? coplanar_front : coplanar_back ).push( polygon );

      } else if ( classification === FRONT ) {

         front.push( polygon );

      } else if ( classification === BACK ) {

         back.push( polygon );

      } else {

         let vertice_count = polygon.vertices.length,
             nnx = this.normal.x,
             nny = this.normal.y,
             nnz = this.normal.z,
             i, j, ti, tj, vi, vj,
             t, v,
             f = [], b = [];

         for ( i = 0; i < vertice_count; ++i ) {

            j = (i + 1) % vertice_count;
            vi = polygon.vertices[i];
            vj = polygon.vertices[j];
            ti = this.classifyVertex( vi );
            tj = this.classifyVertex( vj );

            if ( ti != BACK ) f.push( vi );
            if ( ti != FRONT ) b.push( vi );
            if ( (ti | tj) === SPANNING ) {
               // t = ( this.w - this.normal.dot( vi ) ) / this.normal.dot( vj.clone().subtract( vi ) );
               //v = vi.clone().lerp( vj, t );

               t = (this.w - (nnx*vi.x + nny*vi.y + nnz*vi.z)) / (nnx*(vj.x-vi.x) + nny*(vj.y-vi.y) + nnz*(vj.z-vi.z));

               v = vi.interpolate( vj, t );
               f.push( v );
               b.push( v );
            }
         }

         //if ( f.length >= 3 ) front.push( new Polygon( f ).calculateProperties() );
         //if ( b.length >= 3 ) back.push( new Polygon( b ).calculateProperties() );
         if ( f.length >= 3 ) front.push( new Polygon( f ).copyProperties(polygon, true) );
         if ( b.length >= 3 ) back.push( new Polygon( b ).copyProperties(polygon, true) );
      }
   }

   // ========================================================================================

   function Node( polygons, nodeid ) {
      this.polygons = [];
      this.front = this.back = undefined;

      if (!polygons) return;

      this.divider = polygons[0].clone();

      let polygon_count = polygons.length,
          front = [], back = [];

      for (let i = 0; i < polygon_count; ++i) {
         if (nodeid!==undefined) {
            polygons[i].id = nodeid++;
            delete polygons[i].parent;
         }

         this.divider.splitPolygon( polygons[i], this.polygons, this.polygons, front, back );
      }

      if (nodeid !== undefined) this.maxnodeid = nodeid;

      if ( front.length > 0 )
         this.front = new Node( front );

      if ( back.length > 0 )
         this.back = new Node( back );
   }

   Node.isConvex = function( polygons ) {
      let i, j, len = polygons.length;
      for ( i = 0; i < len; ++i )
         for ( j = 0; j < len; ++j )
            if ( i !== j && polygons[i].classifySide( polygons[j] ) !== BACK ) return false;
      return true;
   }

   Node.prototype.build = function( polygons ) {
      let polygon_count = polygons.length,
          front = [], back = [];

      if ( !this.divider )
         this.divider = polygons[0].clone();

      for (let i = 0; i < polygon_count; ++i )
         this.divider.splitPolygon( polygons[i], this.polygons, this.polygons, front, back );

      if ( front.length > 0 ) {
         if ( !this.front ) this.front = new Node();
         this.front.build( front );
      }

      if ( back.length > 0 ) {
         if ( !this.back ) this.back = new Node();
         this.back.build( back );
      }
   }

   Node.prototype.collectPolygons = function(arr) {
      let len = this.polygons.length;
      for (let i = 0; i < len; ++i) arr.push(this.polygons[i]);
      if ( this.front ) this.front.collectPolygons(arr);
      if ( this.back ) this.back.collectPolygons(arr);
      return arr;
   }

   Node.prototype.allPolygons = function() {
      let polygons = this.polygons.slice();
      if ( this.front ) polygons = polygons.concat( this.front.allPolygons() );
      if ( this.back ) polygons = polygons.concat( this.back.allPolygons() );
      return polygons;
   }

   Node.prototype.numPolygons = function() {
      let res = this.polygons.length;
      if ( this.front ) res += this.front.numPolygons();
      if ( this.back ) res += this.back.numPolygons();
      return res;
   }

   Node.prototype.clone = function() {
      let node = new Node();

      node.divider = this.divider.clone();
      node.polygons = this.polygons.map( function( polygon ) { return polygon.clone(); } );
      node.front = this.front && this.front.clone();
      node.back = this.back && this.back.clone();

      return node;
   }

   Node.prototype.invert = function() {
      let polygon_count = this.polygons.length;

      for (let i = 0; i < polygon_count; ++i )
         this.polygons[i].flip();

      this.divider.flip();
      if ( this.front ) this.front.invert();
      if ( this.back ) this.back.invert();

      let temp = this.front;
      this.front = this.back;
      this.back = temp;

      return this;
   }

   Node.prototype.clipPolygons = function( polygons ) {

      if ( !this.divider ) return polygons.slice();

      let polygon_count = polygons.length, front = [], back = [];

      for (let i = 0; i < polygon_count; ++i )
         this.divider.splitPolygon( polygons[i], front, back, front, back );

      if ( this.front ) front = this.front.clipPolygons( front );
      if ( this.back ) back = this.back.clipPolygons( back );
      else back = [];

      return front.concat( back );
   }

   Node.prototype.clipTo = function( node ) {
      this.polygons = node.clipPolygons( this.polygons );
      if ( this.front ) this.front.clipTo( node );
      if ( this.back ) this.back.clipTo( node );
   }

   // ========================================================================================

   let CreateBufferGeometry = polygons => {
      let i, j, polygon_count = polygons.length, buf_size = 0;

      for ( i = 0; i < polygon_count; ++i )
         buf_size += (polygons[i].vertices.length - 2) * 9;

      let positions_buf = new Float32Array(buf_size),
          normals_buf = new Float32Array(buf_size),
          iii = 0, polygon;

      function CopyVertex(vertex) {

         positions_buf[iii] = vertex.x;
         positions_buf[iii+1] = vertex.y;
         positions_buf[iii+2] = vertex.z;

         normals_buf[iii] = polygon.nsign * vertex.nx;
         normals_buf[iii+1] = polygon.nsign * vertex.ny;
         normals_buf[iii+2] = polygon.nsign * vertex.nz;
         iii+=3;
      }

      for ( i = 0; i < polygon_count; ++i ) {
         polygon = polygons[i];
         for ( j = 2; j < polygon.vertices.length; ++j ) {
            CopyVertex(polygon.vertices[0]);
            CopyVertex(polygon.vertices[j-1]);
            CopyVertex(polygon.vertices[j]);
         }
      }

      let geometry = new THREE.BufferGeometry();
      geometry.setAttribute( 'position', new THREE.BufferAttribute( positions_buf, 3 ) );
      geometry.setAttribute( 'normal', new THREE.BufferAttribute( normals_buf, 3 ) );

      // geometry.computeVertexNormals();
      return geometry;
   }

   function Geometry(geometry, transfer_matrix, nodeid, flippedMesh ) {
      // Convert THREE.Geometry to ThreeBSP

      if ( geometry instanceof THREE.Geometry ) {
         this.matrix = null; // new THREE.Matrix4; not create matrix when do not needed
      } else if ( geometry instanceof THREE.Mesh ) {
         // #todo: add hierarchy support
         geometry.updateMatrix();
         transfer_matrix = this.matrix = geometry.matrix.clone();
         geometry = geometry.geometry;
      } else if ( geometry instanceof Node ) {
         this.tree = geometry;
         this.matrix = null; // new THREE.Matrix4;
         return this;
      } else if ( geometry instanceof THREE.BufferGeometry ) {
         let pos_buf = geometry.getAttribute('position').array,
             norm_buf = geometry.getAttribute('normal').array,
             polygons = [], polygon, vert1, vert2, vert3;

         for (let i=0; i < pos_buf.length; i+=9) {
            polygon = new Polygon;

            vert1 = new Vertex( pos_buf[i], pos_buf[i+1], pos_buf[i+2], norm_buf[i], norm_buf[i+1], norm_buf[i+2]);
            if (transfer_matrix) vert1.applyMatrix4(transfer_matrix);

            vert2 = new Vertex( pos_buf[i+3], pos_buf[i+4], pos_buf[i+5], norm_buf[i+3], norm_buf[i+4], norm_buf[i+5]);
            if (transfer_matrix) vert2.applyMatrix4(transfer_matrix);

            vert3 = new Vertex( pos_buf[i+6], pos_buf[i+7], pos_buf[i+8], norm_buf[i+6], norm_buf[i+7], norm_buf[i+8]);
            if (transfer_matrix) vert3.applyMatrix4(transfer_matrix);

            if (flippedMesh) polygon.vertices.push( vert1, vert3, vert2 );
                        else polygon.vertices.push( vert1, vert2, vert3 );

            polygon.calculateProperties();
            polygons.push( polygon );
         }

         this.tree = new Node( polygons, nodeid );
         if (nodeid!==undefined) this.maxid = this.tree.maxnodeid;
         return this;

      } else if (geometry.polygons && (geometry.polygons[0] instanceof Polygon)) {
         let polygons = geometry.polygons;

         for (let i=0;i<polygons.length;++i) {
            let polygon = polygons[i];
            if (transfer_matrix) {
               for (let n=0;n<polygon.vertices.length;++n)
                  polygon.vertices[n].applyMatrix4(transfer_matrix);
            }

            polygon.calculateProperties();
         }

         this.tree = new Node( polygons, nodeid );
         if (nodeid!==undefined) this.maxid = this.tree.maxnodeid;
         return this;

      } else {
         throw 'ThreeBSP: Given geometry is unsupported';
      }

      let polygons = [],
          nfaces = geometry.faces.length,
          face, polygon, vertex, normal, useVertexNormals;

      for (let i = 0; i < nfaces; ++i ) {
         face = geometry.faces[i];
         normal = face.normal;
         // faceVertexUvs = geometry.faceVertexUvs[0][i];
         polygon = new Polygon;

         useVertexNormals = face.vertexNormals && (face.vertexNormals.length==3);

         vertex = geometry.vertices[ face.a ];
         if (useVertexNormals) normal = face.vertexNormals[0];
         // uvs = faceVertexUvs ? new THREE.Vector2( faceVertexUvs[0].x, faceVertexUvs[0].y ) : null;
         vertex = new Vertex( vertex.x, vertex.y, vertex.z, normal.x, normal.y, normal.z /*face.normal , uvs */ );
         if (transfer_matrix) vertex.applyMatrix4(transfer_matrix);
         polygon.vertices.push( vertex );

         vertex = geometry.vertices[ face.b ];
         if (useVertexNormals) normal = face.vertexNormals[1];
         //uvs = faceVertexUvs ? new THREE.Vector2( faceVertexUvs[1].x, faceVertexUvs[1].y ) : null;
         vertex = new Vertex( vertex.x, vertex.y, vertex.z, normal.x, normal.y, normal.z /*face.normal , uvs */ );
         if (transfer_matrix) vertex.applyMatrix4(transfer_matrix);
         polygon.vertices.push( vertex );

         vertex = geometry.vertices[ face.c ];
         if (useVertexNormals) normal = face.vertexNormals[2];
         // uvs = faceVertexUvs ? new THREE.Vector2( faceVertexUvs[2].x, faceVertexUvs[2].y ) : null;
         vertex = new Vertex( vertex.x, vertex.y, vertex.z, normal.x, normal.y, normal.z /*face.normal, uvs */ );
         if (transfer_matrix) vertex.applyMatrix4(transfer_matrix);
         polygon.vertices.push( vertex );

         polygon.calculateProperties();
         polygons.push( polygon );
      }

      this.tree = new Node( polygons, nodeid );
      if (nodeid!==undefined) this.maxid = this.tree.maxnodeid;
   }

   Geometry.prototype.subtract = function( other_tree ) {
      let a = this.tree.clone(),
          b = other_tree.tree.clone();

      a.invert();
      a.clipTo( b );
      b.clipTo( a );
      b.invert();
      b.clipTo( a );
      b.invert();
      a.build( b.allPolygons() );
      a.invert();
      a = new Geometry( a );
      a.matrix = this.matrix;
      return a;
   }

   Geometry.prototype.union = function( other_tree ) {
      let a = this.tree.clone(),
         b = other_tree.tree.clone();

      a.clipTo( b );
      b.clipTo( a );
      b.invert();
      b.clipTo( a );
      b.invert();
      a.build( b.allPolygons() );
      a = new Geometry( a );
      a.matrix = this.matrix;
      return a;
   }

   Geometry.prototype.intersect = function( other_tree ) {
      let a = this.tree.clone(),
         b = other_tree.tree.clone();

      a.invert();
      b.clipTo( a );
      b.invert();
      a.clipTo( b );
      b.clipTo( a );
      a.build( b.allPolygons() );
      a.invert();
      a = new Geometry( a );
      a.matrix = this.matrix;
      return a;
   }

   Geometry.prototype.tryToCompress = function(polygons) {

      if (this.maxid === undefined) return;

      let arr = [], parts, foundpair,
          nreduce = 0, n, len = polygons.length,
          p, p1, p2, i1, i2;

      // sort out polygons
      for (n=0;n<len;++n) {
         p = polygons[n];
         if (p.id === undefined) continue;
         if (arr[p.id] === undefined) arr[p.id] = [];

         arr[p.id].push(p);
      }

      for(n=0; n<arr.length; ++n) {
         parts = arr[n];
         if (parts===undefined) continue;

         len = parts.length;

         foundpair = (len > 1);

         while (foundpair) {
            foundpair = false;

            for (i1 = 0; i1<len-1; ++i1) {
               p1 = parts[i1];
               if (!p1 || !p1.parent) continue;
               for (i2 = i1+1; i2 < len; ++i2) {
                  p2 = parts[i2];
                  if (p2 && (p1.parent === p2.parent) && (p1.nsign === p2.nsign)) {

                     if (p1.nsign !== p1.parent.nsign) p1.parent.flip();

                     nreduce++;
                     parts[i1] = p1.parent;
                     parts[i2] = null;
                     if (p1.parent.vertices.length < 3) console.log('something wrong with parent');
                     foundpair = true;
                     break;
                  }
               }
            }
         }
      }

      if (nreduce>0) {
         polygons.splice(0, polygons.length);

         for(n=0;n<arr.length;++n) {
            parts = arr[n];
            if (parts !== undefined)
               for (i1=0,len=parts.length; i1<len;++i1)
                  if (parts[i1]) polygons.push(parts[i1]);
         }

      }
   }

   Geometry.prototype.direct_subtract = function( other_tree ) {
      let a = this.tree,
          b = other_tree.tree;
      a.invert();
      a.clipTo( b );
      b.clipTo( a );
      b.invert();
      b.clipTo( a );
      b.invert();
      a.build( b.collectPolygons([]) );
      a.invert();
      return this;
   }

   Geometry.prototype.direct_union = function( other_tree ) {
      let a = this.tree,
          b = other_tree.tree;

      a.clipTo( b );
      b.clipTo( a );
      b.invert();
      b.clipTo( a );
      b.invert();
      a.build( b.collectPolygons([]) );
      return this;
   }

   Geometry.prototype.direct_intersect = function( other_tree ) {
      let a = this.tree,
          b = other_tree.tree;

      a.invert();
      b.clipTo( a );
      b.invert();
      a.clipTo( b );
      b.clipTo( a );
      a.build( b.collectPolygons([]) );
      a.invert();
      return this;
   }

   Geometry.prototype.cut_from_plane = function( other_tree) {
      // just cut peaces from second geometry, which just simple plane

      let a = this.tree,
          b = other_tree.tree;

      a.invert();
      b.clipTo( a );

      return this;
   }

   Geometry.prototype.toGeometry = function() {
      let i, j,
         matrix = this.matrix ? new THREE.Matrix4().getInverse( this.matrix ) : null,
         geometry = new THREE.Geometry(),
         polygons = this.tree.collectPolygons([]),
         polygon_count = polygons.length,
         polygon, polygon_vertice_count,
         vertice_dict = {},
         vertex_idx_a, vertex_idx_b, vertex_idx_c,
         vertex, face;

      for ( i = 0; i < polygon_count; ++i ) {
         polygon = polygons[i];
         polygon_vertice_count = polygon.vertices.length;

         for ( j = 2; j < polygon_vertice_count; ++j ) {
            // verticeUvs = [];

            vertex = polygon.vertices[0];
            // verticeUvs.push( new THREE.Vector2( vertex.uv.x, vertex.uv.y ) );
            vertex = new THREE.Vector3( vertex.x, vertex.y, vertex.z );
            if (matrix) vertex.applyMatrix4(matrix);

            if ( typeof vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ] !== 'undefined' ) {
               vertex_idx_a = vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ];
            } else {
               geometry.vertices.push( vertex );
               vertex_idx_a = vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ] = geometry.vertices.length - 1;
            }

            vertex = polygon.vertices[j-1];
            // verticeUvs.push( new THREE.Vector2( vertex.uv.x, vertex.uv.y ) );
            vertex = new THREE.Vector3( vertex.x, vertex.y, vertex.z );
            if (matrix) vertex.applyMatrix4(matrix);
            if ( typeof vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ] !== 'undefined' ) {
               vertex_idx_b = vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ];
            } else {
               geometry.vertices.push( vertex );
               vertex_idx_b = vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ] = geometry.vertices.length - 1;
            }

            vertex = polygon.vertices[j];
            // verticeUvs.push( new THREE.Vector2( vertex.uv.x, vertex.uv.y ) );
            vertex = new THREE.Vector3( vertex.x, vertex.y, vertex.z );
            if (matrix) vertex.applyMatrix4(matrix);
            if ( typeof vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ] !== 'undefined' ) {
               vertex_idx_c = vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ];
            } else {
               geometry.vertices.push( vertex );
               vertex_idx_c = vertice_dict[ vertex.x + ',' + vertex.y + ',' + vertex.z ] = geometry.vertices.length - 1;
            }

            face = new THREE.Face3(
               vertex_idx_a,
               vertex_idx_b,
               vertex_idx_c,
               new THREE.Vector3( polygon.normal.x, polygon.normal.y, polygon.normal.z )
            );

            geometry.faces.push( face );
            // geometry.faceVertexUvs[0].push( verticeUvs );
         }

      }
      return geometry;
   }

   Geometry.prototype.scale = function(x,y,z) {
      // try to scale as THREE.BufferGeometry
      let polygons = this.tree.collectPolygons([]);

      for (let i = 0; i < polygons.length; ++i) {
         let polygon = polygons[i];
         for (let k=0; k < polygon.vertices.length; ++k) {
            let v = polygon.vertices[k];
            v.x *= x;
            v.y *= y;
            v.z *= z;
         }
         delete polygon.normal;
         polygon.calculateProperties();
      }
   }

   Geometry.prototype.toPolygons = function() {
      let polygons = this.tree.collectPolygons([]);

      this.tryToCompress(polygons);

      for (let i = 0; i < polygons.length; ++i ) {
         delete polygons[i].id;
         delete polygons[i].parent;
      }

      return polygons;
   }

   Geometry.prototype.toBufferGeometry = function() {
      return CreateBufferGeometry(this.toPolygons());
   }

   Geometry.prototype.toMesh = function( material ) {
      let geometry = this.toGeometry(),
         mesh = new THREE.Mesh( geometry, material );

      if (this.matrix) {
         mesh.position.setFromMatrixPosition( this.matrix );
         mesh.rotation.setFromRotationMatrix( this.matrix );
      }

      return mesh;
   }

   // ================================================================================================

   let ThreeBSP = {};

   ThreeBSP.CreateNormal = function(axis_name, pos, size) {
      // create geometry to make cut on specified axis

      let vert1, vert2, vert3;

      if (!size || (size<10000)) size = 10000;

      switch(axis_name) {
         case "x":
            vert1 = new Vertex(pos, -3*size,    size, 1, 0, 0),
            vert3 = new Vertex(pos,    size,    size, 1, 0, 0),
            vert2 = new Vertex(pos,    size, -3*size, 1, 0, 0);
            break;
         case "y":
            vert1 = new Vertex(-3*size,  pos,    size, 0, 1, 0),
            vert2 = new Vertex(   size,  pos,    size, 0, 1, 0),
            vert3 = new Vertex(   size,  pos, -3*size, 0, 1, 0);
            break;
         case "z":
            vert1 = new Vertex(-3*size,    size, pos, 0, 0, 1),
            vert3 = new Vertex(   size,    size, pos, 0, 0, 1),
            vert2 = new Vertex(   size, -3*size, pos, 0, 0, 1);
            break;
      }

      let polygon = new Polygon([vert1, vert2, vert3]);
      polygon.calculateProperties();

      let node = new Node([polygon]);

      return new Geometry(node);
   }

   ThreeBSP.Vertex = Vertex;
   ThreeBSP.Geometry = Geometry;
   ThreeBSP.Polygon = Polygon;

   ThreeBSP.CreateBufferGeometry = CreateBufferGeometry;

   if (JSROOT.nodejs) module.exports = ThreeBSP;
   return ThreeBSP;

});
