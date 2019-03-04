/** @file EveElements.js 
 * used only together with OpenUI5 */

// TODO: add dependency from JSROOT components

sap.ui.define(['rootui5/eve7/lib/EveManager'], function(EveManager) {

   "use strict";

   function EveElemControl(mesh) {
      // JSROOT.Painter.GeoDrawingControl.call(this);
      this.mesh = mesh;
   }

   EveElemControl.prototype = Object.create(JSROOT.Painter.GeoDrawingControl.prototype);

   EveElemControl.prototype.invokeSceneMethod = function(fname, arg1, arg2)
   {
      if (!this.mesh) return false;

      var s = this.mesh.scene;
      if (s && (typeof s[fname] == "function"))
         return s[fname](this.mesh, arg1, arg2, this.evnt);
      return false;
   }

   EveElemControl.prototype.separateDraw = false;

   EveElemControl.prototype.setHighlight = function(col, indx)
   {
      // special hook here
      if (this.invokeSceneMethod("processElementHighlighted", col, indx))
         return true;
      // should not be used, keep as fallback
      return this.drawSpecial(col || this.mesh.select_col, indx);
   }

   EveElemControl.prototype.setSelected = function(col, indx)
   {
      if (this.invokeSceneMethod("processElementSelected", col, indx))
         return true;
      // should not be used, keep as fallback
      var m = this.mesh;
      if ((m.select_col == col) && (m.select_indx == indx)) { col = null; indx = undefined; }
      m.select_col  = col;
      m.select_indx = indx;
      // special hook here
      return this.drawSpecial(col, indx);
   }

   ////////////////////////////////////////////////


   var GL = { POINTS: 0, LINES: 1, LINE_LOOP: 2, LINE_STRIP: 3, TRIANGLES: 4 };


   function EveElements()
   {
   }

   /** Test if render data has vertex buffer. Make logging if not. Only for debug purposes */
   EveElements.prototype.TestRnr = function(name, obj, rnrData) {

      if (obj && rnrData && rnrData.vtxBuff) return false;

      var cnt = this[name] || 0;

      if (cnt++ < 5) console.log(name, obj, rnrData);

      this[name] = cnt;

      return true;
   }

   EveElements.prototype.makeHit = function(hit, rnrData) {
      if (this.TestRnr("hit", hit, rnrData)) return null;

      var hit_size = 8 * rnrData.fMarkerSize;
      var size     = rnrData.vtxBuff.length / 3;
      var pnts     = new JSROOT.Painter.PointsCreator(size, true, hit_size);

      for (var i=0; i<size; i++)
         pnts.AddPoint(rnrData.vtxBuff[i*3],rnrData.vtxBuff[i*3+1],rnrData.vtxBuff[i*3+2]);

      var mesh = pnts.CreatePoints(JSROOT.Painter.root_colors[hit.fMarkerColor]);

      // use points control to toggle highlight and selection
      // mesh.get_ctrl = function() { return new JSROOT.Painter.PointsControl(this); }

      mesh.get_ctrl = function() { return new EveElemControl(this); }

      mesh.highlightScale = 2;

      mesh.object = hit;
      mesh.geo_name = hit.fName;
      mesh.geo_object = hit.fMasterId || hit.fElementId;

      mesh.material.sizeAttenuation = false;
      mesh.material.size = hit.fMarkerSize;
      return mesh;
   }

   EveElements.prototype.makeTrack = function(track, rnrData)
   {
      if (this.TestRnr("track", track, rnrData)) return null;

      var N = rnrData.vtxBuff.length/3;
      var track_width = track.fLineWidth || 1;
      var track_color = JSROOT.Painter.root_colors[track.fLineColor] || "rgb(255,0,255)";

      if (JSROOT.browser.isWin) track_width = 1;  // not supported on windows

      var buf = new Float32Array((N-1) * 6), pos = 0;
      for (var k=0;k<(N-1);++k) {
         buf[pos]   = rnrData.vtxBuff[k*3];
         buf[pos+1] = rnrData.vtxBuff[k*3+1];
         buf[pos+2] = rnrData.vtxBuff[k*3+2];

         var breakTrack = false;
         if (rnrData.idxBuff)
            for (var b = 0; b < rnrData.idxBuff.length; b++) {
               if ( (k+1) == rnrData.idxBuff[b]) {
                  breakTrack = true;
                  break;
               }
            }

         if (breakTrack) {
            buf[pos+3] = rnrData.vtxBuff[k*3];
            buf[pos+4] = rnrData.vtxBuff[k*3+1];
            buf[pos+5] = rnrData.vtxBuff[k*3+2];
         } else {
            buf[pos+3] = rnrData.vtxBuff[k*3+3];
            buf[pos+4] = rnrData.vtxBuff[k*3+4];
            buf[pos+5] = rnrData.vtxBuff[k*3+5];
         }

         // console.log(" vertex ", buf[pos],buf[pos+1], buf[pos+2],buf[pos+3], buf[pos+4],  buf[pos+5]);
         pos+=6;
      }

      var lineMaterial;
      if (track.fLineStyle == 1) {
         lineMaterial = new THREE.LineBasicMaterial({ color: track_color, linewidth: track_width });
      } else {
         //lineMaterial = new THREE.LineDashedMaterial({ color: track_color, linewidth: track_width, gapSize: parseInt(track.fLineStyle) });
         lineMaterial = new THREE.LineDashedMaterial({ color: track_color, linewidth: track_width, dashSize:3, gapSize: 1 });
      }

      var geom = new THREE.BufferGeometry();
      geom.addAttribute( 'position', new THREE.BufferAttribute( buf, 3 )  );
      var line = new THREE.LineSegments(geom, lineMaterial);

      line.object = track;
      line.geo_name = track.fName;
      line.geo_object = track.fMasterId || track.fElementId;
      line.hightlightWidthScale = 2;

      line.get_ctrl = function() { return new EveElemControl(this); }

      return line;
   }

   EveElements.prototype.makeJet = function(jet, rnrData)
   {
      if (this.TestRnr("jet", jet, rnrData)) return null;

      // console.log("make jet ", jet);
      // var jet_ro = new THREE.Object3D();
      var pos_ba = new THREE.BufferAttribute( rnrData.vtxBuff, 3 );
      var N      = rnrData.vtxBuff.length / 3;

      var geo_body = new THREE.BufferGeometry();
      geo_body.addAttribute('position', pos_ba);
      var idcs = [0, N-1, 1];
      for (var i = 1; i < N - 1; ++i)
         idcs.push( 0, i, i + 1 );
      geo_body.setIndex( idcs );
      geo_body.computeVertexNormals();

      var geo_rim = new THREE.BufferGeometry();
      geo_rim.addAttribute('position', pos_ba);
      idcs = new Uint16Array(N-1);
      for (var i = 1; i < N; ++i) idcs[i-1] = i;
      geo_rim.setIndex(new THREE.BufferAttribute( idcs, 1 ));

      var geo_rays = new THREE.BufferGeometry();
      geo_rays.addAttribute('position', pos_ba);
      idcs = [];
      for (var i = 1; i < N; i += 4)
         idcs.push( 0, i );
      geo_rays.setIndex( idcs );

      var mcol = JSROOT.Painter.root_colors[jet.fMainColor];
      var lcol = JSROOT.Painter.root_colors[jet.fLineColor];

      var mesh = new THREE.Mesh(geo_body, new THREE.MeshPhongMaterial({ depthWrite: false, color: mcol, transparent: true, opacity: 0.5, side: THREE.DoubleSide }));
      var line1 = new THREE.LineLoop(geo_rim,  new THREE.LineBasicMaterial({ linewidth: 2,   color: lcol, transparent: true, opacity: 0.5 }))
      var line2 = new THREE.LineSegments(geo_rays, new THREE.LineBasicMaterial({ linewidth: 0.5, color: lcol, transparent: true, opacity: 0.5 }));

      // jet_ro.add( mesh  );
      mesh.add( line1 );
      mesh.add( line2 );

      mesh.object = jet;
      mesh.geo_name = jet.fName;
      mesh.geo_object = jet.fMasterId || jet.fElementId;

      mesh.get_ctrl = function() { return new EveElemControl(this); }

      return mesh;
   }

   EveElements.prototype.makeJetProjected = function(jet, rnrData)
   {
      // JetProjected has 3 or 4 points. 0-th is apex, others are rim.
      // Fourth point is only present in RhoZ when jet hits barrel/endcap transition.

      // console.log("makeJetProjected ", jet);

      if (this.TestRnr("jetp", jet, rnrData)) return null;


      var pos_ba = new THREE.BufferAttribute( rnrData.vtxBuff, 3 );
      var N      = rnrData.vtxBuff.length / 3;

      var geo_body = new THREE.BufferGeometry();
      geo_body.addAttribute('position', pos_ba);
      var idcs = [0, 2, 1];
      if (N > 3)
         idcs.push( 0, 3, 2 );
      geo_body.setIndex( idcs );
      geo_body.computeVertexNormals();

      var geo_rim = new THREE.BufferGeometry();
      geo_rim.addAttribute('position', pos_ba);
      idcs = new Uint16Array(N-1);
      for (var i = 1; i < N; ++i) idcs[i-1] = i;
      geo_rim.setIndex(new THREE.BufferAttribute( idcs, 1 ));

      var geo_rays = new THREE.BufferGeometry();
      geo_rays.addAttribute('position', pos_ba);
      idcs = [ 0, 1, 0, N-1 ];
      geo_rays.setIndex( idcs );

      var fcol = JSROOT.Painter.root_colors[jet.fFillColor];
      var lcol = JSROOT.Painter.root_colors[jet.fLineColor];
      // Process transparency !!!
      // console.log("cols", fcol, lcol);

      // double-side material required for correct tracing of colors - otherwise points sequence should be changed
      var mesh = new THREE.Mesh(geo_body, new THREE.MeshBasicMaterial({ depthWrite: false, color: fcol, transparent: true, opacity: 0.5, side: THREE.DoubleSide }));
      var line1 = new THREE.Line(geo_rim,  new THREE.LineBasicMaterial({ linewidth: 2, color: lcol, transparent: true, opacity: 0.5 }));
      var line2 = new THREE.LineSegments(geo_rays, new THREE.LineBasicMaterial({ linewidth: 1, color: lcol, transparent: true, opacity: 0.5 }));

      // jet_ro.add( mesh  );
      mesh.add( line1 );
      mesh.add( line2 );

      mesh.object = jet;
      mesh.geo_name = jet.fName;
      mesh.geo_object = jet.fMasterId || jet.fElementId;

      mesh.get_ctrl = function() { return new EveElemControl(this); }

      return mesh;
   }

   EveElements.prototype.makeEveGeometry = function(rnr_data, force)
   {
      var nVert = rnr_data.idxBuff[1]*3;

      if (rnr_data.idxBuff[0] != GL.TRIANGLES)  throw "Expect triangles first.";
      if (2 + nVert != rnr_data.idxBuff.length) throw "Expect single list of triangles in index buffer.";

      if (this.useIndexAsIs) {
         var body = new THREE.BufferGeometry();
         body.addAttribute('position', new THREE.BufferAttribute( rnr_data.vtxBuff, 3 ));
         body.setIndex(new THREE.BufferAttribute( rnr_data.idxBuff, 1 ));
         body.setDrawRange(2, nVert);
         // this does not work correctly - draw range ignored when calculating normals
         // even worse - shift 2 makes complete logic wrong while wrong triangle are extracted
         // Let see if it will be fixed https://github.com/mrdoob/three.js/issues/15560
         body.computeVertexNormals();
         return body;
      }

      var vBuf = new Float32Array(nVert*3); // plain buffer with all vertices
      var nBuf = null;                      // plaint buffer with normals per vertex

      if (rnr_data.nrmBuff) {
         if (rnr_data.nrmBuff.length !== nVert) throw "Expect normals per face";
         nBuf = new Float32Array(nVert*3);
      }

      for (var i=0;i<nVert;++i) {
         var pos = rnr_data.idxBuff[i+2];
         vBuf[i*3] = rnr_data.vtxBuff[pos*3];
         vBuf[i*3+1] = rnr_data.vtxBuff[pos*3+1];
         vBuf[i*3+2] = rnr_data.vtxBuff[pos*3+2];
         if (nBuf) {
            pos = i - i%3;
            nBuf[i*3] = rnr_data.nrmBuff[pos];
            nBuf[i*3+1] = rnr_data.nrmBuff[pos+1];
            nBuf[i*3+2] = rnr_data.nrmBuff[pos+2];
         }
      }

      var body = new THREE.BufferGeometry();

      body.addAttribute('position', new THREE.BufferAttribute( vBuf, 3 ));

      if (nBuf)
         body.addAttribute('normal', new THREE.BufferAttribute( nBuf, 3 ));
      else
         body.computeVertexNormals();

      // XXXX Fix this. It seems we could have flat shading with usage of simple shaders.
      // XXXX Also, we could do edge detect on the server for outlines.
      // XXXX a) 3d objects - angle between triangles >= 85 degrees (or something);
      // XXXX b) 2d objects - segment only has one triangle.
      // XXXX Somewhat orthogonal - when we do tesselation, conversion from quads to
      // XXXX triangles is trivial, we could do it before invoking the big guns (if they are even needed).
      // XXXX Oh, and once triangulated, we really don't need to store 3 as number of verts in a poly each time.
      // XXXX Or do we? We might need it for projection stuff.

      return body;
   }

   EveElements.prototype.makeEveGeoShape = function(egs, rnr_data)
   {
      var egs_ro = new THREE.Object3D();

      var geom = this.makeEveGeometry(rnr_data);

      var fcol = JSROOT.Painter.root_colors[egs.fFillColor];

      var material = new THREE.MeshPhongMaterial({// side: THREE.DoubleSide,
                          depthWrite: false, color:fcol, transparent: true, opacity: 0.2 });

      var mesh = new THREE.Mesh(geom, material);

      egs_ro.add(mesh);

      return egs_ro;
   }

   /** Keep this old code for reference, arbitrary referencing via index does not work */
   EveElements.prototype.makePolygonSetProjectedOld = function(psp, rnr_data)
   {
      var psp_ro = new THREE.Object3D();
      var pos_ba = new THREE.BufferAttribute( rnr_data.vtxBuff, 3 );
      var idx_ba = new THREE.BufferAttribute( rnr_data.idxBuff, 1 );

      var ib_len = rnr_data.idxBuff.length;

      var fcol = JSROOT.Painter.root_colors[psp.fMainColor];
      var line_mat = new THREE.LineBasicMaterial({color:fcol });

      for (var ib_pos = 0; ib_pos < ib_len; )
      {
         if (rnr_data.idxBuff[ib_pos] == GL.TRIANGLES)
         {
            // Sergey: make check, for now here many wrong values
            var is_ok = true, maxindx = rnr_data.vtxBuff.length/3;
            for (var k=0;is_ok && (k < 3*rnr_data.idxBuff[ib_pos + 1]); ++k)
               if (rnr_data.idxBuff[ib_pos+2+k] > maxindx) is_ok = false;

            if (is_ok) {
               var body = new THREE.BufferGeometry();
               body.addAttribute('position', pos_ba);
               body.setIndex(idx_ba);
               body.setDrawRange(ib_pos + 2, 3 * rnr_data.idxBuff[ib_pos + 1]);
               body.computeVertexNormals();
               var material = new THREE.MeshBasicMaterial({ side: THREE.DoubleSide, depthWrite: false,
                                               color:fcol, transparent: true, opacity: 0.4 });

               psp_ro.add( new THREE.Mesh(body, material) );
            } else {
               console.log('Error in makePolygonSetProjected - wrong GL.TRIANGLES indexes');
            }

            ib_pos += 2 + 3 * rnr_data.idxBuff[ib_pos + 1];
         }
         else if (rnr_data.idxBuff[ib_pos] == GL.LINE_LOOP)
         {
            var body = new THREE.BufferGeometry();
            body.addAttribute('position', pos_ba);
            body.setIndex(idx_ba);
            body.setDrawRange(ib_pos + 2, rnr_data.idxBuff[ib_pos + 1]);

            psp_ro.add( new THREE.LineLoop(body, line_mat) );

            ib_pos += 2 + rnr_data.idxBuff[ib_pos + 1];
         }
         else
         {
            console.error("Unexpected primitive type " + rnr_data.idxBuff[ib_pos]);
            break;
         }

      }

      return psp_ro;
   }

   EveElements.prototype.makePolygonSetProjected = function(psp, rnr_data)
   {
      if (this.useIndexAsIs)
         return this.makePolygonSetProjectedOld(psp, rnr_data);

      var psp_ro = new THREE.Object3D(),
          ib_len = rnr_data.idxBuff.length,
          fcol = JSROOT.Painter.root_colors[psp.fMainColor];

      for (var ib_pos = 0; ib_pos < ib_len; )
      {
         if (rnr_data.idxBuff[ib_pos] == GL.TRIANGLES) {

            var nVert = rnr_data.idxBuff[ib_pos + 1] * 3,
                vBuf = new Float32Array(nVert*3); // plain buffer with all vertices

            for (var k=0;k<nVert;++k) {
               var pos = rnr_data.idxBuff[ib_pos+2+k];
               if (pos*3 > rnr_data.vtxBuff.length) { vBuf = null; break; }
               vBuf[k*3] = rnr_data.vtxBuff[pos*3];
               vBuf[k*3+1] = rnr_data.vtxBuff[pos*3+1];
               vBuf[k*3+2] = rnr_data.vtxBuff[pos*3+2];
            }

            if (vBuf) {
               var body = new THREE.BufferGeometry();
               body.addAttribute('position', new THREE.BufferAttribute( vBuf, 3 ));
               body.computeVertexNormals();
               var material = new THREE.MeshBasicMaterial({ side: THREE.DoubleSide, depthWrite: false,
                                  color:fcol, transparent: true, opacity: 0.4 });
               psp_ro.add( new THREE.Mesh(body, material) );
            } else {
               console.log('Error in makePolygonSetProjected - wrong GL.TRIANGLES indexes');
            }

            ib_pos += 2 + nVert;
         }
         else if (rnr_data.idxBuff[ib_pos] == GL.LINE_LOOP)
         {
            var nVert = rnr_data.idxBuff[ib_pos + 1],
                vBuf = new Float32Array(nVert*3); // plain buffer with all vertices

            for (var k=0;k<nVert;++k) {
               var pos = rnr_data.idxBuff[ib_pos+2+k];
               if (pos*3 > rnr_data.vtxBuff.length) { vBuf = null; break; }
               vBuf[k*3] = rnr_data.vtxBuff[pos*3];
               vBuf[k*3+1] = rnr_data.vtxBuff[pos*3+1];
               vBuf[k*3+2] = rnr_data.vtxBuff[pos*3+2];
            }

            if (vBuf) {
               var body = new THREE.BufferGeometry();
               body.addAttribute('position', new THREE.BufferAttribute( vBuf, 3 ));
               var line_mat = new THREE.LineBasicMaterial({color:fcol });
               psp_ro.add( new THREE.LineLoop(body, line_mat) );
            } else {
               console.log('Error in makePolygonSetProjected - wrong GL.LINE_LOOP indexes');
            }

            ib_pos += 2 + nVert;
         }
         else
         {
            console.error("Unexpected primitive type " + rnr_data.idxBuff[ib_pos]);
            break;
         }

      }

      return psp_ro;
   }

   ////////////////////////////////////////////////////////////////////////////////////////////


   function StraightLineSetControl(mesh)
   {
      EveElemControl.call(this, mesh);
   }

   StraightLineSetControl.prototype = Object.create(EveElemControl.prototype);

   StraightLineSetControl.prototype.separateDraw = true;

   StraightLineSetControl.prototype.cleanup = function()
   {
      if (!this.mesh) return;
      this.drawSpecial(null, undefined, "h");
      this.drawSpecial(null, undefined, "s");
      delete this.mesh;
   }

   StraightLineSetControl.prototype.extractIndex = function(intersect)
   {
      if (!intersect || intersect.index===undefined) return undefined;
      return intersect.index/2; // return segment id - not a point
   }

   StraightLineSetControl.prototype.setSelected = function(col, indx)
   {
      if (this.invokeSceneMethod("processElementSelected", col, indx))
         return true;

      // old code - should not be used now

      console.error("never come here");

      var m = this.mesh;
      if (!non_recurive) {
         if (this.evnt && this.evnt.ctrlKey) {
            var curr = m.select_indx;
            if (typeof curr === "number") curr = [curr]; else if (!curr) curr = [];
            var pos = curr.indexOf(indx);
            if (col && (pos < 0)) curr.push(indx); else
            if (!col && (pos >= 0)) curr.splice(pos, 1);
            indx = curr.length ? curr : undefined;
            if (!indx) col = null;
         } else if ((m.select_col == col) && (m.select_indx == indx)) {
            col = null;
            indx = undefined;
         }
      }

      m.select_col = col;
      m.select_indx = indx;

      return this.drawSpecial(col, indx, "s");
   }

   StraightLineSetControl.prototype.setHighlight = function(col, indx)
   {
      if (this.invokeSceneMethod("processElementHighlighted", col, indx))
         return true;

      this.mesh.h_index = indx;
      return this.drawSpecial(col, indx, "h");
   }

   StraightLineSetControl.prototype.checkHighlightIndex = function(indx)
   {
      if (this.mesh && this.mesh.scene)
         return this.invokeSceneMethod("processCheckHighlight", indx);

      return true; // means index is different
   }

   StraightLineSetControl.prototype.drawSpecial = function(color, index, prefix)
   {
      if ( ! prefix) prefix = "s";

      var did_change = false;

      var m  = this.mesh;
      var ll = prefix + "l_special";
      var mm = prefix + "m_special";

      if (m[ll])
      {
         m.remove(m[ll]);
         JSROOT.Painter.DisposeThreejsObject(m[ll]);
         delete m[ll];
         did_change = true;
      }
      if (m[mm])
      {
         m.remove(m[mm]);
         JSROOT.Painter.DisposeThreejsObject(m[mm]);
         delete m[mm];
         did_change = true;
      }

      if ( ! color)
         return did_change;

      if (typeof index == "number") index = [ index ]; else
      if ( ! index) index = [];

      var geom = new THREE.BufferGeometry();
      geom.addAttribute( 'position', m.children[0].geometry.getAttribute("position") );
      if (index.length == 1)
      {
         geom.setDrawRange(index[0]*2, 2);
      } else if (index.length > 1)
      {
         var idcs = [];
         for (var i = 0; i < index.length; ++i)
            idcs.push(index[i]*2, index[i]*2+1);
         geom.setIndex( idcs );
      }
      var lineMaterial = new THREE.LineBasicMaterial({ color: color, linewidth: 4 });
      var line         = new THREE.LineSegments(geom, lineMaterial);
      line.jsroot_special = true; // special object, exclude from intersections
      m.add(line);
      m[ll] = line;

      var el = m.eve_el, mindx = []

      for (var i = 0; i < index.length; ++i)
      {
         if (index[i] < el.fLinePlexSize)
         {
            var lineid = m.eve_indx[index[i]];

            for (var k = 0; k < el.fMarkerPlexSize; ++k )
            {
               if (m.eve_indx[ k + el.fLinePlexSize] == lineid) mindx.push(k);
            }
         }
      }

      if (mindx.length > 0)
      {
         var pnts = new JSROOT.Painter.PointsCreator(mindx.length, true, 5);

         var arr = m.children[1].geometry.getAttribute("position").array;

         for (var i = 0; i < mindx.length; ++i)
         {
            var p = mindx[i]*3;
            pnts.AddPoint(arr[p], arr[p+1], arr[p+2] );
         }
         var mark = pnts.CreatePoints(color);
         mark.jsroot_special = true; // special object, exclude from intersections
         m.add(mark);
         m[mm] = mark;
      }

      return true;
   }

   EveElements.prototype.makeStraightLineSet = function(el, rnr_data)
   {
      var obj3d = new THREE.Object3D();

      var mainColor = JSROOT.Painter.root_colors[el.fMainColor];

     // mainColor = "lightgreen";

      var buf = new Float32Array(el.fLinePlexSize * 6);
      for (var i = 0; i < el.fLinePlexSize * 6; ++i)
         buf[i] = rnr_data.vtxBuff[i];
      var lineMaterial = new THREE.LineBasicMaterial({ color: mainColor, linewidth: 1 });

      var geom = new THREE.BufferGeometry();
      geom.addAttribute( 'position', new THREE.BufferAttribute( buf, 3 ) );
      var line = new THREE.LineSegments(geom, lineMaterial);
      obj3d.add(line);

      line.get_ctrl = function() { return new StraightLineSetControl(this.parent, true); }

      line.geo_object = el.fElementId;
      line.geo_name = el.fName || "linesset";

      var msize = el.fMarkerPlexSize;
      var pnts = new JSROOT.Painter.PointsCreator(msize, true, 3);

      var startIdx =el.fLinePlexSize * 6;
      var endIdx = startIdx + msize * 3;
      for (var i = startIdx; i < endIdx; i+=3) {
         pnts.AddPoint(rnr_data.vtxBuff[i], rnr_data.vtxBuff[i+1], rnr_data.vtxBuff[i+2] );
      }
      var marker_mesh = pnts.CreatePoints(mainColor);

      obj3d.add(marker_mesh);

      obj3d.eve_el = el;
      obj3d.eve_indx = rnr_data.idxBuff;
      obj3d.get_ctrl = function() { return new StraightLineSetControl(this); }

      return obj3d;
   }

   //==============================================================================

   JSROOT.EVE.EveElements = EveElements;

   return EveElements;

});
