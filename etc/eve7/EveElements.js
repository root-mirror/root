/** @file EveElements.js */

(function( factory ) {
   if ( typeof define === "function" && define.amd ) {
      define( [ 'JSRootCore', 'threejs' ], factory );
   } else if (typeof exports === 'object' && typeof module !== 'undefined') {
      factory(require("./JSRootCore.js"), require("./three.min.js"));
   } else {

       if (typeof JSROOT == 'undefined')
          throw new Error('JSROOT is not defined', 'EveElements.js');

       if (typeof JSROOT.EVE == 'undefined')
          throw new Error('JSROOT.EVE is not defined', 'EveElements.js');

       if (typeof THREE == 'undefined')
          throw new Error('THREE is not defined', 'EveElements.js');

       factory(JSROOT, THREE);
    }
} (function( JSROOT, THREE ) {

   "use strict";

   var GL = { POINTS: 0, LINES: 1, LINE_LOOP: 2, LINE_STRIP: 3, TRIANGLES: 4 };

   // var flat_material = new THREE.ShaderMaterial( {
   //     uniforms: { time: { value: 1.0 }, resolution: { value: new THREE.Vector2() } },
   //     vertexShader: document.getElementById( 'vertexShader' ).textContent,
   //     fragmentShader: document.getElementById( 'fragmentShader' ).textContent
   // } );

   function EveElements()
   {
   }

   EveElements.prototype.makeHit = function(hit, rnrData) {
      // console.log("drawiHt ", hit, "this type ", this.viewType);
      // console.log("marker size ", hit.fMarkerSize)
      var hit_size = 8*rnrData.fMarkerSize,
          size = rnrData.vtxBuff.length/3,
          pnts = new JSROOT.Painter.PointsCreator(size, true, hit_size);

      for (var i=0;i<size;i++) {
         pnts.AddPoint(rnrData.vtxBuff[i*3],rnrData.vtxBuff[i*3+1],rnrData.vtxBuff[i*3+2]);
         // console.log("add vertex ", rnrData.vtxBuff[i*3],rnrData.vtxBuff[i*3+1],rnrData.vtxBuff[i*3+2]);
      }
      var mesh = pnts.CreatePoints(JSROOT.Painter.root_colors[hit.fMarkerColor]);

      mesh.highlightScale = 2;
      
      mesh.object = hit;
      mesh.geo_name = hit.fName;
      mesh.geo_object = hit.fMasterId || hit.fElementId;

      mesh.visible = hit.fRnrSelf;
      mesh.material.sizeAttenuation = false;
      mesh.material.size = hit.fMarkerSize;
      return mesh;
   }

   EveElements.prototype.makeTrack = function(track, rnrData) {
      var N = rnrData.vtxBuff.length/3;
      var track_width = track.fLineWidth || 1,
          track_color = JSROOT.Painter.root_colors[track.fLineColor] || "rgb(255,0,255)";
      if (JSROOT.browser.isWin) track_width = 1;  // not supported on windows

      var buf = new Float32Array((N-1)*6), pos = 0;
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
      line.visible = track.fRnrSelf;
      line.hightlightWidthScale = 2;

      // console.log("make track ", track, line.visible);
      return line;
   }

   EveElements.prototype.makeJet = function(jet, rnrData)
   {
      // console.log("make jet ", jet);
      var jet_ro = new THREE.Object3D();
      //var geo = new EveJetConeGeometry(jet.geoBuff);
      var pos_ba = new THREE.BufferAttribute( rnrData.vtxBuff, 3 );
      var N      = rnrData.vtxBuff.length / 3;

      var geo_body = new THREE.BufferGeometry();
      geo_body.addAttribute('position', pos_ba);
      var idcs = [];
      idcs.push( 0, N-1, 1 );
      for (var i = 1; i < N - 1; ++i)
         idcs.push( 0, i, i + 1 );
      geo_body.setIndex( idcs );
      geo_body.computeVertexNormals();
      
      var geo_rim = new THREE.BufferGeometry();
      geo_rim.addAttribute('position', pos_ba);
      idcs = [];
      for (var i = 1; i < N; ++i)
         idcs.push( i );
      geo_rim.setIndex( idcs );

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
      
      jet_ro.add( mesh  );
      jet_ro.add( line1 );
      jet_ro.add( line2 );

      mesh.object = jet_ro.object = jet;
      mesh.geo_name = jet_ro.geo_name = jet.fName;
      mesh.geo_object = jet_ro.geo_object = jet.fMasterId || jet.fElementId;
      jet_ro.visible = jet.fRnrSelf;
      
      // redirect highlight to the mesh
      jet_ro.geo_highlight = line1.geo_highlight = line2.geo_highlight = mesh;

      return jet_ro;
   }

   EveElements.prototype.makeJetProjected = function(jet, rnrData)
   {
      // JetProjected has 3 or 4 points. 0-th is apex, others are rim.
      // Fourth point is only present in RhoZ when jet hits barrel/endcap transition.

      // console.log("makeJetProjected ", jet);

      var jet_ro = new THREE.Object3D();
      var pos_ba = new THREE.BufferAttribute( rnrData.vtxBuff, 3 );
      var N      = rnrData.vtxBuff.length / 3;

      var geo_body = new THREE.BufferGeometry();
      geo_body.addAttribute('position', pos_ba);
      var idcs = [];
      idcs.push( 0, 2, 1 );
      if (N > 3) 
         idcs.push( 0, 3, 2 );
      geo_body.setIndex( idcs );
      geo_body.computeVertexNormals();
      
      var geo_rim = new THREE.BufferGeometry();
      geo_rim.addAttribute('position', pos_ba);
      idcs = [];
      for (var i = 1; i < N; ++i) idcs.push( i );
      geo_rim.setIndex( idcs );
      
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
      
      jet_ro.add( mesh  );
      jet_ro.add( line1 );
      jet_ro.add( line2 );

      mesh.object = jet_ro.object = jet;
      mesh.geo_name = jet_ro.geo_name = jet.fName;
      mesh.geo_object = jet_ro.geo_object = jet.fMasterId || jet.fElementId;
      
      // redirect highlight to the line1
      jet_ro.geo_highlight = line1.geo_highlight = line2.geo_highlight = mesh;
      
      jet_ro.visible = jet.fRnrSelf;

      return jet_ro;
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

   EveElements.prototype.makePolygonSetProjected = function(psp, rnr_data)
   {
      // console.log("makePolygonSetProjected ", psp);

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

   JSROOT.EVE.EveElements = EveElements;

   return JSROOT;

}));
