sap.ui.define(['rootui5/eve7/lib/EveManager'], function(EveManager) {

   "use strict";

   // See also EveScene.js makeGLRepresentation(), there several members are
   // set for the top-level Object3D.

   //==============================================================================
   // EveElemControl
   //==============================================================================

   function EveElemControl(o3d)
   {
      // JSROOT.Painter.GeoDrawingControl.call(this);
      this.obj3d = o3d;
   }

   EveElemControl.prototype = Object.create(JSROOT.Painter.GeoDrawingControl.prototype);

   EveElemControl.prototype.invokeSceneMethod = function(fname, arg)
   {
      if ( ! this.obj3d) return false;

      var s = this.obj3d.scene;
      if (s && (typeof s[fname] == "function"))
         return s[fname](this.obj3d, arg, this.event);
      return false;
   }

   EveElemControl.prototype.separateDraw = false;

   EveElemControl.prototype.elementHighlighted = function(indx)
   {
      // default is simple selection, we ignore the indx
      this.invokeSceneMethod("processElementHighlighted"); // , indx);
   }

   EveElemControl.prototype.elementSelected = function(indx)
   {
      // default is simple selection, we ignore the indx
      this.invokeSceneMethod("processElementSelected"); //, indx);
   }


   //==============================================================================
   // EveElements
   //==============================================================================

   var GL = { POINTS: 0, LINES: 1, LINE_LOOP: 2, LINE_STRIP: 3, TRIANGLES: 4 };

   function EveElements(glc)
   {
      console.log("EveElements -- RCore", glc);

      this.glc = glc;
   }

   EveElements.prototype.RC = function()
   {
      return this.glc.RCore;
   }

   //==============================================================================
   // makeEveGeometry / makeEveGeoShape
   //==============================================================================

   EveElements.prototype.makeEveGeometry = function(rnr_data, force)
   {
      let RC = this.glc.RCore;

      var nVert = rnr_data.idxBuff[1]*3;

      if (rnr_data.idxBuff[0] != GL.TRIANGLES)  throw "Expect triangles first.";
      if (2 + nVert != rnr_data.idxBuff.length) throw "Expect single list of triangles in index buffer.";

      if (this.useIndexAsIs) {
         var body = new RC.BufferGeometry();
         body.addAttribute('position', new RC.BufferAttribute( rnr_data.vtxBuff, 3 ));
         body.setIndex(new RC.BufferAttribute( rnr_data.idxBuff, 1 ));
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

      var body = new RC.Geometry();

      body.vertices = new RC.BufferAttribute( vBuf, 3 );

      if (nBuf)
         body.normals = new RC.BufferAttribute( nBuf, 3 );
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
      let RC = this.glc.RCore;

      // var egs_ro = new RC.Object3D();
      var egs_ro = new RC.Group();

      var geom = this.makeEveGeometry(rnr_data);

      var fcol = new RC.Color(JSROOT.Painter.getColor(egs.fFillColor));

      // var material = new RC.MeshPhongMaterial({// side: THREE.DoubleSide,
      //                     depthWrite: false, color:fcol, transparent: true, opacity: 0.2 });
      var material = new RC.MeshPhongMaterial;
      material.side = 2;
      material.depthWrite = false;
      material.color = fcol;

      var mesh = new RC.Mesh(geom, material);

      egs_ro.add(mesh);

      return egs_ro;
   }

   //==============================================================================

   EveElements.prototype.makePolygonSetProjected = function(psp, rnr_data)
   {
      let RC = this.glc.RCore;

      if (this.useIndexAsIs)
         return this.makePolygonSetProjectedOld(psp, rnr_data);

      var psp_ro = new RC.Group(),
          ib_len = rnr_data.idxBuff.length,
          fcol   = new RC.Color(JSROOT.Painter.getColor(psp.fMainColor));

      for (var ib_pos = 0; ib_pos < ib_len; )
      {
         if (rnr_data.idxBuff[ib_pos] == GL.TRIANGLES)
         {
            var nVert = rnr_data.idxBuff[ib_pos + 1] * 3,
                vBuf  = new Float32Array(nVert*3); // plain buffer with all vertices

            for (var k=0;k<nVert;++k)
            {
               var pos = rnr_data.idxBuff[ib_pos+2+k];
               if (pos*3 > rnr_data.vtxBuff.length) { vBuf = null; break; }
               vBuf[k*3] = rnr_data.vtxBuff[pos*3];
               vBuf[k*3+1] = rnr_data.vtxBuff[pos*3+1];
               vBuf[k*3+2] = rnr_data.vtxBuff[pos*3+2];
            }

            if (vBuf)
            {
               var body = new RC.Geometry();
               body.vertices = new RC.BufferAttribute( vBuf, 3 );
               body.computeVertexNormals();

               var material = new RC.MeshBasicMaterial;
               material.side = 2;
               material.depthWrite = false;
               material.color = fcol;
               material.transparent = true;
               material.opacity = 0.4;

               psp_ro.add( new RC.Mesh(body, material) );
            }
            else
            {
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

            if (vBuf)
            {
               let body = new RC.Geometry();
               body.vertices = new RC.BufferAttribute( vBuf, 3 );
               let line_mat = new RC.MeshBasicMaterial;
               line_mat.color = fcol;
               // XXXX var line_mat = new RC.LineBasicMaterial({color:fcol });
               // XXXX psp_ro.add( new RC.LineLoop(body, line_mat) );
               psp_ro.add( new RC.Line(body, line_mat) );
            } else
            {
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

   //==============================================================================

   return EveElements;

});
