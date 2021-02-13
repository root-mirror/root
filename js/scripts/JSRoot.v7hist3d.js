/// @file JSRoot.v7hist3d.js
/// ROOT7 histogram 3D graphics

JSROOT.define(['d3', 'base3d', 'painter', 'v7hist'], (d3, THREE, jsrp) => {

   "use strict";

   /** @summary Text 3d axis visibility
     * @private */
   function testAxisVisibility(camera, toplevel, fb, bb) {
      let top;
      if (toplevel && toplevel.children)
         for (let n=0;n<toplevel.children.length;++n) {
            top = toplevel.children[n];
            if (top.axis_draw) break;
            top = undefined;
         }

      if (!top) return;

      if (!camera) {
         // this is case when axis drawing want to be removed
         toplevel.remove(top);
         return;
      }

      fb = fb ? true : false;
      bb = bb ? true : false;

      let qudrant = 1, pos = camera.position;
      if ((pos.x < 0) && (pos.y >= 0)) qudrant = 2;
      if ((pos.x >= 0) && (pos.y >= 0)) qudrant = 3;
      if ((pos.x >= 0) && (pos.y < 0)) qudrant = 4;

      function testvisible(id, range) {
         if (id <= qudrant) id+=4;
         return (id > qudrant) && (id < qudrant+range);
      }

      for (let n=0;n<top.children.length;++n) {
         let chld = top.children[n];
         if (chld.grid) chld.visible = bb && testvisible(chld.grid, 3); else
         if (chld.zid) chld.visible = testvisible(chld.zid, 2); else
         if (chld.xyid) chld.visible = testvisible(chld.xyid, 3); else
         if (chld.xyboxid) {
            let range = 5, shift = 0;
            if (bb && !fb) { range = 3; shift = -2; } else
            if (fb && !bb) range = 3; else
            if (!fb && !bb) range = (chld.bottom ? 3 : 0);
            chld.visible = testvisible(chld.xyboxid + shift, range);
            if (!chld.visible && chld.bottom && bb)
               chld.visible = testvisible(chld.xyboxid, 3);
         } else
         if (chld.zboxid) {
            let range = 2, shift = 0;
            if (fb && bb) range = 5; else
            if (bb && !fb) range = 4; else
            if (!bb && fb) { shift = -2; range = 4; }
            chld.visible = testvisible(chld.zboxid + shift, range);
         }
      }
   }

   /** @summary Set default camera position
     * @private */
   function setCameraPosition(fp) {
      let max3d = Math.max(0.75*fp.size_xy3d, fp.size_z3d);

      fp.camera.position.set(-1.6*max3d, -3.5*max3d, 1.4*fp.size_z3d);
      fp.camera.lookAt(fp.lookat);
   }

   /** @summary Create all necessary components for 3D drawings
     * @private */
   JSROOT.v7.RFramePainter.prototype.create3DScene = function(render3d) {

      if (render3d === -1) {

         if (!this.mode3d) return;

         if (!this.clear3dCanvas) {
            console.error('Strange, why mode3d is configured!!!!', this.mode3d);
            return;
         }

         testAxisVisibility(null, this.toplevel);

         this.clear3dCanvas();

         jsrp.disposeThreejsObject(this.scene);
         if (this.control) this.control.cleanup();

         if (this.renderer) {
            if (this.renderer.dispose) this.renderer.dispose();
            if (this.renderer.forceContextLoss) this.renderer.forceContextLoss();
         }

         delete this.size_xy3d;
         delete this.size_z3d;
         delete this.tooltip_mesh;
         delete this.scene;
         delete this.toplevel;
         delete this.camera;
         delete this.pointLight;
         delete this.renderer;
         delete this.control;
         if (this.render_tmout) {
            clearTimeout(this.render_tmout);
            delete this.render_tmout;
         }

         this.mode3d = false;

         return;
      }

      this.mode3d = true; // indicate 3d mode as hist painter does

      if ('toplevel' in this) {
         // it is indication that all 3D object created, just replace it with empty
         this.scene.remove(this.toplevel);
         jsrp.disposeThreejsObject(this.toplevel);
         delete this.tooltip_mesh;
         delete this.toplevel;
         if (this.control) this.control.HideTooltip();

         let newtop = new THREE.Object3D();
         this.scene.add(newtop);
         this.toplevel = newtop;

         this.resize3D(); // set actual sizes
         return;
      }

      jsrp.assign3DHandler(this);

      render3d = jsrp.getRender3DKind(render3d);
      let sz = this.getSizeFor3d(undefined, render3d);

      this.size_z3d = 100;
      this.size_xy3d = (sz.height > 10) && (sz.width > 10) ? Math.round(sz.width/sz.height*this.size_z3d) : this.size_z3d;

      // three.js 3D drawing
      this.scene = new THREE.Scene();
      //scene.fog = new THREE.Fog(0xffffff, 500, 3000);

      this.toplevel = new THREE.Object3D();
      this.scene.add(this.toplevel);
      this.scene_width = sz.width;
      this.scene_height = sz.height;

      this.camera = new THREE.PerspectiveCamera(45, this.scene_width / this.scene_height, 1, 40*this.size_z3d);

      this.camera_Phi = 30;
      this.camera_Theta = 30;

      this.pointLight = new THREE.PointLight(0xffffff,1);
      this.camera.add(this.pointLight);
      this.pointLight.position.set(this.size_xy3d/2, this.size_xy3d/2, this.size_z3d/2);
      this.lookat = new THREE.Vector3(0,0,0.8*this.size_z3d);
      this.camera.up = new THREE.Vector3(0,0,1);
      this.scene.add( this.camera );

      setCameraPosition(this);

      this.renderer = jsrp.createRender3D(this.scene_width, this.scene_height, render3d);

      this.webgl = (render3d === JSROOT.constants.Render3D.WebGL);
      this.add3dCanvas(sz, this.renderer.jsroot_dom, this.webgl);

      this.first_render_tm = 0;
      this.enable_highlight = false;

      if (JSROOT.batch_mode || !this.webgl) return;

      this.control = jsrp.createOrbitControl(this, this.camera, this.scene, this.renderer, this.lookat);

      let axis_painter = this, obj_painter = this.getMainPainter();

      this.control.ProcessMouseMove = function(intersects) {
         let tip = null, mesh = null, zoom_mesh = null;

         for (let i = 0; i < intersects.length; ++i) {
            if (intersects[i].object.tooltip) {
               tip = intersects[i].object.tooltip(intersects[i]);
               if (tip) { mesh = intersects[i].object; break; }
            } else if (intersects[i].object.zoom && !zoom_mesh) {
               zoom_mesh = intersects[i].object;
            }
         }

         if (tip && !tip.use_itself) {
            let delta_xy = 1e-4*axis_painter.size_xy3d, delta_z = 1e-4*axis_painter.size_z3d;
            if ((tip.x1 > tip.x2) || (tip.y1 > tip.y2) || (tip.z1 > tip.z2)) console.warn('check 3D hints coordinates');
            tip.x1 -= delta_xy; tip.x2 += delta_xy;
            tip.y1 -= delta_xy; tip.y2 += delta_xy;
            tip.z1 -= delta_z; tip.z2 += delta_z;
         }

         axis_painter.highlightBin3D(tip, mesh);

         if (!tip && zoom_mesh && axis_painter.Get3DZoomCoord) {
            let pnt = zoom_mesh.GlobalIntersect(this.raycaster),
                axis_name = zoom_mesh.zoom,
                axis_value = axis_painter.Get3DZoomCoord(pnt, axis_name);

            if ((axis_name==="z") && zoom_mesh.use_y_for_z) axis_name = "y";

            let hint = { name: axis_name,
                         title: "RAxis",
                         line: "any info",
                         only_status: true };

            hint.line = axis_name + " : " + axis_painter.axisAsText(axis_name, axis_value);

            return hint;
         }

         return (tip && tip.lines) ? tip : "";
      }

      this.control.ProcessMouseLeave = function() {
         axis_painter.highlightBin3D(null);
      }

      this.control.contextMenu = function(pos, intersects) {
         let kind = "painter", p = obj_painter;
         if (intersects)
            for (let n=0;n<intersects.length;++n) {
               let mesh = intersects[n].object;
               if (mesh.zoom) { kind = mesh.zoom; p = null; break; }
               if (mesh.painter && typeof mesh.painter.fillContextMenu === 'function') {
                  p = mesh.painter; break;
               }
            }

         let fp = obj_painter.getFramePainter();
         if (fp && fp.showContextMenu)
            fp.showContextMenu(kind, pos, p);
      }
   }

   /** @summary call 3D rendering of the histogram drawing
     * @desc tmout specified delay, after which actual rendering will be invoked
     * Timeout used to avoid multiple rendering of the picture when several 3D drawings
     * superimposed with each other.
     * If tmeout <= 0, rendering performed immediately
     * If (tmout == -1111), immediate rendering with SVG renderer is performed
     * @private */
   JSROOT.v7.RFramePainter.prototype.render3D = function(tmout) {

      if (tmout === -1111) {
         // special handling for direct SVG renderer
         // probably, here one can use canvas renderer - after modifications
         // let rrr = new THREE.SVGRenderer({ precision: 0, astext: true });
         let doc = JSROOT._.get_document(),
             rrr = THREE.CreateSVGRenderer(false, 0, doc);
         rrr.setSize(this.scene_width, this.scene_height);
         rrr.render(this.scene, this.camera);
         if (rrr.makeOuterHTML) {
            // use text mode, it is faster
            let d = doc.createElement('div');
            d.innerHTML = rrr.makeOuterHTML();
            return d.childNodes[0];
         }
         return rrr.domElement;
      }

      if (tmout === undefined) tmout = 5; // by default, rendering happens with timeout

      if ((tmout > 0) && this.webgl && JSROOT.batch_mode) {
          if (!this.render_tmout)
             this.render_tmout = setTimeout(() => this.render3D(0), tmout);
          return;
       }

      if (this.render_tmout) {
         clearTimeout(this.render_tmout);
         delete this.render_tmout;
      }

      if (!this.renderer) return;

      jsrp.beforeRender3D(this.renderer);

      let tm1 = new Date();

      if (!this.opt3d) this.opt3d = { FrontBox: true, BackBox: true };

      testAxisVisibility(this.camera, this.toplevel, this.opt3d.FrontBox, this.opt3d.BackBox);

      // do rendering, most consuming time
      this.renderer.render(this.scene, this.camera);

      jsrp.afterRender3D(this.renderer);

      let tm2 = new Date();

      if (this.first_render_tm === 0) {
         this.first_render_tm = tm2.getTime() - tm1.getTime();
         this.enable_highlight = (this.first_render_tm < 1200) && this.isTooltipAllowed();
         console.log('three.js r' + THREE.REVISION + ', first render tm = ' + this.first_render_tm);
      }
   }

   /** @summary Check is 3D drawing need to be resized
     * @private */
   JSROOT.v7.RFramePainter.prototype.resize3D = function() {

      let sz = this.getSizeFor3d(this.access3dKind());

      this.apply3dSize(sz);

      if ((this.scene_width === sz.width) && (this.scene_height === sz.height)) return false;

      if ((sz.width<10) || (sz.height<10)) return false;

      // TODO: change xy/z ratio after canvas resize
      // this.size_xy3d = Math.round(sz.width/sz.height*this.size_z3d);

      this.scene_width = sz.width;
      this.scene_height = sz.height;

      this.camera.aspect = this.scene_width / this.scene_height;
      this.camera.updateProjectionMatrix();

      this.renderer.setSize( this.scene_width, this.scene_height );
      if (this.renderer.setJSROOTSize)
         this.renderer.setJSROOTSize(this.scene_width, this.scene_height);

      return true;
   }

   /** @summary Hilight bin in 3D drawing
     * @private */
   JSROOT.v7.RFramePainter.prototype.highlightBin3D = function(tip, selfmesh) {

      let changed = false, tooltip_mesh = null, changed_self = true,
          want_remove = !tip || (tip.x1===undefined) || !this.enable_highlight,
          mainp = this.getMainPainter();

      if (mainp && (!mainp.provideUserTooltip || !mainp.hasUserTooltip())) mainp = null;

      if (this.tooltip_selfmesh) {
         changed_self = (this.tooltip_selfmesh !== selfmesh)
         this.tooltip_selfmesh.material.color = this.tooltip_selfmesh.save_color;
         delete this.tooltip_selfmesh;
         changed = true;
      }

      if (this.tooltip_mesh) {
         tooltip_mesh = this.tooltip_mesh;
         this.toplevel.remove(this.tooltip_mesh);
         delete this.tooltip_mesh;
         changed = true;
      }

      if (want_remove) {
         if (changed) this.render3D();
         if (changed && mainp) mainp.provideUserTooltip(null);
         return;
      }

      if (tip.use_itself) {
         selfmesh.save_color = selfmesh.material.color;
         selfmesh.material.color = new THREE.Color(tip.color);
         this.tooltip_selfmesh = selfmesh;
         changed = changed_self;
      } else {
         changed = true;

         let indicies = jsrp.Box3D.Indexes,
             normals = jsrp.Box3D.Normals,
             vertices = jsrp.Box3D.Vertices,
             pos, norm,
             color = new THREE.Color(tip.color ? tip.color : 0xFF0000),
             opacity = tip.opacity || 1;

         if (!tooltip_mesh) {
            pos = new Float32Array(indicies.length*3);
            norm = new Float32Array(indicies.length*3);
            let geom = new THREE.BufferGeometry();
            geom.setAttribute( 'position', new THREE.BufferAttribute( pos, 3 ) );
            geom.setAttribute( 'normal', new THREE.BufferAttribute( norm, 3 ) );
            let mater = new THREE.MeshBasicMaterial({ color: color, opacity: opacity, flatShading: true });
            tooltip_mesh = new THREE.Mesh(geom, mater);
         } else {
            pos = tooltip_mesh.geometry.attributes.position.array;
            tooltip_mesh.geometry.attributes.position.needsUpdate = true;
            tooltip_mesh.material.color = color;
            tooltip_mesh.material.opacity = opacity;
         }

         if (tip.x1 === tip.x2) console.warn('same tip X', tip.x1, tip.x2);
         if (tip.y1 === tip.y2) console.warn('same tip Y', tip.y1, tip.y2);
         if (tip.z1 === tip.z2) { tip.z2 = tip.z1 + 0.0001; } // avoid zero faces

         for (let k=0,nn=-3;k<indicies.length;++k) {
            let vert = vertices[indicies[k]];
            pos[k*3]   = tip.x1 + vert.x * (tip.x2 - tip.x1);
            pos[k*3+1] = tip.y1 + vert.y * (tip.y2 - tip.y1);
            pos[k*3+2] = tip.z1 + vert.z * (tip.z2 - tip.z1);

            if (norm) {
               if (k%6===0) nn+=3;
               norm[k*3] = normals[nn];
               norm[k*3+1] = normals[nn+1];
               norm[k*3+2] = normals[nn+2];
            }
         }
         this.tooltip_mesh = tooltip_mesh;
         this.toplevel.add(tooltip_mesh);
      }

      if (changed) this.render3D();

      if (changed && tip.$painter && (typeof tip.$painter.redrawProjection == 'function'))
         tip.$painter.redrawProjection(tip.ix-1, tip.ix, tip.iy-1, tip.iy);

      if (changed && mainp && mainp.getObject())
         mainp.provideUserTooltip({ obj: mainp.getObject(),  name: mainp.getObject().fName,
                                    bin: tip.bin, cont: tip.value,
                                    binx: tip.ix, biny: tip.iy, binz: tip.iz,
                                    grx: (tip.x1+tip.x2)/2, gry: (tip.y1+tip.y2)/2, grz: (tip.z1+tip.z2)/2 });
   }

   /** @summary Set options used for 3D drawings
     * @private */
   JSROOT.v7.RFramePainter.prototype.set3DOptions = function(hopt) {
      this.opt3d = hopt;
   }

   /** @summary Draw axes in 3D mode
     * @private */
   JSROOT.v7.RFramePainter.prototype.drawXYZ = function(toplevel, opts) {
      if (!opts) opts = {};

      let grminx = -this.size_xy3d, grmaxx = this.size_xy3d,
          grminy = -this.size_xy3d, grmaxy = this.size_xy3d,
          grminz = 0, grmaxz = 2*this.size_z3d,
          textsize = Math.round(this.size_z3d * 0.05),
          xmin = this.xmin, xmax = this.xmax,
          ymin = this.ymin, ymax = this.ymax,
          zmin = this.zmin, zmax = this.zmax,
          y_zoomed = false, z_zoomed = false;

      if (!this.size_z3d) {
         grminx = this.xmin; grmaxx = this.xmax;
         grminy = this.ymin; grmaxy = this.ymax;
         grminz = this.zmin; grmaxz = this.zmax;
         textsize = (grmaxz - grminz) * 0.05;
      }

      if (('zoom_xmin' in this) && ('zoom_xmax' in this) && (this.zoom_xmin !== this.zoom_xmax)) {
         xmin = this.zoom_xmin; xmax = this.zoom_xmax;
      }
      if (('zoom_ymin' in this) && ('zoom_ymax' in this) && (this.zoom_ymin !== this.zoom_ymax)) {
         ymin = this.zoom_ymin; ymax = this.zoom_ymax; y_zoomed = true;
      }
      if (('zoom_zmin' in this) && ('zoom_zmax' in this) && (this.zoom_zmin !== this.zoom_zmax)) {
         zmin = this.zoom_zmin; zmax = this.zoom_zmax; z_zoomed = true;
      }

      if (opts.use_y_for_z) {
         this.zmin = this.ymin; this.zmax = this.ymax;
         zmin = ymin; zmax = ymax; z_zoomed = y_zoomed;
         // if (!z_zoomed && (this.hmin!==this.hmax)) { zmin = this.hmin; zmax = this.hmax; }
         ymin = 0; ymax = 1;
      }

      // z axis range used for lego plot
      this.lego_zmin = zmin; this.lego_zmax = zmax;

      // factor 1.1 used in ROOT for lego plots
      if ((opts.zmult !== undefined) && !z_zoomed) zmax *= opts.zmult;

      this.x_handle = new JSROOT.v7.RAxisPainter(this.getDom(), this, this.xaxis, "x_");
      this.x_handle.setPadName(this.getPadName());
      this.x_handle.snapid = this.snapid;
      this.x_handle.configureAxis("xaxis", this.xmin, this.xmax, xmin, xmax, false, [grminx, grmaxx]);
      this.x_handle.assignFrameMembers(this,"x");

      this.y_handle = new JSROOT.v7.RAxisPainter(this.getDom(), this, this.yaxis, "y_");
      this.y_handle.setPadName(this.getPadName());
      this.y_handle.snapid = this.snapid;
      this.y_handle.configureAxis("yaxis", this.ymin, this.ymax, ymin, ymax, false, [grminy, grmaxy]);
      this.y_handle.assignFrameMembers(this,"y");

      // this.setRootPadRange(pad, true); // set some coordinates typical for 3D projections in ROOT

      this.z_handle = new JSROOT.v7.RAxisPainter(this.getDom(), this, this.zaxis, "z_");
      this.z_handle.setPadName(this.getPadName());
      this.z_handle.snapid = this.snapid;
      this.z_handle.configureAxis("zaxis", this.zmin, this.zmax, zmin, zmax, false, [grminz, grmaxz]);
      this.z_handle.assignFrameMembers(this,"z");

      let textMaterial = new THREE.MeshBasicMaterial({ color: 0x000000 }),
          lineMaterial = new THREE.LineBasicMaterial({ color: 0x000000 }),
          ticklen = textsize*0.5, lbls = [], text_scale = 1,
          xticks = this.x_handle.createTicks(false, true),
          yticks = this.y_handle.createTicks(false, true),
          zticks = this.z_handle.createTicks(false, true);

      // main element, where all axis elements are placed
      let top = new THREE.Object3D();
      top.axis_draw = true; // mark element as axis drawing
      toplevel.add(top);

      let ticks = [], maxtextheight = 0, xaxis = this.xaxis;

      while (xticks.next()) {
         let grx = xticks.grpos,
            is_major = (xticks.kind===1),
            lbl = this.x_handle.format(xticks.tick, 2);

         if (xticks.last_major()) {
            if (!xaxis || !xaxis.fTitle) lbl = "x";
         } else if (lbl === null) {
            is_major = false; lbl = "";
         }

         if (is_major && lbl && (lbl.length>0)) {
            let text3d = new THREE.TextGeometry(lbl, { font: JSROOT.threejs_font_helvetiker_regular, size: textsize, height: 0, curveSegments: 5 });
            text3d.computeBoundingBox();
            let draw_width = text3d.boundingBox.max.x - text3d.boundingBox.min.x,
                draw_height = text3d.boundingBox.max.y - text3d.boundingBox.min.y;
            text3d.center = true; // place central

            // text3d.translate(-draw_width/2, 0, 0);

            maxtextheight = Math.max(maxtextheight, draw_height);

            text3d.grx = grx;
            lbls.push(text3d);

            if (!xticks.last_major()) {
               let space = (xticks.next_major_grpos() - grx);
               if (draw_width > 0)
                  text_scale = Math.min(text_scale, 0.9*space/draw_width)
               if (this.x_handle.isCenteredLabels()) text3d.grx += space/2;
            }
         }

         ticks.push(grx, 0, 0, grx, (is_major ? -ticklen : -ticklen * 0.6), 0);
      }

      if (xaxis && xaxis.fTitle) {
         let text3d = new THREE.TextGeometry(xaxis.fTitle, { font: JSROOT.threejs_font_helvetiker_regular, size: textsize, height: 0, curveSegments: 5 });
         text3d.computeBoundingBox();
         text3d.center = false; // xaxis.TestBit(JSROOT.EAxisBits.kCenterTitle);
         text3d.gry = 2; // factor 2 shift
         text3d.grx = (grminx + grmaxx)/2; // default position for centered title
         lbls.push(text3d);
      }

      this.Get3DZoomCoord = function(point, kind) {
         // return axis coordinate from intersection point with axis geometry
         let pos = point[kind], min = this['scale_'+kind+'min'], max = this['scale_'+kind+'max'];

         if (kind==="z") pos = pos/2/this.size_z3d;
                   else  pos = (pos+this.size_xy3d)/2/this.size_xy3d;

         if (this["log"+kind]) {
            pos = Math.exp(Math.log(min) + pos*(Math.log(max)-Math.log(min)));
         } else {
            pos = min + pos*(max-min);
         }
         return pos;
      }

      function CreateZoomMesh(kind, size_3d, use_y_for_z) {
         let geom = new THREE.Geometry();

         if (kind==="z")
            geom.vertices.push(
                  new THREE.Vector3(0,0,0),
                  new THREE.Vector3(ticklen*4, 0, 0),
                  new THREE.Vector3(ticklen*4, 0, 2*size_3d),
                  new THREE.Vector3(0, 0, 2*size_3d));
         else
            geom.vertices.push(
                  new THREE.Vector3(-size_3d,0,0),
                  new THREE.Vector3(size_3d,0,0),
                  new THREE.Vector3(size_3d,-ticklen*4,0),
                  new THREE.Vector3(-size_3d,-ticklen*4,0));

         geom.faces.push(new THREE.Face3(0, 2, 1));
         geom.faces.push(new THREE.Face3(0, 3, 2));
         geom.computeFaceNormals();

         let material = new THREE.MeshBasicMaterial({ transparent: true,
                                   vertexColors: THREE.NoColors, //   THREE.FaceColors,
                                   side: THREE.DoubleSide,
                                   opacity: 0 });

         let mesh = new THREE.Mesh(geom, material);
         mesh.zoom = kind;
         mesh.size_3d = size_3d;
         mesh.use_y_for_z = use_y_for_z;
         if (kind=="y") mesh.rotateZ(Math.PI/2).rotateX(Math.PI);

         mesh.GlobalIntersect = function(raycaster) {
            let plane = new THREE.Plane(),
                geom = this.geometry;

            if (!geom || !geom.vertices) return undefined;

            plane.setFromCoplanarPoints(geom.vertices[0], geom.vertices[1], geom.vertices[2]);
            plane.applyMatrix4(this.matrixWorld);

            let v1 = raycaster.ray.origin.clone(),
                v2 = v1.clone().addScaledVector(raycaster.ray.direction, 1e10);

            let pnt = plane.intersectLine(new THREE.Line3(v1,v2), new THREE.Vector3());

            if (!pnt) return undefined;

            let min = -this.size_3d, max = this.size_3d;
            if (this.zoom==="z") { min = 0; max = 2*this.size_3d; }

            if (pnt[this.zoom] < min) pnt[this.zoom] = min; else
            if (pnt[this.zoom] > max) pnt[this.zoom] = max;

            return pnt;
         }

         mesh.ShowSelection = function(pnt1,pnt2) {
            // used to show selection

            let tgtmesh = this.children ? this.children[0] : null, gg, kind = this.zoom;
            if (!pnt1 || !pnt2) {
               if (tgtmesh) {
                  this.remove(tgtmesh)
                  jsrp.disposeThreejsObject(tgtmesh);
               }
               return tgtmesh;
            }

            if (!tgtmesh) {
               gg = this.geometry.clone();
               if (kind==="z") gg.vertices[1].x = gg.vertices[2].x = ticklen;
                          else gg.vertices[2].y = gg.vertices[3].y = -ticklen;
               tgtmesh = new THREE.Mesh(gg, new THREE.MeshBasicMaterial({ color: 0xFF00, side: THREE.DoubleSide }));
               this.add(tgtmesh);
            } else {
               gg = tgtmesh.geometry;
            }

            if (kind=="z") {
               gg.vertices[0].z = gg.vertices[1].z = pnt1[kind];
               gg.vertices[2].z = gg.vertices[3].z = pnt2[kind];
            } else {
               gg.vertices[0].x = gg.vertices[3].x = pnt1[kind];
               gg.vertices[1].x = gg.vertices[2].x = pnt2[kind];
            }

            gg.computeFaceNormals();

            gg.verticesNeedUpdate = true;
            gg.normalsNeedUpdate = true;

            return true;
         }

         return mesh;
      }

      let xcont = new THREE.Object3D();
      xcont.position.set(0, grminy, grminz)
      xcont.rotation.x = 1/4*Math.PI;
      xcont.xyid = 2;
      let xtickslines = jsrp.createLineSegments( ticks, lineMaterial );
      xcont.add(xtickslines);

      lbls.forEach(lbl => {
         let w = lbl.boundingBox.max.x - lbl.boundingBox.min.x,
             posx = lbl.center ? lbl.grx - w/2 : grmaxx - w,
             m = new THREE.Matrix4();
         // matrix to swap y and z scales and shift along z to its position
         m.set(text_scale, 0,           0,  posx,
               0,          text_scale,  0,  (-maxtextheight*text_scale - 1.5*ticklen) * (lbl.gry || 1),
               0,          0,           1,  0,
               0,          0,           0,  1);

         let mesh = new THREE.Mesh(lbl, textMaterial);
         mesh.applyMatrix4(m);
         xcont.add(mesh);
      });

      if (opts.zoom) xcont.add(CreateZoomMesh("x", this.size_xy3d));
      top.add(xcont);

      xcont = new THREE.Object3D();
      xcont.position.set(0, grmaxy, grminz);
      xcont.rotation.x = 3/4*Math.PI;
      xcont.add(new THREE.LineSegments(xtickslines.geometry, lineMaterial));
      lbls.forEach(lbl => {

         let w = lbl.boundingBox.max.x - lbl.boundingBox.min.x,
             posx = lbl.center ? lbl.grx + w/2 : grmaxx,
             m = new THREE.Matrix4();
         // matrix to swap y and z scales and shift along z to its position
         m.set(-text_scale, 0,           0, posx,
               0,           text_scale,  0, (-maxtextheight*text_scale - 1.5*ticklen) * (lbl.gry || 1),
               0,           0,           -1, 0,
               0,            0,           0, 1);
         let mesh = new THREE.Mesh(lbl, textMaterial);
         mesh.applyMatrix4(m);
         xcont.add(mesh);
      });

      //xcont.add(new THREE.Mesh(ggg2, textMaterial));
      xcont.xyid = 4;
      if (opts.zoom) xcont.add(CreateZoomMesh("x", this.size_xy3d));
      top.add(xcont);

      lbls = []; text_scale = 1; maxtextheight = 0; ticks = [];

      let yaxis = this.yaxis;

      while (yticks.next()) {
         let gry = yticks.grpos,
             is_major = (yticks.kind===1),
             lbl = this.y_handle.format(yticks.tick, 2);

         if (yticks.last_major()) {
            if (!yaxis || !yaxis.fTitle) lbl = "y";
         }  else if (lbl === null) {
            is_major = false; lbl = "";
         }

         if (is_major) {
            let text3d = new THREE.TextGeometry(lbl, { font: JSROOT.threejs_font_helvetiker_regular, size: textsize, height: 0, curveSegments: 5 });
            text3d.computeBoundingBox();
            let draw_width = text3d.boundingBox.max.x - text3d.boundingBox.min.x,
                draw_height = text3d.boundingBox.max.y - text3d.boundingBox.min.y;
            // text3d.translate(-draw_width/2, 0, 0);
            text3d.center = true;

            maxtextheight = Math.max(maxtextheight, draw_height);

            text3d.gry = gry;
            lbls.push(text3d);

            if (!yticks.last_major()) {
               let space = (yticks.next_major_grpos() - gry);
               if (draw_width > 0)
                  text_scale = Math.min(text_scale, 0.9*space/draw_width)
               if (this.y_handle.isCenteredLabels()) text3d.gry += space/2;
            }
         }
         ticks.push(0,gry,0, (is_major ? -ticklen : -ticklen*0.6), gry, 0);
      }

      if (yaxis && yaxis.fTitle) {
         let text3d = new THREE.TextGeometry(yaxis.fTitle, { font: JSROOT.threejs_font_helvetiker_regular, size: textsize, height: 0, curveSegments: 5 });
         text3d.computeBoundingBox();
         text3d.center = false; // yaxis.TestBit(JSROOT.EAxisBits.kCenterTitle);
         text3d.grx = 2; // factor 2 shift
         text3d.gry = (grminy + grmaxy)/2; // default position for centered title
         lbls.push(text3d);
      }

      if (!opts.use_y_for_z) {
         let yticksline = jsrp.createLineSegments(ticks, lineMaterial),
             ycont = new THREE.Object3D();
         ycont.position.set(grminx, 0, grminz);
         ycont.rotation.y = -1/4*Math.PI;
         ycont.add(yticksline);
         //ycont.add(new THREE.Mesh(ggg1, textMaterial));

         lbls.forEach(lbl => {

            let w = lbl.boundingBox.max.x - lbl.boundingBox.min.x,
                posy = lbl.center ? lbl.gry + w/2 : grmaxy,
                m = new THREE.Matrix4();
            // matrix to swap y and z scales and shift along z to its position
            m.set(0, text_scale,  0, (-maxtextheight*text_scale - 1.5*ticklen)*(lbl.grx || 1),
                  -text_scale,  0, 0, posy,
                  0, 0,  1, 0,
                  0, 0,  0, 1);

            let mesh = new THREE.Mesh(lbl, textMaterial);
            mesh.applyMatrix4(m);
            ycont.add(mesh);
         });

         ycont.xyid = 3;
         if (opts.zoom) ycont.add(CreateZoomMesh("y", this.size_xy3d));
         top.add(ycont);

         ycont = new THREE.Object3D();
         ycont.position.set(grmaxx, 0, grminz);
         ycont.rotation.y = -3/4*Math.PI;
         ycont.add(new THREE.LineSegments(yticksline.geometry, lineMaterial));
         //ycont.add(new THREE.Mesh(ggg2, textMaterial));
         lbls.forEach(lbl => {
            let w = lbl.boundingBox.max.x - lbl.boundingBox.min.x,
                posy = lbl.center ? lbl.gry - w/2 : grmaxy - w,
                m = new THREE.Matrix4();
            m.set(0, text_scale, 0,  (-maxtextheight*text_scale - 1.5*ticklen)*(lbl.grx || 1),
                  text_scale, 0, 0,  posy,
                  0,         0, -1,  0,
                  0, 0, 0, 1);

            let mesh = new THREE.Mesh(lbl, textMaterial);
            mesh.applyMatrix4(m);
            ycont.add(mesh);
         });
         ycont.xyid = 1;
         if (opts.zoom) ycont.add(CreateZoomMesh("y", this.size_xy3d));
         top.add(ycont);
      }


      lbls = []; text_scale = 1;

      ticks = []; // just array, will be used for the buffer geometry

      let zgridx = null, zgridy = null, lastmajorz = null,
          zaxis = this.zaxis, maxzlblwidth = 0;

      if (this.size_z3d) {
         zgridx = []; zgridy = [];
      }

      while (zticks.next()) {
         let grz = zticks.grpos,
            is_major = (zticks.kind == 1),
            lbl = this.z_handle.format(zticks.tick, 2);

         if (lbl === null) { is_major = false; lbl = ""; }

         if (is_major && lbl) {
            let text3d = new THREE.TextGeometry(lbl, { font: JSROOT.threejs_font_helvetiker_regular, size: textsize, height: 0, curveSegments: 5 });
            text3d.computeBoundingBox();
            let draw_width = text3d.boundingBox.max.x - text3d.boundingBox.min.x,
                draw_height = text3d.boundingBox.max.y - text3d.boundingBox.min.y;
            text3d.translate(-draw_width, -draw_height/2, 0);
            text3d.grz = grz;
            lbls.push(text3d);

            if ((lastmajorz !== null) && (draw_height>0))
               text_scale = Math.min(text_scale, 0.9*(grz - lastmajorz)/draw_height);

            maxzlblwidth = Math.max(maxzlblwidth, draw_width);

            lastmajorz = grz;
         }

         // create grid
         if (zgridx && is_major)
            zgridx.push(grminx,0,grz, grmaxx,0,grz);

         if (zgridy && is_major)
            zgridy.push(0,grminy,grz, 0,grmaxy,grz);

         ticks.push(0, 0, grz, (is_major ? ticklen : ticklen * 0.6), 0, grz);
      }

      if (zgridx && (zgridx.length > 0)) {

         let material = new THREE.LineDashedMaterial({ color: 0x0, dashSize: 2, gapSize: 2 });

         let lines1 = jsrp.createLineSegments(zgridx, material);
         lines1.position.set(0,grmaxy,0);
         lines1.grid = 2; // mark as grid
         lines1.visible = false;
         top.add(lines1);

         let lines2 = new THREE.LineSegments(lines1.geometry, material);
         lines2.position.set(0,grminy,0);
         lines2.grid = 4; // mark as grid
         lines2.visible = false;
         top.add(lines2);
      }

      if (zgridy && (zgridy.length > 0)) {

         let material = new THREE.LineDashedMaterial({ color: 0x0, dashSize: 2, gapSize: 2 });

         let lines1 = jsrp.createLineSegments(zgridy, material);
         lines1.position.set(grmaxx,0, 0);
         lines1.grid = 3; // mark as grid
         lines1.visible = false;
         top.add(lines1);

         let lines2 = new THREE.LineSegments(lines1.geometry, material);
         lines2.position.set(grminx, 0, 0);
         lines2.grid = 1; // mark as grid
         lines2.visible = false;
         top.add(lines2);
      }

      let zcont = [], zticksline = jsrp.createLineSegments( ticks, lineMaterial );
      for (let n=0;n<4;++n) {
         zcont.push(new THREE.Object3D());

         lbls.forEach(lbl => {
            let m = new THREE.Matrix4();
            // matrix to swap y and z scales and shift along z to its position
            m.set(-text_scale,          0,  0, 2*ticklen,
                            0,          0,  1, 0,
                            0, text_scale,  0, lbl.grz);
            let mesh = new THREE.Mesh(lbl, textMaterial);
            mesh.applyMatrix4(m);
            zcont[n].add(mesh);
         });

         if (zaxis && zaxis.fTitle) {
            let text3d = new THREE.TextGeometry(zaxis.fTitle, { font: JSROOT.threejs_font_helvetiker_regular, size: textsize, height: 0, curveSegments: 5 });
            text3d.computeBoundingBox();
            let draw_width = text3d.boundingBox.max.x - text3d.boundingBox.min.x,
                // draw_height = text3d.boundingBox.max.y - text3d.boundingBox.min.y,
                center_title = false, // zaxis.TestBit(JSROOT.EAxisBits.kCenterTitle)
                posz = center_title ? (grmaxz + grminz - draw_width)/2 : grmaxz - draw_width;

            text3d.rotateZ(Math.PI/2);

            let m = new THREE.Matrix4();
            m.set(-text_scale,          0,  0, 3*ticklen + maxzlblwidth,
                            0,          0,  1, 0,
                            0, text_scale,  0, posz);
            let mesh = new THREE.Mesh(text3d, textMaterial);
            mesh.applyMatrix4(m);
            zcont[n].add(mesh);
         }

         zcont[n].add(n==0 ? zticksline : new THREE.LineSegments(zticksline.geometry, lineMaterial));
         if (opts.zoom) zcont[n].add(CreateZoomMesh("z", this.size_z3d, opts.use_y_for_z));

         zcont[n].zid = n + 2;
         top.add(zcont[n]);
      }

      zcont[0].position.set(grminx,grmaxy,0);
      zcont[0].rotation.z = 3/4*Math.PI;

      zcont[1].position.set(grmaxx,grmaxy,0);
      zcont[1].rotation.z = 1/4*Math.PI;

      zcont[2].position.set(grmaxx,grminy,0);
      zcont[2].rotation.z = -1/4*Math.PI;

      zcont[3].position.set(grminx,grminy,0);
      zcont[3].rotation.z = -3/4*Math.PI;


      // for TAxis3D do not show final cube
      if (this.size_z3d === 0) return;

      let linex_geom = jsrp.createLineSegments([grminx,0,0, grmaxx,0,0], lineMaterial, null, true);
      for(let n=0;n<2;++n) {
         let line = new THREE.LineSegments(linex_geom, lineMaterial);
         line.position.set(0, grminy, (n===0) ? grminz : grmaxz);
         line.xyboxid = 2; line.bottom = (n == 0);
         top.add(line);

         line = new THREE.LineSegments(linex_geom, lineMaterial);
         line.position.set(0, grmaxy, (n===0) ? grminz : grmaxz);
         line.xyboxid = 4; line.bottom = (n == 0);
         top.add(line);
      }

      let liney_geom = jsrp.createLineSegments([0,grminy,0, 0,grmaxy,0], lineMaterial, null, true);
      for(let n=0;n<2;++n) {
         let line = new THREE.LineSegments(liney_geom, lineMaterial);
         line.position.set(grminx, 0, (n===0) ? grminz : grmaxz);
         line.xyboxid = 3; line.bottom = (n == 0);
         top.add(line);

         line = new THREE.LineSegments(liney_geom, lineMaterial);
         line.position.set(grmaxx, 0, (n===0) ? grminz : grmaxz);
         line.xyboxid = 1; line.bottom = (n == 0);
         top.add(line);
      }

      let linez_geom = jsrp.createLineSegments([0,0,grminz, 0,0,grmaxz], lineMaterial, null, true);
      for(let n=0;n<4;++n) {
         let line = new THREE.LineSegments(linez_geom, lineMaterial);
         line.zboxid = zcont[n].zid;
         line.position.copy(zcont[n].position);
         top.add(line);
      }
   }

   /** @summary Draw 1D/2D histograms in Lego mode
     * @private */
   JSROOT.v7.RHistPainter.prototype.drawLego = function() {

      if (!this.draw_content) return;

      // Perform RH1/RH2 lego plot with BufferGeometry

      let vertices = jsrp.Box3D.Vertices,
          indicies = jsrp.Box3D.Indexes,
          vnormals = jsrp.Box3D.Normals,
          segments = jsrp.Box3D.Segments,
          // reduced line segments
          rsegments = [0, 1, 1, 2, 2, 3, 3, 0],
          // reduced vertices
          rvertices = [ new THREE.Vector3(0, 0, 0), new THREE.Vector3(0, 1, 0), new THREE.Vector3(1, 1, 0), new THREE.Vector3(1, 0, 0) ],
          main = this.getFramePainter(),
          axis_zmin = main.z_handle.gr.domain()[0],
          axis_zmax = main.z_handle.gr.domain()[1],
          zmin, zmax,
          handle = this.prepareDraw({ rounding: false, use3d: true, extra: 1 }),
          i1 = handle.i1, i2 = handle.i2, j1 = handle.j1, j2 = handle.j2, di = handle.stepi, dj = handle.stepj,
          i, j, k, vert, x1, x2, y1, y2, binz1, binz2, reduced, nobottom, notop,
          histo = this.getHisto(),
          basehisto = histo ? histo.$baseh : null,
          split_faces = (this.options.Lego === 11) || (this.options.Lego === 13); // split each layer on two parts

      if ((i1 >= i2) || (j1 >= j2)) return;

      let GetBinContent = (ii,jj, level) => {
         // return bin content in binz1, binz2, reduced flags
         // return true if bin should be displayed

         binz2 = histo.getBinContent(ii+1, jj+1);
         if (basehisto)
            binz1 = basehisto.getBinContent(ii+1, jj+1);
         else if (this.options.BaseLine !== false)
            binz1 = this.options.BaseLine;
         else
            binz1 = this.options.Zero ? axis_zmin : 0;
         if (binz2 < binz1) { let d = binz1; binz1 = binz2; binz2 = d; }

         if ((binz1 >= zmax) || (binz2 < zmin)) return false;

         reduced = (binz2 === zmin) || (binz1 >= binz2);

         if (!reduced || (level>0)) return true;

         if (histo['$baseh']) return false; // do not draw empty bins on top of other bins

         if (this.options.Zero || (axis_zmin > 0)) return true;

         return this._show_empty_bins;
      }

      // if bin ID fit into 16 bit, use smaller arrays for intersect indexes
      let use16indx = (histo.getBin(i2, j2) < 0xFFFF),
          levels = [ axis_zmin, axis_zmax ], palette = null;

      // DRAW ALL CUBES

      if ((this.options.Lego === 12) || (this.options.Lego === 14)) {
         // drawing colors levels, axis can not exceed palette

         palette = main.getHistPalette();
         this.createContour(main, palette, { full_z_range: true });
         levels = palette.getContour();
         axis_zmin = levels[0];
         axis_zmax = levels[levels.length-1];
      }

      for (let nlevel=0; nlevel<levels.length-1;++nlevel) {

         zmin = levels[nlevel];
         zmax = levels[nlevel+1];

         // artificially extend last level of color palette to maximal visible value
         if (palette && (nlevel==levels.length-2) && (zmax < axis_zmax)) zmax = axis_zmax;

         let z1 = 0, z2 = 0, numvertices = 0, num2vertices = 0,
             grzmin = main.grz(zmin), grzmax = main.grz(zmax);

         // now calculate size of buffer geometry for boxes

         for (i = i1; i < i2; i += di)
            for (j = j1; j < j2; j += dj) {

               if (!GetBinContent(i,j,nlevel)) continue;

               nobottom = !reduced && (nlevel>0);
               notop = !reduced && (binz2 > zmax) && (nlevel < levels.length-2);

               numvertices += (reduced ? 12 : indicies.length);
               if (nobottom) numvertices -= 6;
               if (notop) numvertices -= 6;

               if (split_faces && !reduced) {
                  numvertices -= 12;
                  num2vertices += 12;
               }
            }

         let positions = new Float32Array(numvertices*3),
             normals = new Float32Array(numvertices*3),
             face_to_bins_index = use16indx ? new Uint16Array(numvertices/3) : new Uint32Array(numvertices/3),
             pos2 = null, norm2 = null, face_to_bins_indx2 = null,
             v = 0, v2 = 0, vert, k, nn;

         if (num2vertices > 0) {
            pos2 = new Float32Array(num2vertices*3);
            norm2 = new Float32Array(num2vertices*3);
            face_to_bins_indx2 = use16indx ? new Uint16Array(num2vertices/3) : new Uint32Array(num2vertices/3);
         }

         for (i = i1; i < i2; i += di) {
            x1 = handle.grx[i] + handle.xbar1*(handle.grx[i+di]-handle.grx[i]);
            x2 = handle.grx[i] + handle.xbar2*(handle.grx[i+di]-handle.grx[i]);
            for (j = j1; j < j2; j += dj) {

               if (!GetBinContent(i,j,nlevel)) continue;

               nobottom = !reduced && (nlevel>0);
               notop = !reduced && (binz2 > zmax) && (nlevel < levels.length-2);

               y1 = handle.gry[j] + handle.ybar1*(handle.gry[j+dj] - handle.gry[j]);
               y2 = handle.gry[j] + handle.ybar2*(handle.gry[j+dj] - handle.gry[j]);

               z1 = (binz1 <= zmin) ? grzmin : main.grz(binz1);
               z2 = (binz2 > zmax) ? grzmax : main.grz(binz2);

               nn = 0; // counter over the normals, each normals correspond to 6 vertices
               k = 0; // counter over vertices

               if (reduced) {
                  // we skip all side faces, keep only top and bottom
                  nn += 12;
                  k += 24;
               }

               let size = indicies.length, bin_index = histo.getBin(i+1, j+1);
               if (nobottom) size -= 6;

               // array over all vertices of the single bin
               while(k < size) {

                  vert = vertices[indicies[k]];

                  if (split_faces && (k < 12)) {
                     pos2[v2]   = x1 + vert.x * (x2 - x1);
                     pos2[v2+1] = y1 + vert.y * (y2 - y1);
                     pos2[v2+2] = z1 + vert.z * (z2 - z1);

                     norm2[v2] = vnormals[nn];
                     norm2[v2+1] = vnormals[nn+1];
                     norm2[v2+2] = vnormals[nn+2];
                     if (v2%9===0) face_to_bins_indx2[v2/9] = bin_index; // remember which bin corresponds to the face
                     v2+=3;
                  } else {
                     positions[v]   = x1 + vert.x * (x2 - x1);
                     positions[v+1] = y1 + vert.y * (y2 - y1);
                     positions[v+2] = z1 + vert.z * (z2 - z1);

                     normals[v] = vnormals[nn];
                     normals[v+1] = vnormals[nn+1];
                     normals[v+2] = vnormals[nn+2];
                     if (v%9===0) face_to_bins_index[v/9] = bin_index; // remember which bin corresponds to the face
                     v+=3;
                  }

                  ++k;

                  if (k%6 === 0) {
                     nn+=3;
                     if (notop && (k === indicies.length - 12)) {
                        k+=6; nn+=3; // jump over notop indexes
                     }
                  }
               }
            }
         }

         let geometry = new THREE.BufferGeometry();
         geometry.setAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
         geometry.setAttribute( 'normal', new THREE.BufferAttribute( normals, 3 ) );
         // geometry.computeVertexNormals();

         let rootcolor = 3, fcolor = this.getColor(rootcolor);

         if (palette) {
            fcolor = palette.getColor(nlevel); // calcColor in v6
         } else if ((this.options.Lego === 1) || (rootcolor < 2)) {
            rootcolor = 1;
            fcolor = 'white';
         }

         //var material = new THREE.MeshLambertMaterial( { color: fcolor } );
         let material = new THREE.MeshBasicMaterial( { color: fcolor, flatShading: true } );

         let mesh = new THREE.Mesh(geometry, material);

         mesh.face_to_bins_index = face_to_bins_index;
         mesh.painter = this;
         mesh.zmin = axis_zmin;
         mesh.zmax = axis_zmax;
         mesh.baseline = (this.options.BaseLine!==false) ? this.options.BaseLine : (this.options.Zero ? axis_zmin : 0);
         mesh.tip_color = (rootcolor===3) ? 0xFF0000 : 0x00FF00;
         mesh.handle = handle;

         mesh.tooltip = function(intersect) {
            if (isNaN(intersect.faceIndex)) {
               console.error('faceIndex not provided, check three.js version', THREE.REVISION, 'expected r102');
               return null;
            }

            if ((intersect.faceIndex < 0) || (intersect.faceIndex >= this.face_to_bins_index.length)) return null;

            let p = this.painter,
                handle = this.handle,
                main = p.getFramePainter(),
                histo = p.getHisto(),
                tip = p.get3DToolTip( this.face_to_bins_index[intersect.faceIndex] ),
                i1 = tip.ix - 1, i2 = i1 + handle.stepi,
                j1 = tip.iy - 1, j2 = j1 + handle.stepj;

            tip.x1 = Math.max(-main.size_xy3d, handle.grx[i1] + handle.xbar1*(handle.grx[i2]-handle.grx[i1]));
            tip.x2 = Math.min(main.size_xy3d, handle.grx[i1] + handle.xbar2*(handle.grx[i2]-handle.grx[i1]));

            tip.y1 = Math.max(-main.size_xy3d, handle.gry[j1] + handle.ybar1*(handle.gry[j2] - handle.gry[j1]));
            tip.y2 = Math.min(main.size_xy3d, handle.gry[j1] + handle.ybar2*(handle.gry[j2] - handle.gry[j1]));

            let binz1 = this.baseline, binz2 = tip.value;
            if (histo['$baseh']) binz1 = histo['$baseh'].getBinContent(i1+1, j1+1);
            if (binz2<binz1) { let v = binz1; binz1 = binz2; binz2 = v; }

            tip.z1 = main.grz(Math.max(this.zmin,binz1));
            tip.z2 = main.grz(Math.min(this.zmax,binz2));

            tip.color = this.tip_color;

            if (p.is_projection && (p.getDimension()==2)) tip.$painter = p; // used only for projections

            return tip;
         }

         main.toplevel.add(mesh);

         if (num2vertices > 0) {
            let geom2 = new THREE.BufferGeometry();
            geom2.setAttribute( 'position', new THREE.BufferAttribute( pos2, 3 ) );
            geom2.setAttribute( 'normal', new THREE.BufferAttribute( norm2, 3 ) );
            //geom2.computeVertexNormals();

            //var material2 = new THREE.MeshLambertMaterial( { color: 0xFF0000 } );

            let color2 = (rootcolor<2) ? new THREE.Color(0xFF0000) :
                            new THREE.Color(d3.rgb(fcolor).darker(0.5).toString());

            let material2 = new THREE.MeshBasicMaterial( { color: color2, flatShading: true } );

            let mesh2 = new THREE.Mesh(geom2, material2);
            mesh2.face_to_bins_index = face_to_bins_indx2;
            mesh2.painter = this;
            mesh2.handle = mesh.handle;
            mesh2.tooltip = mesh.tooltip;
            mesh2.zmin = mesh.zmin;
            mesh2.zmax = mesh.zmax;
            mesh2.baseline = mesh.baseline;
            mesh2.tip_color = mesh.tip_color;

            main.toplevel.add(mesh2);
         }
      }

      // lego3 or lego4 do not draw border lines
      if (this.options.Lego > 12) return;

      // DRAW LINE BOXES

      let numlinevertices = 0, numsegments = 0, uselineindx = true;

      zmax = axis_zmax; zmin = axis_zmin;

      for (i = i1; i < i2; i += di)
         for (j = j1; j < j2; j += dj) {
            if (!GetBinContent(i,j,0)) continue;

            // calculate required buffer size for line segments
            numlinevertices += (reduced ? rvertices.length : vertices.length);
            numsegments += (reduced ? rsegments.length : segments.length);
         }

      // On some platforms vertex index required to be Uint16 array
      // While we cannot use index for large vertex list
      // skip index usage at all. It happens for relatively large histograms (100x100 bins)
      if (numlinevertices > 0xFFF0) uselineindx = false;

      if (!uselineindx) numlinevertices = numsegments*3;

      let lpositions = new Float32Array( numlinevertices * 3 ),
          lindicies = uselineindx ? new Uint16Array( numsegments ) : null,
          grzmin = main.grz(axis_zmin),
          grzmax = main.grz(axis_zmax),
          z1 = 0, z2 = 0, ll = 0, ii = 0;

      for (i = i1; i < i2; i += di) {
         x1 = handle.grx[i] + handle.xbar1*(handle.grx[i+di]-handle.grx[i]);
         x2 = handle.grx[i] + handle.xbar2*(handle.grx[i+di]-handle.grx[i]);
         for (j = j1; j < j2; j += dj) {

            if (!GetBinContent(i,j,0)) continue;

            y1 = handle.gry[j] + handle.ybar1*(handle.gry[j+dj] - handle.gry[j]);
            y2 = handle.gry[j] + handle.ybar2*(handle.gry[j+dj] - handle.gry[j]);

            z1 = (binz1 <= axis_zmin) ? grzmin : main.grz(binz1);
            z2 = (binz2 > axis_zmax) ? grzmax : main.grz(binz2);

            let seg = reduced ? rsegments : segments,
                vvv = reduced ? rvertices : vertices;

            if (uselineindx) {
               // array of indicies for the lines, to avoid duplication of points
               for (k = 0; k < seg.length; ++k)
                  lindicies[ii++] = ll/3 + seg[k];

               for (k = 0; k < vvv.length; ++k) {
                  vert = vvv[k];
                  lpositions[ll]   = x1 + vert.x * (x2 - x1);
                  lpositions[ll+1] = y1 + vert.y * (y2 - y1);
                  lpositions[ll+2] = z1 + vert.z * (z2 - z1);
                  ll+=3;
               }
            } else {
               // copy only vertex positions
               for (k = 0; k < seg.length; ++k) {
                  vert = vvv[seg[k]];
                  lpositions[ll]   = x1 + vert.x * (x2 - x1);
                  lpositions[ll+1] = y1 + vert.y * (y2 - y1);
                  lpositions[ll+2] = z1 + vert.z * (z2 - z1);
                  ll+=3;
               }
            }
         }
      }

      // create boxes
      let lcolor = this.v7EvalColor("line_color", "lightblue");
      let material = new THREE.LineBasicMaterial({ color: new THREE.Color(lcolor), linewidth: this.v7EvalAttr("line_width", 1) });

      let line = jsrp.createLineSegments(lpositions, material, uselineindx ? lindicies : null );

      /*
      line.painter = this;
      line.intersect_index = intersect_index;
      line.tooltip = function(intersect) {
         if ((intersect.index<0) || (intersect.index >= this.intersect_index.length)) return null;
         return this.painter.get3DToolTip(this.intersect_index[intersect.index]);
      }
      */

      main.toplevel.add(line);
   }

   // ==========================================================================================

   /** @summary Draw 1-D histogram in 3D mode
     * @private */
   JSROOT.v7.RH1Painter.prototype.draw3D = function(reason) {

      this.mode3d = true;

      let main = this.getFramePainter(), // who makes axis drawing
          is_main = this.isMainPainter(); // is main histogram

      if (reason == "resize")  {
         if (is_main && main.resize3D()) main.render3D();
         return Promise.resolve(this);
      }

      this.deleteAttr();

      this.scanContent(true); // may be required for axis drawings

      if (is_main) {
         main.create3DScene(this.options.Render3D);
         main.setAxesRanges(this.getAxis("x"), this.xmin, this.xmax, null, this.ymin, this.ymax, null, 0, 0);
         main.set3DOptions(this.options);
         main.drawXYZ(main.toplevel, { use_y_for_z: true, zmult: 1.1, zoom: JSROOT.settings.Zooming, ndim: 1 });
      }

      if (!main.mode3d)
         return Promise.resolve(this);

      return this.drawingBins(reason).then(() => {
         // called when bins received from server, must be reentrant
         let main = this.getFramePainter();

         this.drawLego();
         this.updatePaletteDraw();
         main.render3D();
         main.addKeysHandler();
         return this;
      });

   }

   // ==========================================================================================

   JSROOT.v7.RH2Painter.prototype.draw3D = function(reason) {

      this.mode3d = true;

      let main = this.getFramePainter(), // who makes axis drawing
          is_main = this.isMainPainter(); // is main histogram

      if (reason == "resize") {
         if (is_main && main.resize3D()) main.render3D();

         return Promise.resolve(this);
      }

      let zmult = 1.1;

      this.zmin = main.logz ? this.gminposbin * 0.3 : this.gminbin;
      this.zmax = this.gmaxbin;
      if (this.options.minimum !== -1111) this.zmin = this.options.minimum;
      if (this.options.maximum !== -1111) { this.zmax = this.options.maximum; zmult = 1; }
      if (main.logz && (this.zmin<=0)) this.zmin = this.zmax * 1e-5;

      this.deleteAttr();

      if (is_main) {
         main.create3DScene(this.options.Render3D);
         main.setAxesRanges(this.getAxis("x"), this.xmin, this.xmax, this.getAxis("y"), this.ymin, this.ymax, null, this.zmin, this.zmax);
         main.set3DOptions(this.options);
         main.drawXYZ(main.toplevel, { zmult: zmult, zoom: JSROOT.settings.Zooming, ndim: 2 });
      }

      if (!main.mode3d)
         return Promise.resolve(this);

      return this.drawingBins(reason).then(() => {
         // called when bins received from server, must be reentrant
         let main = this.getFramePainter();

         this.draw3DBins();
         main.render3D();
         main.addKeysHandler();
         return this;
      });
   }

   /** Draw histogram bins in 3D, using provided draw options
     * @private */
   JSROOT.v7.RH2Painter.prototype.draw3DBins = function() {

      if (!this.draw_content) return;

      if (this.isRH2Poly())
         return this.drawPolyLego();

      if (this.options.Surf)
         return this.drawSurf();

      if (this.options.Error)
         return this.drawError();

      if (this.options.Contour)
         return this.drawContour3D(true);

      this.drawLego();
      this.updatePaletteDraw();
   }

   // ==============================================================================

   function create3DLineMaterialv7(painter, prefix) {
      if (!painter) return null;
      if (!prefix) prefix = "line_"

      let lcolor = painter.v7EvalColor(prefix+"color", "black"),
          lstyle = painter.v7EvalAttr(prefix+"style", 0),
          lwidth = painter.v7EvalAttr(prefix+"width", 1),
          material = null,
          style = lstyle ? jsrp.root_line_styles[parseInt(lstyle)] : "",
          dash = style ? style.split(",") : [];

      if (dash && dash.length>=2)
         material = new THREE.LineDashedMaterial({ color: lcolor, dashSize: parseInt(dash[0]), gapSize: parseInt(dash[1]) });
      else
         material = new THREE.LineBasicMaterial({ color: lcolor });

      if ((lwidth !== undefined) && (lwidth > 1)) material.linewidth = parseInt(lwidth);

      return material;
   }

   /** @summary Draw RH2 as 3D contour plot
     * @private */
   JSROOT.v7.RH2Painter.prototype.drawContour3D = function(realz) {
      // for contour plots one requires handle with full range
      let main = this.getFramePainter(),
          handle = this.prepareDraw({rounding: false, use3d: true, extra: 100, middle: 0.0 }),
          palette = main.getHistPalette(),
          layerz = 2*main.size_z3d, pnts = [];

      this.createContour(main, palette, { full_z_range: true });

      let levels = palette.getContour();

      this.buildContour(handle, levels, palette,
         function(colindx,xp,yp,iminus,iplus,ilevel) {
             // ignore less than three points
             if (iplus - iminus < 3) return;

             if (realz) {
                layerz = main.grz(levels[ilevel]);
                if ((layerz < 0) || (layerz > 2*main.size_z3d)) return;
             }

             for (let i=iminus;i<iplus;++i) {
                pnts.push(xp[i], yp[i], layerz);
                pnts.push(xp[i+1], yp[i+1], layerz);
             }
         }
      )

      let lines = jsrp.createLineSegments(pnts, create3DLineMaterialv7(this));
      main.toplevel.add(lines);
   }

   /** @summary Draw RH2 histograms in surf mode
     * @private */
   JSROOT.v7.RH2Painter.prototype.drawSurf = function() {
      let histo = this.getHisto(),
          main = this.getFramePainter(),
          handle = this.prepareDraw({rounding: false, use3d: true, extra: 1, middle: 0.5 }),
          i, j, x1, y1, x2, y2, z11, z12, z21, z22,
          di = handle.stepi, dj = handle.stepj,
          numstepi = handle.i2 - handle.i1, numstepj = handle.j2 - handle.j1,
          axis_zmin = main.z_handle.gr.domain()[0];

      if (di > 1) {
         numstepi = Math.round(numstepi/di);
         if (numstepi * di < (handle.i2 - handle.i1)) numstepi++;
      }

      if (dj > 1) {
         numstepj = Math.round(numstepj/dj);
         if (numstepj * dj < (handle.j2 - handle.j1)) numstepj++;
      }

      // first adjust ranges

      let main_grz = !main.logz ? main.grz : function(value) { return value < axis_zmin ? -0.1 : main.grz(value); }

      if ((handle.i2 - handle.i1 < 2) || (handle.j2 - handle.j1 < 2)) return;

      let ilevels = null, levels = null, dolines = true, dogrid = false,
          donormals = false, palette = null, need_palette = 0;

      switch(this.options.Surf) {
         case 11: need_palette = 2; break;
         case 12:
         case 15: // make surf5 same as surf2
         case 17: need_palette = 2; dolines = false; break;
         case 14: dolines = false; donormals = true; break;
         case 16: need_palette = 1; dogrid = true; dolines = false; break;
         default: ilevels = main.z_handle.createTicks(true); dogrid = true; break;
      }

      if (need_palette > 0) {
         palette = main.getHistPalette();
         if (need_palette == 2)
            this.createContour(main, palette, { full_z_range: true });
         ilevels = palette.getContour();
      }

      if (ilevels) {
         // recalculate levels into graphical coordinates
         levels = new Float32Array(ilevels.length);
         for (let ll=0;ll<ilevels.length;++ll)
            levels[ll] = main_grz(ilevels[ll]);
      } else {
         levels = [0, 2*main.size_z3d]; // just cut top/bottom parts
      }

      let loop, nfaces = [], pos = [], indx = [],    // buffers for faces
          nsegments = 0, lpos = null, lindx = 0,     // buffer for lines
          ngridsegments = 0, grid = null, gindx = 0, // buffer for grid lines segments
          normindx = null;                           // buffer to remember place of vertex for each bin

      function CheckSide(z,level1, level2) {
         if (z<level1) return -1;
         if (z>level2) return 1;
         return 0;
      }

      function AddLineSegment(x1,y1,z1, x2,y2,z2) {
         if (!dolines) return;
         let side1 = CheckSide(z1,0,2*main.size_z3d),
             side2 = CheckSide(z2,0,2*main.size_z3d);
         if ((side1===side2) && (side1!==0)) return;
         if (!loop) return ++nsegments;

         if (side1!==0) {
            let diff = z2-z1;
            z1 = (side1<0) ? 0 : 2*main.size_z3d;
            x1 = x2 - (x2-x1)/diff*(z2-z1);
            y1 = y2 - (y2-y1)/diff*(z2-z1);
         }
         if (side2!==0) {
            let diff = z1-z2;
            z2 = (side2<0) ? 0 : 2*main.size_z3d;
            x2 = x1 - (x1-x2)/diff*(z1-z2);
            y2 = y1 - (y1-y2)/diff*(z1-z2);
         }

         lpos[lindx] = x1; lpos[lindx+1] = y1; lpos[lindx+2] = z1; lindx+=3;
         lpos[lindx] = x2; lpos[lindx+1] = y2; lpos[lindx+2] = z2; lindx+=3;
      }

      let pntbuf = new Float32Array(6*3), k = 0, lastpart = 0; // maximal 6 points
      let gridpnts = new Float32Array(2*3), gridcnt = 0;

      function AddCrossingPoint(xx1,yy1,zz1, xx2,yy2,zz2, crossz, with_grid) {
         if (k>=pntbuf.length) console.log('more than 6 points???');

         let part = (crossz - zz1) / (zz2 - zz1), shift = 3;
         if ((lastpart!==0) && (Math.abs(part) < Math.abs(lastpart))) {
            // while second crossing point closer than first to original, move it in memory
            pntbuf[k] = pntbuf[k-3];
            pntbuf[k+1] = pntbuf[k-2];
            pntbuf[k+2] = pntbuf[k-1];
            k-=3; shift = 6;
         }

         pntbuf[k] = xx1 + part*(xx2-xx1);
         pntbuf[k+1] = yy1 + part*(yy2-yy1);
         pntbuf[k+2] = crossz;

         if (with_grid && grid) {
            gridpnts[gridcnt] = pntbuf[k];
            gridpnts[gridcnt+1] = pntbuf[k+1];
            gridpnts[gridcnt+2] = pntbuf[k+2];
            gridcnt+=3;
         }

         k += shift;
         lastpart = part;
      }

      function RememberVertex(indx, ii,jj) {
         let bin = ((ii-handle.i1)/di * numstepj + (jj-handle.j1)/dj)*8;

         if (normindx[bin]>=0)
            return console.error('More than 8 vertexes for the bin');

         let pos = bin+8+normindx[bin]; // position where write index
         normindx[bin]--;
         normindx[pos] = indx; // at this moment index can be overwritten, means all 8 position are there
      }

      function RecalculateNormals(arr) {
         for (let ii=handle.i1; ii<handle.i2; ii += di) {
            for (let jj=handle.j1; jj<handle.j2; jj += dj) {
               let bin = ((ii-handle.i1)/di * numstepj + (jj-handle.j1)/dj) * 8;

               if (normindx[bin] === -1) continue; // nothing there

               let beg = (normindx[bin] >=0) ? bin : bin+9+normindx[bin],
                   end = bin+8, sumx=0, sumy = 0, sumz = 0;

               for (let kk=beg;kk<end;++kk) {
                  let indx = normindx[kk];
                  if (indx<0) return console.error('FAILURE in NORMALS RECALCULATIONS');
                  sumx+=arr[indx];
                  sumy+=arr[indx+1];
                  sumz+=arr[indx+2];
               }

               sumx = sumx/(end-beg); sumy = sumy/(end-beg); sumz = sumz/(end-beg);

               for (let kk=beg;kk<end;++kk) {
                  let indx = normindx[kk];
                  arr[indx] = sumx;
                  arr[indx+1] = sumy;
                  arr[indx+2] = sumz;
               }
            }
         }
      }

      function AddMainTriangle(x1,y1,z1, x2,y2,z2, x3,y3,z3, is_first) {

         for (let lvl=1;lvl<levels.length;++lvl) {

            let side1 = CheckSide(z1, levels[lvl-1], levels[lvl]),
                side2 = CheckSide(z2, levels[lvl-1], levels[lvl]),
                side3 = CheckSide(z3, levels[lvl-1], levels[lvl]),
                side_sum = side1 + side2 + side3;

            if (side_sum === 3) continue;
            if (side_sum === -3) return;

            if (!loop) {
               let npnts = Math.abs(side2-side1) + Math.abs(side3-side2) + Math.abs(side1-side3);
               if (side1===0) ++npnts;
               if (side2===0) ++npnts;
               if (side3===0) ++npnts;

               if ((npnts===1) || (npnts===2)) console.error('FOND npnts', npnts);

               if (npnts>2) {
                  if (nfaces[lvl]===undefined) nfaces[lvl] = 0;
                  nfaces[lvl] += npnts-2;
               }

               // check if any(contours for given level exists
               if (((side1>0) || (side2>0) || (side3>0)) &&
                   ((side1!==side2) || (side2!==side3) || (side3!==side1))) ++ngridsegments;

               continue;
            }

            gridcnt = 0;

            k = 0;
            if (side1 === 0) { pntbuf[k] = x1; pntbuf[k+1] = y1; pntbuf[k+2] = z1; k+=3; }

            if (side1!==side2) {
               // order is important, should move from 1->2 point, checked via lastpart
               lastpart = 0;
               if ((side1<0) || (side2<0)) AddCrossingPoint(x1,y1,z1, x2,y2,z2, levels[lvl-1]);
               if ((side1>0) || (side2>0)) AddCrossingPoint(x1,y1,z1, x2,y2,z2, levels[lvl], true);
            }

            if (side2 === 0) { pntbuf[k] = x2; pntbuf[k+1] = y2; pntbuf[k+2] = z2; k+=3; }

            if (side2!==side3) {
               // order is important, should move from 2->3 point, checked via lastpart
               lastpart = 0;
               if ((side2<0) || (side3<0)) AddCrossingPoint(x2,y2,z2, x3,y3,z3, levels[lvl-1]);
               if ((side2>0) || (side3>0)) AddCrossingPoint(x2,y2,z2, x3,y3,z3, levels[lvl], true);
            }

            if (side3 === 0) { pntbuf[k] = x3; pntbuf[k+1] = y3; pntbuf[k+2] = z3; k+=3; }

            if (side3!==side1) {
               // order is important, should move from 3->1 point, checked via lastpart
               lastpart = 0;
               if ((side3<0) || (side1<0)) AddCrossingPoint(x3,y3,z3, x1,y1,z1, levels[lvl-1]);
               if ((side3>0) || (side1>0)) AddCrossingPoint(x3,y3,z3, x1,y1,z1, levels[lvl], true);
            }

            if (k===0) continue;
            if (k<9) { console.log('found less than 3 points', k/3); continue; }

            if (grid && (gridcnt === 6)) {
               for (let jj=0;jj < 6; ++jj)
                  grid[gindx+jj] = gridpnts[jj];
               gindx+=6;
            }


            // if three points and surf==14, remember vertex for each point

            let buf = pos[lvl], s = indx[lvl];
            if (donormals && (k===9)) {
               RememberVertex(s, i, j);
               RememberVertex(s+3, i+1, is_first ? j+dj : j);
               RememberVertex(s+6, is_first ? i : i+di, j+dj);
            }

            for (let k1=3;k1<k-3;k1+=3) {
               buf[s] = pntbuf[0]; buf[s+1] = pntbuf[1]; buf[s+2] = pntbuf[2]; s+=3;
               buf[s] = pntbuf[k1]; buf[s+1] = pntbuf[k1+1]; buf[s+2] = pntbuf[k1+2]; s+=3;
               buf[s] = pntbuf[k1+3]; buf[s+1] = pntbuf[k1+4]; buf[s+2] = pntbuf[k1+5]; s+=3;
            }
            indx[lvl] = s;

         }
      }

      if (donormals) {
         // for each bin maximal 8 points reserved
         normindx = new Int32Array(numstepi * numstepj * 8);
         for (let n=0;n<normindx.length;++n) normindx[n] = -1;
      }

      for (loop=0;loop<2;++loop) {
         if (loop) {
            for (let lvl=1;lvl<levels.length;++lvl)
               if (nfaces[lvl]) {
                  pos[lvl] = new Float32Array(nfaces[lvl] * 9);
                  indx[lvl] = 0;
               }
            if (dolines && (nsegments > 0))
               lpos = new Float32Array(nsegments * 6);
            if (dogrid && (ngridsegments>0))
               grid = new Float32Array(ngridsegments * 6);
         }
         for (i = handle.i1;i < handle.i2-di; i += di) {
            x1 = handle.grx[i];
            x2 = handle.grx[i+di];
            for (j = handle.j1; j < handle.j2-dj; j += dj) {
               y1 = handle.gry[j];
               y2 = handle.gry[j+dj];
               z11 = main_grz(histo.getBinContent(i+1, j+1));
               z12 = main_grz(histo.getBinContent(i+1, j+1+dj));
               z21 = main_grz(histo.getBinContent(i+1+dj, j+1));
               z22 = main_grz(histo.getBinContent(i+1+dj, j+1+dj));

               AddMainTriangle(x1,y1,z11, x2,y2,z22, x1,y2,z12, true);

               AddMainTriangle(x1,y1,z11, x2,y1,z21, x2,y2,z22, false);

               AddLineSegment(x1,y2,z12, x1,y1,z11);
               AddLineSegment(x1,y1,z11, x2,y1,z21);

               if (i===handle.i2-2) AddLineSegment(x2,y1,z21, x2,y2,z22);
               if (j===handle.j2-2) AddLineSegment(x1,y2,z12, x2,y2,z22);
            }
         }
      }

      for (let lvl=1;lvl<levels.length;++lvl)
         if (pos[lvl]) {
            if (indx[lvl] !== nfaces[lvl]*9)
                 console.error('SURF faces missmatch lvl', lvl, 'faces', nfaces[lvl], 'index', indx[lvl], 'check', nfaces[lvl]*9 - indx[lvl]);
            let geometry = new THREE.BufferGeometry();
            geometry.setAttribute( 'position', new THREE.BufferAttribute( pos[lvl], 3 ) );
            geometry.computeVertexNormals();
            if (donormals && (lvl===1)) RecalculateNormals(geometry.getAttribute('normal').array);

            let fcolor, material, fFillColor = 5;
            if (palette) {
               fcolor = palette.getColor(lvl-1);
            } else {
               fcolor = fFillColor > 1 ? this.getColor(fFillColor) : 'white';
               if ((this.options.Surf === 14) && (fFillColor<2)) fcolor = this.getColor(48);
            }

            if (this.options.Surf === 14)
               material = new THREE.MeshLambertMaterial( { color: fcolor, side: THREE.DoubleSide  } );
            else
               material = new THREE.MeshBasicMaterial( { color: fcolor, side: THREE.DoubleSide  } );

            let mesh = new THREE.Mesh(geometry, material);

            main.toplevel.add(mesh);

            mesh.painter = this; // to let use it with context menu
         }


      if (lpos) {
         if (nsegments*6 !== lindx)
            console.error('SURF lines mismmatch nsegm', nsegments, ' lindx', lindx, 'difference', nsegments*6 - lindx);

         let lcolor = this.getColor(7),
             material = new THREE.LineBasicMaterial({ color: new THREE.Color(lcolor), linewidth: this.v7EvalAttr("line_width", 1) }),
             line = jsrp.createLineSegments(lpos, material);
         line.painter = this;
         main.toplevel.add(line);
      }

      if (grid) {
         if (ngridsegments*6 !== gindx)
            console.error('SURF grid draw mismatch ngridsegm', ngridsegments, 'gindx', gindx, 'diff', ngridsegments*6 - gindx);

         let material;

         if (this.options.Surf === 1)
            material = new THREE.LineDashedMaterial( { color: 0x0, dashSize: 2, gapSize: 2 } );
         else
            material = new THREE.LineBasicMaterial({ color: new THREE.Color(this.v7EvalColor("line_color", "lightblue")) });

         let line = jsrp.createLineSegments(grid, material);
         line.painter = this;
         main.toplevel.add(line);
      }

      if (this.options.Surf === 17)
         this.drawContour3D();

      if (this.options.Surf === 13) {

         handle = this.prepareDraw({rounding: false, use3d: true, extra: 100, middle: 0.0 });

         // get levels
         let levels = this.getContour(), // init contour
             palette = this.getHistPalette(),
             lastcolindx = -1, layerz = 2*main.size_z3d;

         this.buildContour(handle, levels, palette,
            function(colindx,xp,yp,iminus,iplus) {
                // no need for duplicated point
                if ((xp[iplus] === xp[iminus]) && (yp[iplus] === yp[iminus])) iplus--;

                // ignore less than three points
                if (iplus - iminus < 3) return;

                let pnts = [];

                for (let i = iminus; i <= iplus; ++i)
                   if ((i === iminus) || (xp[i] !== xp[i-1]) || (yp[i] !== yp[i-1]))
                      pnts.push(new THREE.Vector2(xp[i], yp[i]));

                if (pnts.length < 3) return;

                let faces = THREE.ShapeUtils.triangulateShape(pnts , []);

                if (!faces || (faces.length === 0)) return;

                if ((lastcolindx < 0) || (lastcolindx !== colindx)) {
                   lastcolindx = colindx;
                   layerz+=0.0001*main.size_z3d; // change layers Z
                }

                let pos = new Float32Array(faces.length*9),
                    norm = new Float32Array(faces.length*9),
                    indx = 0;

                for (let n=0;n<faces.length;++n) {
                   let face = faces[n];
                   for (let v=0;v<3;++v) {
                      let pnt = pnts[face[v]];
                      pos[indx] = pnt.x;
                      pos[indx+1] = pnt.y;
                      pos[indx+2] = layerz;
                      norm[indx] = 0;
                      norm[indx+1] = 0;
                      norm[indx+2] = 1;

                      indx+=3;
                   }
                }

                let geometry = new THREE.BufferGeometry();
                geometry.setAttribute( 'position', new THREE.BufferAttribute( pos, 3 ) );
                geometry.setAttribute( 'normal', new THREE.BufferAttribute( norm, 3 ) );

                let fcolor = palette.getColor(colindx);
                let material = new THREE.MeshBasicMaterial( { color: fcolor, flatShading: true, side: THREE.DoubleSide, opacity: 0.5  } );
                let mesh = new THREE.Mesh(geometry, material);
                mesh.painter = this;
                main.toplevel.add(mesh);
            }
         );
      }

      this.updatePaletteDraw();
   }

   /** @summary Draw RH2 histogram in error mode
     * @private */
   JSROOT.v7.RH2Painter.prototype.drawError = function() {
      let main = this.getFramePainter(),
          histo = this.getHisto(),
          handle = this.prepareDraw({ rounding: false, use3d: true, extra: 1 }),
          zmin = main.z_handle.gr.domain()[0],
          zmax = main.z_handle.gr.domain()[1],
          i, j, binz, binerr, x1, y1, x2, y2, z1, z2,
          nsegments = 0, lpos = null, binindx = null, lindx = 0,
          di = handle.stepi, dj = handle.stepj;

       let check_skip_min = () => {
          // return true if minimal histogram value should be skipped
          if (this.options.Zero || (zmin>0)) return false;
          return !this._show_empty_bins;
       }

       // loop over the points - first loop counts points, second fill arrays
       for (let loop=0;loop<2;++loop) {

          for (i = handle.i1; i < handle.i2; i += di) {
             x1 = handle.grx[i];
             x2 = handle.grx[i+di];
             for (j = handle.j1; j < handle.j2; j += dj) {
                binz = histo.getBinContent(i+1, j+1);
                if ((binz < zmin) || (binz > zmax)) continue;
                if ((binz===zmin) && check_skip_min()) continue;

                // just count number of segments
                if (loop===0) { nsegments+=3; continue; }

                binerr = histo.getBinError(i+1,j+1);
                binindx[lindx/18] = histo.getBin(i+1,j+1);

                y1 = handle.gry[j];
                y2 = handle.gry[j+dj];

                z1 = main.grz((binz - binerr < zmin) ? zmin : binz-binerr);
                z2 = main.grz((binz + binerr > zmax) ? zmax : binz+binerr);

                lpos[lindx] = x1; lpos[lindx+3] = x2;
                lpos[lindx+1] = lpos[lindx+4] = (y1+y2)/2;
                lpos[lindx+2] = lpos[lindx+5] = (z1+z2)/2;
                lindx+=6;

                lpos[lindx] = lpos[lindx+3] = (x1+x2)/2;
                lpos[lindx+1] = y1; lpos[lindx+4] = y2;
                lpos[lindx+2] = lpos[lindx+5] = (z1+z2)/2;
                lindx+=6;

                lpos[lindx] = lpos[lindx+3] = (x1+x2)/2;
                lpos[lindx+1] = lpos[lindx+4] = (y1+y2)/2;
                lpos[lindx+2] = z1; lpos[lindx+5] = z2;
                lindx+=6;
             }
          }

          if (loop===0) {
             if (nsegments===0) return;
             lpos = new Float32Array(nsegments*6);
             binindx = new Int32Array(nsegments/3);
          }
       }

       // create lines
       let lcolor = new THREE.Color(this.v7EvalColor("line_color", "lightblue")),
           material = new THREE.LineBasicMaterial({ color: lcolor, linewidth: this.v7EvalAttr("line_width", 1) }),
           line = jsrp.createLineSegments(lpos, material);

       line.painter = this;
       line.intersect_index = binindx;
       line.zmin = zmin;
       line.zmax = zmax;
       line.tip_color = (lcolor.g < 0.5) ? 0x00FF00 : 0xFF0000;
       line.handle = handle;

       line.tooltip = function(intersect) {
          if (isNaN(intersect.index)) {
             console.error('segment index not provided, check three.js version', THREE.REVISION, 'expected r102');
             return null;
          }

          let pos = Math.floor(intersect.index / 6);
          if ((pos<0) || (pos >= this.intersect_index.length)) return null;
          let p = this.painter,
              main = p.getFramePainter(),
              tip = p.get3DToolTip(this.intersect_index[pos]),
              handle = this.handle,
              i1 = tip.ix - 1, i2 = i1 + handle.stepi,
              j1 = tip.iy - 1, j2 = j1 + handle.stepj;


          tip.x1 = Math.max(-main.size_xy3d, handle.grx[i1]);
          tip.x2 = Math.min(main.size_xy3d, handle.grx[i2]);
          tip.y1 = Math.max(-main.size_xy3d, handle.gry[j1]);
          tip.y2 = Math.min(main.size_xy3d, handle.gry[j2]);

          tip.z1 = main.grz(tip.value-tip.error < this.zmin ? this.zmin : tip.value-tip.error);
          tip.z2 = main.grz(tip.value+tip.error > this.zmax ? this.zmax : tip.value+tip.error);

          tip.color = this.tip_color;

          return tip;
       }

       main.toplevel.add(line);
   }

   /** @summary Draw RH2Poly histogram as lego
     * @private */
   JSROOT.v7.RH2Painter.prototype.drawPolyLego = function() {
      let histo = this.getHisto(),
          palette = this.getHistPalette(),
          pmain = this.getFramePainter(),
          axis_zmin = pmain.z_handle.gr.domain()[0],
          axis_zmax = pmain.z_handle.gr.domain()[1],
          colindx, bin, i, len = histo.fBins.arr.length,
          z0 = pmain.grz(axis_zmin), z1 = z0;

      // use global coordinates
      this.maxbin = this.gmaxbin;
      this.minbin = this.gminbin;
      this.minposbin = this.gminposbin;

      for (i = 0; i < len; ++ i) {
         bin = histo.fBins.arr[i];
         if (bin.fContent < axis_zmin) continue;

         colindx = palette.getContourIndex(bin.fContent);
         if (colindx === null) continue;

         // check if bin outside visible range
         if ((bin.fXmin > pmain.scale_xmax) || (bin.fXmax < pmain.scale_xmin) ||
             (bin.fYmin > pmain.scale_ymax) || (bin.fYmax < pmain.scale_ymin)) continue;

         z1 = pmain.grz((bin.fContent > axis_zmax) ? axis_zmax : bin.fContent);

         let all_pnts = [], all_faces = [],
             ngraphs = 1, gr = bin.fPoly, nfaces = 0;

         if (gr._typename=='TMultiGraph') {
            ngraphs = bin.fPoly.fGraphs.arr.length;
            gr = null;
         }

         for (let ngr = 0; ngr < ngraphs; ++ngr) {
            if (!gr || (ngr>0)) gr = bin.fPoly.fGraphs.arr[ngr];

            let npnts = gr.fNpoints, x = gr.fX, y = gr.fY;
            while ((npnts>2) && (x[0]===x[npnts-1]) && (y[0]===y[npnts-1])) --npnts;

            let pnts, faces;

            for (let ntry=0;ntry<2;++ntry) {
               // run two loops - on the first try to compress data, on second - run as is (removing duplication)

               let lastx, lasty, currx, curry,
                   dist2 = pmain.size_xy3d*pmain.size_z3d,
                   dist2limit = (ntry>0) ? 0 : dist2/1e6;

               pnts = []; faces = null;

               for (let vert = 0; vert < npnts; ++vert) {
                  currx = pmain.grx(x[vert]);
                  curry = pmain.gry(y[vert]);
                  if (vert>0)
                     dist2 = (currx-lastx)*(currx-lastx) + (curry-lasty)*(curry-lasty);
                  if (dist2 > dist2limit) {
                     pnts.push(new THREE.Vector2(currx, curry));
                     lastx = currx;
                     lasty = curry;
                  }
               }

               try {
                  if (pnts.length > 2)
                     faces = THREE.ShapeUtils.triangulateShape(pnts , []);
               } catch(e) {
                  faces = null;
               }

               if (faces && (faces.length>pnts.length-3)) break;
            }

            if (faces && faces.length && pnts) {
               all_pnts.push(pnts);
               all_faces.push(faces);

               nfaces += faces.length * 2;
               if (z1>z0) nfaces += pnts.length*2;
            }
         }

         let pos = new Float32Array(nfaces*9), indx = 0;

         for (let ngr=0;ngr<all_pnts.length;++ngr) {
            let pnts = all_pnts[ngr], faces = all_faces[ngr];

            for (let layer=0;layer<2;++layer) {
               for (let n=0;n<faces.length;++n) {
                  let face = faces[n],
                      pnt1 = pnts[face[0]],
                      pnt2 = pnts[face[(layer===0) ? 2 : 1]],
                      pnt3 = pnts[face[(layer===0) ? 1 : 2]];

                  pos[indx] = pnt1.x;
                  pos[indx+1] = pnt1.y;
                  pos[indx+2] = layer ? z1 : z0;
                  indx+=3;

                  pos[indx] = pnt2.x;
                  pos[indx+1] = pnt2.y;
                  pos[indx+2] = layer ? z1 : z0;
                  indx+=3;

                  pos[indx] = pnt3.x;
                  pos[indx+1] = pnt3.y;
                  pos[indx+2] = layer ? z1 : z0;
                  indx+=3;
               }
            }

            if (z1>z0) {
               for (let n=0;n<pnts.length;++n) {
                  let pnt1 = pnts[n],
                      pnt2 = pnts[(n>0) ? n-1 : pnts.length-1];

                  pos[indx] = pnt1.x;
                  pos[indx+1] = pnt1.y;
                  pos[indx+2] = z0;
                  indx+=3;

                  pos[indx] = pnt2.x;
                  pos[indx+1] = pnt2.y;
                  pos[indx+2] = z0;
                  indx+=3;

                  pos[indx] = pnt2.x;
                  pos[indx+1] = pnt2.y;
                  pos[indx+2] = z1;
                  indx+=3;

                  pos[indx] = pnt1.x;
                  pos[indx+1] = pnt1.y;
                  pos[indx+2] = z0;
                  indx+=3;

                  pos[indx] = pnt2.x;
                  pos[indx+1] = pnt2.y;
                  pos[indx+2] = z1;
                  indx+=3;

                  pos[indx] = pnt1.x;
                  pos[indx+1] = pnt1.y;
                  pos[indx+2] = z1;
                  indx+=3;
               }
            }
         }

         let geometry = new THREE.BufferGeometry();
         geometry.setAttribute( 'position', new THREE.BufferAttribute( pos, 3 ) );
         geometry.computeVertexNormals();

         let fcolor = this.fPalette.getColor(colindx);
         let material = new THREE.MeshBasicMaterial( { color: fcolor, flatShading: true } );
         let mesh = new THREE.Mesh(geometry, material);

         pmain.toplevel.add(mesh);

         mesh.painter = this;
         mesh.bins_index = i;
         mesh.draw_z0 = z0;
         mesh.draw_z1 = z1;
         mesh.tip_color = 0x00FF00;

         mesh.tooltip = function(/*intersects*/) {

            let p = this.painter, main = p.getFramePainter(),
                bin = p.getObject().fBins.arr[this.bins_index];

            let tip = {
              use_itself: true, // indicate that use mesh itself for highlighting
              x1: main.grx(bin.fXmin),
              x2: main.grx(bin.fXmax),
              y1: main.gry(bin.fYmin),
              y2: main.gry(bin.fYmax),
              z1: this.draw_z0,
              z2: this.draw_z1,
              bin: this.bins_index,
              value: bin.fContent,
              color: this.tip_color,
              lines: p.getPolyBinTooltips(this.bins_index)
            };

            return tip;
         };
      }
   }

   // ==============================================================================

   /**
    * @summary Painter for RH3 classes
    *
    * @class
    * @memberof JSROOT.v7
    * @extends JSROOT.v7.RHistPainter
    * @param {object|string} dom - DOM element or id
    * @param {object} histo - histogram object
    * @private
    */

   function RH3Painter(dom, histo) {
      JSROOT.v7.RHistPainter.call(this, dom, histo);

      this.mode3d = true;
   }

   RH3Painter.prototype = Object.create(JSROOT.v7.RHistPainter.prototype);

   RH3Painter.prototype.getDimension = function() {
      return 3;
   }

   RH3Painter.prototype.scanContent = function(when_axis_changed) {

      // no need to rescan histogram while result does not depend from axis selection
      if (when_axis_changed && this.nbinsx && this.nbinsy && this.nbinsz) return;

      let histo = this.getHisto();
      if (!histo) return;

      this.extractAxesProperties(3);

      // global min/max, used at the moment in 3D drawing

      if (this.isDisplayItem()) {
         // take min/max values from the display item
         this.gminbin = histo.fContMin;
         this.gminposbin = histo.fContMinPos > 0 ? histo.fContMinPos : null;
         this.gmaxbin = histo.fContMax;
      } else {
         this.gminbin = this.gmaxbin = histo.getBinContent(1,1,1);

         for (let i = 0; i < this.nbinsx; ++i)
            for (let j = 0; j < this.nbinsy; ++j)
               for (let k = 0; k < this.nbinsz; ++k) {
                  let bin_content = histo.getBinContent(i+1, j+1, k+1);
                  if (bin_content < this.gminbin) this.gminbin = bin_content; else
                  if (bin_content > this.gmaxbin) this.gmaxbin = bin_content;
               }
      }

      this.draw_content = this.gmaxbin > 0;
   }

   RH3Painter.prototype.countStat = function() {
      let histo = this.getHisto(),
          xaxis = this.getAxis("x"),
          yaxis = this.getAxis("y"),
          zaxis = this.getAxis("z"),
          stat_sum0 = 0, stat_sumx1 = 0, stat_sumy1 = 0,
          stat_sumz1 = 0, stat_sumx2 = 0, stat_sumy2 = 0, stat_sumz2 = 0,
          i1 = this.getSelectIndex("x", "left"),
          i2 = this.getSelectIndex("x", "right"),
          j1 = this.getSelectIndex("y", "left"),
          j2 = this.getSelectIndex("y", "right"),
          k1 = this.getSelectIndex("z", "left"),
          k2 = this.getSelectIndex("z", "right"),
          fp = this.getFramePainter(),
          res = { name: histo.fName, entries: 0, integral: 0, meanx: 0, meany: 0, meanz: 0, rmsx: 0, rmsy: 0, rmsz: 0 },
          xi, yi, zi, xx, xside, yy, yside, zz, zside, cont;

      for (xi = 1; xi <= this.nbinsx; ++xi) {

         xx = xaxis.GetBinCoord(xi - 0.5);
         xside = (xi <= i1+1) ? 0 : (xi > i2+1 ? 2 : 1);

         for (yi = 1; yi <= this.nbinsy; ++yi) {

            yy = yaxis.GetBinCoord(yi - 0.5);
            yside = (yi <= j1+1) ? 0 : (yi > j2+1 ? 2 : 1);

            for (zi = 1; zi <= this.nbinsz; ++zi) {

               zz = zaxis.GetBinCoord(zi - 0.5);
               zside = (zi <= k1+1) ? 0 : (zi > k2+1 ? 2 : 1);

               cont = histo.getBinContent(xi, yi, zi);
               res.entries += cont;

               if ((xside==1) && (yside==1) && (zside==1)) {
                  stat_sum0 += cont;
                  stat_sumx1 += xx * cont;
                  stat_sumy1 += yy * cont;
                  stat_sumz1 += zz * cont;
                  stat_sumx2 += xx * xx * cont;
                  stat_sumy2 += yy * yy * cont;
                  stat_sumz2 += zz * zz * cont;
               }
            }
         }
      }

      if ((histo.fTsumw > 0) && !fp.isAxisZoomed("x") && !fp.isAxisZoomed("y") && !fp.isAxisZoomed("z")) {
         stat_sum0  = histo.fTsumw;
         stat_sumx1 = histo.fTsumwx;
         stat_sumx2 = histo.fTsumwx2;
         stat_sumy1 = histo.fTsumwy;
         stat_sumy2 = histo.fTsumwy2;
         stat_sumz1 = histo.fTsumwz;
         stat_sumz2 = histo.fTsumwz2;
      }

      if (stat_sum0 > 0) {
         res.meanx = stat_sumx1 / stat_sum0;
         res.meany = stat_sumy1 / stat_sum0;
         res.meanz = stat_sumz1 / stat_sum0;
         res.rmsx = Math.sqrt(Math.abs(stat_sumx2 / stat_sum0 - res.meanx * res.meanx));
         res.rmsy = Math.sqrt(Math.abs(stat_sumy2 / stat_sum0 - res.meany * res.meany));
         res.rmsz = Math.sqrt(Math.abs(stat_sumz2 / stat_sum0 - res.meanz * res.meanz));
      }

      res.integral = stat_sum0;

      if (histo.fEntries > 1) res.entries = histo.fEntries;

      return res;
   }

   RH3Painter.prototype.fillStatistic = function(stat, dostat /*, dofit */) {

      let data = this.countStat(),
          print_name = dostat % 10,
          print_entries = Math.floor(dostat / 10) % 10,
          print_mean = Math.floor(dostat / 100) % 10,
          print_rms = Math.floor(dostat / 1000) % 10,
          // print_under = Math.floor(dostat / 10000) % 10,
          // print_over = Math.floor(dostat / 100000) % 10,
          print_integral = Math.floor(dostat / 1000000) % 10;
          // var print_skew = Math.floor(dostat / 10000000) % 10;
          // var print_kurt = Math.floor(dostat / 100000000) % 10;

      stat.clearStat();

      if (print_name > 0)
         stat.addText(data.name);

      if (print_entries > 0)
         stat.addText("Entries = " + stat.format(data.entries,"entries"));

      if (print_mean > 0) {
         stat.addText("Mean x = " + stat.format(data.meanx));
         stat.addText("Mean y = " + stat.format(data.meany));
         stat.addText("Mean z = " + stat.format(data.meanz));
      }

      if (print_rms > 0) {
         stat.addText("Std Dev x = " + stat.format(data.rmsx));
         stat.addText("Std Dev y = " + stat.format(data.rmsy));
         stat.addText("Std Dev z = " + stat.format(data.rmsz));
      }

      if (print_integral > 0) {
         stat.addText("Integral = " + stat.format(data.integral,"entries"));
      }

      return true;
   }

   /** @summary Provide text information (tooltips) for histogram bin
     * @private */
   RH3Painter.prototype.getBinTooltips = function (ix, iy, iz) {
      let lines = [], histo = this.getHisto(),
          dx = 1, dy = 1, dz = 1;

      if (this.isDisplayItem()) {
         dx = histo.stepx || 1;
         dy = histo.stepy || 1;
         dz = histo.stepz || 1;
      }

      lines.push(this.getObjectHint());

      lines.push("x = " + this.getAxisBinTip("x", ix, dx) + "  xbin=" + ix);
      lines.push("y = " + this.getAxisBinTip("y", iy, dy) + "  ybin=" + iy);
      lines.push("z = " + this.getAxisBinTip("z", iz, dz) + "  zbin=" + iz);

      let binz = histo.getBinContent(ix+1, iy+1, iz+1),
          lbl = "entries = "+ ((dx>1) || (dy>1) || (dz>1) ? "~" : "");
      if (binz === Math.round(binz))
         lines.push(lbl + binz);
      else
         lines.push(lbl + jsrp.floatToString(binz, JSROOT.gStyle.fStatFormat));

      return lines;
   }

   /** @summary Try to draw 3D histogram as scatter plot
     * @desc If there are too many points, box will be displayed
     * @private */
   RH3Painter.prototype.draw3DScatter = function(handle) {

      let histo = this.getHisto(),
          main = this.getFramePainter(),
          i1 = handle.i1, i2 = handle.i2, di = handle.stepi,
          j1 = handle.j1, j2 = handle.j2, dj = handle.stepj,
          k1 = handle.k1, k2 = handle.k2, dk = handle.stepk,
          i, j, k, bin_content;

      if ((i2<=i1) || (j2<=j1) || (k2<=k1)) return true;

      // scale down factor if too large values
      let coef = (this.gmaxbin > 1000) ? 1000/this.gmaxbin : 1,
          numpixels = 0, sumz = 0, content_lmt = Math.max(0, this.gminbin);

      for (i = i1; i < i2; i += di) {
         for (j = j1; j < j2; j += dj) {
            for (k = k1; k < k2; k += dk) {
               bin_content = histo.getBinContent(i+1, j+1, k+1);
               sumz += bin_content;
               if (bin_content <= content_lmt) continue;
               numpixels += Math.round(bin_content*coef);
            }
         }
      }

      // too many pixels - use box drawing
      if (numpixels > (main.webgl ? 100000 : 30000)) return false;

      JSROOT.seed(sumz);

      let pnts = new jsrp.PointsCreator(numpixels, main.webgl, main.size_xy3d/200),
          bins = new Int32Array(numpixels), nbin = 0,
          xaxis = this.getAxis("x"), yaxis = this.getAxis("y"), zaxis = this.getAxis("z");

      for (i = i1; i < i2; i += di) {
         for (j = j1; j < j2; j += dj) {
            for (k = k1; k < k2; k += dk) {
               bin_content = histo.getBinContent(i+1, j+1, k+1);
               if (bin_content <= content_lmt) continue;
               let num = Math.round(bin_content*coef);

               for (let n=0;n<num;++n) {
                  let binx = xaxis.GetBinCoord(i+JSROOT.random()),
                      biny = yaxis.GetBinCoord(j+JSROOT.random()),
                      binz = zaxis.GetBinCoord(k+JSROOT.random());

                  // remember bin index for tooltip
                  bins[nbin++] = histo.getBin(i+1, j+1, k+1);

                  pnts.addPoint(main.grx(binx), main.gry(biny), main.grz(binz));
               }
            }
         }
      }

      let color = this.v7EvalColor("fill_color", "red");
      let mesh = pnts.createPoints(color);
      main.toplevel.add(mesh);

      mesh.bins = bins;
      mesh.painter = this;
      mesh.tip_color = 0x00FF00;

      mesh.tooltip = function(intersect) {
         if (isNaN(intersect.index)) {
            console.error('intersect.index not provided, check three.js version', THREE.REVISION, 'expected r102');
            return null;
         }

         let indx = Math.floor(intersect.index / this.nvertex);
         if ((indx<0) || (indx >= this.bins.length)) return null;

         let p = this.painter,
             main = p.getFramePainter(),
             tip = p.get3DToolTip(this.bins[indx]);

         tip.x1 = main.grx(p.getAxis("x").GetBinLowEdge(tip.ix));
         tip.x2 = main.grx(p.getAxis("x").GetBinLowEdge(tip.ix+di));
         tip.y1 = main.gry(p.getAxis("y").GetBinLowEdge(tip.iy));
         tip.y2 = main.gry(p.getAxis("y").GetBinLowEdge(tip.iy+dj));
         tip.z1 = main.grz(p.getAxis("z").GetBinLowEdge(tip.iz));
         tip.z2 = main.grz(p.getAxis("z").GetBinLowEdge(tip.iz+dk));
         tip.color = this.tip_color;
         tip.opacity = 0.3;

         return tip;
      }

      return true;
   }

   /** @summary Drawing of 3D histogram
     * @private */
   RH3Painter.prototype.draw3DBins = function() {

      if (!this.draw_content) return;

      let handle = this.prepareDraw({ only_indexes: true, extra: -0.5, right_extra: -1 });

      if (this.options.Scatter)
         if (this.draw3DScatter(handle)) return;

      let fillcolor = this.v7EvalColor("fill_color", "red"),
          main = this.getFramePainter(),
          buffer_size = 0, use_lambert = false,
          use_helper = false, use_colors = false, use_opacity = 1, use_scale = true,
          single_bin_verts, single_bin_norms,
          tipscale = 0.5;

      if (this.options.Sphere) {

         // drawing spheres
         tipscale = 0.4;
         use_lambert = true;
         if (this.options.Sphere === 11) use_colors = true;

         let geom = main.webgl ? new THREE.SphereGeometry(0.5, 16, 12) : new THREE.SphereGeometry(0.5, 8, 6);
         geom.applyMatrix4( new THREE.Matrix4().makeRotationX( Math.PI / 2 ) );

         buffer_size = geom.faces.length*9;
         single_bin_verts = new Float32Array(buffer_size);
         single_bin_norms = new Float32Array(buffer_size);

         // Fill a typed array with cube geometry that will be shared by all
         // (This technically could be put into an InstancedBufferGeometry but
         // performance gain is likely not huge )
         for (let face = 0; face < geom.faces.length; ++face) {
            single_bin_verts[9*face  ] = geom.vertices[geom.faces[face].a].x;
            single_bin_verts[9*face+1] = geom.vertices[geom.faces[face].a].y;
            single_bin_verts[9*face+2] = geom.vertices[geom.faces[face].a].z;
            single_bin_verts[9*face+3] = geom.vertices[geom.faces[face].b].x;
            single_bin_verts[9*face+4] = geom.vertices[geom.faces[face].b].y;
            single_bin_verts[9*face+5] = geom.vertices[geom.faces[face].b].z;
            single_bin_verts[9*face+6] = geom.vertices[geom.faces[face].c].x;
            single_bin_verts[9*face+7] = geom.vertices[geom.faces[face].c].y;
            single_bin_verts[9*face+8] = geom.vertices[geom.faces[face].c].z;

            single_bin_norms[9*face  ] = geom.faces[face].vertexNormals[0].x;
            single_bin_norms[9*face+1] = geom.faces[face].vertexNormals[0].y;
            single_bin_norms[9*face+2] = geom.faces[face].vertexNormals[0].z;
            single_bin_norms[9*face+3] = geom.faces[face].vertexNormals[1].x;
            single_bin_norms[9*face+4] = geom.faces[face].vertexNormals[1].y;
            single_bin_norms[9*face+5] = geom.faces[face].vertexNormals[1].z;
            single_bin_norms[9*face+6] = geom.faces[face].vertexNormals[2].x;
            single_bin_norms[9*face+7] = geom.faces[face].vertexNormals[2].y;
            single_bin_norms[9*face+8] = geom.faces[face].vertexNormals[2].z;
         }

      } else {

         let indicies = jsrp.Box3D.Indexes,
             normals = jsrp.Box3D.Normals,
             vertices = jsrp.Box3D.Vertices;

         buffer_size = indicies.length*3;
         single_bin_verts = new Float32Array(buffer_size);
         single_bin_norms = new Float32Array(buffer_size);

         for (let k=0,nn=-3;k<indicies.length;++k) {
            let vert = vertices[indicies[k]];
            single_bin_verts[k*3]   = vert.x-0.5;
            single_bin_verts[k*3+1] = vert.y-0.5;
            single_bin_verts[k*3+2] = vert.z-0.5;

            if (k%6===0) nn+=3;
            single_bin_norms[k*3]   = normals[nn];
            single_bin_norms[k*3+1] = normals[nn+1];
            single_bin_norms[k*3+2] = normals[nn+2];
         }
         use_helper = true;

         if (this.options.Box == 11) { use_colors = true; } else
         if (this.options.Box == 12) { use_colors = true; use_helper = false; }  else
         if (this.options.Color) { use_colors = true; use_opacity = 0.5; use_scale = false; use_helper = false; use_lambert = true; }
      }

      if (use_scale)
         use_scale = (this.gminbin || this.gmaxbin) ? 1 / Math.max(Math.abs(this.gminbin), Math.abs(this.gmaxbin)) : 1;

      let histo = this.getHisto(),
          i1 = handle.i1, i2 = handle.i2, di = handle.stepi,
          j1 = handle.j1, j2 = handle.j2, dj = handle.stepj,
          k1 = handle.k1, k2 = handle.k2, dk = handle.stepk,
          palette = null;

      if (use_colors) {
         palette = main.getHistPalette();
         this.createContour(main, palette);
      }

      if ((i2<=i1) || (j2<=j1) || (k2<=k1)) return;

      let xaxis = this.getAxis("x"), yaxis = this.getAxis("y"), zaxis = this.getAxis("z"),
          scalex = (main.grx(xaxis.GetBinCoord(i2)) - main.grx(xaxis.GetBinCoord(i1))) / (i2 - i1) * di,
          scaley = (main.gry(yaxis.GetBinCoord(j2)) - main.gry(yaxis.GetBinCoord(j1))) / (j2 - j1) * dj,
          scalez = (main.grz(zaxis.GetBinCoord(k2)) - main.grz(zaxis.GetBinCoord(k1))) / (k2 - k1) * dk;

      let nbins = 0, i, j, k, wei, bin_content, cols_size = [], num_colors = 0, cols_sequence = [];

      for (i = i1; i < i2; i += di) {
         for (j = j1; j < j2; j += dj) {
            for (k = k1; k < k2; k += dk) {
               bin_content = histo.getBinContent(i+1, j+1, k+1);
               if (!this.options.Color && ((bin_content===0) || (bin_content < this.gminbin))) continue;
               wei = use_scale ? Math.pow(Math.abs(bin_content*use_scale), 0.3333) : 1;
               if (wei < 1e-3) continue; // do not draw empty or very small bins

               nbins++;

               if (!use_colors) continue;

               let colindx = palette.getContourIndex(bin_content);
               if (colindx >= 0) {
                  if (cols_size[colindx] === undefined) {
                     cols_size[colindx] = 0;
                     cols_sequence[colindx] = num_colors++;
                  }
                  cols_size[colindx]+=1;
               } else {
                  console.error('not found color for', bin_content);
               }
            }
         }
      }

      if (!use_colors) {
         cols_size.push(nbins);
         num_colors = 1;
         cols_sequence = [0];
      }

      let cols_nbins = new Array(num_colors),
          bin_verts = new Array(num_colors),
          bin_norms = new Array(num_colors),
          bin_tooltips = new Array(num_colors),
          helper_kind = new Array(num_colors),
          helper_indexes = new Array(num_colors),  // helper_kind == 1, use original vertices
          helper_positions = new Array(num_colors);  // helper_kind == 2, all vertices copied into separate buffer

      for(let ncol = 0; ncol < cols_size.length; ++ncol) {
         if (!cols_size[ncol]) continue; // ignore dummy colors

         nbins = cols_size[ncol]; // how many bins with specified color
         let nseq = cols_sequence[ncol];

         cols_nbins[nseq] = 0; // counter for the filled bins

         helper_kind[nseq] = 0;

         // 1 - use same vertices to create helper, one can use maximal 64K vertices
         // 2 - all vertices copied into separate buffer
         if (use_helper)
            helper_kind[nseq] = (nbins * buffer_size / 3 > 0xFFF0) ? 2 : 1;

         bin_verts[nseq] = new Float32Array(nbins * buffer_size);
         bin_norms[nseq] = new Float32Array(nbins * buffer_size);
         bin_tooltips[nseq] = new Int32Array(nbins);

         if (helper_kind[nseq]===1)
            helper_indexes[nseq] = new Uint16Array(nbins * jsrp.Box3D.MeshSegments.length);

         if (helper_kind[nseq]===2)
            helper_positions[nseq] = new Float32Array(nbins * jsrp.Box3D.Segments.length * 3);
      }

      let binx, grx, biny, gry, binz, grz;
      xaxis = this.getAxis("x"),
      yaxis = this.getAxis("y"),
      zaxis = this.getAxis("z");

      for (i = i1; i < i2; i += di) {
         binx = xaxis.GetBinCenter(i+1); grx = main.grx(binx);
         for (j = j1; j < j2; j += dj) {
            biny = yaxis.GetBinCenter(j+1); gry = main.gry(biny);
            for (k = k1; k < k2; k +=dk) {
               bin_content = histo.getBinContent(i+1, j+1, k+1);
               if (!this.options.Color && ((bin_content===0) || (bin_content < this.gminbin))) continue;

               wei = use_scale ? Math.pow(Math.abs(bin_content*use_scale), 0.3333) : 1;
               if (wei < 1e-3) continue; // do not show very small bins

               let nseq = 0;
               if (use_colors) {
                  let colindx = palette.getContourIndex(bin_content);
                  if (colindx < 0) continue;
                  nseq = cols_sequence[colindx];
               }

               nbins = cols_nbins[nseq];

               binz = zaxis.GetBinCenter(k+1); grz = main.grz(binz);

               // remember bin index for tooltip
               bin_tooltips[nseq][nbins] = histo.getBin(i+1, j+1, k+1);

               let vvv = nbins * buffer_size, bin_v = bin_verts[nseq], bin_n = bin_norms[nseq];

               // Grab the coordinates and scale that are being assigned to each bin
               for (let vi = 0; vi < buffer_size; vi+=3, vvv+=3) {
                  bin_v[vvv]   = grx + single_bin_verts[vi]*scalex*wei;
                  bin_v[vvv+1] = gry + single_bin_verts[vi+1]*scaley*wei;
                  bin_v[vvv+2] = grz + single_bin_verts[vi+2]*scalez*wei;

                  bin_n[vvv]   = single_bin_norms[vi];
                  bin_n[vvv+1] = single_bin_norms[vi+1];
                  bin_n[vvv+2] = single_bin_norms[vi+2];
               }

               if (helper_kind[nseq]===1) {
                  // reuse vertices created for the mesh
                  let helper_segments = jsrp.Box3D.MeshSegments;
                  vvv = nbins * helper_segments.length;
                  let shift = Math.round(nbins * buffer_size/3),
                      helper_i = helper_indexes[nseq];
                  for (let n=0;n<helper_segments.length;++n)
                     helper_i[vvv+n] = shift + helper_segments[n];
               }

               if (helper_kind[nseq]===2) {
                  let helper_segments = jsrp.Box3D.Segments,
                      helper_p = helper_positions[nseq];
                  vvv = nbins * helper_segments.length * 3;
                  for (let n=0;n<helper_segments.length;++n, vvv+=3) {
                     let vert = jsrp.Box3D.Vertices[helper_segments[n]];
                     helper_p[vvv]   = grx + (vert.x-0.5)*scalex*wei;
                     helper_p[vvv+1] = gry + (vert.y-0.5)*scaley*wei;
                     helper_p[vvv+2] = grz + (vert.z-0.5)*scalez*wei;
                  }
               }

               cols_nbins[nseq] = nbins+1;
            }
         }
      }

      for(let ncol=0;ncol<cols_size.length;++ncol) {
         if (!cols_size[ncol]) continue; // ignore dummy colors

         nbins = cols_size[ncol]; // how many bins with specified color
         let nseq = cols_sequence[ncol];

         // BufferGeometries that store geometry of all bins
         let all_bins_buffgeom = new THREE.BufferGeometry();

         // Create mesh from bin buffergeometry
         all_bins_buffgeom.setAttribute('position', new THREE.BufferAttribute( bin_verts[nseq], 3 ) );
         all_bins_buffgeom.setAttribute('normal', new THREE.BufferAttribute( bin_norms[nseq], 3 ) );

         if (use_colors) fillcolor = palette.getColor(ncol);

         let material = use_lambert ? new THREE.MeshLambertMaterial({ color: fillcolor, opacity: use_opacity, transparent: (use_opacity<1) })
                                    : new THREE.MeshBasicMaterial({ color: fillcolor, opacity: use_opacity });

         let combined_bins = new THREE.Mesh(all_bins_buffgeom, material);

         combined_bins.bins = bin_tooltips[nseq];
         combined_bins.bins_faces = buffer_size/9;
         combined_bins.painter = this;

         combined_bins.scalex = tipscale*scalex;
         combined_bins.scaley = tipscale*scaley;
         combined_bins.scalez = tipscale*scalez;
         combined_bins.tip_color = 0x00FF00;
         combined_bins.use_scale = use_scale;

         combined_bins.tooltip = function(intersect) {
            if (isNaN(intersect.faceIndex)) {
               console.error('intersect.faceIndex not provided, check three.js version', THREE.REVISION, 'expected r102');
               return null;
            }
            let indx = Math.floor(intersect.faceIndex / this.bins_faces);
            if ((indx<0) || (indx >= this.bins.length)) return null;

            let p = this.painter,
                main = p.getFramePainter(),
                tip = p.get3DToolTip(this.bins[indx]),
                grx = main.grx(p.getAxis("x").GetBinCoord(tip.ix-0.5)),
                gry = main.gry(p.getAxis("y").GetBinCoord(tip.iy-0.5)),
                grz = main.grz(p.getAxis("z").GetBinCoord(tip.iz-0.5)),
                wei = this.use_scale ? Math.pow(Math.abs(tip.value*this.use_scale), 0.3333) : 1;

            tip.x1 = grx - this.scalex*wei; tip.x2 = grx + this.scalex*wei;
            tip.y1 = gry - this.scaley*wei; tip.y2 = gry + this.scaley*wei;
            tip.z1 = grz - this.scalez*wei; tip.z2 = grz + this.scalez*wei;

            tip.color = this.tip_color;

            return tip;
         }

         main.toplevel.add(combined_bins);

         if (helper_kind[nseq] > 0) {
            let lcolor = this.v7EvalColor("line_color", "lightblue"),
                helper_material = new THREE.LineBasicMaterial( { color: lcolor } ),
                lines = null;

            if (helper_kind[nseq] === 1) {
               // reuse positions from the mesh - only special index was created
               lines = jsrp.createLineSegments( bin_verts[nseq], helper_material, helper_indexes[nseq] );
            } else {
               lines = jsrp.createLineSegments( helper_positions[nseq], helper_material );
            }

            main.toplevel.add(lines);
         }
      }

      if (use_colors)
         this.updatePaletteDraw();
   }

   /** @summary Redraw of 3D histogram
     * @private */
   RH3Painter.prototype.redraw = function(reason) {

      let main = this.getFramePainter(); // who makes axis and 3D drawing

      if (reason == "resize") {
         if (main.resize3D()) main.render3D();
         return;
      }

      main.create3DScene(this.options.Render3D);
      main.setAxesRanges(this.getAxis("x"), this.xmin, this.xmax, this.getAxis("y"), this.ymin, this.ymax, this.getAxis("z"), this.zmin, this.zmax);
      main.set3DOptions(this.options);
      main.drawXYZ(main.toplevel, { zoom: JSROOT.settings.Zooming, ndim: 3 });

      this.drawingBins(reason).then(() => {
         // called when bins received from server, must be reentrant

         let main = this.getFramePainter();

         this.draw3DBins();
         main.render3D();
         main.addKeysHandler();
      });
   }

   /** @summary Fill pad toolbar with RH3-related functions
     * @private */
   RH3Painter.prototype.fillToolbar = function() {
      let pp = this.getPadPainter();
      if (!pp) return;

      pp.addPadButton("auto_zoom", 'Unzoom all axes', 'ToggleZoom', "Ctrl *");
      if (this.draw_content)
         pp.addPadButton("statbox", 'Toggle stat box', "ToggleStatBox");
      pp.showPadButtons();
   }

   /** @summary Checks if it makes sense to zoom inside specified axis range */
   RH3Painter.prototype.canZoomInside = function(axis,min,max) {
      let obj = this.getHisto();
      if (obj) obj = obj["f"+axis.toUpperCase()+"axis"];
      return !obj || (obj.FindBin(max,0.5) - obj.FindBin(min,0) > 1);
   }

   /** @summary Perform automatic zoom inside non-zero region of histogram
     * @private */
   RH3Painter.prototype.autoZoom = function() {
      let i1 = this.getSelectIndex("x", "left"),
          i2 = this.getSelectIndex("x", "right"),
          j1 = this.getSelectIndex("y", "left"),
          j2 = this.getSelectIndex("y", "right"),
          k1 = this.getSelectIndex("z", "left"),
          k2 = this.getSelectIndex("z", "right"),
          i, j, k, histo = this.getHisto();

      if ((i1 === i2) || (j1 === j2) || (k1 === k2)) return;

      // first find minimum
      let min = histo.getBinContent(i1 + 1, j1 + 1, k1+1);
      for (i = i1; i < i2; ++i)
         for (j = j1; j < j2; ++j)
            for (k = k1; k < k2; ++k)
               min = Math.min(min, histo.getBinContent(i+1, j+1, k+1));

      if (min>0) return; // if all points positive, no chance for autoscale

      let ileft = i2, iright = i1, jleft = j2, jright = j1, kleft = k2, kright = k1;

      for (i = i1; i < i2; ++i)
         for (j = j1; j < j2; ++j)
            for (k = k1; k < k2; ++k)
               if (histo.getBinContent(i+1, j+1, k+1) > min) {
                  if (i < ileft) ileft = i;
                  if (i >= iright) iright = i + 1;
                  if (j < jleft) jleft = j;
                  if (j >= jright) jright = j + 1;
                  if (k < kleft) kleft = k;
                  if (k >= kright) kright = k + 1;
               }

      let xmin, xmax, ymin, ymax, zmin, zmax, isany = false;

      if ((ileft === iright-1) && (ileft > i1+1) && (iright < i2-1)) { ileft--; iright++; }
      if ((jleft === jright-1) && (jleft > j1+1) && (jright < j2-1)) { jleft--; jright++; }
      if ((kleft === kright-1) && (kleft > k1+1) && (kright < k2-1)) { kleft--; kright++; }

      if ((ileft > i1 || iright < i2) && (ileft < iright - 1)) {
         xmin = this.getAxis("x").GetBinLowEdge(ileft+1);
         xmax = this.getAxis("x").GetBinLowEdge(iright+1);
         isany = true;
      }

      if ((jleft > j1 || jright < j2) && (jleft < jright - 1)) {
         ymin = this.getAxis("y").GetBinLowEdge(jleft+1);
         ymax = this.getAxis("y").GetBinLowEdge(jright+1);
         isany = true;
      }

      if ((kleft > k1 || kright < k2) && (kleft < kright - 1)) {
         zmin = this.getAxis("z").GetBinLowEdge(kleft+1);
         zmax = this.getAxis("z").GetBinLowEdge(kright+1);
         isany = true;
      }

      if (isany) this.getFramePainter().zoom(xmin, xmax, ymin, ymax, zmin, zmax);
   }

   /** @summary Fill histogram context menu
     * @private */
   RH3Painter.prototype.fillHistContextMenu = function(menu) {

      let sett = jsrp.getDrawSettings("ROOT." + this.getObject()._typename, 'nosame');

      menu.addDrawMenu("Draw with", sett.opts, arg => {
         if (arg==='inspect')
            return this.showInspector();

         this.decodeOptions(arg);

         this.interactiveRedraw(true, "drawopt");
      });
   }

   let drawHist3 = (divid, histo /*, opt*/) => {
      // create painter and add it to canvas
      let painter = new RH3Painter(divid, histo);

      return jsrp.ensureRCanvas(painter, "3d").then(() => {

         painter.setAsMainPainter();

         painter.options = { Box: 0, Scatter: false, Sphere: 0, Color: false, minimum: -1111, maximum: -1111 };

         let kind = painter.v7EvalAttr("kind", ""),
             sub = painter.v7EvalAttr("sub", 0),
             o = painter.options;

         switch(kind) {
            case "box": o.Box = 10 + sub; break;
            case "sphere": o.Sphere = 10 + sub; break;
            case "col": o.Color = true; break;
            case "scat": o.Scatter = true;  break;
            default: o.Box = 10;
         }

         painter.scanContent();
         painter.redraw();
         return painter;
      });
   }

   JSROOT.v7.RH3Painter = RH3Painter;

   JSROOT.v7.drawHist3 = drawHist3;

   return JSROOT;
});
