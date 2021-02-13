/// @file JSRoot.geom.js
/// JavaScript ROOT 3D geometry painter

JSROOT.define(['d3', 'three', 'geobase', 'painter', 'base3d'], (d3, THREE, geo, jsrp) => {

   "use strict";

   JSROOT.loadScript('$$$style/JSRoot.geom');

   const _ENTIRE_SCENE = 0, _BLOOM_SCENE = 1;

   // ============================================================================================

   function Toolbar(container, bright) {
      this.bright = bright;
      this.element = container.append("div").attr('class','geo_toolbar_group');
   }

   Toolbar.prototype.addButtons = function(buttons) {
      this.buttonsNames = [];

      JSROOT.require(["interactive"]).then(inter => {
         buttons.forEach(buttonConfig => {
            let buttonName = buttonConfig.name;
            if (!buttonName) {
               throw new Error('must provide button \'name\' in button config');
            }
            if (this.buttonsNames.indexOf(buttonName) !== -1) {
               throw new Error('button name \'' + buttonName + '\' is taken');
            }

            this.buttonsNames.push(buttonName);

            let title = buttonConfig.title || buttonConfig.name;

            if (typeof buttonConfig.click !== 'function')
               throw new Error('must provide button \'click\' function in button config');

            let button = this.element.append('a')
                              .attr('class', this.bright ? 'geo_toolbar_btn_bright' : 'geo_toolbar_btn')
                              .attr('rel', 'tooltip')
                              .attr('data-title', title)
                              .on('click', buttonConfig.click);

            inter.ToolbarIcons.createSVG(button, inter.ToolbarIcons[buttonConfig.icon], 16, title);
         });

      });
   }

   Toolbar.prototype.changeBrightness = function(bright) {
      this.bright = bright;
      if (this.element)
         this.element.selectAll(bright ? '.geo_toolbar_btn' : ".geo_toolbar_btn_bright")
                     .attr("class", !bright ? 'geo_toolbar_btn' : "geo_toolbar_btn_bright");
   }

   Toolbar.prototype.cleanup = function() {
      if (this.element) {
         this.element.remove();
         delete this.element;
      }
   }

   //////////////////////

   function GeoDrawingControl(mesh, bloom) {
      jsrp.InteractiveControl.call(this);
      this.mesh = (mesh && mesh.material) ? mesh : null;
      this.bloom = bloom;
   }

   GeoDrawingControl.prototype = Object.create(jsrp.InteractiveControl.prototype);

   GeoDrawingControl.prototype.setHighlight = function(col, indx) {
      return this.drawSpecial(col, indx);
   }

   GeoDrawingControl.prototype.drawSpecial = function(col /*, indx*/) {
      let c = this.mesh;
      if (!c || !c.material) return;

      if (col) {
         if (!c.origin)
            c.origin = {
              color: c.material.color,
              emissive: c.material.emissive,
              opacity: c.material.opacity,
              width: c.material.linewidth,
              size: c.material.size
           };
         if (this.bloom) {
            c.layers.enable(_BLOOM_SCENE);
            c.material.emissive = new THREE.Color(0x00ff00);
         } else {
            c.material.color = new THREE.Color( col );
            c.material.opacity = 1.;
         }

         if (c.hightlightWidthScale && !JSROOT.browser.isWin)
            c.material.linewidth = c.origin.width * c.hightlightWidthScale;
         if (c.highlightScale)
            c.material.size = c.origin.size * c.highlightScale;
         return true;
      } else if (c.origin) {
         if (this.bloom) {
            c.material.emissive = c.origin.emissive;
            c.layers.enable(_ENTIRE_SCENE);
         } else {
            c.material.color = c.origin.color;
            c.material.opacity = c.origin.opacity;
         }
         if (c.hightlightWidthScale)
            c.material.linewidth = c.origin.width;
         if (c.highlightScale)
            c.material.size = c.origin.size;
         return true;
      }
   }

   ///////////////////////////////////////////////////////////////////////////////

   /**
    * @summary Painter class for geometries drawing
    *
    * @class
    * @memberof JSROOT
    * @summary Painter for TGeo object.
    * @param {object|string} dom - DOM element for drawing or element id
    * @param {object} obj - supported TGeo object
    * @private
    */

   function TGeoPainter(dom, obj) {

      if (obj && (obj._typename === "TGeoManager")) {
         this.geo_manager = obj;
         obj = obj.fMasterVolume;
      }

      if (obj && (obj._typename.indexOf('TGeoVolume') === 0))
         obj = { _typename:"TGeoNode", fVolume: obj, fName: obj.fName, $geoh: obj.$geoh, _proxy: true };

      JSROOT.ObjectPainter.call(this, dom, obj);

      this.no_default_title = true; // do not set title to main DIV
      this.mode3d = true; // indication of 3D mode
      this.drawing_stage = 0; //
      this.ctrl = {
         clipIntersect: true,
         clip: [{ name:"x", enabled: false, value: 0, min: -100, max: 100}, { name:"y", enabled: false, value: 0, min: -100, max: 100}, { name:"z", enabled: false, value: 0, min: -100, max: 100}],
         ssao: { enabled: false, output: THREE.SSAOPass.OUTPUT.Default, kernelRadius: 0, minDistance: 0.001, maxDistance: 0.1 },
         bloom: { enabled: true, strength: 1.5 },
         info: { num_meshes: 0, num_faces: 0, num_shapes: 0 },
         highlight: false,
         highlight_scene: false,
         depthTest: true,
         depthMethod: "dflt",
         select_in_view: false,
         update_browser: true,
         light: { kind: "points", top: false, bottom: false, left: false, right: false, front: false, specular: true, power: 1 },
         trans_radial: 0,
         trans_z: 0
      };

      this.ctrl.depthMethodItems = [
         {name: 'Default', value: "dflt"},
         {name: 'Raytraicing', value: "ray"},
         {name: 'Boundary box', value: "box"},
         {name: 'Mesh size', value: "size"},
         {name: 'Central point', value: "pnt" }
       ];

      this.ctrl.ssao.outputItems = [
         {name: 'Default', value: THREE.SSAOPass.OUTPUT.Default},
         {name: 'SSAO Only', value: THREE.SSAOPass.OUTPUT.SSAO},
         {name: 'SSAO Only + Blur', value: THREE.SSAOPass.OUTPUT.Blur},
         {name: 'Beauty', value: THREE.SSAOPass.OUTPUT.Beauty},
         {name: 'Depth', value: THREE.SSAOPass.OUTPUT.Depth},
         {name: 'Normal', value: THREE.SSAOPass.OUTPUT.Normal}
      ];

      this.cleanup(true);
   }

   TGeoPainter.prototype = Object.create(JSROOT.ObjectPainter.prototype);

   /** @summary Create toolbar */
   TGeoPainter.prototype.createToolbar = function() {
      if (this._toolbar || !this._webgl || this.ctrl.notoolbar || JSROOT.batch_mode) return;
      let buttonList = [{
         name: 'toImage',
         title: 'Save as PNG',
         icon: "camera",
         click: () => this.createSnapshot()
      }, {
         name: 'control',
         title: 'Toggle control UI',
         icon: "rect",
         click: () => this.showControlOptions('toggle')
      }, {
         name: 'enlarge',
         title: 'Enlarge geometry drawing',
         icon: "circle",
         click: () => this.toggleEnlarge()
      }];

      // Only show VR icon if WebVR API available.
      if (navigator.getVRDisplays) {
         buttonList.push({
            name: 'entervr',
            title: 'Enter VR (It requires a VR Headset connected)',
            icon: "vrgoggles",
            click: () => this.toggleVRMode()
         });
         this.initVRMode();
      }

      if (JSROOT.settings.ContextMenu)
      buttonList.push({
         name: 'menu',
         title: 'Show context menu',
         icon: "question",
         click: evnt => {

            evnt.preventDefault();
            evnt.stopPropagation();

            if (jsrp.closeMenu && jsrp.closeMenu()) return;

            jsrp.createMenu(evnt, this).then(menu => {
                menu.painter.fillContextMenu(menu);
                menu.show();
            });
         }
      });

      let bkgr = new THREE.Color(this.ctrl.background);

      this._toolbar = new Toolbar(this.selectDom(), (bkgr.r + bkgr.g + bkgr.b) < 1);

      this._toolbar.addButtons(buttonList);
   }

   /** @summary Initialize VR mode */
   TGeoPainter.prototype.initVRMode = function() {
      // Dolly contains camera and controllers in VR Mode
      // Allows moving the user in the scene
      this._dolly = new THREE.Group();
      this._scene.add(this._dolly);
      this._standingMatrix = new THREE.Matrix4();

      // Raycaster temp variables to avoid one per frame allocation.
      this._raycasterEnd = new THREE.Vector3();
      this._raycasterOrigin = new THREE.Vector3();

      navigator.getVRDisplays().then(displays => {
         let vrDisplay = displays[0];
         if (!vrDisplay) return;
         this._renderer.vr.setDevice(vrDisplay);
         this._vrDisplay = vrDisplay;
         if (vrDisplay.stageParameters) {
            this._standingMatrix.fromArray(vrDisplay.stageParameters.sittingToStandingTransform);
         }
         this.initVRControllersGeometry();
      });
   }

   /** @summary Init VR controllers geometry
     * @private */
   TGeoPainter.prototype.initVRControllersGeometry = function() {
      let geometry = new THREE.SphereGeometry(0.025, 18, 36);
      let material = new THREE.MeshBasicMaterial({color: 'grey'});
      let rayMaterial = new THREE.MeshBasicMaterial({color: 'fuchsia'});
      let rayGeometry = new THREE.BoxBufferGeometry(0.001, 0.001, 2);
      let ray1Mesh = new THREE.Mesh(rayGeometry, rayMaterial);
      let ray2Mesh = new THREE.Mesh(rayGeometry, rayMaterial);
      let sphere1 = new THREE.Mesh(geometry, material);
      let sphere2 = new THREE.Mesh(geometry, material);

      this._controllersMeshes = [];
      this._controllersMeshes.push(sphere1);
      this._controllersMeshes.push(sphere2);
      ray1Mesh.position.z -= 1;
      ray2Mesh.position.z -= 1;
      sphere1.add(ray1Mesh);
      sphere2.add(ray2Mesh);
      this._dolly.add(sphere1);
      this._dolly.add(sphere2);
      // Controller mesh hidden by default
      sphere1.visible = false;
      sphere2.visible = false;
   }

   /** @summary Update VR controllers list
     * @private */
   TGeoPainter.prototype.updateVRControllersList = function() {
      let gamepads = navigator.getGamepads && navigator.getGamepads();
      // Has controller list changed?
      if (this.vrControllers && (gamepads.length === this.vrControllers.length)) { return; }
      // Hide meshes.
      this._controllersMeshes.forEach(mesh => { mesh.visible = false; });
      this._vrControllers = [];
      for (let i = 0; i < gamepads.length; ++i) {
         if (!gamepads[i] || !gamepads[i].pose) { continue; }
         this._vrControllers.push({
            gamepad: gamepads[i],
            mesh: this._controllersMeshes[i]
         });
         this._controllersMeshes[i].visible = true;
      }
   }

   /** @summary Process VR controller intersection
     * @private */
   TGeoPainter.prototype.processVRControllerIntersections = function() {
      let intersects = []
      for (let i = 0; i < this._vrControllers.length; ++i) {
         let controller = this._vrControllers[i].mesh;
         let end = controller.localToWorld(this._raycasterEnd.set(0, 0, -1));
         let origin = controller.localToWorld(this._raycasterOrigin.set(0, 0, 0));
         end.sub(origin).normalize();
         intersects = intersects.concat(this._controls.getOriginDirectionIntersects(origin, end));
      }
      // Remove duplicates.
      intersects = intersects.filter(function (item, pos) {return intersects.indexOf(item) === pos});
      this._controls.ProcessMouseMove(intersects);
   }

   /** @summary Update VR controllers
     * @private */
   TGeoPainter.prototype.updateVRControllers = function() {
      this.updateVRControllersList();
      // Update pose.
      for (let i = 0; i < this._vrControllers.length; ++i) {
         let controller = this._vrControllers[i];
         let orientation = controller.gamepad.pose.orientation;
         let position = controller.gamepad.pose.position;
         let controllerMesh = controller.mesh;
         if (orientation) { controllerMesh.quaternion.fromArray(orientation); }
         if (position) { controllerMesh.position.fromArray(position); }
         controllerMesh.updateMatrix();
         controllerMesh.applyMatrix4(this._standingMatrix);
         controllerMesh.matrixWorldNeedsUpdate = true;
      }
      this.processVRControllerIntersections();
   }

   /** @summary Toggle VR mode
     * @private */
   TGeoPainter.prototype.toggleVRMode = function() {
      if (!this._vrDisplay) return;
      // Toggle VR mode off
      if (this._vrDisplay.isPresenting) {
         this.exitVRMode();
         return;
      }
      this._previousCameraPosition = this._camera.position.clone();
      this._previousCameraRotation = this._camera.rotation.clone();
      this._vrDisplay.requestPresent([{ source: this._renderer.domElement }]).then(() => {
         this._previousCameraNear = this._camera.near;
         this._dolly.position.set(this._camera.position.x/4, - this._camera.position.y/8, - this._camera.position.z/4);
         this._camera.position.set(0,0,0);
         this._dolly.add(this._camera);
         this._camera.near = 0.1;
         this._camera.updateProjectionMatrix();
         this._renderer.vr.enabled = true;
         this._renderer.setAnimationLoop(() => {
            this.updateVRControllers();
            this.render3D(0);
         });
      });
      this._renderer.vr.enabled = true;

      window.addEventListener( 'keydown', event => {
         // Esc Key turns VR mode off
         if (event.keyCode === 27) this.exitVRMode();
      });
   }

   /** @summary Exit VR mode
     * @private */
   TGeoPainter.prototype.exitVRMode = function() {
      if (!this._vrDisplay.isPresenting) return;
      this._renderer.vr.enabled = false;
      this._dolly.remove(this._camera);
      this._scene.add(this._camera);
      // Restore Camera pose
      this._camera.position.copy(this._previousCameraPosition);
      this._previousCameraPosition = undefined;
      this._camera.rotation.copy(this._previousCameraRotation);
      this._previousCameraRotation = undefined;
      this._camera.near = this._previousCameraNear;
      this._camera.updateProjectionMatrix();
      this._vrDisplay.exitPresent();
   }

   /** @summary Returns main geometry object */
   TGeoPainter.prototype.getGeometry = function() {
      return this.getObject();
   }

   /** @summary Modify visibility of provided node by name */
   TGeoPainter.prototype.modifyVisisbility = function(name, sign) {
      if (geo.getNodeKind(this.getGeometry()) !== 0) return;

      if (name == "")
         return geo.SetBit(this.getGeometry().fVolume, geo.BITS.kVisThis, (sign === "+"));

      let regexp, exact = false;

      //arg.node.fVolume
      if (name.indexOf("*") < 0) {
         regexp = new RegExp("^"+name+"$");
         exact = true;
      } else {
         regexp = new RegExp("^" + name.split("*").join(".*") + "$");
         exact = false;
      }

      this.findNodeWithVolume(regexp, function(arg) {
         geo.setInvisibleAll(arg.node.fVolume, (sign !== "+"));
         return exact ? arg : null; // continue search if not exact expression provided
      });
   }

   /** @summary Decode drawing options */
   TGeoPainter.prototype.decodeOptions = function(opt) {
      if (typeof opt != "string") opt = "";

      let res = { _grid: false, _bound: false, _debug: false,
                  _full: false, _axis: 0,
                  _count: false, wireframe: false,
                   scale: new THREE.Vector3(1,1,1), zoom: 1.0, rotatey: 0, rotatez: 0,
                   more: 1, maxlimit: 100000,
                   vislevel: undefined, maxnodes: undefined, dflt_colors: false,
                   use_worker: false, show_controls: false,
                   highlight: false, highlight_scene: false, no_screen: false,
                   project: '', is_main: false, tracks: false, showtop: false, can_rotate: true, ortho_camera: false,
                   clipx: false, clipy: false, clipz: false, usessao: false, usebloom: true, outline: false,
                   script_name: "", transparency: 0, rotate: false, background: '#FFFFFF',
                   depthMethod: "dflt", mouse_tmout: 50, trans_radial: 0, trans_z: 0 };

      let dd = JSROOT.decodeUrl();
      if (dd.get('_grid') == "true") res._grid = true;
      let _opt = dd.get('_debug');
      if (_opt == "true") { res._debug = true; res._grid = true; }
      if (_opt == "bound") { res._debug = true; res._grid = true; res._bound = true; }
      if (_opt == "full") { res._debug = true; res._grid = true; res._full = true; res._bound = true; }

      let macro = opt.indexOf("macro:");
      if (macro >= 0) {
         let separ = opt.indexOf(";", macro+6);
         if (separ<0) separ = opt.length;
         res.script_name = opt.substr(macro+6,separ-macro-6);
         opt = opt.substr(0, macro) + opt.substr(separ+1);
         console.log('script', res.script_name, 'rest', opt);
      }

      while (true) {
         let pp = opt.indexOf("+"), pm = opt.indexOf("-");
         if ((pp<0) && (pm<0)) break;
         let p1 = pp, sign = "+";
         if ((p1<0) || ((pm>=0) && (pm<pp))) { p1 = pm; sign = "-"; }

         let p2 = p1+1, regexp = new RegExp('[,; .]');
         while ((p2<opt.length) && !regexp.test(opt[p2]) && (opt[p2]!='+') && (opt[p2]!='-')) p2++;

         let name = opt.substring(p1+1, p2);
         opt = opt.substr(0,p1) + opt.substr(p2);
         // console.log("Modify visibility", sign,':',name);

         this.modifyVisisbility(name, sign);
      }

      let d = new JSROOT.DrawOptions(opt);

      if (d.check("MAIN")) res.is_main = true;

      if (d.check("TRACKS")) res.tracks = true; // only for TGeoManager
      if (d.check("SHOWTOP")) res.showtop = true; // only for TGeoManager
      if (d.check("NO_SCREEN")) res.no_screen = true; // use kVisOnScreen bits as visibility

      if (d.check("ORTHO_CAMERA_ROTATE")) { res.ortho_camera = true; res.can_rotate = true; }
      if (d.check("ORTHO_CAMERA")) { res.ortho_camera = true; res.can_rotate = false; }
      if (d.check("MOUSE_CLICK")) res.mouse_click = true;

      if (d.check("DEPTHRAY") || d.check("DRAY")) res.depthMethod = "ray";
      if (d.check("DEPTHBOX") || d.check("DBOX")) res.depthMethod = "box";
      if (d.check("DEPTHPNT") || d.check("DPNT")) res.depthMethod = "pnt";
      if (d.check("DEPTHSIZE") || d.check("DSIZE")) res.depthMethod = "size";
      if (d.check("DEPTHDFLT") || d.check("DDFLT")) res.depthMethod = "dflt";

      if (d.check("ZOOM", true)) res.zoom = d.partAsFloat(0, 100) / 100;
      if (d.check("ROTY", true)) res.rotatey = d.partAsFloat();
      if (d.check("ROTZ", true)) res.rotatez = d.partAsFloat();
      if (d.check("VISLVL", true)) res.vislevel = d.partAsInt();

      if (d.check('BLACK')) res.background = "#000000";
      if (d.check('WHITE')) res.background = "#FFFFFF";

      if (d.check('BKGR_', true)) {
         let bckgr = null;
         if (d.partAsInt(1) > 0) {
           bckgr = jsrp.getColor(d.partAsInt());
         } else {
            for (let col=0;col<8;++col)
               if (jsrp.getColor(col).toUpperCase() === d.part)
                  bckgr = jsrp.getColor(col);
         }
         if (bckgr) res.background = "#" + new THREE.Color(bckgr).getHexString();
      }

      if (d.check('R3D_', true))
         res.Render3D = JSROOT.constants.Render3D.fromString(d.part.toLowerCase());

      if (d.check("MORE3")) res.more = 3;
      if (d.check("MORE")) res.more = 2;
      if (d.check("ALL")) { res.more = 10; res.vislevel = 9; }

      if (d.check("CONTROLS") || d.check("CTRL")) res.show_controls = true;

      if (d.check("CLIPXYZ")) res.clipx = res.clipy = res.clipz = true;
      if (d.check("CLIPX")) res.clipx = true;
      if (d.check("CLIPY")) res.clipy = true;
      if (d.check("CLIPZ")) res.clipz = true;
      if (d.check("CLIP")) res.clipx = res.clipy = res.clipz = true;

      if (d.check("PROJX", true)) { res.project = 'x'; if (d.partAsInt(1)>0) res.projectPos = d.partAsInt(); res.can_rotate = false; }
      if (d.check("PROJY", true)) { res.project = 'y'; if (d.partAsInt(1)>0) res.projectPos = d.partAsInt(); res.can_rotate = false; }
      if (d.check("PROJZ", true)) { res.project = 'z'; if (d.partAsInt(1)>0) res.projectPos = d.partAsInt(); res.can_rotate = false; }

      if (d.check("DFLT_COLORS") || d.check("DFLT")) res.dflt_colors = true;
      if (d.check("SSAO")) res.usessao = true;
      if (d.check("NOBLOOM")) res.usebloom = false;
      if (d.check("BLOOM")) res.usebloom = true;
      if (d.check("OUTLINE")) res.outline = true;

      if (d.check("NOWORKER")) res.use_worker = -1;
      if (d.check("WORKER")) res.use_worker = 1;

      if (d.check("NOHIGHLIGHT") || d.check("NOHIGH")) res.highlight_scene = res.highlight = 0;
      if (d.check("HIGHLIGHT")) res.highlight_scene = res.highlight = true;
      if (d.check("HSCENEONLY")) { res.highlight_scene = true; res.highlight = 0; }
      if (d.check("NOHSCENE")) res.highlight_scene = 0;
      if (d.check("HSCENE")) res.highlight_scene = true;

      if (d.check("WIREFRAME") || d.check("WIRE")) res.wireframe = true;
      if (d.check("ROTATE")) res.rotate = true;

      if (d.check("INVX") || d.check("INVERTX")) res.scale.x = -1;
      if (d.check("INVY") || d.check("INVERTY")) res.scale.y = -1;
      if (d.check("INVZ") || d.check("INVERTZ")) res.scale.z = -1;

      if (d.check("COUNT")) res._count = true;

      if (d.check('TRANSP',true))
         res.transparency = d.partAsInt(0,100)/100;

      if (d.check('OPACITY',true))
         res.transparency = 1 - d.partAsInt(0,100)/100;

      if (d.check("AXISCENTER") || d.check("AC")) res._axis = 2;

      if (d.check('TRR',true)) res.trans_radial = d.partAsInt()/100;
      if (d.check('TRZ',true)) res.trans_z = d.partAsInt()/100;

      if (d.check("AXIS") || d.check("A")) res._axis = true;

      if (d.check("D")) res._debug = true;
      if (d.check("G")) res._grid = true;
      if (d.check("B")) res._bound = true;
      if (d.check("W")) res.wireframe = true;
      if (d.check("F")) res._full = true;
      if (d.check("Y")) res._yup = true;
      if (d.check("Z")) res._yup = false;

      // when drawing geometry without TCanvas, yup = true by default
      if (res._yup === undefined)
         res._yup = this.getCanvSvg().empty();

      return res;
   }

   /** @summary Activate specified items in the browser */
   TGeoPainter.prototype.activateInBrowser = function(names, force) {
      // if (this.getItemName() === null) return;

      if (typeof names == 'string') names = [ names ];

      if (this._hpainter) {
         // show browser if it not visible

         this._hpainter.activateItems(names, force);

         // if highlight in the browser disabled, suppress in few seconds
         if (!this.ctrl.update_browser)
            setTimeout(() => this._hpainter.activateItems([]), 2000);
      }
   }

   /** @summary  method used to check matrix calculations performance with current three.js model */
   TGeoPainter.prototype.testMatrixes = function() {

      let errcnt = 0, totalcnt = 0, totalmax = 0;

      let arg = {
            domatrix: true,
            func: (/*node*/) => {

               let m2 = this.getmatrix();

               let entry = this.CopyStack();

               let mesh = this._clones.createObject3D(entry.stack, this._toplevel, 'mesh');

               if (!mesh) return true;

               totalcnt++;

               let m1 = mesh.matrixWorld, flip;

               if (m1.equals(m2)) return true;
               if ((m1.determinant()>0) && (m2.determinant()<-0.9)) {
                  flip = THREE.Vector3(1,1,-1);
                  m2 = m2.clone().scale(flip);
                  if (m1.equals(m2)) return true;
               }

               let max = 0;
               for (let k=0;k<16;++k)
                  max = Math.max(max, Math.abs(m1.elements[k] - m2.elements[k]));

               totalmax = Math.max(max, totalmax);

               if (max < 1e-4) return true;

               console.log(this._clones.resolveStack(entry.stack).name, 'maxdiff', max, 'determ', m1.determinant(), m2.determinant());

               errcnt++;

               return false;
            }
         };


      let tm1 = new Date().getTime();

      /* let cnt = */ this._clones.scanVisible(arg);

      let tm2 = new Date().getTime();

      console.log('Compare matrixes total', totalcnt, 'errors', errcnt, 'takes', tm2-tm1, 'maxdiff', totalmax);
   }

   /** @summary Fills context menu */
   TGeoPainter.prototype.fillContextMenu = function(menu) {
      menu.add("header: Draw options");

      menu.addchk(this.ctrl.update_browser, "Browser update", () => {
         this.ctrl.update_browser = !this.ctrl.update_browser;
         if (!this.ctrl.update_browser) this.activateInBrowser([]);
      });
      menu.addchk(this.ctrl.show_controls, "Show Controls", () => this.showControlOptions('toggle'));

      menu.addchk(this.ctrl._axis, "Show axes", () => this.setAxesDraw('toggle'));

      if (this.geo_manager)
         menu.addchk(this.ctrl.showtop, "Show top volume", () => this.setShowTop(!this.ctrl.showtop));

      menu.addchk(this.ctrl.wireframe, "Wire frame", () => this.toggleWireFrame());

      menu.addchk(this.ctrl.highlight, "Highlight volumes", () => {
         this.ctrl.highlight = !this.ctrl.highlight;
      });
      menu.addchk(this.ctrl.highlight_scene, "Highlight scene", () => {
         this.ctrl.highlight_scene = !this.ctrl.highlight_scene;
      });
      menu.add("Reset camera position", () => this.focusCamera());

      if (!this._geom_viewer)
         menu.add("Get camera position", () => alert("Position (as url): &opt=" + this.produceCameraUrl()));

      if (!this.ctrl.project)
         menu.addchk(this.ctrl.rotate, "Autorotate", () => this.setAutoRotate(!this.ctrl.rotate));
      menu.addchk(this.ctrl.select_in_view, "Select in view", () => {
         this.ctrl.select_in_view = !this.ctrl.select_in_view;
         if (this.ctrl.select_in_view) this.startDrawGeometry();
      });
   }

   /** @summary Method used to set transparency for all geometrical shapes
    * @param {number|Function} transparency - one could provide function
    * @param {boolean} [skip_render] - if specified, do not perform rendering */
   TGeoPainter.prototype.changedGlobalTransparency = function(transparency, skip_render) {
      let func = (typeof transparency == 'function') ? transparency : null;
      if (func || (transparency === undefined)) transparency = this.ctrl.transparency;
      this._toplevel.traverse( function (node) {
         if (node && node.material && (node.material.inherentOpacity !== undefined)) {
            let t = func ? func(node) : undefined;
            if (t !== undefined)
               node.material.opacity = 1 - t;
            else
               node.material.opacity = Math.min(1 - (transparency || 0), node.material.inherentOpacity);
            node.material.transparent = node.material.opacity < 1;
         }
      });
      if (!skip_render) this.render3D(-1);
   }

   /** @summary Reset transformation */
   TGeoPainter.prototype.resetTransformation = function() {
      this.changedTransformation("reset");
   }

   /** @summary Method should be called when transformation parameters were changed */
   TGeoPainter.prototype.changedTransformation = function(arg) {
      if (!this._toplevel) return;

      let ctrl = this.ctrl,
          translation = new THREE.Matrix4(),
          vect2 = new THREE.Vector3();

      if (arg == "reset")
         ctrl.trans_z = ctrl.trans_radial = 0;

      this._toplevel.traverse(function(mesh) {
         if (mesh.stack === undefined) return;

         let node = mesh.parent;

         if (arg == "reset") {
            if (node.matrix0) {
               node.matrix.copy(node.matrix0);
               node.matrix.decompose( node.position, node.quaternion, node.scale );
               node.matrixWorldNeedsUpdate = true;
            }
            delete node.matrix0;
            delete node.vect0;
            delete node.vect1;
            delete node.minvert;
            return;
         }

         if (node.vect0 === undefined) {
            node.matrix0 = node.matrix.clone();
            node.minvert = new THREE.Matrix4().getInverse( node.matrixWorld );

            let box3 = geo.getBoundingBox(mesh, null, true),
                signz = mesh._flippedMesh ? -1 : 1;

            // real center of mesh in local coordinates
            node.vect0 = new THREE.Vector3((box3.max.x  + box3.min.x) / 2, (box3.max.y  + box3.min.y) / 2, signz * (box3.max.z  + box3.min.z) / 2).applyMatrix4(node.matrixWorld);
            node.vect1 = new THREE.Vector3(0,0,0).applyMatrix4(node.minvert);
         }

         vect2.set(ctrl.trans_radial * node.vect0.x, ctrl.trans_radial * node.vect0.y, ctrl.trans_z * node.vect0.z).applyMatrix4(node.minvert).sub(node.vect1);

         node.matrix.multiplyMatrices(node.matrix0, translation.makeTranslation(vect2.x, vect2.y, vect2.z));
         node.matrix.decompose( node.position, node.quaternion, node.scale );
         node.matrixWorldNeedsUpdate = true;
      });

      this._toplevel.updateMatrixWorld();

      // axes drawing always triggers rendering
      if (arg != "norender")
         this.drawSimpleAxis();
   }

   /** @summary Should be called when autorotate property changed */
   TGeoPainter.prototype.changedAutoRotate = function() {
      this.autorotate(2.5);
   }

   /** @summary Method should be called when changing axes drawing */
   TGeoPainter.prototype.changedAxes = function() {
      if (typeof this.ctrl._axis == 'string')
         this.ctrl._axis = parseInt(this.ctrl._axis);

      this.drawSimpleAxis();
   }

   /** @summary Method should be called to change background color */
   TGeoPainter.prototype.changedBackground = function(val) {
      if (val !== undefined) this.ctrl.background = val;
      this._renderer.setClearColor(this.ctrl.background, 1);
      this.render3D(0);

      if (this._toolbar) {
         let bkgr = new THREE.Color(this.ctrl.background);
         this._toolbar.changeBrightness((bkgr.r + bkgr.g + bkgr.b) < 1);
      }
   }

   /** @summary Method called when SSAO configuration changed via GUI */
   TGeoPainter.prototype.changedSSAO = function() {
      if (!this.ctrl.ssao.enabled) {
         this.removeSSAO();
      } else {
         this.createSSAO();

         this._ssaoPass.output = parseInt(this.ctrl.ssao.output);
         this._ssaoPass.kernelRadius = this.ctrl.ssao.kernelRadius;
         this._ssaoPass.minDistance = this.ctrl.ssao.minDistance;
         this._ssaoPass.maxDistance = this.ctrl.ssao.maxDistance;
      }

      this.updateClipping();

      if (this._slave_painters)
         this._slave_painters.forEach(p => {
            JSROOT.extend(p.ctrl.ssao, this.ctrl.ssao);
            p.changedSSAO();
         });
   }

   /** @summary Display control GUI */
   TGeoPainter.prototype.showControlOptions = function(on) {
      // while complete geo drawing can be removed until dat is loaded - just check and ignore callback
      if (!this.ctrl) return;

      if (on === 'toggle') {
         on = !this._datgui;
      } else if (on === undefined) {
         on = this.ctrl.show_controls;
      }

      this.ctrl.show_controls = on;

      if (this._datgui) {
         if (!on) {
            d3.select(this._datgui.domElement).remove();
            this._datgui.destroy();
            delete this._datgui;
         }
         return;
      }

      if (on)
         JSROOT.require('dat.gui')
               .then(dat => this.buildDatGui(dat));
   }

   /** @summary build dat.gui elements
     * @private */
   TGeoPainter.prototype.buildDatGui = function(dat) {
      // can happen when dat gui loaded after drawing is already cleaned
      if (!this._renderer) return;

      if (!dat)
         throw Error('Fail to load dat.gui');

      this._datgui = new dat.GUI({ autoPlace: false, width: Math.min(650, this._renderer.domElement.width / 2) });

      let main = this.selectDom();
      if (main.style('position')=='static') main.style('position','relative');

      d3.select(this._datgui.domElement)
               .style('position','absolute')
               .style('top',0).style('right',0);

      main.node().appendChild(this._datgui.domElement);

      this._datgui.painter = this;

      if (this.ctrl.project) {

         let bound = this.getGeomBoundingBox(this.getProjectionSource(), 0.01);

         let axis = this.ctrl.project;

         if (this.ctrl.projectPos === undefined)
            this.ctrl.projectPos = (bound.min[axis] + bound.max[axis])/2;

         this._datgui.add(this.ctrl, 'projectPos', bound.min[axis], bound.max[axis])
             .name(axis.toUpperCase() + ' projection')
             .onChange(this.startDrawGeometry.bind(this));

      } else {
         // Clipping Options

         let clipFolder = this._datgui.addFolder('Clipping'),
             clip_handler = this.changedClipping.bind(this, -1);

         for (let naxis=0;naxis<3;++naxis) {
            let cc = this.ctrl.clip[naxis],
                axisC = cc.name.toUpperCase();

            clipFolder.add(cc, 'enabled')
                .name('Enable ' + axisC)
                .listen() // react if option changed outside
                .onChange(clip_handler);

            clipFolder.add(cc, "value", cc.min, cc.max)
                .name(axisC + ' position')
                .onChange(this.changedClipping.bind(this, naxis));
         }

         clipFolder.add(this.ctrl, 'clipIntersect').name("Clip intersection")
                   .listen().onChange(clip_handler);

      }

      // Appearance Options

      let appearance = this._datgui.addFolder('Appearance');

      appearance.add(this.ctrl, 'highlight').name('Highlight Selection')
                .listen().onChange(this.changedHighlight.bind(this));

      appearance.add(this.ctrl, 'transparency', 0.0, 1.0, 0.001)
                     .listen().onChange(this.changedGlobalTransparency.bind(this));

      appearance.addColor(this.ctrl, 'background').name('Background')
                .onChange(this.changedBackground.bind(this));

      appearance.add(this.ctrl, 'wireframe').name('Wireframe')
                     .listen().onChange(this.changedWireFrame.bind(this));

      this.ctrl._axis_cfg = 0;
      appearance.add(this.ctrl, '_axis', { "none" : 0, "show": 1, "center": 2}).name('Axes')
                    .onChange(this.changedAxes.bind(this));

      if (!this.ctrl.project)
         appearance.add(this.ctrl, 'rotate').name("Autorotate")
                      .listen().onChange(this.changedAutoRotate.bind(this));

      appearance.add(this, 'focusCamera').name('Reset camera position');

      // Advanced Options

      if (this._webgl) {
         let advanced = this._datgui.addFolder('Advanced'), depthcfg = {};
         this.ctrl.depthMethodItems.forEach(i => { depthcfg[i.name] = i.value; });

         advanced.add(this.ctrl, 'depthTest').name("Depth test")
            .listen().onChange(this.changedDepthTest.bind(this));

         advanced.add( this.ctrl, 'depthMethod', depthcfg)
             .name("Rendering order")
             .onChange(this.changedDepthMethod.bind(this));

        advanced.add(this, 'resetAdvanced').name('Reset');
      }

      // Transformation Options
      if (!this.ctrl.project) {
         let transform = this._datgui.addFolder('Transform');
         transform.add(this.ctrl, 'trans_z', 0., 3., 0.01)
                     .name('Z axis')
                     .listen().onChange(this.changedTransformation.bind(this));
         transform.add(this.ctrl, 'trans_radial', 0., 3., 0.01)
                  .name('Radial')
                  .listen().onChange(this.changedTransformation.bind(this));

         transform.add(this, 'resetTransformation').name('Reset');

         if (this.ctrl.trans_z || this.ctrl.trans_radial) transform.open();
      }

      // no SSAO folder if outline is enabled
      if (this.ctrl.outline) return;

      let ssaofolder = this._datgui.addFolder('Smooth Lighting (SSAO)'),
          ssao_handler = this.changedSSAO.bind(this), ssaocfg = {};

      this.ctrl.ssao.outputItems.forEach(i => { ssaocfg[i.name] = i.value; });

      ssaofolder.add(this.ctrl.ssao, 'enabled').name('Enable SSAO')
                .listen().onChange(ssao_handler);

      ssaofolder.add( this.ctrl.ssao, 'output', ssaocfg)
                .listen().onChange(ssao_handler);

      ssaofolder.add( this.ctrl.ssao, 'kernelRadius', 0, 32)
                .listen().onChange(ssao_handler);

      ssaofolder.add( this.ctrl.ssao, 'minDistance', 0.001, 0.02)
                .listen().onChange(ssao_handler);

      ssaofolder.add( this.ctrl.ssao, 'maxDistance', 0.01, 0.3)
                .listen().onChange(ssao_handler);

      let blooming = this._datgui.addFolder('Unreal Bloom');
      let bloom_handler = this.changedBloomSettings.bind(this);

      blooming.add(this.ctrl.bloom, 'enabled').name('Enable Blooming')
                .listen().onChange(bloom_handler);

      blooming.add( this.ctrl.bloom, 'strength', 0.0, 3.0).name("Strength")
            .listen().onChange(bloom_handler);
   }

   /** @summary Method called when bloom configuration changed via GUI */
   TGeoPainter.prototype.changedBloomSettings = function() {
      if (this.ctrl.bloom.enabled) {
         this.createBloom();
         this._bloomPass.strength = this.ctrl.bloom.strength;
      } else {
         this.removeBloom();
      }

      if (this._slave_painters)
         this._slave_painters.forEach(p => {
            JSROOT.extend(p.ctrl.bloom, this.ctrl.bloom);
            p.changedBloomSettings();
         });
   }

   TGeoPainter.prototype.createBloom = function() {
      if (this._bloomPass) return;

      this._camera.layers.enable( _BLOOM_SCENE );
      this._bloomComposer = new THREE.EffectComposer( this._renderer );
      this._bloomComposer.addPass( new THREE.RenderPass( this._scene, this._camera ) );
      this._bloomPass = new THREE.UnrealBloomPass( new THREE.Vector2( window.innerWidth, window.innerHeight ), 1.5, 0.4, 0.85 );
      this._bloomPass.threshold = 0;
      this._bloomPass.strength = this.ctrl.bloom.strength;
      this._bloomPass.radius = 0;
      this._bloomPass.renderToScreen = true;
      this._bloomComposer.addPass( this._bloomPass );
      this._renderer.autoClear = false;
   }

   TGeoPainter.prototype.removeBloom = function() {
      if (!this._bloomPass) return;
      delete this._bloomPass;
      delete this._bloomComposer;
      this._renderer.autoClear = true;
      this._camera.layers.disable( _BLOOM_SCENE );
   }

   TGeoPainter.prototype.removeSSAO = function() {
      // we cannot remove pass from composer - just disable it
      delete this._ssaoPass;
      delete this._effectComposer;
   }

   TGeoPainter.prototype.createSSAO = function() {
      if (!this._webgl) return;

      // this._depthRenderTarget = new THREE.WebGLRenderTarget( this._scene_width, this._scene_height, { minFilter: THREE.LinearFilter, magFilter: THREE.LinearFilter } );
      // Setup SSAO pass
      if (!this._ssaoPass) {
         if (!this._effectComposer) {
            this._effectComposer = new THREE.EffectComposer( this._renderer );
            this._effectComposer.addPass( new THREE.RenderPass( this._scene, this._camera ) );
         }

         this._ssaoPass = new THREE.SSAOPass( this._scene, this._camera, this._scene_width, this._scene_height );
         this._ssaoPass.kernelRadius = 16;
         this._ssaoPass.renderToScreen = true;

         // Add pass to effect composer
         this._effectComposer.addPass( this._ssaoPass );
      }
   }

   /** @summary Show context menu for orbit control
     * @private */
   TGeoPainter.prototype.orbitContext = function(evnt, intersects) {

      jsrp.createMenu(evnt, this).then(menu => {
         let numitems = 0, numnodes = 0, cnt = 0;
         if (intersects)
            for (let n=0;n<intersects.length;++n) {
               if (intersects[n].object.stack) numnodes++;
               if (intersects[n].object.geo_name) numitems++;
            }

         if (numnodes + numitems === 0) {
            this.fillContextMenu(menu);
         } else {
            let many = (numnodes + numitems) > 1;

            if (many) menu.add("header:" + ((numitems > 0) ? "Items" : "Nodes"));

            for (let n=0;n<intersects.length;++n) {
               let obj = intersects[n].object,
                   name, itemname, hdr;

               if (obj.geo_name) {
                  itemname = obj.geo_name;
                  if (itemname.indexOf("<prnt>") == 0)
                     itemname = (this.getItemName() || "top") + itemname.substr(6);
                  name = itemname.substr(itemname.lastIndexOf("/")+1);
                  if (!name) name = itemname;
                  hdr = name;
               } else if (obj.stack) {
                  name = this._clones.resolveStack(obj.stack).name;
                  itemname = this.getStackFullName(obj.stack);
                  hdr = this.getItemName();
                  if (name.indexOf("Nodes/") === 0) hdr = name.substr(6); else
                  if (name.length > 0) hdr = name; else
                  if (!hdr) hdr = "header";

               } else
                  continue;

               menu.add((many ? "sub:" : "header:") + hdr, itemname, arg => this.activateInBrowser([arg], true));

               menu.add("Browse", itemname, arg => this.activateInBrowser([arg], true));

               if (this._hpainter)
                  menu.add("Inspect", itemname, arg => this._hpainter.display(arg, "inspect"));

               if (obj.geo_name) {
                  menu.add("Hide", n, indx => {
                     let mesh = intersects[indx].object;
                     mesh.visible = false; // just disable mesh
                     if (mesh.geo_object) mesh.geo_object.$hidden_via_menu = true; // and hide object for further redraw
                     menu.painter.render3D();
                  });

                  if (many) menu.add("endsub:");

                  continue;
               }

               let wireframe = this.accessObjectWireFrame(obj);

               if (wireframe!==undefined)
                  menu.addchk(wireframe, "Wireframe", n, function(indx) {
                     let m = intersects[indx].object.material;
                     m.wireframe = !m.wireframe;
                     this.render3D();
                  });

               if (++cnt>1)
                  menu.add("Manifest", n, function(indx) {

                     if (this._last_manifest)
                        this._last_manifest.wireframe = !this._last_manifest.wireframe;

                     if (this._last_hidden)
                        this._last_hidden.forEach(obj => { obj.visible = true; });

                     this._last_hidden = [];

                     for (let i=0;i<indx;++i)
                        this._last_hidden.push(intersects[i].object);

                     this._last_hidden.forEach(obj => { obj.visible = false; });

                     this._last_manifest = intersects[indx].object.material;

                     this._last_manifest.wireframe = !this._last_manifest.wireframe;

                     this.render3D();
                  });


               menu.add("Focus", n, function(indx) {
                  this.focusCamera(intersects[indx].object);
               });

               if (!this._geom_viewer)
               menu.add("Hide", n, function(indx) {
                  let resolve = menu.painter._clones.resolveStack(intersects[indx].object.stack);
                  const kindGeo = 0, kindEve = 1;
                  if (resolve.obj && (resolve.node.kind === kindGeo) && resolve.obj.fVolume) {
                     geo.SetBit(resolve.obj.fVolume, geo.BITS.kVisThis, false);
                     geo.updateBrowserIcons(resolve.obj.fVolume, this._hpainter);
                  } else if (resolve.obj && (resolve.node.kind === kindEve)) {
                     resolve.obj.fRnrSelf = false;
                     geo.updateBrowserIcons(resolve.obj, this._hpainter);
                  }

                  this.testGeomChanges();// while many volumes may disappear, recheck all of them
               });

               if (many) menu.add("endsub:");
            }
         }
         menu.show();
      });
   }

   /** @summary Filter some objects from three.js intersects array */
   TGeoPainter.prototype.filterIntersects = function(intersects) {

      if (!intersects.length) return intersects;

      // check redirections
      for (let n=0;n<intersects.length;++n)
         if (intersects[n].object.geo_highlight)
            intersects[n].object = intersects[n].object.geo_highlight;

      // remove all elements without stack - indicator that this is geometry object
      // also remove all objects which are mostly transparent
      for (let n=intersects.length-1; n>=0; --n) {

         let obj = intersects[n].object;

         let unique = (obj.stack !== undefined) || (obj.geo_name !== undefined);

         if (unique && obj.material && (obj.material.opacity !== undefined))
            unique = (obj.material.opacity >= 0.1);

         if (obj.jsroot_special) unique = false;

         for (let k=0;(k<n) && unique;++k)
            if (intersects[k].object === obj) unique = false;

         if (!unique) intersects.splice(n,1);
      }

      let clip = this.ctrl.clip;

      if (clip[0].enabled || clip[1].enabled || clip[2].enabled) {
         let clippedIntersects = [];

         for (let i = 0; i < intersects.length; ++i) {
            let point = intersects[i].point, special = (intersects[i].object.type == "Points"), clipped = true;

            if (clip[0].enabled && ((this._clipPlanes[0].normal.dot(point) > this._clipPlanes[0].constant) ^ special)) clipped = false;
            if (clip[1].enabled && ((this._clipPlanes[1].normal.dot(point) > this._clipPlanes[1].constant) ^ special)) clipped = false;
            if (clip[2].enabled && (this._clipPlanes[2].normal.dot(point) > this._clipPlanes[2].constant)) clipped = false;

            if (!clipped) clippedIntersects.push(intersects[i]);
         }

         intersects = clippedIntersects;
      }

      return intersects;
   }

    /** @summary test camera position
      * @desc function analyzes camera position and start redraw of geometry
        if objects in view may be changed */
   TGeoPainter.prototype.testCameraPositionChange = function() {

      if (!this.ctrl.select_in_view || this._draw_all_nodes) return;

      let matrix = geo.createProjectionMatrix(this._camera);

      let frustum = geo.createFrustum(matrix);

      // check if overall bounding box seen
      if (!frustum.CheckBox(this.getGeomBoundingBox(this._toplevel)))
         this.startDrawGeometry();
   }

   /** @summary Resolve stack */
   TGeoPainter.prototype.resolveStack = function(stack) {
      return this._clones && stack ? this._clones.resolveStack(stack) : null;
   }

   /** @summary Returns stack full name
     * @desc Includes item name of top geo object */
   TGeoPainter.prototype.getStackFullName = function(stack) {
      let mainitemname = this.getItemName(),
          sub = this.resolveStack(stack);
      if (!sub || !sub.name) return mainitemname;
      return mainitemname ? (mainitemname + "/" + sub.name) : sub.name;
   }

   /** @summary Add handler which will be called when element is highlighted in geometry drawing
    * @desc Handler should have highlightMesh function with same arguments as TGeoPainter  */
   TGeoPainter.prototype.addHighlightHandler = function(handler) {
      if (!handler || typeof handler.highlightMesh != 'function') return;
      if (!this._highlight_handlers) this._highlight_handlers = [];
      this._highlight_handlers.push(handler);
   }

   /** @summary perform mesh highlight */
   TGeoPainter.prototype.highlightMesh = function(active_mesh, color, geo_object, geo_index, geo_stack, no_recursive) {

      if (geo_object) {
         active_mesh = active_mesh ? [ active_mesh ] : [];
         let extras = this.getExtrasContainer();
         if (extras)
            extras.traverse(obj3d => {
               if ((obj3d.geo_object === geo_object) && (active_mesh.indexOf(obj3d)<0)) active_mesh.push(obj3d);
            });
      } else if (geo_stack && this._toplevel) {
         active_mesh = [];
         this._toplevel.traverse(mesh => {
            if ((mesh instanceof THREE.Mesh) && geo.isSameStack(mesh.stack, geo_stack)) active_mesh.push(mesh);
         });
      } else {
         active_mesh = active_mesh ? [ active_mesh ] : [];
      }

      if (!active_mesh.length) active_mesh = null;

      if (active_mesh) {
         // check if highlight is disabled for correspondent objects kinds
         if (active_mesh[0].geo_object) {
            if (!this.ctrl.highlight_scene) active_mesh = null;
         } else {
            if (!this.ctrl.highlight) active_mesh = null;
         }
      }

      if (!no_recursive) {
         // check all other painters

         if (active_mesh) {
            if (!geo_object) geo_object = active_mesh[0].geo_object;
            if (!geo_stack) geo_stack = active_mesh[0].stack;
         }

         let lst = this._highlight_handlers || (!this._main_painter ? this._slave_painters : this._main_painter._slave_painters.concat([this._main_painter]));

         for (let k=0;k<lst.length;++k)
            if (lst[k] !== this) lst[k].highlightMesh(null, color, geo_object, geo_index, geo_stack, true);
      }

      let curr_mesh = this._selected_mesh;

      let get_ctrl = mesh => {
         return mesh.get_ctrl ? mesh.get_ctrl() : new GeoDrawingControl(mesh, this.ctrl.bloom.enabled);
      }

      // check if selections are the same
      if (!curr_mesh && !active_mesh) return false;
      let same = false;
      if (curr_mesh && active_mesh && (curr_mesh.length == active_mesh.length)) {
         same = true;
         for (let k=0;(k<curr_mesh.length) && same;++k) {
            if ((curr_mesh[k] !== active_mesh[k]) || get_ctrl(curr_mesh[k]).checkHighlightIndex(geo_index)) same = false;
         }
      }
      if (same) return !!curr_mesh;

      if (curr_mesh)
         for (let k = 0; k < curr_mesh.length; ++k)
            get_ctrl(curr_mesh[k]).setHighlight();

      this._selected_mesh = active_mesh;

      if (active_mesh)
         for (let k = 0; k < active_mesh.length; ++k)
            get_ctrl(active_mesh[k]).setHighlight(color || 0x00ff00, geo_index);

      this.render3D(0);

      return !!active_mesh;
   }

   /** @summary handle mouse click event */
   TGeoPainter.prototype.processMouseClick = function(pnt, intersects, evnt) {
      if (!intersects.length) return;

      let mesh = intersects[0].object;
      if (!mesh.get_ctrl) return;

      let ctrl = mesh.get_ctrl();

      let click_indx = ctrl.extractIndex(intersects[0]);

      ctrl.evnt = evnt;

      if (ctrl.setSelected("blue", click_indx))
         this.render3D();

      ctrl.evnt = null;
   }

   /** @summary Configure mouse delay, required for complex geometries */
   TGeoPainter.prototype.setMouseTmout = function(val) {
      if (this.ctrl)
         this.ctrl.mouse_tmout = val;

      if (this._controls)
         this._controls.mouse_tmout = val;
   }

   /** @summary Configure depth method, used for render order production.
    * @param {string} method - Allowed values: "ray", "box","pnt", "size", "dflt" */
   TGeoPainter.prototype.setDepthMethod = function(method) {
      if (this.ctrl)
         this.ctrl.depthMethod = method;
   }

   /** @summary Add orbit control */
   TGeoPainter.prototype.addOrbitControls = function() {

      if (this._controls || !this._webgl || JSROOT.batch_mode) return;

      let painter = this;

      this.setTooltipAllowed(JSROOT.settings.Tooltip);

      this._controls = jsrp.createOrbitControl(this, this._camera, this._scene, this._renderer, this._lookat);

      this._controls.mouse_tmout = this.ctrl.mouse_tmout; // set larger timeout for geometry processing

      if (!this.ctrl.can_rotate) this._controls.enableRotate = false;

      this._controls.contextMenu = this.orbitContext.bind(this);

      this._controls.ProcessMouseMove = function(intersects) {

         // painter already cleaned up, ignore any incoming events
         if (!painter.ctrl || !painter._controls) return;

         let active_mesh = null, tooltip = null, resolve = null, names = [], geo_object, geo_index;

         // try to find mesh from intersections
         for (let k = 0; k < intersects.length; ++k) {
            let obj = intersects[k].object, info = null;
            if (!obj) continue;
            if (obj.geo_object) info = obj.geo_name; else
            if (obj.stack) info = painter.getStackFullName(obj.stack);
            if (!info) continue;

            if (info.indexOf("<prnt>")==0)
               info = painter.getItemName() + info.substr(6);

            names.push(info);

            if (!active_mesh) {
               active_mesh = obj;
               tooltip = info;
               geo_object = obj.geo_object;
               if (obj.get_ctrl) {
                  geo_index = obj.get_ctrl().extractIndex(intersects[k]);
                  if ((geo_index !== undefined) && (typeof tooltip == "string")) tooltip += " indx:" + JSON.stringify(geo_index);
               }
               if (active_mesh.stack) resolve = painter.resolveStack(active_mesh.stack);
            }
         }

         painter.highlightMesh(active_mesh, undefined, geo_object, geo_index);

         if (painter.ctrl.update_browser) {
            if (painter.ctrl.highlight && tooltip) names = [ tooltip ];
            painter.activateInBrowser(names);
         }

         if (!resolve || !resolve.obj) return tooltip;

         let lines = geo.provideObjectInfo(resolve.obj);
         lines.unshift(tooltip);

         return { name: resolve.obj.fName, title: resolve.obj.fTitle || resolve.obj._typename, lines: lines };
      }

      this._controls.ProcessMouseLeave = function() {
         this.ProcessMouseMove([]); // to disable highlight and reset browser
      }

      this._controls.processDblClick = function() {
         // painter already cleaned up, ignore any incoming events
         if (!painter.ctrl || !painter._controls) return;

         if (painter._last_manifest) {
            painter._last_manifest.wireframe = !painter._last_manifest.wireframe;
            if (painter._last_hidden)
               painter._last_hidden.forEach(obj => { obj.visible = true; });
            delete painter._last_hidden;
            delete painter._last_manifest;
            painter.render3D();
         } else {
            painter.adjustCameraPosition();
         }
      }
   }

   /** @summary add transformation control */
   TGeoPainter.prototype.addTransformControl = function() {
      if (this._tcontrols) return;

      if ( !this.ctrl._debug && !this.ctrl._grid ) return;

      this._tcontrols = new THREE.TransformControls( this._camera, this._renderer.domElement );
      this._scene.add( this._tcontrols );
      this._tcontrols.attach( this._toplevel );
      //this._tcontrols.setSize( 1.1 );

      window.addEventListener( 'keydown', event => {
         switch ( event.keyCode ) {
         case 81: // Q
            this._tcontrols.setSpace( this._tcontrols.space === "local" ? "world" : "local" );
            break;
         case 17: // Ctrl
            this._tcontrols.setTranslationSnap( Math.ceil( this._overall_size ) / 50 );
            this._tcontrols.setRotationSnap( THREE.Math.degToRad( 15 ) );
            break;
         case 84: // T (Translate)
            this._tcontrols.setMode( "translate" );
            break;
         case 82: // R (Rotate)
            this._tcontrols.setMode( "rotate" );
            break;
         case 83: // S (Scale)
            this._tcontrols.setMode( "scale" );
            break;
         case 187:
         case 107: // +, =, num+
            this._tcontrols.setSize( this._tcontrols.size + 0.1 );
            break;
         case 189:
         case 109: // -, _, num-
            this._tcontrols.setSize( Math.max( this._tcontrols.size - 0.1, 0.1 ) );
            break;
         }
      });
      window.addEventListener( 'keyup', function ( event ) {
         switch ( event.keyCode ) {
         case 17: // Ctrl
            this._tcontrols.setTranslationSnap( null );
            this._tcontrols.setRotationSnap( null );
            break;
         }
      });

      this._tcontrols.addEventListener( 'change', () => this.render3D(0));
   }

   /** @summary Main function in geometry creation loop
    * @desc Return false when nothing todo
    * return true if one could perform next action immediately
    * return 1 when call after short timeout required
    * return 2 when call must be done from processWorkerReply
    * @returns {number} next operation kind, see desc */
   TGeoPainter.prototype.nextDrawAction = function() {

      if (!this._clones || (this.drawing_stage == 0)) return false;

      if (this.drawing_stage == 1) {

         if (this._geom_viewer) {
            this._draw_all_nodes = false;
            this.drawing_stage = 3;
            return true;
         }

         // wait until worker is really started
         if (this.ctrl.use_worker > 0) {
            if (!this._worker) { this.startWorker(); return 1; }
            if (!this._worker_ready) return 1;
         }

         // first copy visibility flags and check how many unique visible nodes exists
         let numvis = this._clones.countVisibles() || this._clones.markVisibles(),
             matrix = null, frustum = null;

         if (this.ctrl.select_in_view && !this._first_drawing) {
            // extract camera projection matrix for selection

            matrix = geo.createProjectionMatrix(this._camera);

            frustum = geo.createFrustum(matrix);

            // check if overall bounding box seen
            if (frustum.CheckBox(this.getGeomBoundingBox(this._toplevel))) {
               matrix = null; // not use camera for the moment
               frustum = null;
            }
         }

         this._current_face_limit = this.ctrl.maxlimit;
         if (matrix) this._current_face_limit*=1.25;

         // here we decide if we need worker for the drawings
         // main reason - too large geometry and large time to scan all camera positions
         let need_worker = !JSROOT.batch_mode && ((numvis > 10000) || (matrix && (this._clones.scanVisible() > 1e5)));

         // worker does not work when starting from file system
         if (need_worker && JSROOT.source_dir.indexOf("file://")==0) {
            console.log('disable worker for jsroot from file system');
            need_worker = false;
         }

         if (need_worker && !this._worker && (this.ctrl.use_worker >= 0))
            this.startWorker(); // we starting worker, but it may not be ready so fast

         if (!need_worker || !this._worker_ready) {
            // let tm1 = new Date().getTime();
            let res = this._clones.collectVisibles(this._current_face_limit, frustum);
            this._new_draw_nodes = res.lst;
            this._draw_all_nodes = res.complete;
            // let tm2 = new Date().getTime();
            // console.log('Collect visibles', this._new_draw_nodes.length, 'takes', tm2-tm1);
            this.drawing_stage = 3;
            return true;
         }

         let job = {
            collect: this._current_face_limit,   // indicator for the command
            flags: this._clones.getVisibleFlags(),
            matrix: matrix ? matrix.elements : null
         };

         this.submitToWorker(job);

         this.drawing_stage = 2;

         this.drawing_log = "Worker select visibles";

         return 2; // we now waiting for the worker reply
      }

      if (this.drawing_stage == 2) {
         // do nothing, we are waiting for worker reply

         this.drawing_log = "Worker select visibles";

         return 2;
      }

      if (this.drawing_stage == 3) {
         // here we merge new and old list of nodes for drawing,
         // normally operation is fast and can be implemented with one call

         this.drawing_log = "Analyse visibles";

         if (this._new_append_nodes) {

            this._new_draw_nodes = this._draw_nodes.concat(this._new_append_nodes);

            delete this._new_append_nodes;

         } else if (this._draw_nodes) {

            let del;
            if (this._geom_viewer)
               del = this._draw_nodes;
            else
               del = this._clones.mergeVisibles(this._new_draw_nodes, this._draw_nodes);

            // remove should be fast, do it here
            for (let n=0;n<del.length;++n)
               this._clones.createObject3D(del[n].stack, this._toplevel, 'delete_mesh');

            if (del.length > 0)
               this.drawing_log = "Delete " + del.length + " nodes";
         }

         this._draw_nodes = this._new_draw_nodes;
         delete this._new_draw_nodes;
         this.drawing_stage = 4;
         return true;
      }

      if (this.drawing_stage === 4) {

         this.drawing_log = "Collect shapes";

         // collect shapes
         let shapes = this._clones.collectShapes(this._draw_nodes);

         // merge old and new list with produced shapes
         this._build_shapes = this._clones.mergeShapesLists(this._build_shapes, shapes);

         this.drawing_stage = 5;
         return true;
      }


      if (this.drawing_stage === 5) {
         // this is building of geometries,
         // one can ask worker to build them or do it ourself

         if (this.canSubmitToWorker()) {
            let job = { limit: this._current_face_limit, shapes: [] }, cnt = 0;
            for (let n=0;n<this._build_shapes.length;++n) {
               let clone = null, item = this._build_shapes[n];
               // only submit not-done items
               if (item.ready || item.geom) {
                  // this is place holder for existing geometry
                  clone = { id: item.id, ready: true, nfaces: geo.countGeometryFaces(item.geom), refcnt: item.refcnt };
               } else {
                  clone = JSROOT.clone(item, null, true);
                  cnt++;
               }

               job.shapes.push(clone);
            }

            if (cnt > 0) {
               /// only if some geom missing, submit job to the worker
               this.submitToWorker(job);
               this.drawing_log = "Worker build shapes";
               this.drawing_stage = 6;
               return 2;
            }
         }

         this.drawing_stage = 7;
      }

      if (this.drawing_stage === 6) {
         // waiting shapes from the worker, worker should activate our code
         return 2;
      }

      if ((this.drawing_stage === 7) || (this.drawing_stage === 8)) {

         if (this.drawing_stage === 7) {
            // building shapes
            let res = this._clones.buildShapes(this._build_shapes, this._current_face_limit, 500);
            if (res.done) {
               this.ctrl.info.num_shapes = this._build_shapes.length;
               this.drawing_stage = 8;
            } else {
               this.ctrl.info.num_shapes = res.shapes;
               this.drawing_log = "Creating: " + res.shapes + " / " + this._build_shapes.length + " shapes,  "  + res.faces + " faces";
               if (res.notusedshapes < 30) return true;
            }
         }

         // final stage, create all meshes

         let tm0 = new Date().getTime(), ready = true,
             toplevel = this.ctrl.project ? this._full_geom : this._toplevel;

         for (let n = 0; n < this._draw_nodes.length; ++n) {
            let entry = this._draw_nodes[n];
            if (entry.done) continue;

            /// shape can be provided with entry itself
            let shape = entry.server_shape || this._build_shapes[entry.shapeid];
            if (!shape.ready) {
               if (this.drawing_stage === 8) console.warn('shape marked as not ready when should');
               ready = false;
               continue;
            }

            entry.done = true;
            shape.used = true; // indicate that shape was used in building

            if (this.createEntryMesh(entry, shape, toplevel)) {
               this.ctrl.info.num_meshes++;
               this.ctrl.info.num_faces += shape.nfaces;
            }

            let tm1 = new Date().getTime();
            if (tm1 - tm0 > 500) { ready = false; break; }
         }

         if (ready) {
            if (this.ctrl.project) {
               this.drawing_log = "Build projection";
               this.drawing_stage = 10;
               return true;
            }

            this.drawing_log = "Building done";
            this.drawing_stage = 0;
            return false;
         }

         if (this.drawing_stage > 7)
            this.drawing_log = "Building meshes " + this.ctrl.info.num_meshes + " / " + this.ctrl.info.num_faces;
         return true;
      }

      if (this.drawing_stage === 9) {
         // wait for main painter to be ready

         if (!this._main_painter) {
            console.warn('MAIN PAINTER DISAPPER');
            this.drawing_stage = 0;
            return false;
         }
         if (!this._main_painter._drawing_ready) return 1;

         this.drawing_stage = 10; // just do projection
      }

      if (this.drawing_stage === 10) {
         this.doProjection();
         this.drawing_log = "Building done";
         this.drawing_stage = 0;
         return false;
      }

      console.error('never come here stage = ' + this.drawing_stage);

      return false;
   }

   /** @summary Insert appropriate mesh for given entry */
   TGeoPainter.prototype.createEntryMesh = function(entry, shape, toplevel) {

      if (!shape.geom || (shape.nfaces === 0)) {
         // node is visible, but shape does not created
         this._clones.createObject3D(entry.stack, toplevel, 'delete_mesh');
         return false;
      }

      // workaround for the TGeoOverlap, where two branches should get predefined color
      if (this._splitColors && entry.stack) {
         if (entry.stack[0]===0) entry.custom_color = "green"; else
         if (entry.stack[0]===1) entry.custom_color = "blue";
      }

      let mesh, prop = this._clones.getDrawEntryProperties(entry),
          obj3d = this._clones.createObject3D(entry.stack, toplevel, this.ctrl);

      prop.material.wireframe = this.ctrl.wireframe;

      prop.material.side = this.ctrl.bothSides ? THREE.DoubleSide : THREE.FrontSide;

      let matrix = obj3d.absMatrix || obj3d.matrixWorld;

      if (matrix.determinant() > -0.9) {
         mesh = new THREE.Mesh( shape.geom, prop.material );
      } else {
         mesh = geo.createFlippedMesh(shape, prop.material);
      }

      obj3d.add(mesh);

      if (obj3d.absMatrix) {
         mesh.matrix.copy(obj3d.absMatrix);
         mesh.matrix.decompose( mesh.position, mesh.quaternion, mesh.scale );
         mesh.updateMatrixWorld();
      }

      // keep full stack of nodes
      mesh.stack = entry.stack;
      mesh.renderOrder = this._clones.maxdepth - entry.stack.length; // order of transparency handling

      // keep hierarchy level
      mesh.$jsroot_order = obj3d.$jsroot_depth;

      // set initial render order, when camera moves, one must refine it
      //mesh.$jsroot_order = mesh.renderOrder =
      //   this._clones.maxdepth - ((obj3d.$jsroot_depth !== undefined) ? obj3d.$jsroot_depth : entry.stack.length);

      if (this.ctrl._debug || this.ctrl._full) {
         let wfg = new THREE.WireframeGeometry( mesh.geometry ),
             wfm = new THREE.LineBasicMaterial( { color: prop.fillcolor, linewidth: prop.linewidth || 1 } ),
             helper = new THREE.LineSegments(wfg, wfm);
         obj3d.add(helper);
      }

      if (this.ctrl._bound || this.ctrl._full) {
         let boxHelper = new THREE.BoxHelper( mesh );
         obj3d.add( boxHelper );
      }

      return true;
   }

   /** @summary used by geometry viewer to show more nodes
    * @desc These nodes excluded from selection logic and always inserted into the model
    * Shape already should be created and assigned to the node */
   TGeoPainter.prototype.appendMoreNodes = function(nodes, from_drawing) {
      if (this.drawing_stage && !from_drawing) {
         this._provided_more_nodes = nodes;
         return;
      }

      // delete old nodes
      if (this._more_nodes)
         for (let n=0;n<this._more_nodes.length;++n) {
            let entry = this._more_nodes[n];
            let obj3d = this._clones.createObject3D(entry.stack, this._toplevel, 'delete_mesh');
            jsrp.disposeThreejsObject(obj3d);
            geo.cleanupShape(entry.server_shape);
            delete entry.server_shape;
         }

      delete this._more_nodes;

      if (!nodes) return;

      let real_nodes = [];

      for (let k=0;k<nodes.length;++k) {
         let entry = nodes[k];
         let shape = entry.server_shape;
         if (!shape || !shape.ready) continue;

         entry.done = true;
         shape.used = true; // indicate that shape was used in building

         if (this.createEntryMesh(entry, shape, this._toplevel))
            real_nodes.push(entry);
      }

      // remember additional nodes only if they include shape - otherwise one can ignore them
      if (real_nodes) this._more_nodes = real_nodes;

      if (!from_drawing) this.render3D();
   }

   /** @summary Returns hierarchy of 3D objects used to produce projection.
    * @desc Typically external master painter is used, but also internal data can be used */
   TGeoPainter.prototype.getProjectionSource = function() {
      if (this._clones_owner)
         return this._full_geom;
      if (!this._main_painter) {
         console.warn('MAIN PAINTER DISAPPER');
         return null;
      }
      if (!this._main_painter._drawing_ready) {
         console.warn('MAIN PAINTER NOT READY WHEN DO PROJECTION');
         return null;
      }
      return this._main_painter._toplevel;
   }

   /** @summary Calculate geometry bounding box */
   TGeoPainter.prototype.getGeomBoundingBox = function(topitem, scalar) {
      let box3 = new THREE.Box3(), check_any = !this._clones;

      if (!topitem) {
         box3.min.x = box3.min.y = box3.min.z = -1;
         box3.max.x = box3.max.y = box3.max.z = 1;
         return box3;
      }

      box3.makeEmpty();

      topitem.traverse(function(mesh) {
         if (check_any || (mesh.stack && (mesh instanceof THREE.Mesh)) ||
             (mesh.main_track && (mesh instanceof THREE.LineSegments)))
            geo.getBoundingBox(mesh, box3);
      });

      if (scalar !== undefined) box3.expandByVector(box3.getSize(new THREE.Vector3()).multiplyScalar(scalar));

      return box3;
   }


   /** @summary Create geometry projection */
   TGeoPainter.prototype.doProjection = function() {
      let toplevel = this.getProjectionSource();

      if (!toplevel) return false;

      jsrp.disposeThreejsObject(this._toplevel, true);

      // let axis = this.ctrl.project;

      if (this.ctrl.projectPos === undefined) {

         let bound = this.getGeomBoundingBox(toplevel),
             min = bound.min[this.ctrl.project], max = bound.max[this.ctrl.project],
             mean = (min+max)/2;

         if ((min<0) && (max>0) && (Math.abs(mean) < 0.2*Math.max(-min,max))) mean = 0; // if middle is around 0, use 0

         this.ctrl.projectPos = mean;
      }

      toplevel.traverse(mesh => {
         if (!(mesh instanceof THREE.Mesh) || !mesh.stack) return;

         let geom2 = geo.projectGeometry(mesh.geometry, mesh.parent.absMatrix || mesh.parent.matrixWorld, this.ctrl.project, this.ctrl.projectPos, mesh._flippedMesh);

         if (!geom2) return;

         let mesh2 = new THREE.Mesh( geom2, mesh.material.clone() );

         this._toplevel.add(mesh2);

         mesh2.stack = mesh.stack;
      });

      return true;
   }

   /** @summary Should be invoked when light configuration changed */
   TGeoPainter.prototype.changedLight = function(box) {
      if (!this._camera) return;

      let need_render = !box;

      if (!box) box = this.getGeomBoundingBox(this._toplevel);

      let sizex = box.max.x - box.min.x,
          sizey = box.max.y - box.min.y,
          sizez = box.max.z - box.min.z,
          plights = [], p = this.ctrl.light.power;

      if (p === undefined) p = 1;

      if (this._camera._lights != this.ctrl.light.kind) {
         // remove all childs and recreate only necessary lights
         jsrp.disposeThreejsObject(this._camera, true);

         this._camera._lights = this.ctrl.light.kind;

         switch (this._camera._lights) {
            case "ambient" : this._camera.add(new THREE.AmbientLight(0xefefef, p)); break;
            case "hemisphere" : this._camera.add(new THREE.HemisphereLight(0xffffbb, 0x080820, p)); break;
            default: // 6 point lights
               for (let n=0;n<6;++n)
                  this._camera.add( new THREE.PointLight(0xefefef, p) );
         }
      }

      for (let k=0;k<this._camera.children.length;++k) {
         let light = this._camera.children[k], enabled = false;
         if (light.isAmbientLight || light.isHemisphereLight) {
            light.intensity = p;
            continue;
         }

         if (!light.isPointLight) continue;
         switch (k) {
            case 0: light.position.set(sizex/5, sizey/5, sizez/5); enabled = this.ctrl.light.specular; break;
            case 1: light.position.set(0, 0, sizez/2); enabled = this.ctrl.light.front; break;
            case 2: light.position.set(0, 2*sizey, 0); enabled = this.ctrl.light.top; break;
            case 3: light.position.set(0, -2*sizey, 0); enabled = this.ctrl.light.bottom; break;
            case 4: light.position.set(-2*sizex, 0, 0); enabled = this.ctrl.light.left; break;
            case 5: light.position.set(2*sizex, 0, 0); enabled = this.ctrl.light.right; break;
         }
         light.power = enabled ? p*Math.PI*4 : 0;
         if (enabled) plights.push(light);
      }

      // keep light power of all soources constant
      plights.forEach(ll => { ll.power = p*4*Math.PI/plights.length; })

      if (need_render) this.render3D();
   }

   /** @summary Initial scene creation */
   TGeoPainter.prototype.createScene = function(w, h) {
      // three.js 3D drawing
      this._scene = new THREE.Scene();
      this._scene.fog = new THREE.Fog(0xffffff, 1, 10000);
      this._scene.overrideMaterial = new THREE.MeshLambertMaterial( { color: 0x7000ff, transparent: true, opacity: 0.2, depthTest: false } );

      this._scene_width = w;
      this._scene_height = h;

      if (this.ctrl.ortho_camera) {
         this._camera =  new THREE.OrthographicCamera(-600, 600, -600, 600, 1, 10000);
      } else {
         this._camera = new THREE.PerspectiveCamera(25, w / h, 1, 10000);
         this._camera.up = this.ctrl._yup ? new THREE.Vector3(0,1,0) : new THREE.Vector3(0,0,1);
      }

      this._scene.add( this._camera );

      this._selected_mesh = null;

      this._overall_size = 10;

      this._toplevel = new THREE.Object3D();

      this._scene.add(this._toplevel);

      this._renderer = jsrp.createRender3D(w, h, this.options.Render3D, { antialias: true, logarithmicDepthBuffer: false, preserveDrawingBuffer: true });

      this._webgl = (this._renderer.jsroot_render3d === JSROOT.constants.Render3D.WebGL);

      if (this._renderer.setPixelRatio && !JSROOT.nodejs)
         this._renderer.setPixelRatio(window.devicePixelRatio);
      this._renderer.setSize(w, h, !this._fit_main_area);
      this._renderer.localClippingEnabled = true;

      this._renderer.setClearColor(this.ctrl.background, 1);

      if (this._fit_main_area && this._webgl) {
         this._renderer.domElement.style.width = "100%";
         this._renderer.domElement.style.height = "100%";
         let main = this.selectDom();
         if (main.style('position')=='static') main.style('position','relative');
      }

      this._animating = false;

      // Clipping Planes

      this.ctrl.bothSides = false; // which material kind should be used
      this._clipPlanes = [ new THREE.Plane(new THREE.Vector3(1, 0, 0), 0),
                           new THREE.Plane(new THREE.Vector3(0, this.ctrl._yup ? -1 : 1, 0), 0),
                           new THREE.Plane(new THREE.Vector3(0, 0, this.ctrl._yup ? 1 : -1), 0) ];

      // Light - add default point light, adjust later

      let light = new THREE.PointLight(0xefefef, 1);
      light.position.set(10, 10, 10);
      this._camera.add( light );

      // Smooth Lighting Shader (Screen Space Ambient Occlusion)
      // http://threejs.org/examples/webgl_postprocessing_ssao.html

      if (this._webgl && (this.ctrl.ssao.enabled || this.ctrl.outline)) {

         if (this.ctrl.outline && (typeof this.createOutline == "function")) {
            this._effectComposer = new THREE.EffectComposer( this._renderer );
            this._effectComposer.addPass( new THREE.RenderPass( this._scene, this._camera ) );
            this.createOutline(w, h);
         } else if (this.ctrl.ssao.enabled) {
            this.createSSAO();
         }
      }

      if (this._webgl && this.ctrl.bloom.enabled) {
         this.createBloom();
      }

      if (this._fit_main_area && !this._webgl) {
         // create top-most SVG for geomtery drawings
         let doc = JSROOT._.get_document();
         let svg = doc.createElementNS("http://www.w3.org/2000/svg", "svg");
         d3.select(svg).attr("width",w).attr("height",h);
         svg.appendChild(this._renderer.jsroot_dom);
         return svg;
      }

      return this._renderer.jsroot_dom;
   }

   /** @summary Start geometry drawing */
   TGeoPainter.prototype.startDrawGeometry = function(force) {

      if (!force && (this.drawing_stage!==0)) {
         this._draw_nodes_again = true;
         return;
      }

      if (this._clones_owner && this._clones)
         this._clones.setDefaultColors(this.ctrl.dflt_colors);

      this._startm = new Date().getTime();
      this._last_render_tm = this._startm;
      this._last_render_meshes = 0;
      this.drawing_stage = 1;
      this._drawing_ready = false;
      this.drawing_log = "collect visible";
      this.ctrl.info.num_meshes = 0;
      this.ctrl.info.num_faces = 0;
      this.ctrl.info.num_shapes = 0;
      this._selected_mesh = null;

      if (this.ctrl.project) {
         if (this._clones_owner) {
            if (this._full_geom) {
               this.drawing_stage = 10;
               this.drawing_log = "build projection";
            } else {
               this._full_geom = new THREE.Object3D();
            }

         } else {
            this.drawing_stage = 9;
            this.drawing_log = "wait for main painter";
         }
      }

      delete this._last_manifest;
      delete this._last_hidden; // clear list of hidden objects

      delete this._draw_nodes_again; // forget about such flag

      this.continueDraw();
   }

   /** @summary reset all kind of advanced features like SSAO or depth test changes */
   TGeoPainter.prototype.resetAdvanced = function() {
      this.ctrl.ssao.kernelRadius = 16;
      this.ctrl.ssao.output = THREE.SSAOPass.OUTPUT.Default;

      this.ctrl.depthTest = true;
      this.ctrl.clipIntersect = true;
      this.ctrl.depthMethod = "ray";

      this.changedDepthMethod("norender");
      this.changedDepthTest();
   }

   /** @summary returns maximal dimension */
   TGeoPainter.prototype.getOverallSize = function(force) {
      if (!this._overall_size || force) {
         let box = this.getGeomBoundingBox(this._toplevel);

         // if detect of coordinates fails - ignore
         if (isNaN(box.min.x)) return 1000;

         let sizex = box.max.x - box.min.x,
             sizey = box.max.y - box.min.y,
             sizez = box.max.z - box.min.z;

         this._overall_size = 2 * Math.max(sizex, sizey, sizez);
      }

      return this._overall_size;
   }

   /** @summary Create png image with drawing snapshot. */
   TGeoPainter.prototype.createSnapshot = function(filename) {
      if (!this._renderer) return;
      this.render3D(0);
      let dataUrl = this._renderer.domElement.toDataURL("image/png");
      if (filename==="asis") return dataUrl;
      dataUrl.replace("image/png", "image/octet-stream");
      let doc = JSROOT._.get_document();
      let link = doc.createElement('a');
      if (typeof link.download === 'string') {
         doc.body.appendChild(link); //Firefox requires the link to be in the body
         link.download = filename || "geometry.png";
         link.href = dataUrl;
         link.click();
         doc.body.removeChild(link); //remove the link when done
      }
   }

   /** @summary Returns url parameters defining camera position.
    * @desc It is zoom, roty, rotz parameters
    * These parameters applied from default position which is shift along X axis */
   TGeoPainter.prototype.produceCameraUrl = function(prec) {

      if (!this._lookat || !this._camera0pos || !this._camera || !this.ctrl) return;

      let pos1 = new THREE.Vector3().add(this._camera0pos).sub(this._lookat),
          pos2 = new THREE.Vector3().add(this._camera.position).sub(this._lookat),
          len1 = pos1.length(), len2 = pos2.length(),
          zoom = this.ctrl.zoom * len2 / len1 * 100;

      if (zoom < 1) zoom = 1; else if (zoom>10000) zoom = 10000;

      pos1.normalize();
      pos2.normalize();

      let quat = new THREE.Quaternion();
      quat.setFromUnitVectors(pos1, pos2);

      let euler = new THREE.Euler();
      euler.setFromQuaternion(quat, "YZX");

      let roty = euler.y / Math.PI * 180,
          rotz = euler.z / Math.PI * 180;

      if (roty<0) roty += 360;
      if (rotz<0) rotz += 360;

      return "roty" + roty.toFixed(prec || 0) + ",rotz" + rotz.toFixed(prec || 0) + ",zoom" + zoom.toFixed(prec || 0);
   }

   /** @summary Calculates current zoom factor */
   TGeoPainter.prototype.calculateZoom = function() {
      if (this._camera0pos && this._camera && this._lookat) {
         let pos1 = new THREE.Vector3().add(this._camera0pos).sub(this._lookat),
             pos2 = new THREE.Vector3().add(this._camera.position).sub(this._lookat);
         return pos2.length() / pos1.length();
      }

      return 0;
   }


   /** @summary Place camera to default position */
   TGeoPainter.prototype.adjustCameraPosition = function(first_time, keep_zoom) {
      if (!this._toplevel) return;

      let box = this.getGeomBoundingBox(this._toplevel);

      // if detect of coordinates fails - ignore
      if (isNaN(box.min.x)) return;

      let sizex = box.max.x - box.min.x,
          sizey = box.max.y - box.min.y,
          sizez = box.max.z - box.min.z,
          midx = (box.max.x + box.min.x)/2,
          midy = (box.max.y + box.min.y)/2,
          midz = (box.max.z + box.min.z)/2;

      this._overall_size = 2 * Math.max(sizex, sizey, sizez);

      this._camera.near = this._overall_size / 350;
      this._camera.far = this._overall_size * 12;
      this._scene.fog.near = this._overall_size * 2;
      this._scene.fog.far = this._overall_size * 12;

      if (first_time)
         for (let naxis=0;naxis<3;++naxis) {
            let cc = this.ctrl.clip[naxis];
            cc.min = box.min[cc.name];
            cc.max = box.max[cc.name];
            let sz = cc.max - cc.min;
            cc.max += sz*0.01;
            cc.min -= sz*0.01;
            if (!cc.value)
               cc.value = (cc.min + cc.max) / 2;
            else if (cc.value < cc.min)
               cc.value = cc.min;
            else if (cc.value > cc.max)
               cc.value = cc.max;
         }

      if (this.ctrl.ortho_camera) {
         this._camera.left = box.min.x;
         this._camera.right = box.max.x;
         this._camera.top = box.max.y;
         this._camera.bottom = box.min.y;
      }

      // this._camera.far = 100000000000;

      this._camera.updateProjectionMatrix();

      let k = 2*this.ctrl.zoom,
          max_all = Math.max(sizex,sizey,sizez);

      if ((this.ctrl.rotatey || this.ctrl.rotatez) && this.ctrl.can_rotate) {

         let prev_zoom = this.calculateZoom();
         if (keep_zoom && prev_zoom) k = 2*prev_zoom;

         let euler = new THREE.Euler( 0, this.ctrl.rotatey/180.*Math.PI, this.ctrl.rotatez/180.*Math.PI, 'YZX' );

         this._camera.position.set(-k*max_all, 0, 0);
         this._camera.position.applyEuler(euler);
         this._camera.position.add(new THREE.Vector3(midx,midy,midz));

         if (keep_zoom && prev_zoom) {
            let actual_zoom = this.calculateZoom();
            k *= prev_zoom/actual_zoom;

            this._camera.position.set(-k*max_all, 0, 0);
            this._camera.position.applyEuler(euler);
            this._camera.position.add(new THREE.Vector3(midx,midy,midz));
         }

      } else if (this.ctrl.ortho_camera) {
         this._camera.position.set(midx, midy, Math.max(sizex,sizey));
      } else if (this.ctrl.project) {
         switch (this.ctrl.project) {
            case 'x': this._camera.position.set(k*1.5*Math.max(sizey,sizez), 0, 0); break;
            case 'y': this._camera.position.set(0, k*1.5*Math.max(sizex,sizez), 0); break;
            case 'z': this._camera.position.set(0, 0, k*1.5*Math.max(sizex,sizey)); break;
         }
      } else if (this.ctrl._yup) {
         this._camera.position.set(midx-k*Math.max(sizex,sizez), midy+k*sizey, midz-k*Math.max(sizex,sizez));
      } else {
         this._camera.position.set(midx-k*Math.max(sizex,sizey), midy-k*Math.max(sizex,sizey), midz+k*sizez);
      }

      this._lookat = new THREE.Vector3(midx, midy, midz);
      this._camera0pos = new THREE.Vector3(-2*max_all, 0, 0); // virtual 0 position, where rotation starts
      this._camera.lookAt(this._lookat);

      this.changedLight(box);

      if (this._controls) {
         this._controls.target.copy(this._lookat);
         this._controls.update();
      }

      // recheck which elements to draw
      if (this.ctrl.select_in_view)
         this.startDrawGeometry();
   }

   /** @summary Specifies camera position */
   TGeoPainter.prototype.setCameraPosition = function(rotatey, rotatez, zoom) {
      if (!this.ctrl) return;
      this.ctrl.rotatey = rotatey || 0;
      this.ctrl.rotatez = rotatez || 0;
      let preserve_zoom = false;
      if (zoom && !isNaN(zoom)) {
         this.ctrl.zoom = zoom;
      } else {
         preserve_zoom = true;
      }
      this.adjustCameraPosition(false, preserve_zoom);
   }

   TGeoPainter.prototype.focusOnItem = function(itemname) {

      if (!itemname || !this._clones) return;

      let stack = this._clones.findStackByName(itemname);

      if (!stack) return;

      let info = this._clones.resolveStack(stack, true);

      this.focusCamera( info, false );
   }

   /** @summary focus camera on speicifed position */
   TGeoPainter.prototype.focusCamera = function( focus, autoClip ) {

      if (this.ctrl.project || this.ctrl.ortho_camera)
         return this.adjustCameraPosition();

      let box = new THREE.Box3();
      if (focus === undefined) {
         box = this.getGeomBoundingBox(this._toplevel);
      } else if (focus instanceof THREE.Mesh) {
         box.setFromObject(focus);
      } else {
         let center = new THREE.Vector3().setFromMatrixPosition(focus.matrix),
             node = focus.node,
             halfDelta = new THREE.Vector3( node.fDX, node.fDY, node.fDZ ).multiplyScalar(0.5);
         box.min = center.clone().sub(halfDelta);
         box.max = center.clone().add(halfDelta);
      }

      let sizex = box.max.x - box.min.x,
          sizey = box.max.y - box.min.y,
          sizez = box.max.z - box.min.z,
          midx = (box.max.x + box.min.x)/2,
          midy = (box.max.y + box.min.y)/2,
          midz = (box.max.z + box.min.z)/2;

      let position;
      if (this.ctrl._yup)
         position = new THREE.Vector3(midx-2*Math.max(sizex,sizez), midy+2*sizey, midz-2*Math.max(sizex,sizez));
      else
         position = new THREE.Vector3(midx-2*Math.max(sizex,sizey), midy-2*Math.max(sizex,sizey), midz+2*sizez);

      let target = new THREE.Vector3(midx, midy, midz);

      // Find to points to animate "lookAt" between
      // let dist = this._camera.position.distanceTo(target);
      let oldTarget = this._controls.target;

      // probably, reduce number of frames
      let frames = 50, step = 0,
         // Amount to change camera position at each step
         posIncrement = position.sub(this._camera.position).divideScalar(frames),
         // Amount to change "lookAt" so it will end pointed at target
         targetIncrement = target.sub(oldTarget).divideScalar(frames);

      autoClip = autoClip && this._webgl;

      // Automatic Clipping
      if (autoClip) {
         for (let axis = 0; axis < 3; ++axis) {
            let cc = this.ctrl.clip[axis];
            if (!cc.enabled) { cc.value = cc.min; cc.enabled = true; }
            cc.inc = ((cc.min + cc.max) / 2 - cc.value) / frames;
         }
         this.updateClipping();
      }

      this._animating = true;

      // Interpolate //

      let animate = () => {
         if (this._animating === undefined) return;

         if (this._animating) {
            requestAnimationFrame( animate );
         } else {
            if (!this._geom_viewer)
               this.startDrawGeometry();
         }
         let smoothFactor = -Math.cos( ( 2.0 * Math.PI * step ) / frames ) + 1.0;
         this._camera.position.add( posIncrement.clone().multiplyScalar( smoothFactor ) );
         oldTarget.add( targetIncrement.clone().multiplyScalar( smoothFactor ) );
         this._lookat = oldTarget;
         this._camera.lookAt( this._lookat );
         this._camera.updateProjectionMatrix();

         let tm1 = new Date().getTime();
         if (autoClip) {
            for (let axis = 0; axis < 3; ++axis)
               this.ctrl.clip[axis].value += this.ctrl.clip[axis].inc * smoothFactor;
            this.updateClipping();
         } else {
            this.render3D(0);
         }
         let tm2 = new Date().getTime();
         if ((step==0) && (tm2-tm1 > 200)) frames = 20;
         step++;
         this._animating = step < frames;
      }

      animate();

   //   this._controls.update();
   }

   /** @summary actiavte auto rotate */
   TGeoPainter.prototype.autorotate = function(speed) {

      let rotSpeed = (speed === undefined) ? 2.0 : speed,
          last = new Date();

      let animate = () => {
         if (!this._renderer || !this.ctrl) return;

         let current = new Date();

         if ( this.ctrl.rotate ) requestAnimationFrame( animate );

         if (this._controls) {
            this._controls.autoRotate = this.ctrl.rotate;
            this._controls.autoRotateSpeed = rotSpeed * ( current.getTime() - last.getTime() ) / 16.6666;
            this._controls.update();
         }
         last = new Date();
         this.render3D(0);
      }

      if (this._webgl) animate();
   }

   /** @summary called at the end of scene drawing */
   TGeoPainter.prototype.completeScene = function() {

      if ( this.ctrl._debug || this.ctrl._grid ) {
         if ( this.ctrl._full ) {
            let boxHelper = new THREE.BoxHelper(this._toplevel);
            this._scene.add( boxHelper );
         }
         this._scene.add( new THREE.AxesHelper( 2 * this._overall_size ) );
         this._scene.add( new THREE.GridHelper( Math.ceil( this._overall_size), Math.ceil( this._overall_size ) / 50 ) );
         this.helpText("<font face='verdana' size='1' color='red'><center>Transform Controls<br>" +
               "'T' translate | 'R' rotate | 'S' scale<br>" +
               "'+' increase size | '-' decrease size<br>" +
               "'W' toggle wireframe/solid display<br>"+
         "keep 'Ctrl' down to snap to grid</center></font>");
      }
   }

   /** @summary Drawing with "count" option
     * @desc Scans hieararchy and check for unique nodes
     * @returns {Promise} with object drawing ready */
   TGeoPainter.prototype.drawCount = function(unqievis, clonetm) {

      function makeTime(tm) {
         return JSROOT.batch_mode ? "anytime ms" : tm.toString() + " ms";
      }

      let res = [ 'Unique nodes: ' + this._clones.nodes.length,
                  'Unique visible: ' + unqievis,
                  'Time to clone: ' + makeTime(clonetm) ];

      // need to fill cached value line numvischld
      this._clones.scanVisible();

      let nshapes = 0;

      let arg = {
         clones: this._clones,
         cnt: [],
         func: function(node) {
            if (this.cnt[this.last]===undefined)
               this.cnt[this.last] = 1;
            else
               this.cnt[this.last]++;

            nshapes += geo.countNumShapes(this.clones.getNodeShape(node.id));

            // for debugging - search if there some TGeoHalfSpace
            //if (geo.HalfSpace) {
            //    let entry = this.CopyStack();
            //    let res = painter._clones.resolveStack(entry.stack);
            //    console.log('SAW HALF SPACE', res.name);
            //    geo.HalfSpace = false;
            //}
            return true;
         }
      };

      let tm1 = new Date().getTime();
      let numvis = this._clones.scanVisible(arg);
      let tm2 = new Date().getTime();

      res.push('Total visible nodes: ' + numvis, 'Total shapes: ' + nshapes);

      for (let lvl=0;lvl<arg.cnt.length;++lvl) {
         if (arg.cnt[lvl] !== undefined)
            res.push('  lvl' + lvl + ': ' + arg.cnt[lvl]);
      }

      res.push("Time to scan: " + makeTime(tm2-tm1), "", "Check timing for matrix calculations ...");

      let elem = this.selectDom().style('overflow', 'auto');

      if (JSROOT.batch_mode)
         elem.property("_json_object_", res);
      else
         res.forEach(str => elem.append("p").text(str));

      return new Promise(resolveFunc => {
         setTimeout(() => {
            arg.domatrix = true;
            tm1 = new Date().getTime();
            numvis = this._clones.scanVisible(arg);
            tm2 = new Date().getTime();

            let last_str = "Time to scan with matrix: " + makeTime(tm2-tm1);
            if (JSROOT.batch_mode)
               res.push(last_str);
            else
               elem.append("p").text(last_str);
            resolveFunc(this);
         }, 100);
      });
   }

   /** @summary Handle drop operation
     * @desc opt parameter can include function name like opt$func_name
     * Such function should be possible to find via {@link JSROOT.findFunction}
     * Function has to return Promise with objects to draw on geometry
     * By default function with name "extract_geo_tracks" is checked
     * @returns {Promise} handling of drop operation */
   TGeoPainter.prototype.performDrop = function(obj, itemname, hitem, opt) {

      if (obj && (obj.$kind==='TTree')) {
         // drop tree means function call which must extract tracks from provided tree

         let funcname = "extract_geo_tracks";

         if (opt && opt.indexOf("$") > 0) {
            funcname = opt.substr(0, opt.indexOf("$"));
            opt = opt.substr(opt.indexOf("$")+1);
         }

         let func = JSROOT.findFunction(funcname);

         if (!func) return Promise.reject(Error(`Function ${funcname} not found`));

         return func(obj, opt).then(tracks => {
            if (tracks) {
               this.drawExtras(tracks, "", false); // FIXME: probably tracks should be remembered??
               this.updateClipping(true);
               this.render3D(100);
            }
            return this;
         });
      }

      if (this.drawExtras(obj, itemname)) {
         if (hitem) hitem._painter = this; // set for the browser item back pointer
         this.render3D(100);
      }

      return Promise.resolve(this);
   }

   /** @summary function called when mouse is going over the item in the browser */
   TGeoPainter.prototype.mouseOverHierarchy = function(on, itemname, hitem) {
      if (!this.ctrl) return; // protection for cleaned-up painter

      let obj = hitem._obj;
      if (this.ctrl._debug)
         console.log('Mouse over', on, itemname, (obj ? obj._typename : "---"));

      // let's highlight tracks and hits only for the time being
      if (!obj || (obj._typename !== "TEveTrack" && obj._typename !== "TEvePointSet" && obj._typename !== "TPolyMarker3D")) return;

      this.highlightMesh(null, 0x00ff00, on ? obj : null);
   }

   /** @summary clear extra drawn objects like tracks or hits */
   TGeoPainter.prototype.clearExtras = function() {
      this.getExtrasContainer("delete");
      delete this._extraObjects; // workaround, later will be normal function
      this.render3D();
   }

   /** @summary Register extra objects like tracks or hits
    * @desc Rendered after main geometry volumes are created
    * Check if object already exists to prevent duplication */
   TGeoPainter.prototype.addExtra = function(obj, itemname) {
      if (this._extraObjects === undefined)
         this._extraObjects = JSROOT.create("TList");

      if (this._extraObjects.arr.indexOf(obj)>=0) return false;

      this._extraObjects.Add(obj, itemname);

      delete obj.$hidden_via_menu; // remove previous hidden property

      return true;
   }

   /** @summary manipulate visisbility of extra objects, used for HierarhcyPainter
     * @private */
   TGeoPainter.prototype.extraObjectVisible = function(hpainter, hitem, toggle) {
      if (!this._extraObjects) return;

      let itemname = hpainter.itemFullName(hitem),
          indx = this._extraObjects.opt.indexOf(itemname);

      if ((indx<0) && hitem._obj) {
         indx = this._extraObjects.arr.indexOf(hitem._obj);
         // workaround - if object found, replace its name
         if (indx>=0) this._extraObjects.opt[indx] = itemname;
      }

      if (indx < 0) return;

      let obj = this._extraObjects.arr[indx],
          res = obj.$hidden_via_menu ? false : true;

      if (toggle) {
         obj.$hidden_via_menu = res; res = !res;

         let mesh = null;
         // either found painted object or just draw once again
         this._toplevel.traverse(function(node) { if (node.geo_object === obj) mesh = node; });

         if (mesh) mesh.visible = res; else
         if (res) {
            this.drawExtras(obj, "", false);
            this.updateClipping(true);
         }

         if (mesh || res) this.render3D();
      }

      return res;
   }

   /** @summary Draw extra object like tracks
     * @returns {Promise} when ready */
   TGeoPainter.prototype.drawExtras = function(obj, itemname, add_objects) {
      if (!obj || !obj._typename)
         return false;

      // if object was hidden via menu, do not redraw it with next draw call
      if (!add_objects && obj.$hidden_via_menu)
         return false;

      let is_any = false, do_render = false;
      if (add_objects === undefined) {
         add_objects = true;
         do_render = true;
      }

      if ((obj._typename === "TList") || (obj._typename === "TObjArray")) {
         if (!obj.arr) return false;
         for (let n=0;n<obj.arr.length;++n) {
            let sobj = obj.arr[n], sname = obj.opt ? obj.opt[n] : "";
            if (!sname) sname = (itemname || "<prnt>") + "/[" + n + "]";
            if (this.drawExtras(sobj, sname, add_objects)) is_any = true;
         }
      } else if (obj._typename === 'THREE.Mesh') {
         // adding mesh as is
         this.addToExtrasContainer(obj);
         is_any = true;
      } else if (obj._typename === 'TGeoTrack') {
         if (add_objects && !this.addExtra(obj, itemname)) return false;
         is_any = this.drawGeoTrack(obj, itemname);
      } else if ((obj._typename === 'TEveTrack') || (obj._typename === 'ROOT::Experimental::TEveTrack')) {
         if (add_objects && !this.addExtra(obj, itemname)) return false;
         is_any = this.drawEveTrack(obj, itemname);
      } else if ((obj._typename === 'TEvePointSet') || (obj._typename === "ROOT::Experimental::TEvePointSet") || (obj._typename === "TPolyMarker3D")) {
         if (add_objects && !this.addExtra(obj, itemname)) return false;
         is_any = this.drawHit(obj, itemname);
      } else if ((obj._typename === "TEveGeoShapeExtract") || (obj._typename === "ROOT::Experimental::REveGeoShapeExtract")) {
         if (add_objects && !this.addExtra(obj, itemname)) return false;
         is_any = this.drawExtraShape(obj, itemname);
      }

      if (do_render && is_any) {
         this.updateClipping(true);
         this.render3D(100);
      }
      return is_any;
   }

   /** @summary returns container for extra objects */
   TGeoPainter.prototype.getExtrasContainer = function(action, name) {
      if (!this._toplevel) return null;

      if (!name) name = "tracks";

      let extras = null, lst = [];
      for (let n=0;n<this._toplevel.children.length;++n) {
         let chld = this._toplevel.children[n];
         if (!chld._extras) continue;
         if (action==='collect') { lst.push(chld); continue; }
         if (chld._extras === name) { extras = chld; break; }
      }

      if (action==='collect') {
         for (let k=0;k<lst.length;++k) this._toplevel.remove(lst[k]);
         return lst;
      }

      if (action==="delete") {
         if (extras) this._toplevel.remove(extras);
         jsrp.disposeThreejsObject(extras);
         return null;
      }

      if ((action!=="get") && !extras) {
         extras = new THREE.Object3D();
         extras._extras = name;
         this._toplevel.add(extras);
      }

      return extras;
   }

   /** @summary add object to extras container.
     * @desc If fail, dispore object */
   TGeoPainter.prototype.addToExtrasContainer = function(obj, name) {
      let container = this.getExtrasContainer("", name);
      if (container) {
         container.add(obj);
      } else {
         console.warn('Fail to add object to extras');
         jsrp.disposeThreejsObject(obj);
      }
   }

   /** @summary drawing TGeoTrack */
   TGeoPainter.prototype.drawGeoTrack = function(track, itemname) {
      if (!track || !track.fNpoints) return false;

      let track_width = track.fLineWidth || 1,
          track_color = jsrp.getColor(track.fLineColor) || "rgb(255,0,255)";

      if (JSROOT.browser.isWin) track_width = 1; // not supported on windows

      let npoints = Math.round(track.fNpoints/4),
          buf = new Float32Array((npoints-1)*6),
          pos = 0, projv = this.ctrl.projectPos,
          projx = (this.ctrl.project === "x"),
          projy = (this.ctrl.project === "y"),
          projz = (this.ctrl.project === "z");

      for (let k=0;k<npoints-1;++k) {
         buf[pos]   = projx ? projv : track.fPoints[k*4];
         buf[pos+1] = projy ? projv : track.fPoints[k*4+1];
         buf[pos+2] = projz ? projv : track.fPoints[k*4+2];
         buf[pos+3] = projx ? projv : track.fPoints[k*4+4];
         buf[pos+4] = projy ? projv : track.fPoints[k*4+5];
         buf[pos+5] = projz ? projv : track.fPoints[k*4+6];
         pos+=6;
      }

      let lineMaterial = new THREE.LineBasicMaterial({ color: track_color, linewidth: track_width }),
          line = jsrp.createLineSegments(buf, lineMaterial);

      line.renderOrder = 1000000; // to bring line to the front
      line.geo_name = itemname;
      line.geo_object = track;
      line.hightlightWidthScale = 2;

      if (itemname && itemname.indexOf("<prnt>/Tracks")==0)
         line.main_track = true;

      this.addToExtrasContainer(line);

      return true;
   }

   /** @summary Drawing TEveTrack */
   TGeoPainter.prototype.drawEveTrack = function(track, itemname) {
      if (!track || (track.fN <= 0)) return false;

      let track_width = track.fLineWidth || 1,
          track_color = jsrp.getColor(track.fLineColor) || "rgb(255,0,255)";

      if (JSROOT.browser.isWin) track_width = 1; // not supported on windows

      let buf = new Float32Array((track.fN-1)*6), pos = 0,
          projv = this.ctrl.projectPos,
          projx = (this.ctrl.project === "x"),
          projy = (this.ctrl.project === "y"),
          projz = (this.ctrl.project === "z");

      for (let k=0;k<track.fN-1;++k) {
         buf[pos]   = projx ? projv : track.fP[k*3];
         buf[pos+1] = projy ? projv : track.fP[k*3+1];
         buf[pos+2] = projz ? projv : track.fP[k*3+2];
         buf[pos+3] = projx ? projv : track.fP[k*3+3];
         buf[pos+4] = projy ? projv : track.fP[k*3+4];
         buf[pos+5] = projz ? projv : track.fP[k*3+5];
         pos+=6;
      }

      let lineMaterial = new THREE.LineBasicMaterial({ color: track_color, linewidth: track_width }),
          line = jsrp.createLineSegments(buf, lineMaterial);

      line.renderOrder = 1000000; // to bring line to the front
      line.geo_name = itemname;
      line.geo_object = track;
      line.hightlightWidthScale = 2;

      this.addToExtrasContainer(line);

      return true;
   }

   /** @summary Drawing different hits types like TPolyMarker3d */
   TGeoPainter.prototype.drawHit = function(hit, itemname) {
      if (!hit || !hit.fN || (hit.fN < 0)) return false;

      // make hit size scaling factor of overall geometry size
      // otherwise it is not possible to correctly see hits at all
      let hit_size = hit.fMarkerSize * this.getOverallSize() * 0.005;
      if (hit_size <= 0.2) hit_size = 0.2;

      let hit_style = hit.fMarkerStyle;
      // FIXME: marker style 2 does not work why?
      if ((hit_style == 4) || (hit_style == 2)) { hit_style = 7; hit_size *= 1.5; } // style 4 is very bad for hits representation

      let size = hit.fN,
          projv = this.ctrl.projectPos,
          projx = (this.ctrl.project === "x"),
          projy = (this.ctrl.project === "y"),
          projz = (this.ctrl.project === "z"),
          pnts = new jsrp.PointsCreator(size, this._webgl, hit_size);

      for (let i=0;i<size;i++)
         pnts.addPoint(projx ? projv : hit.fP[i*3],
                       projy ? projv : hit.fP[i*3+1],
                       projz ? projv : hit.fP[i*3+2]);

      let mesh = pnts.createPoints({ color: jsrp.getColor(hit.fMarkerColor) || "rgb(0,0,255)", style: hit_style });
      mesh.renderOrder = 1000000; // to bring points to the front
      mesh.highlightScale = 2;
      mesh.geo_name = itemname;
      mesh.geo_object = hit;
      this.addToExtrasContainer(mesh);

      return true; // indicate that rendering should be done
   }

   TGeoPainter.prototype.drawExtraShape = function(obj, itemname) {
      let toplevel = geo.build(obj);
      if (!toplevel) return false;

      toplevel.geo_name = itemname;
      toplevel.geo_object = obj;

      this.addToExtrasContainer(toplevel);
      return true;
   }

   /** @summary Serach for specified node
     * @private */
   TGeoPainter.prototype.findNodeWithVolume = function(name, action, prnt, itemname, volumes) {

      let first_level = false, res = null;

      if (!prnt) {
         prnt = this.getGeometry();
         if (!prnt && (geo.getNodeKind(prnt)!==0)) return null;
         itemname = this.geo_manager ? prnt.fName : "";
         first_level = true;
         volumes = [];
      } else {
         if (itemname.length>0) itemname += "/";
         itemname += prnt.fName;
      }

      if (!prnt.fVolume || prnt.fVolume._searched) return null;

      if (name.test(prnt.fVolume.fName)) {
         res = action({ node: prnt, item: itemname });
         if (res) return res;
      }

      prnt.fVolume._searched = true;
      volumes.push(prnt.fVolume);

      if (prnt.fVolume.fNodes)
         for (let n=0;n<prnt.fVolume.fNodes.arr.length;++n) {
            res = this.findNodeWithVolume(name, action, prnt.fVolume.fNodes.arr[n], itemname, volumes);
            if (res) break;
         }

      if (first_level)
         for (let n=0, len=volumes.length; n<len; ++n)
            delete volumes[n]._searched;

      return res;
   }

   /** @summary Process script option - load and execute some gGeoManager-related calls */
   TGeoPainter.prototype.loadMacro = function(script_name) {

      let result = { obj: this.getGeometry(), prefix: "" };

      if (this.geo_manager) result.prefix = result.obj.fName;

      // let result = { obj: draw_obj, prefix: name_prefix };

      if (!script_name || (script_name.length < 3) || (geo.getNodeKind(result.obj)!==0))
         return Promise.resolve(result);

      let painter = this;

      let mgr = {
            GetVolume: name => {
               let regexp = new RegExp("^"+name+"$");
               let currnode = painter.findNodeWithVolume(regexp, arg => arg);

               if (!currnode) console.log('Did not found '+name + ' volume');

               // return proxy object with several methods, typically used in ROOT geom scripts
               return {
                   found: currnode,
                   fVolume: currnode ? currnode.node.fVolume : null,
                   InvisibleAll: function(flag) {
                      geo.setInvisibleAll(this.fVolume, flag);
                   },
                   Draw: function() {
                      if (!this.found || !this.fVolume) return;
                      result.obj = this.found.node;
                      result.prefix = this.found.item;
                      console.log('Select volume for drawing', this.fVolume.fName, result.prefix);
                   },
                   SetTransparency: function(lvl) {
                     if (this.fVolume && this.fVolume.fMedium && this.fVolume.fMedium.fMaterial)
                        this.fVolume.fMedium.fMaterial.fFillStyle = 3000+lvl;
                   },
                   SetLineColor: function(col) {
                      if (this.fVolume) this.fVolume.fLineColor = col;
                   }
                };
            },

            DefaultColors: () => {
               painter.ctrl.dflt_colors = true;
            },

            SetMaxVisNodes: limit => {
               if (!painter.ctrl.maxnodes)
                  painter.ctrl.maxnodes = pasrseInt(limit) || 0;
            },

            SetVisLevel: limit => {
               if (!painter.ctrl.vislevel)
                  painter.ctrl.vislevel = parseInt(limit) || 0;
            }
          };

      jsrp.showProgress('Loading macro ' + script_name);

      return JSROOT.httpRequest(script_name, "text").then(script => {
         let lines = script.split('\n'), indx = 0;

         while (indx < lines.length) {
            let line = lines[indx++].trim();

            if (line.indexOf('//')==0) continue;

            if (line.indexOf('gGeoManager') < 0) continue;
            line = line.replace('->GetVolume','.GetVolume');
            line = line.replace('->InvisibleAll','.InvisibleAll');
            line = line.replace('->SetMaxVisNodes','.SetMaxVisNodes');
            line = line.replace('->DefaultColors','.DefaultColors');
            line = line.replace('->Draw','.Draw');
            line = line.replace('->SetTransparency','.SetTransparency');
            line = line.replace('->SetLineColor','.SetLineColor');
            line = line.replace('->SetVisLevel','.SetVisLevel');
            if (line.indexOf('->') >= 0) continue;

            try {
               let func = new Function('gGeoManager', line);
               func(mgr);
            } catch(err) {
               console.error('Problem by processing ' + line);
            }
         }

         return result;
      }).catch(() => {
         console.error('Fail to load ' + script_name);
         return result;
      });
   }

   /** @summary Assign clones, created outside.
    * @desc Used by geometry painter, where clones are handled by the server */
   TGeoPainter.prototype.assignClones = function(clones) {
      this._clones_owner = true;
      this._clones = clones;
   }

   /** @summary Prepare drawings
     * @desc Return value used as promise for painter */
   TGeoPainter.prototype.prepareObjectDraw = function(draw_obj, name_prefix) {

      // if did cleanup - ignore all kind of activity
      if (this.did_cleanup)
         return Promise.resolve(null);

      if (name_prefix == "__geom_viewer_append__") {
         this._new_append_nodes = draw_obj;
         this.ctrl.use_worker = 0;
         this._geom_viewer = true; // indicate that working with geom viewer
      } else if ((name_prefix == "__geom_viewer_selection__") && this._clones) {
         // these are selection done from geom viewer
         this._new_draw_nodes = draw_obj;
         this.ctrl.use_worker = 0;
         this._geom_viewer = true; // indicate that working with geom viewer
      } else if (this._main_painter) {

         this._clones_owner = false;

         this._clones = this._main_painter._clones;

         console.log('Reuse clones', this._clones.nodes.length, 'from main painter');

      } else if (!draw_obj) {

         this._clones_owner = false;

         this._clones = null;

      } else {

         this._start_drawing_time = new Date().getTime();

         this._clones_owner = true;

         this._clones = new geo.ClonedNodes(draw_obj);

         let lvl = this.ctrl.vislevel, maxnodes = this.ctrl.maxnodes;
         if (this.geo_manager) {
            if (!lvl && this.geo_manager.fVisLevel)
               lvl = this.geo_manager.fVisLevel;
            if (!maxnodes)
               maxnodes = this.geo_manager.fMaxVisNodes;
         }

         this._clones.setVisLevel(lvl);
         this._clones.setMaxVisNodes(maxnodes);

         this._clones.name_prefix = name_prefix;

         let hide_top_volume = !!this.geo_manager && !this.ctrl.showtop;

         let uniquevis = this.ctrl.no_screen ? 0 : this._clones.markVisibles(true, false, hide_top_volume);

         if (uniquevis <= 0)
            uniquevis = this._clones.markVisibles(false, false, hide_top_volume);
         else
            uniquevis = this._clones.markVisibles(true, true, hide_top_volume); // copy bits once and use normal visibility bits

         this._clones.produceIdShifts();

         let spent = new Date().getTime() - this._start_drawing_time;

         if (!this._scene)
            console.log('Creating clones', this._clones.nodes.length, 'takes', spent, 'uniquevis', uniquevis);

         if (this.options._count)
            return this.drawCount(uniquevis, spent);
      }

      if (!this._scene) {

         // this is limit for the visible faces, number of volumes does not matter
         this.ctrl.maxlimit = (this._webgl ? 200000 : 100000) * this.ctrl.more;

         this._first_drawing = true;

         // activate worker
         if (this.ctrl.use_worker > 0) this.startWorker();

         jsrp.assign3DHandler(this);

         let size = this.getSizeFor3d(this._webgl ? undefined : 3);

         this._fit_main_area = (size.can3d === -1);

         let dom = this.createScene(size.width, size.height);

         this.add3dCanvas(size, dom, this._webgl);

         // set top painter only when first child exists
         this.setAsMainPainter();
      }

      this.createToolbar();

      if (this._clones)
         return new Promise(resolveFunc => {
            this._resolveFunc = resolveFunc;
            this.showDrawInfo("Drawing geometry");
            this.startDrawGeometry(true);
         });

      this.completeDraw();
      return Promise.resolve(this);
   }

   /** @summary methods show info when first geometry drawing is performed */
   TGeoPainter.prototype.showDrawInfo = function(msg) {
      if (!this._first_drawing || !this._start_drawing_time) return;

      let main = this._renderer.domElement.parentNode,
          info = d3.select(main).select(".geo_info");

      if (!msg) {
         info.remove();
      } else {
         let spent = (new Date().getTime() - this._start_drawing_time)*1e-3;
         if (info.empty()) info = d3.select(main).append("p").attr("class","geo_info");
         info.html(msg + ", " + spent.toFixed(1) + "s");
      }
   }

   /** @summary Reentrant method to perform geometry drawing step by step */
   TGeoPainter.prototype.continueDraw = function() {

      // nothing to do - exit
      if (this.drawing_stage === 0) return;

      let tm0 = new Date().getTime(),
          interval = this._first_drawing ? 1000 : 200,
          now = tm0;

      while(true) {

         let res = this.nextDrawAction();

         if (!res) break;

         now = new Date().getTime();

         // stop creation after 100 sec, render as is
         if (now - this._startm > 1e5) {
            this.drawing_stage = 0;
            break;
         }

         // if we are that fast, do next action
         if ((res === true) && (now - tm0 < interval)) continue;

         if ((now - tm0 > interval) || (res === 1) || (res === 2)) {

            jsrp.showProgress(this.drawing_log);

            this.showDrawInfo(this.drawing_log);

            if (this._first_drawing && this._webgl && (this._num_meshes - this._last_render_meshes > 100) && (now - this._last_render_tm > 2.5*interval)) {
               this.adjustCameraPosition();
               this.render3D(-1);
               this._last_render_meshes = this.ctrl.info.num_meshes;
            }
            if (res !== 2) setTimeout(() => this.continueDraw(), (res === 1) ? 100 : 1);

            return;
         }
      }

      let take_time = now - this._startm;

      if (this._first_drawing || this._full_redrawing)
         console.log(`Create tm = ${take_time} meshes ${this.ctrl.info.num_meshes} faces ${this.ctrl.info.num_faces}`);

      if (take_time > 300) {
         jsrp.showProgress('Rendering geometry');
         this.showDrawInfo("Rendering");
         return setTimeout(this.completeDraw.bind(this, true), 10);
      }

      this.completeDraw(true);
   }

   /** @summary Checks camera position and recalculate rendering order if needed
    * @param force - if specified, forces calculations of render order */
   TGeoPainter.prototype.testCameraPosition = function(force) {
      this._camera.updateMatrixWorld();
      let origin = this._camera.position.clone();

      if (!force && this._last_camera_position) {
         // if camera position does not changed a lot, ignore such change
         let dist = this._last_camera_position.distanceTo(origin);
         if (dist < (this._overall_size || 1000)*1e-4) return;
      }

      this._last_camera_position = origin; // remember current camera position

      if (!this.ctrl.project && this._webgl)
         geo.produceRenderOrder(this._toplevel, origin, this.ctrl.depthMethod, this._clones);
   }

   /** @summary Call 3D rendering of the geometry
     * @param tmout - specifies delay, after which actual rendering will be invoked
     * @desc Timeout used to avoid multiple rendering of the picture when several 3D drawings
     * superimposed with each other. If tmeout<=0, rendering performed immediately
     * Several special values are used:
     *   -1    - force recheck of rendering order based on camera position */
   TGeoPainter.prototype.render3D = function(tmout, measure) {

      if (!this._renderer) {
         console.warn('renderer object not exists - check code');
         return;
      }

      if (tmout === undefined) tmout = 5; // by default, rendering happens with timeout

      if ((tmout > 0) && this._webgl && !JSROOT.batch_mode) {
         // do not shoot timeout many times
         if (!this.render_tmout)
            this.render_tmout = setTimeout(() => this.render3D(0,measure), tmout);
         return;
      }

      if (this.render_tmout) {
         clearTimeout(this.render_tmout);
         delete this.render_tmout;
      }

      jsrp.beforeRender3D(this._renderer);

      let tm1 = new Date();

      this.testCameraPosition(tmout === -1);

      // its needed for outlinePass - do rendering, most consuming time
      if (this._webgl && this._effectComposer && (this._effectComposer.passes.length > 0)) {
         this._effectComposer.render();
      } else if (this._webgl && this._bloomComposer && (this._bloomComposer.passes.length > 0)) {
         this._renderer.clear();
         this._camera.layers.set( _BLOOM_SCENE );
         this._bloomComposer.render();
         this._renderer.clearDepth();
         this._camera.layers.set( _ENTIRE_SCENE );
         this._renderer.render(this._scene, this._camera);
      } else {
    //     this._renderer.logarithmicDepthBuffer = true;
         this._renderer.render(this._scene, this._camera);
      }

      let tm2 = new Date();

      this.last_render_tm = tm2.getTime();

      if ((this.first_render_tm === 0) && measure) {
         this.first_render_tm = tm2.getTime() - tm1.getTime();
         console.log(`three.js r${THREE.REVISION}, first render tm = ${this.first_render_tm}`);
      }

      return jsrp.afterRender3D(this._renderer);
   }

   /** @summary Start geo worker */
   TGeoPainter.prototype.startWorker = function() {

      if (this._worker) return;

      this._worker_ready = false;
      this._worker_jobs = 0; // counter how many requests send to worker

      let painter = this;

      this._worker = new Worker(JSROOT.source_dir + "scripts/JSRoot.geoworker.js");

      this._worker.onmessage = function(e) {

         if (typeof e.data !== 'object') return;

         if ('log' in e.data)
            return console.log(`geo: ${e.data.log}`);

         if ('progress' in e.data)
            return jsrp.showProgress(e.data.progress);

         e.data.tm3 = new Date().getTime();

         if ('init' in e.data) {
            painter._worker_ready = true;
            console.log(`Worker ready: ${e.data.tm3 - e.data.tm0}`);
         } else {
            painter.processWorkerReply(e.data);
         }
      };

      // send initialization message with clones
      this._worker.postMessage({
         init: true,   // indicate init command for worker
         browser: JSROOT.browser,
         tm0: new Date().getTime(),
         vislevel: this._clones.getVisLevel(),
         maxvisnodes: this._clones.getMaxVisNodes(),
         clones: this._clones.nodes,
         sortmap: this._clones.sortmap
      });
   }

   TGeoPainter.prototype.canSubmitToWorker = function(force) {
      if (!this._worker) return false;

      return this._worker_ready && ((this._worker_jobs == 0) || force);
   }

   TGeoPainter.prototype.submitToWorker = function(job) {
      if (!this._worker) return false;

      this._worker_jobs++;

      job.tm0 = new Date().getTime();

      this._worker.postMessage(job);
   }

   TGeoPainter.prototype.processWorkerReply = function(job) {
      this._worker_jobs--;

      if ('collect' in job) {
         this._new_draw_nodes = job.new_nodes;
         this._draw_all_nodes = job.complete;
         this.drawing_stage = 3;
         // invoke methods immediately
         return this.continueDraw();
      }

      if ('shapes' in job) {

         for (let n=0;n<job.shapes.length;++n) {
            let item = job.shapes[n],
                origin = this._build_shapes[n];

            // let shape = this._clones.getNodeShape(item.nodeid);

            if (item.buf_pos && item.buf_norm) {
               if (item.buf_pos.length === 0) {
                  origin.geom = null;
               } else if (item.buf_pos.length !== item.buf_norm.length) {
                  console.error('item.buf_pos',item.buf_pos.length, 'item.buf_norm', item.buf_norm.length);
                  origin.geom = null;
               } else {
                  origin.geom = new THREE.BufferGeometry();

                  origin.geom.setAttribute( 'position', new THREE.BufferAttribute( item.buf_pos, 3 ) );
                  origin.geom.setAttribute( 'normal', new THREE.BufferAttribute( item.buf_norm, 3 ) );
               }

               origin.ready = true;
               origin.nfaces = item.nfaces;
            }
         }

         job.tm4 = new Date().getTime();

         // console.log('Get reply from worker', job.tm3-job.tm2, ' decode json in ', job.tm4-job.tm3);

         this.drawing_stage = 7; // first check which shapes are used, than build meshes

         // invoke methods immediately
         return this.continueDraw();
      }
   }

   TGeoPainter.prototype.testGeomChanges = function() {
      if (this._main_painter) {
         console.warn('Get testGeomChanges call for slave painter');
         return this._main_painter.testGeomChanges();
      }
      this.startDrawGeometry();
      for (let k=0;k<this._slave_painters.length;++k)
         this._slave_painters[k].startDrawGeometry();
   }

   /** @summary Draw axes if configured, otherwise just remove completely */
   TGeoPainter.prototype.drawSimpleAxis = function(norender) {
      this.getExtrasContainer('delete', 'axis');

      if (!this.ctrl._axis)
         return norender ? null : this.render3D();

      let box = this.getGeomBoundingBox(this._toplevel);

      let container = this.getExtrasContainer('create', 'axis');

      let text_size = 0.02 * Math.max( (box.max.x - box.min.x), (box.max.y - box.min.y), (box.max.z - box.min.z)),
          center = [0,0,0],
          names = ['x','y','z'],
          labels = ['X','Y','Z'],
          colors = ["red","green","blue"],
          ortho = this.ctrl.ortho_camera,
          yup = [this.ctrl._yup, this.ctrl._yup, this.ctrl._yup],
          numaxis = 3;

      if (this.ctrl._axis == 2)
         for (let naxis=0;naxis<3;++naxis) {
            let name = names[naxis];
            if ((box.min[name]<=0) && (box.max[name]>=0)) continue;
            center[naxis] = (box.min[name] + box.max[name])/2;
         }

      // only two dimensions are seen by ortho camera, X draws Z, can be configured better later
      if (this.ctrl.ortho_camera) {
         numaxis = 2;
         labels[0] = labels[2];
         colors[0] = colors[2];
         yup[0] = yup[2];
         ortho = true;
      }

      for (let naxis=0;naxis<numaxis;++naxis) {

         let buf = new Float32Array(6), axiscol = colors[naxis], name = names[naxis];

         function Convert(value) {
            let range = box.max[name] - box.min[name];
            if (range<2) return value.toFixed(3);
            if (Math.abs(value)>1e5) return value.toExponential(3);
            return Math.round(value).toString();
         }

         let lbl = Convert(box.max[name]);

         buf[0] = box.min.x;
         buf[1] = box.min.y;
         buf[2] = box.min.z;

         buf[3] = box.min.x;
         buf[4] = box.min.y;
         buf[5] = box.min.z;

         switch (naxis) {
           case 0: buf[3] = box.max.x; if (yup[0] && !ortho) lbl = labels[0] + " " + lbl; else lbl += " " + labels[0]; break;
           case 1: buf[4] = box.max.y; if (yup[1]) lbl += " " + labels[1]; else lbl = labels[1] + " " + lbl; break;
           case 2: buf[5] = box.max.z; lbl += " " + labels[2]; break;
         }

         if (this.ctrl._axis == 2)
            for (let k=0;k<6;++k)
               if ((k % 3) !== naxis) buf[k] = center[k%3];

         let lineMaterial = new THREE.LineBasicMaterial({ color: axiscol }),
             mesh = jsrp.createLineSegments(buf, lineMaterial);

         container.add(mesh);

         let textMaterial = new THREE.MeshBasicMaterial({ color: axiscol });

         if ((center[naxis]===0) && (center[naxis]>=box.min[name]) && (center[naxis]<=box.max[name]))
           if ((this.ctrl._axis != 2) || (naxis===0)) {
               let geom = ortho ? new THREE.CircleBufferGeometry(text_size*0.25) :
                                  new THREE.SphereBufferGeometry(text_size*0.25);
               mesh = new THREE.Mesh(geom, textMaterial);
               mesh.translateX((naxis===0) ? center[0] : buf[0]);
               mesh.translateY((naxis===1) ? center[1] : buf[1]);
               mesh.translateZ((naxis===2) ? center[2] : buf[2]);
               container.add(mesh);
           }

         let text3d = new THREE.TextGeometry(lbl, { font: JSROOT.threejs_font_helvetiker_regular, size: text_size, height: 0, curveSegments: 5 });
         mesh = new THREE.Mesh(text3d, textMaterial);
         let textbox = new THREE.Box3().setFromObject(mesh);

         mesh.translateX(buf[3]);
         mesh.translateY(buf[4]);
         mesh.translateZ(buf[5]);

         if (yup[naxis]) {
            switch (naxis) {
               case 0:
                  if (!ortho) {
                     mesh.rotateY(Math.PI);
                     mesh.translateX(-textbox.max.x-text_size*0.5);
                  } else {
                     mesh.translateX(text_size*0.5);
                  }
                  mesh.translateY(-textbox.max.y/2);
                  break;
               case 1:
                  if (!ortho) {
                     mesh.rotateX(-Math.PI/2);
                     mesh.rotateY(-Math.PI/2);
                  } else {
                     mesh.rotateZ(Math.PI/2);
                  }
                  mesh.translateX(text_size*0.5);
                  mesh.translateY(-textbox.max.y/2);
                  break;
               case 2: mesh.rotateY(-Math.PI/2); mesh.translateX(text_size*0.5); mesh.translateY(-textbox.max.y/2); break;
           }
         } else {
            switch (naxis) {
               case 0: mesh.rotateX(Math.PI/2); mesh.translateY(-textbox.max.y/2); mesh.translateX(text_size*0.5); break;
               case 1: mesh.rotateX(Math.PI/2); mesh.rotateY(-Math.PI/2); mesh.translateX(-textbox.max.x-text_size*0.5); mesh.translateY(-textbox.max.y/2); break;
               case 2: mesh.rotateX(Math.PI/2); mesh.rotateZ(Math.PI/2); mesh.translateX(text_size*0.5); mesh.translateY(-textbox.max.y/2); break;
            }
         }

         container.add(mesh);

         text3d = new THREE.TextGeometry(Convert(box.min[name]), { font: JSROOT.threejs_font_helvetiker_regular, size: text_size, height: 0, curveSegments: 5 });

         mesh = new THREE.Mesh(text3d, textMaterial);
         textbox = new THREE.Box3().setFromObject(mesh);

         mesh.translateX(buf[0]);
         mesh.translateY(buf[1]);
         mesh.translateZ(buf[2]);

         if (yup[naxis]) {
            switch (naxis) {
               case 0:
                  if (!ortho) {
                     mesh.rotateY(Math.PI);
                     mesh.translateX(text_size*0.5);
                  } else {
                     mesh.translateX(-textbox.max.x-text_size*0.5);
                  }
                  mesh.translateY(-textbox.max.y/2);
                  break;
               case 1:
                  if (!ortho) {
                     mesh.rotateX(-Math.PI/2);
                     mesh.rotateY(-Math.PI/2);
                  } else {
                     mesh.rotateZ(Math.PI/2);
                  }
                  mesh.translateY(-textbox.max.y/2);
                  mesh.translateX(-textbox.max.x-text_size*0.5);
                  break;
               case 2: mesh.rotateY(-Math.PI/2);  mesh.translateX(-textbox.max.x-text_size*0.5); mesh.translateY(-textbox.max.y/2); break;
            }
         } else {
            switch (naxis) {
               case 0: mesh.rotateX(Math.PI/2); mesh.translateX(-textbox.max.x-text_size*0.5); mesh.translateY(-textbox.max.y/2); break;
               case 1: mesh.rotateX(Math.PI/2); mesh.rotateY(-Math.PI/2); mesh.translateY(-textbox.max.y/2); mesh.translateX(text_size*0.5); break;
               case 2: mesh.rotateX(Math.PI/2); mesh.rotateZ(Math.PI/2);  mesh.translateX(-textbox.max.x-text_size*0.5); mesh.translateY(-textbox.max.y/2); break;
            }
         }

         container.add(mesh);
      }

      // after creating axes trigger rendering and recalculation of depth
      this.changedDepthMethod(norender ? "norender" : undefined);
   }

   /** @summary  Toggle axes visibility */
   TGeoPainter.prototype.toggleAxesDraw = function() {
      this.setAxesDraw("toggle");
   }

   /** @summary Set axes visibility 0 - off, 1 - on, 2 - centered */
   TGeoPainter.prototype.setAxesDraw = function(on) {
      if (on === "toggle")
         this.ctrl._axis = this.ctrl._axis ? 0 : 1;
      else
         this.ctrl._axis = (typeof on == 'number') ? on : (on ? 1 : 0);
      this.drawSimpleAxis();
   }

   /** @summary Set auto rotate mode */
   TGeoPainter.prototype.setAutoRotate = function(on) {
      if (this.ctrl.project) return;
      if (on !== undefined) this.ctrl.rotate = on;
      this.autorotate(2.5);
   }

   /** @summary Toggle wireframe mode */
   TGeoPainter.prototype.toggleWireFrame = function() {
      this.ctrl.wireframe = !this.ctrl.wireframe;
      this.changedWireFrame();
   }

   /** @summary Specify wireframe mode */
   TGeoPainter.prototype.setWireFrame = function(on) {
      this.ctrl.wireframe = on ? true : false;
      this.changedWireFrame();
   }

   /** @summary Specify showtop draw options, relevant only for TGeoManager */
   TGeoPainter.prototype.setShowTop = function(on) {
      this.ctrl.showtop = on ? true : false;
      this.redrawObject('same');
   }

   /** @summary Should be called when configuration of particular axis is changed */
   TGeoPainter.prototype.changedClipping = function(naxis) {
      let clip = this.ctrl.clip;

      if ((naxis !== undefined) && (naxis >= 0)) {
         if (!clip[naxis].enabled) return;
      }

      if (clip[0].enabled || clip[1].enabled || clip[2].enabled) {
         this.ctrl.ssao.enabled = false;
         this.removeSSAO();
      }

      this.updateClipping(false, true);
   }

   /** @summary Should be called when depth test flag is changed */
   TGeoPainter.prototype.changedDepthTest = function() {
      if (!this._toplevel) return;
      let flag = this.ctrl.depthTest;
      this._toplevel.traverse( function (node) {
         if (node instanceof THREE.Mesh) {
            node.material.depthTest = flag;
         }
      });

      this.render3D(0);
   }

   /** @summary Should be called when depth method is changed */
   TGeoPainter.prototype.changedDepthMethod = function(arg) {
      // force recalculatiion of render order
      delete this._last_camera_position;
      if (arg !== "norender") this.render3D();
   }

   /** @summary Should be called when configuration of highlight is changed */
   TGeoPainter.prototype.changedHighlight = function() {
      if (!this.ctrl.highlight)
         this.highlightMesh(null);
   }

   /** @summary Assign clipping attributes to the meshes - supported only for webgl */
   TGeoPainter.prototype.updateClipping = function(without_render, force_traverse) {
      if (!this._webgl) return;

      let clip = this.ctrl.clip, panels = [], changed = false,
          constants = [ clip[0].value, -1 * clip[1].value, (this.ctrl._yup ? -1 : 1) * clip[2].value ],
          clip_cfg = this.ctrl.clipIntersect ? 16 : 0;

      for (let k=0;k<3;++k) {
         if (clip[k].enabled) clip_cfg += 2 << k;
         if (this._clipPlanes[k].constant != constants[k]) {
            changed = true;
            this._clipPlanes[k].constant = constants[k];
         }
      }

      if (!this.ctrl.ssao.enabled) {
         if (clip[0].enabled) panels.push(this._clipPlanes[0]);
         if (clip[1].enabled) panels.push(this._clipPlanes[1]);
         if (clip[2].enabled) panels.push(this._clipPlanes[2]);
         clip_cfg += panels.length*1000;
      }
      if (panels.length == 0) panels = null;

      if (this._clipCfg !== clip_cfg) changed = true;

      this._clipCfg = clip_cfg;

      let any_clipping = !!panels, ci = this.ctrl.clipIntersect,
          material_side = any_clipping ? THREE.DoubleSide : THREE.FrontSide;

      if (force_traverse || changed)
      this._scene.traverse( function (node) {
         if (node.hasOwnProperty("material") && node.material && (node.material.clippingPlanes !== undefined)) {

            if (node.material.clippingPlanes !== panels) {
               node.material.clipIntersection = ci;
               node.material.clippingPlanes = panels;
               node.material.needsUpdate = true;
            }

            if (node.material.emissive !== undefined) {
               if (node.material.side != material_side) {
                  node.material.side = material_side;
                  node.material.needsUpdate = true;
               }
            }
         }
      });

      this.ctrl.bothSides = any_clipping;

      if (!without_render) this.render3D(0);

      return changed;
   }

   /** @summary Assign callback, invoked every time when drawing is completed
     * @desc Used together with web-based geometry viewer
     * @private */
   TGeoPainter.prototype.setCompleteHandler = function(callback) {
      this._complete_handler = callback;
   }

   /** @summary Completes drawing procedure */
   TGeoPainter.prototype.completeDraw = function(close_progress) {

      let first_time = false, full_redraw = false, check_extras = true;

      if (!this.ctrl) {
         console.warn('ctrl object does not exist in completeDraw - something went wrong');
         return;
      }

      if (!this._clones) {
         check_extras = false;
         // if extra object where append, redraw them at the end
         this.getExtrasContainer("delete"); // delete old container
         let extras = (this._main_painter ? this._main_painter._extraObjects : null) || this._extraObjects;
         this.drawExtras(extras, "", false);
      } else if (this._first_drawing || this._full_redrawing) {
         if (this.ctrl.tracks && this.geo_manager)
            this.drawExtras(this.geo_manager.fTracks, "<prnt>/Tracks");
      }

      if (this._full_redrawing) {
         this.adjustCameraPosition(true);
         this._full_redrawing = false;
         full_redraw = true;
         this.changedDepthMethod("norender");
      }

      if (this._first_drawing) {
         this.adjustCameraPosition(true);
         this.showDrawInfo();
         this._first_drawing = false;
         first_time = true;
         full_redraw = true;
      }

      if (this.ctrl.transparency!==0)
         this.changedGlobalTransparency(this.ctrl.transparency, true);

      if (first_time)
         this.completeScene();

      if (full_redraw && (this.ctrl.trans_radial || this.ctrl.trans_z))
         this.changedTransformation("norender");

      if (full_redraw && this.ctrl._axis)
         this.drawSimpleAxis(true);

      this._scene.overrideMaterial = null;

      if (this._provided_more_nodes !== undefined) {
         this.appendMoreNodes(this._provided_more_nodes, true);
         delete this._provided_more_nodes;
      }

      if (check_extras) {
         // if extra object where append, redraw them at the end
         this.getExtrasContainer("delete"); // delete old container
         let extras = (this._main_painter ? this._main_painter._extraObjects : null) || this._extraObjects;
         this.drawExtras(extras, "", false);
      }

      this.updateClipping(true); // do not render

      this.render3D(0, true);

      if (close_progress) jsrp.showProgress();

      this.addOrbitControls();

      this.addTransformControl();

      if (first_time) {

         // after first draw check if highlight can be enabled
         if (this.ctrl.highlight === false)
            this.ctrl.highlight = (this.first_render_tm < 1000);

         // also highlight of scene object can be assigned at the first draw
         if (this.ctrl.highlight_scene === false)
            this.ctrl.highlight_scene = this.ctrl.highlight;

         // if rotation was enabled, do it
         if (this._webgl && this.ctrl.rotate && !this.ctrl.project) this.autorotate(2.5);
         if (this._webgl && this.ctrl.show_controls && !JSROOT.batch_mode) this.showControlOptions(true);
      }

      this.setAsMainPainter();

      if (typeof this._resolveFunc == 'function') {
         this._resolveFunc(this);
         delete this._resolveFunc;
      }

      if (typeof this._complete_handler == 'function')
         this._complete_handler(this);

      if (this._draw_nodes_again)
         return this.startDrawGeometry(); // relaunch drawing

      this._drawing_ready = true; // indicate that drawing is completed
   }

   /** @summary Returns true if geometry drawing is completed */
   TGeoPainter.prototype.isDrawingReady = function() {
      return this._drawing_ready || false;
   }

   /** @summary Remove already drawn node. Used by geom viewer */
   TGeoPainter.prototype.removeDrawnNode = function(nodeid) {
      if (!this._draw_nodes) return;

      let new_nodes = [];

      for (let n = 0; n < this._draw_nodes.length; ++n) {
         let entry = this._draw_nodes[n];
         if ((entry.nodeid === nodeid) || this._clones.isIdInStack(nodeid, entry.stack)) {
            this._clones.createObject3D(entry.stack, this._toplevel, 'delete_mesh');
         } else {
            new_nodes.push(entry);
         }
      }

      if (new_nodes.length < this._draw_nodes.length) {
         this._draw_nodes = new_nodes;
         this.render3D();
      }
   }

   /** @summary Cleanup geometry painter */
   TGeoPainter.prototype.cleanup = function(first_time) {

      if (!first_time) {

         this.removeSSAO();

         this.clearTopPainter(); // remove as pointer

         let can3d = this.clear3dCanvas(); // remove 3d canvas from main HTML element

         if (this._toolbar) this._toolbar.cleanup(); // remove toolbar

         this.helpText();

         jsrp.disposeThreejsObject(this._scene);

         jsrp.disposeThreejsObject(this._full_geom);

         if (this._tcontrols)
            this._tcontrols.dispose();

         if (this._controls)
            this._controls.cleanup();

         if (this._context_menu)
            this._renderer.domElement.removeEventListener( 'contextmenu', this._context_menu, false );

         if (this._datgui)
            this._datgui.destroy();

         if (this._worker) this._worker.terminate();

         delete this._animating;

         let obj = this.getGeometry();
         if (obj && this.ctrl.is_main) {
            if (obj.$geo_painter===this) delete obj.$geo_painter; else
            if (obj.fVolume && obj.fVolume.$geo_painter===this) delete obj.fVolume.$geo_painter;
         }

         if (this._main_painter) {
            let pos = this._main_painter._slave_painters.indexOf(this);
            if (pos>=0) this._main_painter._slave_painters.splice(pos,1);
         }

         for (let k=0;k<this._slave_painters.length;++k) {
            let slave = this._slave_painters[k];
            if (slave && (slave._main_painter===this)) slave._main_painter = null;
         }

         delete this.geo_manager;
         delete this._highlight_handlers;

         JSROOT.ObjectPainter.prototype.cleanup.call(this);

         delete this.ctrl;
         delete this.options;

         this.did_cleanup = true;

         if (can3d < 0) this.selectDom().html("");
      }

      if (this._slave_painters)
         for (let k in this._slave_painters) {
            let slave = this._slave_painters[k];
            slave._main_painter = null;
            if (slave._clones === this._clones) slave._clones = null;
         }

      this._main_painter = null;
      this._slave_painters = [];

      if (this._renderer) {
         if (this._renderer.dispose) this._renderer.dispose();
         if (this._renderer.forceContextLoss) this._renderer.forceContextLoss();
      }

      delete this._scene;
      this._scene_width = 0;
      this._scene_height = 0;
      this._renderer = null;
      this._toplevel = null;
      delete this._full_geom;
      delete this._camera;
      delete this._camera0pos;
      delete this._lookat;
      delete this._selected_mesh;

      if (this._clones && this._clones_owner)
         this._clones.cleanup(this._draw_nodes, this._build_shapes);
      delete this._clones;
      delete this._clones_owner;
      delete this._draw_nodes;
      delete this._drawing_ready;
      delete this._build_shapes;
      delete this._new_draw_nodes;
      delete this._new_append_nodes;
      delete this._last_camera_position;

      this.first_render_tm = 0; // time needed for first rendering
      this.last_render_tm = 0;

      this.drawing_stage = 0;
      delete this.drawing_log;

      delete this._datgui;
      delete this._controls;
      delete this._context_menu;
      delete this._tcontrols;
      delete this._toolbar;

      delete this._worker;
   }

   TGeoPainter.prototype.helpText = function(msg) {
      jsrp.showProgress(msg);
   }

   /** @summary Check if HTML element was resized and drawing need to be adjusted */
   TGeoPainter.prototype.checkResize = function(arg) {
      let cp = this.getCanvPainter();

      // firefox is the only browser which correctly supports resize of embedded canvas,
      // for others we should force canvas redrawing at every step
      if (cp && !cp.checkCanvasResize(arg)) return false;

      let sz = this.getSizeFor3d();

      if ((this._scene_width === sz.width) && (this._scene_height === sz.height)) return false;
      if ((sz.width<10) || (sz.height<10)) return false;

      this._scene_width = sz.width;
      this._scene_height = sz.height;

      if (this._camera && this._renderer) {
         if (this._camera.type == "PerspectiveCamera")
            this._camera.aspect = this._scene_width / this._scene_height;
         this._camera.updateProjectionMatrix();
         this._renderer.setSize( this._scene_width, this._scene_height, !this._fit_main_area );
         if (this._effectComposer)
            this._effectComposer.setSize( this._scene_width, this._scene_height );
         if (this._bloomComposer)
            this._bloomComposer.setSize( this._scene_width, this._scene_height );

         if (!this.drawing_stage) this.render3D();
      }

      return true;
   }

   /** @summary Toggle enlarge state */
   TGeoPainter.prototype.toggleEnlarge = function() {
      if (this.enlargeMain('toggle'))
        this.checkResize();
   }

   TGeoPainter.prototype.ownedByTransformControls = function(child) {
      let obj = child.parent;
      while (obj && !(obj instanceof THREE.TransformControls) ) {
         obj = obj.parent;
      }
      return (obj && (obj instanceof THREE.TransformControls));
   }

   TGeoPainter.prototype.accessObjectWireFrame = function(obj, on) {
      // either change mesh wireframe or return current value
      // return undefined when wireframe cannot be accessed

      if (!obj.hasOwnProperty("material") || (obj instanceof THREE.GridHelper)) return;

      if (this.ownedByTransformControls(obj)) return;

      if ((on !== undefined) && obj.stack)
         obj.material.wireframe = on;

      return obj.material.wireframe;
   }

   TGeoPainter.prototype.changedWireFrame = function() {
      if (!this._scene) return;

      let on = this.ctrl.wireframe;

      this._scene.traverse(obj => this.accessObjectWireFrame(obj, on));

      this.render3D();
   }

   /** @summary Update object in geo painter */
   TGeoPainter.prototype.updateObject = function(obj) {
      if (obj === "same") return true;
      if (!obj || !obj._typename) return false;
      if (obj === this.getObject()) return true;

      if (this.geo_manager && (obj._typename == "TGeoManager")) {
         this.geo_manager = obj;
         this.assignObject({ _typename:"TGeoNode", fVolume: obj.fMasterVolume, fName: obj.fMasterVolume.fName, $geoh: obj.fMasterVolume.$geoh, _proxy: true });
         return true;
      }

      if (!this.matchObjectType(obj._typename)) return false;

      this.assignObject(obj);
      return true;
   }

   /** @summary Cleanup TGeo drawings */
   TGeoPainter.prototype.clearDrawings = function() {
      if (this._clones && this._clones_owner)
         this._clones.cleanup(this._draw_nodes, this._build_shapes);
      delete this._clones;
      delete this._clones_owner;
      delete this._draw_nodes;
      delete this._drawing_ready;
      delete this._build_shapes;

      delete this._extraObjects;
      delete this._clipCfg;

      // only remove all childs from top level object
      jsrp.disposeThreejsObject(this._toplevel, true);

      this._full_redrawing = true;
   }

   /** @summary Redraw TGeo object */
   TGeoPainter.prototype.redrawObject = function(obj) {
      if (!this.updateObject(obj))
         return false;

      this.clearDrawings();

      let draw_obj = this.getGeometry(), name_prefix = "";
      if (this.geo_manager) name_prefix = draw_obj.fName;

      this.prepareObjectDraw(draw_obj, name_prefix);

      return true;
   }

   /** @summary Create geo painter
     * @private */
   jsrp.createGeoPainter = function(divid, obj, opt) {
      geo.GradPerSegm = JSROOT.settings.GeoGradPerSegm;
      geo.CompressComp = JSROOT.settings.GeoCompressComp;

      let painter = new TGeoPainter(divid, obj);

      // one could use TGeoManager setting, but for some example JSROOT does not build composites
      // if (obj && obj._typename=='TGeoManager' && (obj.fNsegments > 3))
      //   geo.GradPerSegm = 360/obj.fNsegments;

      // painter.addToPadPrimitives(); // will add to pad primitives if any

      painter.options = painter.decodeOptions(opt); // indicator of initialization

      // copy all attributes from options to control
      JSROOT.extend(painter.ctrl, painter.options);

      painter.ctrl.ssao.enabled = painter.options.usessao;
      painter.ctrl.bloom.enabled = painter.options.usebloom;

      // special handling for array of clips
      painter.ctrl.clip[0].enabled = painter.options.clipx;
      painter.ctrl.clip[1].enabled = painter.options.clipy;
      painter.ctrl.clip[2].enabled = painter.options.clipz;

      return painter;
   }

   let drawGeoObject = (divid, obj, opt) => {
      if (!obj) return null;

      let shape = null, extras = null, extras_path = "", is_eve = false;

      if (('fShapeBits' in obj) && ('fShapeId' in obj)) {
         shape = obj; obj = null;
      } else if ((obj._typename === 'TGeoVolumeAssembly') || (obj._typename === 'TGeoVolume')) {
         shape = obj.fShape;
      } else if ((obj._typename === "TEveGeoShapeExtract") || (obj._typename === "ROOT::Experimental::REveGeoShapeExtract")) {
         shape = obj.fShape; is_eve = true;
      } else if (obj._typename === 'TGeoManager') {
         shape = obj.fMasterVolume.fShape;
      } else if (obj._typename === 'TGeoOverlap') {
         extras = obj.fMarker; extras_path = "<prnt>/Marker";
         obj = geo.buildOverlapVolume(obj);
         if (!opt) opt = "wire";
      } else if ('fVolume' in obj) {
         if (obj.fVolume) shape = obj.fVolume.fShape;
      } else {
         obj = null;
      }

      if ((typeof opt == "string") && opt.indexOf("comp")==0 && shape && (shape._typename == 'TGeoCompositeShape') && shape.fNode) {
         let maxlvl = 1;
         opt = opt.substr(4);
         if (opt[0] == "x") {  maxlvl = 999; opt = opt.substr(1) + "_vislvl999"; }
         obj = geo.buildCompositeVolume(shape, maxlvl);
      }

      if (!obj && shape)
         obj = JSROOT.extend(JSROOT.create("TEveGeoShapeExtract"),
                   { fTrans: null, fShape: shape, fRGBA: [0, 1, 0, 1], fElements: null, fRnrSelf: true });

      if (!obj) return null;

      let painter = jsrp.createGeoPainter(divid, obj, opt);

      if (painter.ctrl.is_main && !obj.$geo_painter)
         obj.$geo_painter = painter;

      if (!painter.ctrl.is_main && painter.ctrl.project && obj.$geo_painter) {
         painter._main_painter = obj.$geo_painter;
         painter._main_painter._slave_painters.push(painter);
      }

      if (is_eve && !painter.ctrl.vislevel || (painter.ctrl.vislevel < 9))
         painter.ctrl.vislevel = 9;

      if (extras) {
         painter._splitColors = true;
         painter.addExtra(extras, extras_path);
      }

      return painter.loadMacro(painter.ctrl.script_name).then(arg => painter.prepareObjectDraw(arg.obj, arg.prefix));
   }

   // ===============================================================================

   /** @summary Function used to build hierarchy of elements of composite shapes
    * @private */
   geo.buildCompositeVolume = function(comp, maxlvl, side) {

      if (maxlvl === undefined) maxlvl = 1;
      if (!side) {
         this.$comp_col_cnt = 0;
         side = "";
      }

      let vol = JSROOT.create("TGeoVolume");
      geo.SetBit(vol, geo.BITS.kVisThis, true);
      geo.SetBit(vol, geo.BITS.kVisDaughters, true);

      if ((side && (comp._typename!=='TGeoCompositeShape')) || (maxlvl<=0)) {
         vol.fName = side;
         vol.fLineColor = (this.$comp_col_cnt++ % 8) + 2;
         vol.fShape = comp;
         return vol;
      }

      if (side) side += "/";
      vol.$geoh = true; // workaround, let know browser that we are in volumes hierarchy
      vol.fName = "";

      let node1 = JSROOT.create("TGeoNodeMatrix");
      geo.SetBit(node1, geo.BITS.kVisThis, true);
      geo.SetBit(node1, geo.BITS.kVisDaughters, true);
      node1.fName = "Left";
      node1.fMatrix = comp.fNode.fLeftMat;
      node1.fVolume = geo.buildCompositeVolume(comp.fNode.fLeft, maxlvl-1, side + "Left");

      let node2 = JSROOT.create("TGeoNodeMatrix");
      geo.SetBit(node2, geo.BITS.kVisThis, true);
      geo.SetBit(node2, geo.BITS.kVisDaughters, true);
      node2.fName = "Right";
      node2.fMatrix = comp.fNode.fRightMat;
      node2.fVolume = geo.buildCompositeVolume(comp.fNode.fRight, maxlvl-1, side + "Right");

      vol.fNodes = JSROOT.create("TList");
      vol.fNodes.Add(node1);
      vol.fNodes.Add(node2);

      if (!side) delete this.$comp_col_cnt;

      return vol;
   }

   /** @summary Function used to build hierarchy of elements of overlap object
    * @private */
   geo.buildOverlapVolume = function(overlap) {

      let vol = JSROOT.create("TGeoVolume");

      geo.SetBit(vol, geo.BITS.kVisDaughters, true);
      vol.$geoh = true; // workaround, let know browser that we are in volumes hierarchy
      vol.fName = "";

      let node1 = JSROOT.create("TGeoNodeMatrix");
      node1.fName = overlap.fVolume1.fName || "Overlap1";
      node1.fMatrix = overlap.fMatrix1;
      node1.fVolume = overlap.fVolume1;
      // node1.fVolume.fLineColor = 2; // color assigned with _splitColors

      let node2 = JSROOT.create("TGeoNodeMatrix");
      node2.fName = overlap.fVolume2.fName || "Overlap2";
      node2.fMatrix = overlap.fMatrix2;
      node2.fVolume = overlap.fVolume2;
      // node2.fVolume.fLineColor = 3;  // color assigned with _splitColors

      vol.fNodes = JSROOT.create("TList");
      vol.fNodes.Add(node1);
      vol.fNodes.Add(node2);

      return vol;
   }

   geo.provideVisStyle = function(obj) {
      if ((obj._typename === 'TEveGeoShapeExtract') || (obj._typename === 'ROOT::Experimental::REveGeoShapeExtract'))
         return obj.fRnrSelf ? " geovis_this" : "";

      let vis = !geo.TestBit(obj, geo.BITS.kVisNone) &&
                geo.TestBit(obj, geo.BITS.kVisThis),
          chld = geo.TestBit(obj, geo.BITS.kVisDaughters);

      if (chld && (!obj.fNodes || (obj.fNodes.arr.length === 0))) chld = false;

      if (vis && chld) return " geovis_all";
      if (vis) return " geovis_this";
      if (chld) return " geovis_daughters";
      return "";
   }


   geo.createItem = function(node, obj, name) {
      let sub = {
         _kind: "ROOT." + obj._typename,
         _name: name ? name : geo.getObjectName(obj),
         _title: obj.fTitle,
         _parent: node,
         _geoobj: obj,
         _get: function(item /* ,itemname */) {
             // mark object as belong to the hierarchy, require to
             if (item._geoobj) item._geoobj.$geoh = true;
             return Promise.resolve(item._geoobj);
         }
      };

      let volume, shape, subnodes, iseve = false;

      if (obj._typename == "TGeoMaterial") sub._icon = "img_geomaterial"; else
      if (obj._typename == "TGeoMedium") sub._icon = "img_geomedium"; else
      if (obj._typename == "TGeoMixture") sub._icon = "img_geomixture"; else
      if ((obj._typename.indexOf("TGeoNode")===0) && obj.fVolume) {
         sub._title = "node:"  + obj._typename;
         if (obj.fTitle.length > 0) sub._title += " " + obj.fTitle;
         volume = obj.fVolume;
      } else if (obj._typename.indexOf("TGeoVolume")===0) {
         volume = obj;
      } else if ((obj._typename == "TEveGeoShapeExtract") || (obj._typename == "ROOT::Experimental::REveGeoShapeExtract")) {
         iseve = true;
         shape = obj.fShape;
         subnodes = obj.fElements ? obj.fElements.arr : null;
      } else if ((obj.fShapeBits !== undefined) && (obj.fShapeId !== undefined)) {
         shape = obj;
      }

      if (volume) {
         shape = volume.fShape;
         subnodes = volume.fNodes ? volume.fNodes.arr : null;
      }

      if (volume || shape || subnodes) {
         if (volume) sub._volume = volume;

         if (subnodes) {
            sub._more = true;
            sub._expand = geo.expandObject;
         } else
         if (shape && (shape._typename === "TGeoCompositeShape") && shape.fNode) {
            sub._more = true;
            sub._shape = shape;
            sub._expand = function(node /*, obj */) {
               geo.createItem(node, node._shape.fNode.fLeft, 'Left');
               geo.createItem(node, node._shape.fNode.fRight, 'Right');
               return true;
            }
         }

         if (!sub._title && (obj._typename != "TGeoVolume")) sub._title = obj._typename;

         if (shape) {
            if (sub._title == "")
               sub._title = shape._typename;

            sub._icon = geo.getShapeIcon(shape);
         } else {
            sub._icon = sub._more ? "img_geocombi" : "img_geobbox";
         }

         if (volume)
            sub._icon += geo.provideVisStyle(volume);
         else if (iseve)
            sub._icon += geo.provideVisStyle(obj);

         sub._menu = geo.provideMenu;
         sub._icon_click  = geo.browserIconClick;
      }

      if (!node._childs) node._childs = [];

      if (!sub._name)
         if (typeof node._name === 'string') {
            sub._name = node._name;
            if (sub._name.lastIndexOf("s")===sub._name.length-1)
               sub._name = sub._name.substr(0, sub._name.length-1);
            sub._name += "_" + node._childs.length;
         } else {
            sub._name = "item_" + node._childs.length;
         }

      node._childs.push(sub);

      return sub;
   }

   geo.createList = function(parent, lst, name, title) {

      if (!lst || !('arr' in lst) || (lst.arr.length==0)) return;

      let item = {
          _name: name,
          _kind: "ROOT.TList",
          _title: title,
          _more: true,
          _geoobj: lst,
          _parent: parent,
      }

      item._get = function(item /*, itemname */) {
         return Promise.resolve(item._geoobj || null);
      }

      item._expand = function(node, lst) {
         // only childs

         if ('fVolume' in lst)
            lst = lst.fVolume.fNodes;

         if (!('arr' in lst)) return false;

         node._childs = [];

         geo.checkDuplicates(null, lst.arr);

         for (let n in lst.arr)
            geo.createItem(node, lst.arr[n]);

         return true;
      }

      if (!parent._childs) parent._childs = [];
      parent._childs.push(item);
   }

   geo.provideMenu = function(menu, item, hpainter) {

      if (!item._geoobj) return false;

      let obj = item._geoobj, vol = item._volume,
          iseve = ((obj._typename === 'TEveGeoShapeExtract') || (obj._typename === 'ROOT::Experimental::REveGeoShapeExtract'));

      if (!vol && !iseve) return false;

      menu.add("separator");

      function ScanEveVisible(obj, arg, skip_this) {
         if (!arg) arg = { visible: 0, hidden: 0 };

         if (!skip_this) {
            if (arg.assign!==undefined) obj.fRnrSelf = arg.assign; else
            if (obj.fRnrSelf) arg.vis++; else arg.hidden++;
         }

         if (obj.fElements)
            for (let n=0;n<obj.fElements.arr.length;++n)
               ScanEveVisible(obj.fElements.arr[n], arg, false);

         return arg;
      }

      function ToggleEveVisibility(arg) {
         if (arg === 'self') {
            obj.fRnrSelf = !obj.fRnrSelf;
            item._icon = item._icon.split(" ")[0] + geo.provideVisStyle(obj);
            hpainter.updateTreeNode(item);
         } else {
            ScanEveVisible(obj, { assign: (arg === "true") }, true);
            hpainter.forEachItem(m => {
               // update all child items
               if (m._geoobj && m._icon) {
                  m._icon = item._icon.split(" ")[0] + geo.provideVisStyle(m._geoobj);
                  hpainter.updateTreeNode(m);
               }
            }, item);
         }

         geo.findItemWithPainter(item, 'testGeomChanges');
      }

      function ToggleMenuBit(arg) {
         geo.ToggleBit(vol, arg);
         let newname = item._icon.split(" ")[0] + geo.provideVisStyle(vol);
         hpainter.forEachItem(m => {
            // update all items with that volume
            if (item._volume === m._volume) {
               m._icon = newname;
               hpainter.updateTreeNode(m);
            }
         });

         hpainter.updateTreeNode(item);
         geo.findItemWithPainter(item, 'testGeomChanges');
      }

      if ((item._geoobj._typename.indexOf("TGeoNode")===0) && geo.findItemWithPainter(item))
         menu.add("Focus", function() {

           let drawitem = geo.findItemWithPainter(item);

           if (!drawitem) return;

           let fullname = hpainter.itemFullName(item, drawitem);

           if (drawitem._painter && typeof drawitem._painter.focusOnItem == 'function')
              drawitem._painter.focusOnItem(fullname);
         });

      if (iseve) {
         menu.addchk(obj.fRnrSelf, "Visible", "self", ToggleEveVisibility);
         let res = ScanEveVisible(obj, undefined, true);

         if (res.hidden + res.visible > 0)
            menu.addchk((res.hidden==0), "Daughters", (res.hidden!=0) ? "true" : "false", ToggleEveVisibility);

      } else {
         menu.addchk(geo.TestBit(vol, geo.BITS.kVisNone), "Invisible",
               geo.BITS.kVisNone, ToggleMenuBit);
         menu.addchk(geo.TestBit(vol, geo.BITS.kVisThis), "Visible",
               geo.BITS.kVisThis, ToggleMenuBit);
         menu.addchk(geo.TestBit(vol, geo.BITS.kVisDaughters), "Daughters",
               geo.BITS.kVisDaughters, ToggleMenuBit);
      }

      return true;
   }

   geo.findItemWithPainter = function(hitem, funcname) {
      while (hitem) {
         if (hitem._painter && hitem._painter._camera) {
            if (funcname && typeof hitem._painter[funcname] == 'function')
               hitem._painter[funcname]();
            return hitem;
         }
         hitem = hitem._parent;
      }
      return null;
   }

   geo.updateBrowserIcons = function(obj, hpainter) {
      if (!obj || !hpainter) return;

      hpainter.forEachItem(m => {
         // update all items with that volume
         if ((obj === m._volume) || (obj === m._geoobj)) {
            m._icon = m._icon.split(" ")[0] + geo.provideVisStyle(obj);
            hpainter.updateTreeNode(m);
         }
      });
   }

   geo.browserIconClick = function(hitem, hpainter) {
      if (hitem._volume) {
         if (hitem._more && hitem._volume.fNodes && (hitem._volume.fNodes.arr.length>0))
            geo.ToggleBit(hitem._volume, geo.BITS.kVisDaughters);
         else
            geo.ToggleBit(hitem._volume, geo.BITS.kVisThis);

         geo.updateBrowserIcons(hitem._volume, hpainter);

         geo.findItemWithPainter(hitem, 'testGeomChanges');
         return false; // no need to update icon - we did it ourself
      }

      if (hitem._geoobj && (( hitem._geoobj._typename == "TEveGeoShapeExtract") || ( hitem._geoobj._typename == "ROOT::Experimental::REveGeoShapeExtract"))) {
         hitem._geoobj.fRnrSelf = !hitem._geoobj.fRnrSelf;

         geo.updateBrowserIcons(hitem._geoobj, hpainter);
         geo.findItemWithPainter(hitem, 'testGeomChanges');
         return false; // no need to update icon - we did it ourself
      }


      // first check that geo painter assigned with the item
      let drawitem = geo.findItemWithPainter(hitem);
      if (!drawitem) return false;

      let newstate = drawitem._painter.extraObjectVisible(hpainter, hitem, true);

      // return true means browser should update icon for the item
      return (newstate!==undefined) ? true : false;
   }

   geo.getShapeIcon = function(shape) {
      switch (shape._typename) {
         case "TGeoArb8" : return "img_geoarb8";
         case "TGeoCone" : return "img_geocone";
         case "TGeoConeSeg" : return "img_geoconeseg";
         case "TGeoCompositeShape" : return "img_geocomposite";
         case "TGeoTube" : return "img_geotube";
         case "TGeoTubeSeg" : return "img_geotubeseg";
         case "TGeoPara" : return "img_geopara";
         case "TGeoParaboloid" : return "img_geoparab";
         case "TGeoPcon" : return "img_geopcon";
         case "TGeoPgon" : return "img_geopgon";
         case "TGeoShapeAssembly" : return "img_geoassembly";
         case "TGeoSphere" : return "img_geosphere";
         case "TGeoTorus" : return "img_geotorus";
         case "TGeoTrd1" : return "img_geotrd1";
         case "TGeoTrd2" : return "img_geotrd2";
         case "TGeoXtru" : return "img_geoxtru";
         case "TGeoTrap" : return "img_geotrap";
         case "TGeoGtra" : return "img_geogtra";
         case "TGeoEltu" : return "img_geoeltu";
         case "TGeoHype" : return "img_geohype";
         case "TGeoCtub" : return "img_geoctub";
      }
      return "img_geotube";
   }

   /** @summary Get icon for the browser
     * @private */
   geo.getBrowserIcon = function(hitem, hpainter) {
      let icon = "";
      if (hitem._kind == 'ROOT.TEveTrack') icon = 'img_evetrack'; else
      if (hitem._kind == 'ROOT.TEvePointSet') icon = 'img_evepoints'; else
      if (hitem._kind == 'ROOT.TPolyMarker3D') icon = 'img_evepoints';
      if (icon.length>0) {
         let drawitem = geo.findItemWithPainter(hitem);
         if (drawitem)
            if (drawitem._painter.extraObjectVisible(hpainter, hitem))
               icon += " geovis_this";
      }
      return icon;
   }

   /** @summary Expand geo object
     * @private */
   geo.expandObject = function(parent, obj) {
      if (!parent || !obj) return false;

      let isnode = (obj._typename.indexOf('TGeoNode') === 0),
          isvolume = (obj._typename.indexOf('TGeoVolume') === 0),
          ismanager = (obj._typename === 'TGeoManager'),
          iseve = ((obj._typename === 'TEveGeoShapeExtract') || (obj._typename === 'ROOT::Experimental::REveGeoShapeExtract')),
          isoverlap = (obj._typename === 'TGeoOverlap');

      if (!isnode && !isvolume && !ismanager && !iseve && !isoverlap) return false;

      if (parent._childs) return true;

      if (ismanager) {
         geo.createList(parent, obj.fMaterials, "Materials", "list of materials");
         geo.createList(parent, obj.fMedia, "Media", "list of media");
         geo.createList(parent, obj.fTracks, "Tracks", "list of tracks");
         geo.createList(parent, obj.fOverlaps, "Overlaps", "list of detected overlaps");
         geo.createItem(parent, obj.fMasterVolume);
         return true;
      }

      if (isoverlap) {
         geo.createItem(parent, obj.fVolume1);
         geo.createItem(parent, obj.fVolume2);
         geo.createItem(parent, obj.fMarker, 'Marker');
         return true;
      }

      let volume, subnodes, shape;

      if (iseve) {
         subnodes = obj.fElements ? obj.fElements.arr : null;
         shape = obj.fShape;
      } else {
         volume = (isnode ? obj.fVolume : obj);
         subnodes = volume && volume.fNodes ? volume.fNodes.arr : null;
         shape = volume ? volume.fShape : null;
      }

      if (!subnodes && shape && (shape._typename === "TGeoCompositeShape") && shape.fNode) {
         if (!parent._childs) {
            geo.createItem(parent, shape.fNode.fLeft, 'Left');
            geo.createItem(parent, shape.fNode.fRight, 'Right');
         }

         return true;
      }

      if (!subnodes) return false;

      geo.checkDuplicates(obj, subnodes);

      for (let i = 0; i < subnodes.length; ++i)
         geo.createItem(parent, subnodes[i]);

      return true;
   }

   jsrp.addDrawFunc({ name: "TGeoVolumeAssembly", icon: 'img_geoassembly', func: drawGeoObject, expand: geo.expandObject, opt: ";more;all;count" });
   jsrp.addDrawFunc({ name: "TEvePointSet", icon_get: geo.getBrowserIcon, icon_click: geo.browserIconClick });
   jsrp.addDrawFunc({ name: "TEveTrack", icon_get: geo.getBrowserIcon, icon_click: geo.browserIconClick });

   JSROOT.TGeoPainter = TGeoPainter;

   jsrp.GeoDrawingControl = GeoDrawingControl;
   jsrp.drawGeoObject = drawGeoObject;

   return geo;

});
