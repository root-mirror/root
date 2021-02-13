/// @file JSRoot.base3d.js
/// JavaScript ROOT 3D graphics

JSROOT.define(['d3', 'threejs_jsroot', 'painter'], (d3, THREE, jsrp) => {

   "use strict";

   /** @ummary Define rendering kind which will be used for rendering of 3D elements
    * @memberOf JSROOT.Painter
    * @param {value} [render3d] - preconfigured value, will be used if applicable
    * @returns {value} - rendering kind, see JSROOT.constants.Render3D
    * @private */
   jsrp.getRender3DKind = function(render3d) {
      if (!render3d) render3d = JSROOT.batch_mode ? JSROOT.settings.Render3DBatch : JSROOT.settings.Render3D;
      let rc = JSROOT.constants.Render3D;

      if (render3d == rc.Default) render3d = JSROOT.batch_mode ? rc.WebGLImage : rc.WebGL;
      if (JSROOT.batch_mode && (render3d == rc.WebGL)) render3d = rc.WebGLImage;

      return render3d;
   }

   let Handling3DDrawings = {};

   /** @summary Access current 3d mode
     * @param {string} [new_value] - when specified, set new 3d mode
     * @returns current value
     * @private*/
   Handling3DDrawings.access3dKind = function(new_value) {
      let svg = this.getPadSvg();
      if (svg.empty()) return -1;

      // returns kind of currently created 3d canvas
      let kind = svg.property('can3d');
      if (new_value !== undefined) svg.property('can3d', new_value);
      return ((kind === null) || (kind === undefined)) ? -1 : kind;
   }

   /** @summary Returns size which availble for 3D drawing.
     * @desc One uses frame sizes for the 3D drawing - like TH2/TH3 objects
     * @private */
   Handling3DDrawings.getSizeFor3d = function(can3d, render3d) {

      if (can3d === undefined) {
         // analyze which render/embed mode can be used
         can3d = jsrp.getRender3DKind();
         // all non-webgl elements can be embedded into SVG as is
         if (can3d !== JSROOT.constants.Render3D.WebGL)
            can3d = JSROOT.constants.Embed3D.EmbedSVG;
         else if (JSROOT.settings.Embed3D != JSROOT.constants.Embed3D.Default)
            can3d = JSROOT.settings.Embed3D;
         else if (JSROOT.browser.isFirefox)
            can3d = JSROOT.constants.Embed3D.Embed;
         else
            can3d = JSROOT.constants.Embed3D.Overlay;
      }

      let pad = this.getPadSvg(),
          clname = "draw3d_" + (this.getPadName() || 'canvas');

      if (pad.empty()) {
         // this is a case when object drawn without canvas

         let rect = jsrp.getElementRect(this.selectDom());

         if ((rect.height < 10) && (rect.width > 10)) {
            rect.height = Math.round(0.66 * rect.width);
            this.selectDom().style('height', rect.height + "px");
         }
         rect.x = 0; rect.y = 0; rect.clname = clname; rect.can3d = -1;
         return rect;
      }

      let elem = pad, fp = this.getFramePainter();
      if (can3d === 0) elem = this.getCanvSvg();

      let size = { x: 0, y: 0, width: 100, height: 100, clname: clname, can3d: can3d };

      if (fp && !fp.mode3d) {
         elem = this.getFrameSvg();
         size.x = elem.property("draw_x");
         size.y = elem.property("draw_y");
      }

      size.width = elem.property("draw_width");
      size.height = elem.property("draw_height");

      if ((!fp || fp.mode3d) && (can3d > 0)) {
         size.x = Math.round(size.x + size.width * JSROOT.gStyle.fPadLeftMargin);
         size.y = Math.round(size.y + size.height * JSROOT.gStyle.fPadTopMargin);
         size.width = Math.round(size.width * (1 - JSROOT.gStyle.fPadLeftMargin - JSROOT.gStyle.fPadRightMargin));
         size.height = Math.round(size.height * (1 - JSROOT.gStyle.fPadTopMargin - JSROOT.gStyle.fPadBottomMargin));
      }

      let pp = this.getPadPainter(),
          rect = pp ? pp.getPadRect() : { width: 100, height: 100 },
          x2 = rect.width - size.x - size.width,
          y2 = rect.height - size.y - size.height;

      if ((x2 >= 0) && (y2 >= 0)) {
         // while 3D canvas uses area also for the axis labels, extend area relative to normal frame
         size.x = Math.round(size.x * 0.3);
         size.y = Math.round(size.y * 0.9);
         size.width = rect.width - size.x - Math.round(x2 * 0.3);
         size.height = rect.height - size.y - Math.round(y2 * 0.5);
      }

      if (can3d === 1)
         size = jsrp.getAbsPosInCanvas(this.getPadSvg(), size);

      return size;
   }

   /** @summary Clear all 3D drawings
     * @returns can3d value - how webgl canvas was placed
     * @private */
   Handling3DDrawings.clear3dCanvas = function() {
      let can3d = this.access3dKind(null);
      if (can3d < 0) {
         // remove first child from main element - if it is canvas
         let main = this.selectDom().node();
         if (main && main.firstChild && main.firstChild.$jsroot) {
            delete main.firstChild.painter;
            main.removeChild(main.firstChild);
         }
         return can3d;
      }

      let size = this.getSizeFor3d(can3d);

      if (size.can3d === 0) {
         d3.select(this.getCanvSvg().node().nextSibling).remove(); // remove html5 canvas
         this.getCanvSvg().style('display', null); // show SVG canvas
      } else {
         if (this.getPadSvg().empty()) return;

         this.apply3dSize(size).remove();

         this.getFrameSvg().style('display', null);  // clear display property
      }
      return can3d;
   }

   /** @summary Add 3D canvas
    * @private */
   Handling3DDrawings.add3dCanvas = function(size, canv, webgl) {

      if (!canv || (size.can3d < -1)) return;

      if (size.can3d === -1) {
         // case when 3D object drawn without canvas

         let main = this.selectDom().node();
         if (main !== null) {
            main.appendChild(canv);
            canv.painter = this;
            canv.$jsroot = true; // mark canvas as added by jsroot
         }

         return;
      }

      if ((size.can3d > 0) && !webgl)
         size.can3d = JSROOT.constants.Embed3D.EmbedSVG;

      this.access3dKind(size.can3d);

      if (size.can3d === 0) {
         this.getCanvSvg().style('display', 'none'); // hide SVG canvas

         this.getCanvSvg().node().parentNode.appendChild(canv); // add directly
      } else {
         if (this.getPadSvg().empty()) return;

         // first hide normal frame
         this.getFrameSvg().style('display', 'none');

         let elem = this.apply3dSize(size);

         elem.attr('title', '').node().appendChild(canv);
      }
   }

   /** @summary Apply size to 3D elements
    * @private */
   Handling3DDrawings.apply3dSize = function(size, onlyget) {

      if (size.can3d < 0) return d3.select(null);

      let elem;

      if (size.can3d > 1) {

         elem = this.getLayerSvg(size.clname);

         // elem = layer.select("." + size.clname);
         if (onlyget) return elem;

         let svg = this.getPadSvg();

         if (size.can3d === JSROOT.constants.Embed3D.EmbedSVG) {
            // this is SVG mode or image mode - just create group to hold element

            if (elem.empty())
               elem = svg.insert("g", ".primitives_layer").attr("class", size.clname);

            elem.attr("transform", "translate(" + size.x + "," + size.y + ")");

         } else {

            if (elem.empty())
               elem = svg.insert("foreignObject", ".primitives_layer").attr("class", size.clname);

            elem.attr('x', size.x)
               .attr('y', size.y)
               .attr('width', size.width)
               .attr('height', size.height)
               .attr('viewBox', "0 0 " + size.width + " " + size.height)
               .attr('preserveAspectRatio', 'xMidYMid');
         }

      } else {
         let prnt = this.getCanvSvg().node().parentNode;

         elem = d3.select(prnt).select("." + size.clname);
         if (onlyget) return elem;

         // force redraw by resize
         this.getCanvSvg().property('redraw_by_resize', true);

         if (elem.empty())
            elem = d3.select(prnt).append('div').attr("class", size.clname + " jsroot_noselect");

         // our position inside canvas, but to set 'absolute' position we should use
         // canvas element offset relative to first parent with non-static position
         // now try to use getBoundingClientRect - it should be more precise

         let pos0 = prnt.getBoundingClientRect();

         while (prnt) {
            if (prnt === document) { prnt = null; break; }
            try {
               if (getComputedStyle(prnt).position !== 'static') break;
            } catch (err) {
               break;
            }
            prnt = prnt.parentNode;
         }

         let pos1 = prnt ? prnt.getBoundingClientRect() : { top: 0, left: 0 };

         let offx = Math.round(pos0.left - pos1.left),
             offy = Math.round(pos0.top - pos1.top);

         elem.style('position', 'absolute').style('left', (size.x + offx) + 'px').style('top', (size.y + offy) + 'px').style('width', size.width + 'px').style('height', size.height + 'px');
      }

      return elem;
   }

   /** @summary Assigns method to handle 3D drawings inside SVG
     * @private */
   jsrp.assign3DHandler = function(painter) {
      JSROOT.extend(painter, Handling3DDrawings);
   }

   /** @summary Creates renderer for the 3D drawings
     * @memberOf JSROOT.Painter
     * @param {value} width - rendering width
     * @param {value} height - rendering height
     * @param {value} render3d - render type, see {@link JSROOT.constants.Render3D}
     * @param {object} args - different arguments for creating 3D renderer
     * @private */
   jsrp.createRender3D = function(width, height, render3d, args) {

      let rc = JSROOT.constants.Render3D;

      render3d = jsrp.getRender3DKind(render3d);

      if (!args) args = { antialias: true, alpha: true };

      let need_workaround = false, renderer,
          doc = JSROOT._.get_document();

      if (render3d == rc.WebGL) {
         // interactive WebGL Rendering
         renderer = new THREE.WebGLRenderer(args);

      } else if (render3d == rc.SVG) {
         // SVG rendering
         renderer = THREE.CreateSVGRenderer(false, 0, doc);

         if (JSROOT.batch_mode) {
            need_workaround = true;
         } else {
            renderer.jsroot_dom = doc.createElementNS('http://www.w3.org/2000/svg', 'svg');
            // d3.select(renderer.jsroot_dom).attr("width", width).attr("height", height);
         }
      } else if (JSROOT.nodejs) {
         // try to use WebGL inside node.js - need to create headless context
         let gl = require('gl')(width, height, { preserveDrawingBuffer: true });

         const { createCanvas } = require('canvas');

         args.canvas = createCanvas(width, height);
         args.canvas.addEventListener = function() { }; // dummy
         args.canvas.style = {};

         args.context = gl;

         // in node.js too many similar warnings makes it difficult to debug
         let warn = console.warn;
         console.warn = function(msg) {
            if (msg && msg.indexOf("OES_texture_half_float") < 0 && msg.indexOf("EXT_texture_filter_anisotropic") < 0 &&
                       msg.indexOf("WEBGL_depth_texture") < 0 && msg.indexOf("OES_vertex_array_object") < 0)
              console.log("NEW: " + msg);
         }
         renderer = new THREE.WebGLRenderer(args);

         renderer.jsroot_output = new THREE.WebGLRenderTarget(width, height);

         renderer.setRenderTarget(renderer.jsroot_output);

         console.warn = warn;

         need_workaround = true;
      } else {
         // rendering with WebGL directly into svg image
         renderer = new THREE.WebGLRenderer(args);
         renderer.jsroot_dom = doc.createElementNS('http://www.w3.org/2000/svg', 'image');
         d3.select(renderer.jsroot_dom).attr("width", width).attr("height", height);
      }

      if (need_workaround) {
         if (!JSROOT._.svg_3ds) JSROOT._.svg_3ds = [];
         renderer.workaround_id = JSROOT._.svg_3ds.length;
         JSROOT._.svg_3ds[renderer.workaround_id] = "<svg></svg>"; // dummy, provided in afterRender3D

         // replace DOM element in renderer
         renderer.jsroot_dom = doc.createElementNS('http://www.w3.org/2000/svg', 'path');
         renderer.jsroot_dom.setAttribute('jsroot_svg_workaround', renderer.workaround_id);
      } else if (!renderer.jsroot_dom) {
         renderer.jsroot_dom = renderer.domElement;
      }

      // res.renderer.setClearColor("#000000", 1);
      // res.renderer.setClearColor(0x0, 0);
      renderer.setSize(width, height);
      renderer.jsroot_render3d = render3d;

      // apply size to dom element
      renderer.setJSROOTSize = function(width, height) {
         if ((this.jsroot_render3d === JSROOT.constants.Render3D.WebGLImage) && !JSROOT.batch_mode && !JSROOT.nodejs)
            return d3.select(this.jsroot_dom).attr("width", width).attr("height", height);
      }

      return renderer;
   }

   /** @summary Cleanup previous renderings before doing next one
     * @desc used together with SVG
     * @private */
   jsrp.beforeRender3D = function(renderer) {
      // cleanup previous rendering, from SVG renderer
      if (renderer.clearHTML) renderer.clearHTML();
   }

   /** @summary Post-process result of rendering
     * @desc used together with SVG or node.js image rendering
     * @private */
   jsrp.afterRender3D = function(renderer) {

      let rc = JSROOT.constants.Render3D;
      if (renderer.jsroot_render3d == rc.WebGL) return;

      if (renderer.jsroot_render3d == rc.SVG) {
         // case of SVGRenderer
         if (JSROOT.batch_mode) {
            JSROOT._.svg_3ds[renderer.workaround_id] = renderer.makeOuterHTML();
         } else {
            let parent = renderer.jsroot_dom.parentNode;
            if (parent) {
               parent.innerHTML = renderer.makeOuterHTML();
               renderer.jsroot_dom = parent.firstChild;
            }
         }
      } else if (JSROOT.nodejs) {
         // this is WebGL rendering in node.js
         let canvas = renderer.domElement,
            context = canvas.getContext('2d');

         let pixels = new Uint8Array(4 * canvas.width * canvas.height);
         renderer.readRenderTargetPixels(renderer.jsroot_output, 0, 0, canvas.width, canvas.height, pixels);

         // small code to flip Y scale
         let indx1 = 0, indx2 = (canvas.height - 1) * 4 * canvas.width, k, d;
         while (indx1 < indx2) {
            for  (k = 0; k < 4 * canvas.width; ++k) {
               d = pixels[indx1 + k]; pixels[indx1 + k] = pixels[indx2 + k]; pixels[indx2 + k] = d;
            }
            indx1 += 4 * canvas.width;
            indx2 -= 4 * canvas.width;
         }

         let imageData = context.createImageData( canvas.width, canvas.height );
         imageData.data.set( pixels );
         context.putImageData( imageData, 0, 0 );

         let dataUrl = canvas.toDataURL("image/png"),
             svg = '<image width="' + canvas.width + '" height="' + canvas.height + '" xlink:href="' + dataUrl + '"></image>';
         JSROOT._.svg_3ds[renderer.workaround_id] = svg;
      } else {
         let dataUrl = renderer.domElement.toDataURL("image/png");
         d3.select(renderer.jsroot_dom).attr("xlink:href", dataUrl);
      }
   }

   /** @summary Special way to insert WebGL drawing into produced SVG batch code
     * @desc Used only in batch mode for SVG images generation
     * @private */
   jsrp.processSvgWorkarounds = function(svg, keep_workarounds) {
      if (!JSROOT._.svg_3ds) return svg;
      for (let k = 0;  k < JSROOT._.svg_3ds.length; ++k)
         svg = svg.replace('<path jsroot_svg_workaround="' + k + '"></path>', JSROOT._.svg_3ds[k]);
      if (!keep_workarounds)
         JSROOT._.svg_3ds = undefined;
      return svg;
   }

   // ========================================================================================================

   /**
    * @summary Tooltip handler for 3D drawings
    *
    * @class
    * @memberof JSROOT.Painter
    * @param {object} dom - DOM element
    * @param {object} canvas - canvas for 3D rendering
    * @private
    */

   function TooltipFor3D(prnt, canvas) {
      this.tt = null;
      this.cont = null;
      this.lastlbl = '';
      this.parent = prnt ? prnt : document.body;
      this.canvas = canvas; // we need canvas to recalculate mouse events
      this.abspos = !prnt;
   }

   /** @summary check parent */
   TooltipFor3D.prototype.checkParent = function(prnt) {
      if (prnt && (this.parent !== prnt)) {
         this.hide();
         this.parent = prnt;
      }
   }

   /** @summary extract position from event
     * @desc can be used to process it later when event is gone */
   TooltipFor3D.prototype.extract_pos = function(e) {
      if (typeof e == 'object' && (e.u !== undefined) && (e.l !== undefined)) return e;
      let res = { u: 0, l: 0 };
      if (this.abspos) {
         res.l = e.pageX;
         res.u = e.pageY;
      } else {
         res.l = e.offsetX;
         res.u = e.offsetY;
      }

      return res;
   }

   /** @summary Method used to define position of next tooltip
     * @desc event is delivered from canvas,
     * but position should be calculated relative to the element where tooltip is placed */
   TooltipFor3D.prototype.pos = function(e) {

      if (!this.tt) return;

      let pos = this.extract_pos(e);
      if (!this.abspos) {
         let rect1 = this.parent.getBoundingClientRect(),
             rect2 = this.canvas.getBoundingClientRect();

         if ((rect1.left !== undefined) && (rect2.left!== undefined)) pos.l += (rect2.left-rect1.left);

         if ((rect1.top !== undefined) && (rect2.top!== undefined)) pos.u += rect2.top-rect1.top;

         if (pos.l + this.tt.offsetWidth + 3 >= this.parent.offsetWidth)
            pos.l = this.parent.offsetWidth - this.tt.offsetWidth - 3;

         if (pos.u + this.tt.offsetHeight + 15 >= this.parent.offsetHeight)
            pos.u = this.parent.offsetHeight - this.tt.offsetHeight - 15;

         // one should find parent with non-static position,
         // all absolute coordinates calculated relative to such node
         let abs_parent = this.parent;
         while (abs_parent) {
            let style = getComputedStyle(abs_parent);
            if (!style || (style.position !== 'static')) break;
            if (!abs_parent.parentNode || (abs_parent.parentNode.nodeType != 1)) break;
            abs_parent = abs_parent.parentNode;
         }

         if (abs_parent && (abs_parent !== this.parent)) {
            let rect0 = abs_parent.getBoundingClientRect();
            pos.l += (rect1.left - rect0.left);
            pos.u += (rect1.top - rect0.top);
         }
      }

      this.tt.style.top = (pos.u + 15) + 'px';
      this.tt.style.left = (pos.l + 3) + 'px';
   }

   /** @summary Show tooltip */
   TooltipFor3D.prototype.show = function(v /*, mouse_pos, status_func*/) {
      if (!v || (v==="")) return this.hide();

      if (v && (typeof v =='object') && (v.lines || v.line)) {
         if (v.only_status) return this.hide();

         if (v.line) {
            v = v.line;
         } else {
            let res = v.lines[0];
            for (let n=1;n<v.lines.length;++n) res+= "<br/>" + v.lines[n];
            v = res;
         }
      }

      if (this.tt === null) {
         this.tt = document.createElement('div');
         this.tt.setAttribute('class', 'jsroot_tt3d_main');
         this.cont = document.createElement('div');
         this.cont.setAttribute('class', 'jsroot_tt3d_cont');
         this.tt.appendChild(this.cont);
         this.parent.appendChild(this.tt);
      }

      if (this.lastlbl !== v) {
         this.cont.innerHTML = v;
         this.lastlbl = v;
         this.tt.style.width = 'auto'; // let it be automatically resizing...
      }
   }

   /** @summary Hide tooltip */
   TooltipFor3D.prototype.hide = function() {
      if (this.tt !== null)
         this.parent.removeChild(this.tt);

      this.tt = null;
      this.lastlbl = "";
   }


   /** @summary Create THREE.OrbitControl for painter
     * @private */
   jsrp.createOrbitControl = function(painter, camera, scene, renderer, lookat) {

      if (JSROOT.settings.Zooming && JSROOT.settings.ZoomWheel)
         renderer.domElement.addEventListener( 'wheel', control_mousewheel);

      let enable_zoom = JSROOT.settings.Zooming && JSROOT.settings.ZoomMouse,
          enable_select = (typeof painter.processMouseClick == "function"),
          control = null;

      function control_mousedown(evnt) {

         // function used to hide some events from orbit control and redirect them to zooming rect
         if (control.mouse_zoom_mesh) {
            evnt.stopImmediatePropagation();
            evnt.stopPropagation();
            return;
         }

         // only left-button is considered
         if ((evnt.button!==undefined) && (evnt.button !==0)) return;
         if ((evnt.buttons!==undefined) && (evnt.buttons !== 1)) return;

         if (control.enable_zoom) {
            control.mouse_zoom_mesh = control.detectZoomMesh(evnt);
            if (control.mouse_zoom_mesh) {
               // just block orbit control
               evnt.stopImmediatePropagation();
               evnt.stopPropagation();
               return;
            }
         }

         if (control.enable_select)
            control.mouse_select_pnt = control.getMousePos(evnt, {});
      }

      function control_mouseup(evnt) {

         if (control.mouse_zoom_mesh && control.mouse_zoom_mesh.point2 && control.painter.Get3DZoomCoord) {

            let kind = control.mouse_zoom_mesh.object.zoom,
                pos1 = control.painter.Get3DZoomCoord(control.mouse_zoom_mesh.point, kind),
                pos2 = control.painter.Get3DZoomCoord(control.mouse_zoom_mesh.point2, kind);

            if (pos1>pos2) { let v = pos1; pos1 = pos2; pos2 = v; }

            if ((kind==="z") && control.mouse_zoom_mesh.object.use_y_for_z) kind="y";

            if ((kind==="z") && control.mouse_zoom_mesh.object.use_y_for_z) kind="y";

            // try to zoom
            if (pos1 < pos2)
              if (control.painter.zoom(kind, pos1, pos2))
                 control.mouse_zoom_mesh = null;
         }

         // if selection was drawn, it should be removed and picture rendered again
         if (control.enable_zoom)
            control.removeZoomMesh();

         // only left-button is considered
         //if ((evnt.button!==undefined) && (evnt.button !==0)) return;
         //if ((evnt.buttons!==undefined) && (evnt.buttons !== 1)) return;

         if (control.enable_select && control.mouse_select_pnt) {

            let pnt = control.getMousePos(evnt, {});

            let same_pnt = (pnt.x == control.mouse_select_pnt.x) && (pnt.y == control.mouse_select_pnt.y);
            delete control.mouse_select_pnt;

            if (same_pnt) {
               let intersects = control.getMouseIntersects(pnt);
               control.painter.processMouseClick(pnt, intersects, evnt);
            }
         }
      }

      if (enable_zoom || enable_select) {
         renderer.domElement.addEventListener( 'pointerdown', control_mousedown);
         renderer.domElement.addEventListener( 'pointerup', control_mouseup);
      }

      control = new THREE.OrbitControls(camera, renderer.domElement);

      control.enableDamping = false;
      control.dampingFactor = 1.0;
      control.enableZoom = true;
      control.enableKeys = JSROOT.settings.HandleKeys;

      if (lookat) {
         control.target.copy(lookat);
         control.target0.copy(lookat);
         control.update();
      }

      control.tooltip = new TooltipFor3D(painter.selectDom().node(), renderer.domElement);

      control.painter = painter;
      control.camera = camera;
      control.scene = scene;
      control.renderer = renderer;
      control.raycaster = new THREE.Raycaster();
      control.raycaster.params.Line.threshold = 10;
      control.raycaster.params.Points.threshold = 5;
      control.mouse_zoom_mesh = null; // zoom mesh, currently used in the zooming
      control.block_ctxt = false; // require to block context menu command appearing after control ends, required in chrome which inject contextmenu when key released
      control.block_mousemove = false; // when true, tooltip or cursor will not react on mouse move
      control.cursor_changed = false;
      control.control_changed = false;
      control.control_active = false;
      control.mouse_ctxt = { x:0, y: 0, on: false };
      control.enable_zoom = enable_zoom;
      control.enable_select = enable_select;

      control.cleanup = function() {
         if (JSROOT.settings.Zooming && JSROOT.settings.ZoomWheel)
            this.domElement.removeEventListener( 'wheel', control_mousewheel);
         if (this.enable_zoom || this.enable_select) {
            this.domElement.removeEventListener( 'mousedown', control_mousedown);
            this.domElement.removeEventListener( 'mouseup', control_mouseup);
         }

         if (this.lstn_click)
            this.domElement.removeEventListener('click', this.lstn_click);
         this.domElement.removeEventListener('dblclick', this.lstn_dblclick);
         this.domElement.removeEventListener('contextmenu', this.lstn_contextmenu);
         this.domElement.removeEventListener('mousemove', this.lstn_mousemove);
         this.domElement.removeEventListener('mouseleave', this.lstn_mouseleave);

         this.dispose(); // this is from OrbitControl itself

         this.tooltip.hide();
         delete this.tooltip;
         delete this.painter;
         delete this.camera;
         delete this.scene;
         delete this.renderer;
         delete this.raycaster;
         delete this.mouse_zoom_mesh;
      }

      control.HideTooltip = function() {
         this.tooltip.hide();
      }

      control.getMousePos = function(evnt, mouse) {
         mouse.x = ('offsetX' in evnt) ? evnt.offsetX : evnt.layerX;
         mouse.y = ('offsetY' in evnt) ? evnt.offsetY : evnt.layerY;
         mouse.clientX = evnt.clientX;
         mouse.clientY = evnt.clientY;
         return mouse;
      }

      control.getOriginDirectionIntersects = function(origin, direction) {
         this.raycaster.set(origin, direction);
         let intersects = this.raycaster.intersectObjects(this.scene.children, true);
         // painter may want to filter intersects
         if (typeof this.painter.filterIntersects == 'function')
            intersects = this.painter.filterIntersects(intersects);
         return intersects;
      }

      control.getMouseIntersects = function(mouse) {
         // domElement gives correct coordinate with canvas render, but isn't always right for webgl renderer
         if (!this.renderer) return [];

         let sz = (this.renderer instanceof THREE.WebGLRenderer) ?
                     this.renderer.getSize(new THREE.Vector2()) :
                     this.renderer.domElement;

         let pnt = { x: mouse.x / sz.width * 2 - 1, y: -mouse.y / sz.height * 2 + 1 };

         this.camera.updateMatrix();
         this.camera.updateMatrixWorld();
         this.raycaster.setFromCamera( pnt, this.camera );
         let intersects = this.raycaster.intersectObjects(this.scene.children, true);

         // painter may want to filter intersects
         if (typeof this.painter.filterIntersects == 'function')
            intersects = this.painter.filterIntersects(intersects);

         return intersects;
      }

      control.detectZoomMesh = function(evnt) {
         let mouse = this.getMousePos(evnt, {}),
             intersects = this.getMouseIntersects(mouse);
         if (intersects)
            for (let n = 0; n < intersects.length; ++n)
               if (intersects[n].object.zoom)
                  return intersects[n];

         return null;
      }

      control.processDblClick = function(evnt) {
         let intersect = this.detectZoomMesh(evnt);
         if (intersect && this.painter) {
            this.painter.unzoom(intersect.object.use_y_for_z ? "y" : intersect.object.zoom);
         } else {
            this.reset();
         }
         // this.painter.render3D();
      }

      control.changeEvent = function() {
         this.mouse_ctxt.on = false; // disable context menu if any changes where done by orbit control
         this.painter.render3D(0);
         this.control_changed = true;
      }

      control.startEvent = function() {
         this.control_active = true;
         this.block_ctxt = false;
         this.mouse_ctxt.on = false;

         this.tooltip.hide();

         // do not reset here, problem of events sequence in orbitcontrol
         // it issue change/start/stop event when do zooming
         // control.control_changed = false;
      }

      control.endEvent = function() {
         this.control_active = false;
         if (this.mouse_ctxt.on) {
            this.mouse_ctxt.on = false;
            this.contextMenu(this.mouse_ctxt, this.getMouseIntersects(this.mouse_ctxt));
         } /* else if (this.control_changed) {
            // react on camera change when required
         } */
         this.control_changed = false;
      }

      control.mainProcessContextMenu = function(evnt) {
         evnt.preventDefault();
         this.getMousePos(evnt, this.mouse_ctxt);
         if (this.control_active)
            this.mouse_ctxt.on = true;
         else if (this.block_ctxt)
            this.block_ctxt = false;
         else
            this.contextMenu(this.mouse_ctxt, this.getMouseIntersects(this.mouse_ctxt));
      }

      control.contextMenu = function(/* pos, intersects */) {
         // do nothing, function called when context menu want to be activated
      }

      control.setTooltipEnabled = function(on) {
         this.block_mousemove = !on;
         if (on === false) {
            this.tooltip.hide();
            this.removeZoomMesh();
         }
      }

      control.removeZoomMesh = function() {
         if (this.mouse_zoom_mesh && this.mouse_zoom_mesh.object.ShowSelection())
            this.painter.render3D();
         this.mouse_zoom_mesh = null; // in any case clear mesh, enable orbit control again
      }

      control.mainProcessMouseMove = function(evnt) {
         if (!this.painter) return; // protect when cleanup

         if (this.control_active && evnt.buttons && (evnt.buttons & 2))
            this.block_ctxt = true; // if right button in control was active, block next context menu

         if (this.control_active || this.block_mousemove || !this.ProcessMouseMove) return;

         if (this.mouse_zoom_mesh) {
            // when working with zoom mesh, need special handling

            let zoom2 = this.detectZoomMesh(evnt), pnt2 = null;

            if (zoom2 && (zoom2.object === this.mouse_zoom_mesh.object)) {
               pnt2 = zoom2.point;
            } else {
               pnt2 = this.mouse_zoom_mesh.object.GlobalIntersect(this.raycaster);
            }

            if (pnt2) this.mouse_zoom_mesh.point2 = pnt2;

            if (pnt2 && this.painter.enable_highlight)
               if (this.mouse_zoom_mesh.object.ShowSelection(this.mouse_zoom_mesh.point, pnt2))
                  this.painter.render3D(0);

            this.tooltip.hide();
            return;
         }

         evnt.preventDefault();

         // extract mouse position
         this.tmout_mouse = this.getMousePos(evnt, {});
         this.tmout_ttpos =  this.tooltip ? this.tooltip.extract_pos(evnt) : null;

         if (this.tmout_handle) {
            clearTimeout(this.tmout_handle);
            delete this.tmout_handle;
         }

         if (!this.mouse_tmout)
            this.delayedProcessMouseMove();
         else
            this.tmout_handle = setTimeout(() => this.delayedProcessMouseMove(), this.mouse_tmout);
      }

      control.delayedProcessMouseMove = function() {
         // remove handle - allow to trigger new timeout
         delete this.tmout_handle;
         if (!this.painter) return; // protect when cleanup

         let mouse = this.tmout_mouse,
             intersects = this.getMouseIntersects(mouse),
             tip = this.ProcessMouseMove(intersects);

         if (tip) {
            let name = "", title = "", coord = "", info = "";
            if (mouse) coord = mouse.x.toFixed(0)+ "," + mouse.y.toFixed(0);
            if (typeof tip == "string") {
               info = tip;
            } else {
               name = tip.name; title = tip.title;
               if (tip.line) info = tip.line; else
               if (tip.lines) { info = tip.lines.slice(1).join(' '); name = tip.lines[0]; }
            }
            this.painter.showObjectStatus(name, title, info, coord);
         }

         this.cursor_changed = false;
         if (tip && this.painter && this.painter.isTooltipAllowed()) {
            this.tooltip.checkParent(this.painter.selectDom().node());

            this.tooltip.show(tip, mouse);
            this.tooltip.pos(this.tmout_ttpos);
         } else {
            this.tooltip.hide();
            if (intersects)
               for (let n=0;n<intersects.length;++n)
                  if (intersects[n].object.zoom) this.cursor_changed = true;
         }

         document.body.style.cursor = this.cursor_changed ? 'pointer' : 'auto';
      }

      control.mainProcessMouseLeave = function() {
         if (!this.painter) return; // protect when cleanup

         // do not enter main event at all
         if (this.tmout_handle) {
            clearTimeout(this.tmout_handle);
            delete this.tmout_handle;
         }
         this.tooltip.hide();
         if (typeof this.ProcessMouseLeave === 'function')
            this.ProcessMouseLeave();
         if (this.cursor_changed) {
            document.body.style.cursor = 'auto';
            this.cursor_changed = false;
         }
      }

      function render3DFired(painter) {
         if (!painter || painter.renderer === undefined) return false;
         return painter.render_tmout !== undefined; // when timeout configured, object is prepared for rendering
      }

      function control_mousewheel(evnt) {
         // try to handle zoom extra

         if (render3DFired(control.painter) || control.mouse_zoom_mesh) {
            evnt.preventDefault();
            evnt.stopPropagation();
            evnt.stopImmediatePropagation();
            return; // already fired redraw, do not react on the mouse wheel
         }

         let intersect = control.detectZoomMesh(evnt);
         if (!intersect) return;

         evnt.preventDefault();
         evnt.stopPropagation();
         evnt.stopImmediatePropagation();

         if (control.painter && (typeof control.painter.analyzeMouseWheelEvent == 'function')) {
            let kind = intersect.object.zoom,
                position = intersect.point[kind],
                item = { name: kind, ignore: false };

            // z changes from 0..2*size_z3d, others -size_xy3d..+size_xy3d
            if (kind!=="z") position = (position + control.painter.size_xy3d)/2/control.painter.size_xy3d;
                       else position = position/2/control.painter.size_z3d;

            control.painter.analyzeMouseWheelEvent(evnt, item, position, false);

            if ((kind==="z") && intersect.object.use_y_for_z) kind = "y";

            control.painter.zoom(kind, item.min, item.max);
         }
      }

      control.mainProcessDblClick = function(evnt) {
         this.processDblClick(evnt);
      }

      if (painter && painter.options && painter.options.mouse_click) {
         control.processClick = function(mouse) {
            if (typeof this.ProcessSingleClick == 'function') {
               let intersects = this.getMouseIntersects(mouse);
               this.ProcessSingleClick(intersects);
            }
         }

         control.lstn_click = function(evnt) {
            if (this.single_click_tm) {
               clearTimeout(this.single_click_tm);
               delete this.single_click_tm;
            }

            // if normal event, set longer timeout waiting if double click not detected
            if (evnt.detail != 2)
               this.single_click_tmout = setTimeout(this.ProcessClick.bind(this, this.getMousePos(evnt, {})), 300);
         }.bind(control);
      }

      control.addEventListener('change', () => control.changeEvent());
      control.addEventListener('start', () => control.startEvent());
      control.addEventListener('end', () => control.endEvent());

      control.lstn_contextmenu = evnt => control.mainProcessContextMenu(evnt);
      control.lstn_dblclick = evnt => control.mainProcessDblClick(evnt);
      control.lstn_mousemove = evnt => control.mainProcessMouseMove(evnt);
      control.lstn_mouseleave = () => control.mainProcessMouseLeave();

      if (control.lstn_click)
         renderer.domElement.addEventListener('click', control.lstn_click);
      renderer.domElement.addEventListener('dblclick', control.lstn_dblclick);
      renderer.domElement.addEventListener('contextmenu', control.lstn_contextmenu);
      renderer.domElement.addEventListener('mousemove', control.lstn_mousemove);
      renderer.domElement.addEventListener('mouseleave', control.lstn_mouseleave);

      return control;
   }

   /** @summary Method cleanup three.js object as much as possible.
     * @desc Simplify JS engine to remove it from memory
     * @private */
   jsrp.disposeThreejsObject = function(obj, only_childs) {
      if (!obj) return;

      if (obj.children) {
         for (let i = 0; i < obj.children.length; i++)
            jsrp.disposeThreejsObject(obj.children[i]);
      }

      if (only_childs) {
         obj.children = [];
         return;
      }

      obj.children = undefined;

      if (obj.geometry) {
         obj.geometry.dispose();
         obj.geometry = undefined;
      }
      if (obj.material) {
         if (obj.material.map) {
            obj.material.map.dispose();
            obj.material.map = undefined;
         }
         obj.material.dispose();
         obj.material = undefined;
      }

      // cleanup jsroot fields to simplify browser cleanup job
      delete obj.painter;
      delete obj.bins_index;
      delete obj.tooltip;
      delete obj.stack; // used in geom painter
      delete obj.drawn_highlight; // special highlight object

      obj = undefined;
   }


   /** @summary Create THREE.LineSegments mesh (or only geometry)
     * @desc If required, calculates lineDistance attribute for dashed geometries
     * @private */
   jsrp.createLineSegments = function(arr, material, index, only_geometry) {

      let geom = new THREE.BufferGeometry();

      geom.setAttribute( 'position', arr instanceof Float32Array ? new THREE.BufferAttribute( arr, 3 ) : new THREE.Float32BufferAttribute( arr, 3 ) );
      if (index) geom.setIndex(  new THREE.BufferAttribute(index, 1) );

      if (material.isLineDashedMaterial) {

         let v1 = new THREE.Vector3(),
             v2 = new THREE.Vector3(),
             d = 0, distances = null;

         if (index) {
            distances = new Float32Array(index.length);
            for (let n=0; n<index.length; n+=2) {
               let i1 = index[n], i2 = index[n+1];
               v1.set(arr[i1],arr[i1+1],arr[i1+2]);
               v2.set(arr[i2],arr[i2+1],arr[i2+2]);
               distances[n] = d;
               d += v2.distanceTo( v1 );
               distances[n+1] = d;
            }
         } else {
            distances = new Float32Array(arr.length/3);
            for (let n=0; n<arr.length; n+=6) {
               v1.set(arr[n],arr[n+1],arr[n+2]);
               v2.set(arr[n+3],arr[n+4],arr[n+5]);
               distances[n/3] = d;
               d += v2.distanceTo( v1 );
               distances[n/3+1] = d;
            }
         }
         geom.setAttribute( 'lineDistance', new THREE.BufferAttribute(distances, 1) );
      }

      return only_geometry ? geom : new THREE.LineSegments(geom, material);
   }

   /** @summary Help structures for calculating Box mesh
     * @private */
   jsrp.Box3D = {
       Vertices: [ new THREE.Vector3(1, 1, 1), new THREE.Vector3(1, 1, 0),
                   new THREE.Vector3(1, 0, 1), new THREE.Vector3(1, 0, 0),
                   new THREE.Vector3(0, 1, 0), new THREE.Vector3(0, 1, 1),
                   new THREE.Vector3(0, 0, 0), new THREE.Vector3(0, 0, 1) ],
       Indexes: [ 0,2,1, 2,3,1, 4,6,5, 6,7,5, 4,5,1, 5,0,1, 7,6,2, 6,3,2, 5,7,0, 7,2,0, 1,3,4, 3,6,4 ],
       Normals: [ 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1 ],
       Segments: [0, 2, 2, 7, 7, 5, 5, 0, 1, 3, 3, 6, 6, 4, 4, 1, 1, 0, 3, 2, 6, 7, 4, 5],  // segments addresses Vertices
   };

   // these segments address vertices from the mesh, we can use positions from box mesh
   jsrp.Box3D.MeshSegments = (function() {
      let box3d = jsrp.Box3D,
          arr = new Int32Array(box3d.Segments.length);

      for (let n=0;n<arr.length;++n) {
         for (let k=0;k<box3d.Indexes.length;++k)
            if (box3d.Segments[n] === box3d.Indexes[k]) {
               arr[n] = k; break;
            }
      }
      return arr;
   })();

   // ==============================================================================

   function InteractiveControl() {}

   InteractiveControl.prototype.cleanup = function() {}

   InteractiveControl.prototype.extractIndex = function(/*intersect*/) { return undefined; }

   InteractiveControl.prototype.setSelected = function(/*col, indx*/) {}

   InteractiveControl.prototype.setHighlight = function(/*col, indx*/) {}

   InteractiveControl.prototype.checkHighlightIndex = function(/*indx*/) { return undefined; }

   // ==============================================================================

   /**
    * @summary Special class to control highliht and selection of single points, used in geo painter
    *
    * @class
    * @memberof JSROOT
    * @param {object} mesh - draw object
    * @private
    */

   function PointsControl(mesh) {
      InteractiveControl.call(this);
      this.mesh = mesh;
   }

   PointsControl.prototype = Object.create(InteractiveControl.prototype);

   /** @summary cleanup object */
   PointsControl.prototype.cleanup = function() {
      if (!this.mesh) return;
      delete this.mesh.is_selected;
      this.createSpecial(null);
      delete this.mesh;
   }

   /** @summary extract intersect index */
   PointsControl.prototype.extractIndex = function(intersect) {
      return intersect && intersect.index!==undefined ? intersect.index : undefined;
   }

   /** @summary set selection */
   PointsControl.prototype.setSelected = function(col, indx) {
      let m = this.mesh;
      if ((m.select_col == col) && (m.select_indx == indx)) {
         col = null; indx = undefined;
      }
      m.select_col = col;
      m.select_indx = indx;
      this.createSpecial(col, indx);
      return true;
   }

   /** @summary set highlight */
   PointsControl.prototype.setHighlight = function(col, indx) {
      let m = this.mesh;
      m.h_index = indx;
      if (col)
         this.createSpecial(col, indx);
      else
         this.createSpecial(m.select_col, m.select_indx);
      return true;
   }

   /** @summary create special object */
   PointsControl.prototype.createSpecial = function(color, index) {
      let m = this.mesh;
      if (!color) {
         if (m.js_special) {
            m.remove(m.js_special);
            jsrp.disposeThreejsObject(m.js_special);
            delete m.js_special;
         }
         return;
      }

      if (!m.js_special) {
         let geom = new THREE.BufferGeometry();
         geom.setAttribute( 'position', m.geometry.getAttribute("position"));
         let material = new THREE.PointsMaterial( { size: m.material.size*2, color: color } );
         material.sizeAttenuation = m.material.sizeAttenuation;

         m.js_special = new THREE.Points(geom, material);
         m.js_special.jsroot_special = true; // special object, exclude from intersections
         m.add(m.js_special);
      }

      if (color) m.js_special.material.color = new THREE.Color(color);
      if (index !== undefined) m.js_special.geometry.setDrawRange(index, 1);
   }

   // ==============================================================================

   /**
    * @summary Class for creation of 3D points
    *
    * @class
    * @memberof JSROOT
    * @param {number} size - number of points
    * @param {booleand} [iswebgl=true] - if WebGL is used
    * @param {number} [scale=1] - scale factor
    * @private
    */

   function PointsCreator(size, iswebgl, scale) {
      this.webgl = (iswebgl === undefined) ? true : iswebgl;
      this.scale = scale || 1.;

      this.pos = new Float32Array(size*3);
      this.geom = new THREE.BufferGeometry();
      this.geom.setAttribute('position', new THREE.BufferAttribute(this.pos, 3));
      this.indx = 0;
   }

   /** @summary Add point */
   PointsCreator.prototype.addPoint = function(x,y,z) {
      this.pos[this.indx]   = x;
      this.pos[this.indx+1] = y;
      this.pos[this.indx+2] = z;
      this.indx+=3;
   }

   /** @summary Create points */
   PointsCreator.prototype.createPoints = function(args) {

      if (typeof args !== 'object')
         args = { color: args };
      if (!args.color)
         args.color = 'black';

      let k = 1;

      // special dots
      if (!args.style) k = 1.1; else
      if (args.style === 1) k = 0.3; else
      if (args.style === 2) args.style = 3; else // just avoid plot of "+" sign, issue #205
      if (args.style === 6) k = 0.5; else
      if (args.style === 7) k = 0.7;

      let material;

      if (!args.style || (k !== 1) || JSROOT.nodejs) {
         // this is plain creation of points, no texture loading, which does not work in node.js
         material = new THREE.PointsMaterial( { size: (this.webgl ? 3 : 1) * this.scale * k, color: args.color } );

      } else {

         let handler = new JSROOT.TAttMarkerHandler({ style: args.style, color: args.color, size: 8 });

         let plainSVG = '<svg width="70" height="70" xmlns="http://www.w3.org/2000/svg">' +
                        '<path d="' + handler.create(35,35) + '" stroke="' + handler.getStrokeColor() + '" stroke-width="4" fill="' + handler.getFillColor() + '"/>' +
                        '</svg>';

         // let need_replace = JSROOT.nodejs && !globalThis.document;
         // if (need_replace) globalThis.document = JSROOT._.get_document();

         let texture = new THREE.TextureLoader().load( 'data:image/svg+xml;utf8,' + plainSVG);

         // if (need_replace) globalThis.document = undefined;

         material = new THREE.PointsMaterial( { size: (this.webgl ? 3 : 1) * this.scale, map: texture, transparent: true } );
      }

      let pnts = new THREE.Points(this.geom, material);
      pnts.nvertex = 1;

      return pnts;
   }


   // ==============================================================================

   /** @summary Create material for 3D line
     * @desc Takes into account dashed properties
     * @private
     * @memberof JSROOT.Painter */
   function create3DLineMaterial(painter, obj) {
      if (!painter || !obj) return null;

      let lcolor = painter.getColor(obj.fLineColor),
          material = null,
          style = obj.fLineStyle ? jsrp.root_line_styles[obj.fLineStyle] : "",
          dash = style ? style.split(",") : [];

      if (dash && dash.length>=2)
         material = new THREE.LineDashedMaterial( { color: lcolor, dashSize: parseInt(dash[0]), gapSize: parseInt(dash[1]) } );
      else
         material = new THREE.LineBasicMaterial({ color: lcolor });

      if ((obj.fLineWidth !== undefined) && (obj.fLineWidth > 1)) material.linewidth = obj.fLineWidth;

      return material;
   }

   // ============================================================================================================

   /** @summary Draw TPolyLine3D object
     * @desc Takes into account dashed properties
     * @private
     * @memberof JSROOT.Painter */
   function drawPolyLine3D() {
      let line = this.getObject(),
          main = this.getFramePainter();

      if (!main || !main.mode3d || !main.toplevel || !line)
         return null;

      let fN, fP, pnts = [];

      if (line._blob && (line._blob.length==4)) {
         // workaround for custom streamer for JSON, should be resolved
         fN = line._blob[1];
         fP = line._blob[2];
         // fOption = line._blob[3];
      } else {
         fN = line.fN;
         fP = line.fP;
         // fOption = line.fOption;
      }

      for (let n = 3; n < 3*fN; n += 3)
         pnts.push(main.grx(fP[n-3]), main.gry(fP[n-2]), main.grz(fP[n-1]),
                   main.grx(fP[n]), main.gry(fP[n+1]), main.grz(fP[n+2]));

      let lines = jsrp.createLineSegments(pnts, create3DLineMaterial(this, line));

      main.toplevel.add(lines);

      return true;
   }

   // ==============================================================================================


   jsrp.PointsCreator = PointsCreator;
   jsrp.InteractiveControl = InteractiveControl;
   jsrp.PointsControl = PointsControl;

   jsrp.drawPolyLine3D = drawPolyLine3D;
   jsrp.create3DLineMaterial = create3DLineMaterial;

   if (JSROOT.nodejs) module.exports = THREE;
   return THREE;

});
