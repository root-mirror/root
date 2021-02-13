/// @file JSRoot.interactive.js
/// Basic interactive functionality

JSROOT.define(['d3', 'painter'], (d3, jsrp) => {

   "use strict";

   let TooltipHandler = {

      /** @desc only canvas info_layer can be used while other pads can overlay
        * @returns layer where frame tooltips are shown */
      hints_layer: function() {
         let pp = this.getCanvPainter();
         return pp ? pp.getLayerSvg("info_layer") : d3.select(null);
      },

      /** @returns true if tooltip is shown, use to prevent some other action */
      isTooltipShown: function() {
         if (!this.tooltip_enabled || !this.isTooltipAllowed()) return false;
         let hintsg = this.hints_layer().select(".objects_hints");
         return hintsg.empty() ? false : hintsg.property("hints_pad") == this.getPadName();
      },

      setTooltipEnabled: function(enabled) {
         if (enabled !== undefined) this.tooltip_enabled = enabled;
      },

      processFrameTooltipEvent: function(pnt, evnt) {
         // make central function which let show selected hints for the object

         if (pnt && pnt.handler) {
            // special use of interactive handler in the frame painter
            let rect = this.draw_g ? this.draw_g.select(".main_layer") : null;
            if (!rect || rect.empty()) {
               pnt = null; // disable
            } else if (pnt.touch && evnt) {
               let pos = d3.pointers(evnt, rect.node());
               pnt = (pos && pos.length == 1) ? { touch: true, x: pos[0][0], y: pos[0][1] } : null;
            } else if (evnt) {
               let pos = d3.pointer(evnt, rect.node());
               pnt = { touch: false, x: pos[0], y: pos[1] };
            }
         }

         let hints = [], nhints = 0, maxlen = 0, lastcolor1 = 0, usecolor1 = false,
            textheight = 11, hmargin = 3, wmargin = 3, hstep = 1.2,
            frame_rect = this.getFrameRect(),
            pp = this.getPadPainter(),
            pad_width = pp.getPadWidth(),
            font = new JSROOT.FontHandler(160, textheight),
            disable_tootlips = !this.isTooltipAllowed() || !this.tooltip_enabled;

         if (pnt && disable_tootlips) pnt.disabled = true; // indicate that highlighting is not required
         if (pnt) pnt.painters = true; // get also painter

         // collect tooltips from pad painter - it has list of all drawn objects
         if (pp) hints = pp.processPadTooltipEvent(pnt);

         if (pnt && pnt.touch) textheight = 15;

         for (let n = 0; n < hints.length; ++n) {
            let hint = hints[n];
            if (!hint) continue;

            if (hint.painter && (hint.user_info !== undefined))
               hint.painter.provideUserTooltip(hint.user_info);

            if (!hint.lines || (hint.lines.length === 0)) {
               hints[n] = null; continue;
            }

            // check if fully duplicated hint already exists
            for (let k = 0; k < n; ++k) {
               let hprev = hints[k], diff = false;
               if (!hprev || (hprev.lines.length !== hint.lines.length)) continue;
               for (let l = 0; l < hint.lines.length && !diff; ++l)
                  if (hprev.lines[l] !== hint.lines[l]) diff = true;
               if (!diff) { hints[n] = null; break; }
            }
            if (!hints[n]) continue;

            nhints++;

            for (let l = 0; l < hint.lines.length; ++l)
               maxlen = Math.max(maxlen, hint.lines[l].length);

            hint.height = Math.round(hint.lines.length * textheight * hstep + 2 * hmargin - textheight * (hstep - 1));

            if ((hint.color1 !== undefined) && (hint.color1 !== 'none')) {
               if ((lastcolor1 !== 0) && (lastcolor1 !== hint.color1)) usecolor1 = true;
               lastcolor1 = hint.color1;
            }
         }

         let layer = this.hints_layer(),
            hintsg = layer.select(".objects_hints"); // group with all tooltips

         let title = "", name = "", info = "",
            hint = null, best_dist2 = 1e10, best_hint = null,
            coordinates = pnt ? Math.round(pnt.x) + "," + Math.round(pnt.y) : "";
         // try to select hint with exact match of the position when several hints available
         for (let k = 0; k < (hints ? hints.length : 0); ++k) {
            if (!hints[k]) continue;
            if (!hint) hint = hints[k];
            if (hints[k].exact && (!hint || !hint.exact)) { hint = hints[k]; break; }

            if (!pnt || (hints[k].x === undefined) || (hints[k].y === undefined)) continue;

            let dist2 = (pnt.x - hints[k].x) * (pnt.x - hints[k].x) + (pnt.y - hints[k].y) * (pnt.y - hints[k].y);
            if (dist2 < best_dist2) { best_dist2 = dist2; best_hint = hints[k]; }
         }

         if ((!hint || !hint.exact) && (best_dist2 < 400)) hint = best_hint;

         if (hint) {
            name = (hint.lines && hint.lines.length > 1) ? hint.lines[0] : hint.name;
            title = hint.title || "";
            info = hint.line;
            if (!info && hint.lines) info = hint.lines.slice(1).join(' ');
         }

         this.showObjectStatus(name, title, info, coordinates);

         // end of closing tooltips
         if (!pnt || disable_tootlips || (hints.length === 0) || (maxlen === 0) || (nhints > 15)) {
            hintsg.remove();
            return;
         }

         // we need to set pointer-events=none for all elements while hints
         // placed in front of so-called interactive rect in frame, used to catch mouse events

         if (hintsg.empty())
            hintsg = layer.append("svg:g")
               .attr("class", "objects_hints")
               .style("pointer-events", "none");

         let frame_shift = { x: 0, y: 0 }, trans = frame_rect.transform || "";
         if (!pp.iscan) {
            frame_shift = jsrp.getAbsPosInCanvas(this.getPadSvg(), frame_shift);
            trans = "translate(" + frame_shift.x + "," + frame_shift.y + ") " + trans;
         }

         // copy transform attributes from frame itself
         hintsg.attr("transform", trans)
            .property("last_point", pnt)
            .property("hints_pad", this.getPadName());

         let viewmode = hintsg.property('viewmode') || "",
            actualw = 0, posx = pnt.x + frame_rect.hint_delta_x;

         if (nhints > 1) {
            // if there are many hints, place them left or right

            let bleft = 0.5, bright = 0.5;

            if (viewmode == "left") bright = 0.7; else
               if (viewmode == "right") bleft = 0.3;

            if (posx <= bleft * frame_rect.width) {
               viewmode = "left";
               posx = 20;
            } else if (posx >= bright * frame_rect.width) {
               viewmode = "right";
               posx = frame_rect.width - 60;
            } else {
               posx = hintsg.property('startx');
            }
         } else {
            viewmode = "single";
            posx += 15;
         }

         if (viewmode !== hintsg.property('viewmode')) {
            hintsg.property('viewmode', viewmode);
            hintsg.selectAll("*").remove();
         }

         let curry = 10, // normal y coordinate
            gapy = 10,  // y coordinate, taking into account all gaps
            gapminx = -1111, gapmaxx = -1111,
            minhinty = -frame_shift.y,
            maxhinty = this.getCanvPainter().getPadHeight() - frame_rect.y - frame_shift.y;

         function FindPosInGap(y) {
            for (let n = 0; (n < hints.length) && (y < maxhinty); ++n) {
               let hint = hints[n];
               if (!hint) continue;
               if ((hint.y >= y - 5) && (hint.y <= y + hint.height + 5)) {
                  y = hint.y + 10;
                  n = -1;
               }
            }
            return y;
         }

         for (let n = 0; n < hints.length; ++n) {
            let hint = hints[n],
               group = hintsg.select(".painter_hint_" + n);
            if (hint === null) {
               group.remove();
               continue;
            }

            let was_empty = group.empty();

            if (was_empty)
               group = hintsg.append("svg:svg")
                  .attr("class", "painter_hint_" + n)
                  .attr('opacity', 0) // use attribute, not style to make animation with d3.transition()
                  .style('overflow', 'hidden')
                  .style("pointer-events", "none");

            if (viewmode == "single") {
               curry = pnt.touch ? (pnt.y - hint.height - 5) : Math.min(pnt.y + 15, maxhinty - hint.height - 3) + frame_rect.hint_delta_y;
            } else {
               gapy = FindPosInGap(gapy);
               if ((gapminx === -1111) && (gapmaxx === -1111)) gapminx = gapmaxx = hint.x;
               gapminx = Math.min(gapminx, hint.x);
               gapmaxx = Math.min(gapmaxx, hint.x);
            }

            group.attr("x", posx)
               .attr("y", curry)
               .property("curry", curry)
               .property("gapy", gapy);

            curry += hint.height + 5;
            gapy += hint.height + 5;

            if (!was_empty)
               group.selectAll("*").remove();

            group.attr("width", 60)
               .attr("height", hint.height);

            let r = group.append("rect")
               .attr("x", 0)
               .attr("y", 0)
               .attr("width", 60)
               .attr("height", hint.height)
               .attr("fill", "lightgrey")
               .style("pointer-events", "none");

            if (nhints > 1) {
               let col = usecolor1 ? hint.color1 : hint.color2;
               if ((col !== undefined) && (col !== 'none'))
                  r.attr("stroke", col).attr("stroke-width", hint.exact ? 3 : 1);
            }

            for (let l = 0; l < (hint.lines ? hint.lines.length : 0); l++)
               if (hint.lines[l] !== null) {
                  let txt = group.append("svg:text")
                     .attr("text-anchor", "start")
                     .attr("x", wmargin)
                     .attr("y", hmargin + l * textheight * hstep)
                     .attr("dy", ".8em")
                     .attr("fill", "black")
                     .style("pointer-events", "none")
                     .call(font.func)
                     .text(hint.lines[l]);

                  let box = jsrp.getElementRect(txt, 'bbox');

                  actualw = Math.max(actualw, box.width);
               }

            function translateFn() {
               // We only use 'd', but list d,i,a as params just to show can have them as params.
               // Code only really uses d and t.
               return function(/*d, i, a*/) {
                  return function(t) {
                     return t < 0.8 ? "0" : (t - 0.8) * 5;
                  };
               };
            }

            if (was_empty)
               if (JSROOT.settings.TooltipAnimation > 0)
                  group.transition().duration(JSROOT.settings.TooltipAnimation).attrTween("opacity", translateFn());
               else
                  group.attr('opacity', 1);
         }

         actualw += 2 * wmargin;

         let svgs = hintsg.selectAll("svg");

         if ((viewmode == "right") && (posx + actualw > frame_rect.width - 20)) {
            posx = frame_rect.width - actualw - 20;
            svgs.attr("x", posx);
         }

         if ((viewmode == "single") && (posx + actualw > pad_width - frame_rect.x) && (posx > actualw + 20)) {
            posx -= (actualw + 20);
            svgs.attr("x", posx);
         }

         // if gap not very big, apply gapy coordinate to open view on the histogram
         if ((viewmode !== "single") && (gapy < maxhinty) && (gapy !== curry)) {
            if ((gapminx <= posx + actualw + 5) && (gapmaxx >= posx - 5))
               svgs.attr("y", function() { return d3.select(this).property('gapy'); });
         } else if ((viewmode !== 'single') && (curry > maxhinty)) {
            let shift = Math.max((maxhinty - curry - 10), minhinty);
            if (shift < 0)
               svgs.attr("y", function() { return d3.select(this).property('curry') + shift; });
         }

         if (actualw > 10)
            svgs.attr("width", actualw)
               .select('rect').attr("width", actualw);

         hintsg.property('startx', posx);
      },

      /** @summary Assigns tooltip methods */
      assign: function(painter) {
         painter.tooltip_enabled = true;
         painter.hints_layer = this.hints_layer;
         painter.isTooltipShown = this.isTooltipShown;
         painter.setTooltipEnabled = this.setTooltipEnabled;
         painter.processFrameTooltipEvent = this.processFrameTooltipEvent;
      }

   } // TooltipHandler

   let setPainterTooltipEnabled = (painter,on) => {
      if (!painter) return;

      let fp = painter.getFramePainter();
      if (fp && typeof fp.setTooltipEnabled == 'function') {
         fp.setTooltipEnabled(on);
         fp.processFrameTooltipEvent(null);
      }
      // this is 3D control object
      if (this.control && (typeof this.control.setTooltipEnabled == 'function'))
         this.control.setTooltipEnabled(on);
   }


   /** @summary Add drag for interactive rectangular elements for painter */
   function addDragHandler(painter, callback) {
      if (!JSROOT.settings.MoveResize || JSROOT.batch_mode) return;

      let pthis = painter, drag_rect = null, pp = pthis.getPadPainter();
      if (pp && pp._fast_drawing) return;

      function detectRightButton(event) {
         if ('buttons' in event) return event.buttons === 2;
         else if ('which' in event) return event.which === 3;
         else if ('button' in event) return event.button === 2;
         return false;
      }

      function rect_width() { return Number(pthis.draw_g.attr("width")); }
      function rect_height() { return Number(pthis.draw_g.attr("height")); }

      function MakeResizeElements(group, width, height, handler) {
         function make(cursor, d) {
            let clname = "js_" + cursor.replace('-', '_'),
               elem = group.select('.' + clname);
            if (elem.empty()) elem = group.append('path').classed(clname, true);
            elem.style('opacity', 0).style('cursor', cursor).attr('d', d);
            if (handler) elem.call(handler);
         }

         make("nw-resize", "M2,2h15v-5h-20v20h5Z");
         make("ne-resize", "M" + (width - 2) + ",2h-15v-5h20v20h-5 Z");
         make("sw-resize", "M2," + (height - 2) + "h15v5h-20v-20h5Z");
         make("se-resize", "M" + (width - 2) + "," + (height - 2) + "h-15v5h20v-20h-5Z");

         if (!callback.no_change_x) {
            make("w-resize", "M-3,18h5v" + Math.max(0, height - 2 * 18) + "h-5Z");
            make("e-resize", "M" + (width + 3) + ",18h-5v" + Math.max(0, height - 2 * 18) + "h5Z");
         }
         if (!callback.no_change_y) {
            make("n-resize", "M18,-3v5h" + Math.max(0, width - 2 * 18) + "v-5Z");
            make("s-resize", "M18," + (height + 3) + "v-5h" + Math.max(0, width - 2 * 18) + "v5Z");
         }
      }

      function complete_drag() {
         drag_rect.style("cursor", "auto");

         if (!pthis.draw_g) {
            drag_rect.remove();
            drag_rect = null;
            return false;
         }

         let oldx = Number(pthis.draw_g.attr("x")),
            oldy = Number(pthis.draw_g.attr("y")),
            newx = Number(drag_rect.attr("x")),
            newy = Number(drag_rect.attr("y")),
            newwidth = Number(drag_rect.attr("width")),
            newheight = Number(drag_rect.attr("height"));

         if (callback.minwidth && newwidth < callback.minwidth) newwidth = callback.minwidth;
         if (callback.minheight && newheight < callback.minheight) newheight = callback.minheight;

         let change_size = (newwidth !== rect_width()) || (newheight !== rect_height()),
            change_pos = (newx !== oldx) || (newy !== oldy);

         pthis.draw_g.attr('x', newx).attr('y', newy)
            .attr("transform", "translate(" + newx + "," + newy + ")")
            .attr('width', newwidth).attr('height', newheight);

         drag_rect.remove();
         drag_rect = null;

         setPainterTooltipEnabled(pthis, true);

         MakeResizeElements(pthis.draw_g, newwidth, newheight);

         if (change_size || change_pos) {
            if (change_size && ('resize' in callback)) callback.resize(newwidth, newheight);
            if (change_pos && ('move' in callback)) callback.move(newx, newy, newx - oldxx, newy - oldy);

            if (change_size || change_pos) {
               if ('obj' in callback) {
                  let rect = pp.getPadRect();
                  callback.obj.fX1NDC = newx / rect.width;
                  callback.obj.fX2NDC = (newx + newwidth) / rect.width;
                  callback.obj.fY1NDC = 1 - (newy + newheight) / rect.height;
                  callback.obj.fY2NDC = 1 - newy / rect.height;
                  callback.obj.modified_NDC = true; // indicate that NDC was interactively changed, block in updated
               }
               if ('redraw' in callback) callback.redraw();
            }
         }

         return change_size || change_pos;
      }

      let drag_move = d3.drag().subject(Object),
          drag_resize = d3.drag().subject(Object);

      drag_move
         .on("start", function(evnt) {
            if (detectRightButton(evnt.sourceEvent)) return;

            if (jsrp.closeMenu) jsrp.closeMenu(); // close menu

            setPainterTooltipEnabled(pthis, false); // disable tooltip

            evnt.sourceEvent.preventDefault();
            evnt.sourceEvent.stopPropagation();

            let pad_rect = pp.getPadRect();

            let handle = {
               acc_x1: Number(pthis.draw_g.attr("x")),
               acc_y1: Number(pthis.draw_g.attr("y")),
               pad_w: pad_rect.width - rect_width(),
               pad_h: pad_rect.height - rect_height(),
               drag_tm: new Date(),
               path: "v" + rect_height() + "h" + rect_width() + "v" + (-rect_height()) + "z"
            };

            drag_rect = d3.select(pthis.draw_g.node().parentNode).append("path")
               .classed("zoom", true)
               .attr("x", handle.acc_x1)
               .attr("y", handle.acc_y1)
               .attr("width", rect_width())
               .attr("height", rect_height())
               .attr("d", "M" + handle.acc_x1 + "," + handle.acc_y1 + handle.path)
               .style("cursor", "move")
               .style("pointer-events", "none") // let forward double click to underlying elements
               .property('drag_handle', handle);


         }).on("drag", function(evnt) {
            if (!drag_rect) return;

            evnt.sourceEvent.preventDefault();
            evnt.sourceEvent.stopPropagation();

            let handle = drag_rect.property('drag_handle');

            if (!callback.no_change_x)
               handle.acc_x1 += evnt.dx;
            if (!callback.no_change_y)
               handle.acc_y1 += evnt.dy;

            let x = Math.min(Math.max(handle.acc_x1, 0), handle.pad_w),
                y = Math.min(Math.max(handle.acc_y1, 0), handle.pad_h);

            drag_rect.attr("x", x)
                     .attr("y", y)
                     .attr("d", "M" + x + "," + y + handle.path);

         }).on("end", function(evnt) {
            if (!drag_rect) return;

            evnt.sourceEvent.preventDefault();

            let handle = drag_rect.property('drag_handle');

            if (complete_drag() === false) {
               let spent = (new Date()).getTime() - handle.drag_tm.getTime();
               if (callback.ctxmenu && (spent > 600) && pthis.showContextMenu) {
                  let rrr = resize_se.node().getBoundingClientRect();
                  pthis.showContextMenu('main', { clientX: rrr.left, clientY: rrr.top });
               } else if (callback.canselect && (spent <= 600)) {
                  let pp = pthis.getPadPainter();
                  if (pp) pp.selectObjectPainter(pthis);
               }
            }
         });

      drag_resize
         .on("start", function(evnt) {
            if (detectRightButton(evnt.sourceEvent)) return;

            evnt.sourceEvent.stopPropagation();
            evnt.sourceEvent.preventDefault();

            setPainterTooltipEnabled(pthis, false); // disable tooltip

            let pad_rect = pp.getPadRect();

            let handle = {
               acc_x1: Number(pthis.draw_g.attr("x")),
               acc_y1: Number(pthis.draw_g.attr("y")),
               pad_w: pad_rect.width,
               pad_h: pad_rect.height
            };

            handle.acc_x2 = handle.acc_x1 + rect_width();
            handle.acc_y2 = handle.acc_y1 + rect_height();

            drag_rect = d3.select(pthis.draw_g.node().parentNode)
               .append("rect")
               .classed("zoom", true)
               .style("cursor", d3.select(this).style("cursor"))
               .attr("x", handle.acc_x1)
               .attr("y", handle.acc_y1)
               .attr("width", handle.acc_x2 - handle.acc_x1)
               .attr("height", handle.acc_y2 - handle.acc_y1)
               .property('drag_handle', handle);

         }).on("drag", function(evnt) {
            if (!drag_rect) return;

            evnt.sourceEvent.preventDefault();
            evnt.sourceEvent.stopPropagation();

            let handle = drag_rect.property('drag_handle'),
               dx = evnt.dx, dy = evnt.dy, elem = d3.select(this);

            if (callback.no_change_x) dx = 0;
            if (callback.no_change_y) dy = 0;

            if (elem.classed('js_nw_resize')) { handle.acc_x1 += dx; handle.acc_y1 += dy; }
            else if (elem.classed('js_ne_resize')) { handle.acc_x2 += dx; handle.acc_y1 += dy; }
            else if (elem.classed('js_sw_resize')) { handle.acc_x1 += dx; handle.acc_y2 += dy; }
            else if (elem.classed('js_se_resize')) { handle.acc_x2 += dx; handle.acc_y2 += dy; }
            else if (elem.classed('js_w_resize')) { handle.acc_x1 += dx; }
            else if (elem.classed('js_n_resize')) { handle.acc_y1 += dy; }
            else if (elem.classed('js_e_resize')) { handle.acc_x2 += dx; }
            else if (elem.classed('js_s_resize')) { handle.acc_y2 += dy; }

            let x1 = Math.max(0, handle.acc_x1), x2 = Math.min(handle.acc_x2, handle.pad_w),
               y1 = Math.max(0, handle.acc_y1), y2 = Math.min(handle.acc_y2, handle.pad_h);

            drag_rect.attr("x", x1).attr("y", y1).attr("width", Math.max(0, x2 - x1)).attr("height", Math.max(0, y2 - y1));

         }).on("end", function(evnt) {
            if (!drag_rect) return;

            evnt.sourceEvent.preventDefault();

            complete_drag();
         });

      if (!callback.only_resize)
         pthis.draw_g.style("cursor", "move").call(drag_move);

      if (!callback.only_move)
         MakeResizeElements(pthis.draw_g, rect_width(), rect_height(), drag_resize);
   }

   /** @summary Add move handlers for drawn element
     * @private */
   function addMoveHandler(painter) {

      if (!JSROOT.settings.MoveResize || JSROOT.batch_mode ||
         !painter.draw_g || painter.draw_g.property("assigned_move")) return;

      function detectRightButton(event) {
         if ('buttons' in event) return event.buttons === 2;
         else if ('which' in event) return event.which === 3;
         else if ('button' in event) return event.button === 2;
         return false;
      }

      let drag_move = d3.drag().subject(Object),
         not_changed = true;

      drag_move
         .on("start", function(evnt) {
            if (detectRightButton(evnt.sourceEvent)) return;
            evnt.sourceEvent.preventDefault();
            evnt.sourceEvent.stopPropagation();
            let pos = d3.pointer(evnt, this.draw_g.node());
            not_changed = true;
            if (this.moveStart)
               this.moveStart(pos[0], pos[1]);
         }.bind(painter)).on("drag", function(evnt) {
            evnt.sourceEvent.preventDefault();
            evnt.sourceEvent.stopPropagation();
            not_changed = false;
            if (this.moveDrag)
               this.moveDrag(evnt.dx, evnt.dy);
         }.bind(painter)).on("end", function(evnt) {
            evnt.sourceEvent.preventDefault();
            evnt.sourceEvent.stopPropagation();
            if (this.moveEnd)
               this.moveEnd(not_changed);
            let pp = this.getPadPainter();
            if (pp) pp.selectObjectPainter(this);
         }.bind(painter));

      painter.draw_g
             .style("cursor", "move")
             .property("assigned_move", true)
             .call(drag_move);
   }

   // ================================================================================

   let FrameInteractive = {

      addBasicInteractivity: function() {

         TooltipHandler.assign(this);

         this.draw_g.attr("x", this._frame_x)
                    .attr("y", this._frame_y)
                    .attr("width", this.getFrameWidth())
                    .attr("height", this.getFrameHeight());

         if (!this._frame_rotate && !this._frame_fixpos)
            addDragHandler(this, { obj: this, only_resize: true,
                                    minwidth: 20, minheight: 20, redraw: this.sizeChanged.bind(this) });

         let main_svg = this.draw_g.select(".main_layer");

         main_svg.style("pointer-events","visibleFill")
                 .property('handlers_set', 0);

         let pp = this.getPadPainter(),
             handlers_set = (pp && pp._fast_drawing) ? 0 : 1;

         if (main_svg.property('handlers_set') != handlers_set) {
            let close_handler = handlers_set ? this.processFrameTooltipEvent.bind(this, null) : null,
                mouse_handler = handlers_set ? this.processFrameTooltipEvent.bind(this, { handler: true, touch: false }) : null;

            main_svg.property('handlers_set', handlers_set)
                    .on('mouseenter', mouse_handler)
                    .on('mousemove', mouse_handler)
                    .on('mouseleave', close_handler);

            if (JSROOT.browser.touches) {
               let touch_handler = handlers_set ? this.processFrameTooltipEvent.bind(this, { handler: true, touch: true }) : null;

               main_svg.on("touchstart", touch_handler)
                       .on("touchmove", touch_handler)
                       .on("touchend", close_handler)
                       .on("touchcancel", close_handler);
            }
         }

         main_svg.attr("x", 0)
                 .attr("y", 0)
                 .attr("width", this.getFrameWidth())
                 .attr("height", this.getFrameHeight());

         let hintsg = this.hints_layer().select(".objects_hints");
         // if tooltips were visible before, try to reconstruct them after short timeout
         if (!hintsg.empty() && this.isTooltipAllowed() && (hintsg.property("hints_pad") == this.getPadName()))
            setTimeout(this.processFrameTooltipEvent.bind(this, hintsg.property('last_point'), null), 10);
      },

      addInteractivity: function() {

         let pp = this.getPadPainter(),
             svg = this.getFrameSvg();
         if ((pp && pp._fast_drawing) || svg.empty())
            return Promise.resolve(this);

         let svg_x = svg.selectAll(".xaxis_container"),
             svg_y = svg.selectAll(".yaxis_container");

         if (!svg.property('interactive_set')) {
            this.addKeysHandler();

            this.last_touch = new Date(0);
            this.zoom_kind = 0; // 0 - none, 1 - XY, 2 - only X, 3 - only Y, (+100 for touches)
            this.zoom_rect = null;
            this.zoom_origin = null;  // original point where zooming started
            this.zoom_curr = null;    // current point for zooming
            this.touch_cnt = 0;
         }

         if (JSROOT.settings.Zooming && !this.projection) {
            if (JSROOT.settings.ZoomMouse) {
               svg.on("mousedown", this.startRectSel.bind(this));
               svg.on("dblclick", this.mouseDoubleClick.bind(this));
            }
            if (JSROOT.settings.ZoomWheel)
               svg.on("wheel", this.mouseWheel.bind(this));
         }

         if (JSROOT.browser.touches && ((JSROOT.settings.Zooming && JSROOT.settings.ZoomTouch && !this.projection) || JSROOT.settings.ContextMenu))
            svg.on("touchstart", this.startTouchZoom.bind(this));

         if (JSROOT.settings.ContextMenu) {
            if (JSROOT.browser.touches) {
               svg_x.on("touchstart", this.startTouchMenu.bind(this,"x"));
               svg_y.on("touchstart", this.startTouchMenu.bind(this,"y"));
            }
            svg.on("contextmenu", evnt => this.showContextMenu("", evnt));
            svg_x.on("contextmenu", evnt => this.showContextMenu("x", evnt));
            svg_y.on("contextmenu", evnt => this.showContextMenu("y", evnt));
         }

         svg_x.on("mousemove", evnt => this.showAxisStatus("x", evnt));
         svg_y.on("mousemove", evnt => this.showAxisStatus("y", evnt));

         svg.property('interactive_set', true);

         return Promise.resolve(this);
      },

      addKeysHandler: function() {
         if (this.keys_handler || (typeof window == 'undefined')) return;

         this.keys_handler = evnt => this.processKeyPress(evnt);

         window.addEventListener('keydown', this.keys_handler, false);
      },

      processKeyPress: function(evnt) {
         let main = this.selectDom();
         if (!JSROOT.settings.HandleKeys || main.empty() || (this.enabledKeys === false)) return;

         let key = "";
         switch (evnt.keyCode) {
            case 33: key = "PageUp"; break;
            case 34: key = "PageDown"; break;
            case 37: key = "ArrowLeft"; break;
            case 38: key = "ArrowUp"; break;
            case 39: key = "ArrowRight"; break;
            case 40: key = "ArrowDown"; break;
            case 42: key = "PrintScreen"; break;
            case 106: key = "*"; break;
            default: return false;
         }

         let pp = this.getPadPainter();
         if (jsrp.getActivePad() !== pp) return;

         if (evnt.shiftKey) key = "Shift " + key;
         if (evnt.altKey) key = "Alt " + key;
         if (evnt.ctrlKey) key = "Ctrl " + key;

         let zoom = { name: "x", dleft: 0, dright: 0 };

         switch (key) {
            case "ArrowLeft":  zoom.dleft = -1; zoom.dright = 1; break;
            case "ArrowRight":  zoom.dleft = 1; zoom.dright = -1; break;
            case "Ctrl ArrowLeft": zoom.dleft = zoom.dright = -1; break;
            case "Ctrl ArrowRight": zoom.dleft = zoom.dright = 1; break;
            case "ArrowUp":  zoom.name = "y"; zoom.dleft = 1; zoom.dright = -1; break;
            case "ArrowDown":  zoom.name = "y"; zoom.dleft = -1; zoom.dright = 1; break;
            case "Ctrl ArrowUp": zoom.name = "y"; zoom.dleft = zoom.dright = 1; break;
            case "Ctrl ArrowDown": zoom.name = "y"; zoom.dleft = zoom.dright = -1; break;
         }

         if (zoom.dleft || zoom.dright) {
            if (!JSROOT.settings.Zooming) return false;
            // in 3dmode with orbit control ignore simple arrows
            if (this.mode3d && (key.indexOf("Ctrl")!==0)) return false;
            this.analyzeMouseWheelEvent(null, zoom, 0.5);
            this.zoom(zoom.name, zoom.min, zoom.max);
            if (zoom.changed) this.zoomChangedInteractive(zoom.name, true);
            evnt.stopPropagation();
            evnt.preventDefault();
         } else {
            let func = pp && pp.findPadButton ? pp.findPadButton(key) : "";
            if (func) {
               pp.clickPadButton(func);
               evnt.stopPropagation();
               evnt.preventDefault();
            }
         }

         return true; // just process any key press
      },

      /** @summary Function called when frame is clicked and object selection can be performed
        * @desc such event can be used to select */
      processFrameClick: function(pnt, dblckick) {

         let pp = this.getPadPainter();
         if (!pp) return;

         pnt.painters = true; // provide painters reference in the hints
         pnt.disabled = true; // do not invoke graphics

         // collect tooltips from pad painter - it has list of all drawn objects
         let hints = pp.processPadTooltipEvent(pnt), exact = null;
         for (let k=0; (k<hints.length) && !exact; ++k)
            if (hints[k] && hints[k].exact) exact = hints[k];
         //if (exact) console.log('Click exact', pnt, exact.painter.getObjectHint());
         //      else console.log('Click frame', pnt);

         let res;

         if (exact) {
            let handler = dblckick ? this._dblclick_handler : this._click_handler;
            if (handler) res = handler(exact.user_info, pnt);
         }

         if (!dblckick)
            pp.selectObjectPainter(exact ? exact.painter : this,
                  { x: pnt.x + (this._frame_x || 0),  y: pnt.y + (this._frame_y || 0) });

         return res;
      },

      startRectSel: function(evnt) {
         // ignore when touch selection is activated

         if (this.zoom_kind > 100) return;

         // ignore all events from non-left button
         if ((evnt.which || evnt.button) !== 1) return;

         evnt.preventDefault();

         let frame = this.getFrameSvg(),
             pos = d3.pointer(evnt, frame.node());

         this.clearInteractiveElements();

         let w = this.getFrameWidth(), h = this.getFrameHeight();

         this.zoom_lastpos = pos;
         this.zoom_curr = [ Math.max(0, Math.min(w, pos[0])),
                            Math.max(0, Math.min(h, pos[1])) ];

         this.zoom_origin = [0,0];

         if ((pos[0] < 0) || (pos[0] > w)) {
            this.zoom_kind = 3; // only y
            this.zoom_origin[1] = this.zoom_curr[1];
            this.zoom_curr[0] = w;
            this.zoom_curr[1] += 1;
         } else if ((pos[1] < 0) || (pos[1] > h)) {
            this.zoom_kind = 2; // only x
            this.zoom_origin[0] = this.zoom_curr[0];
            this.zoom_curr[0] += 1;
            this.zoom_curr[1] = h;
         } else {
            this.zoom_kind = 1; // x and y
            this.zoom_origin[0] = this.zoom_curr[0];
            this.zoom_origin[1] = this.zoom_curr[1];
         }

         d3.select(window).on("mousemove.zoomRect", this.moveRectSel.bind(this))
                          .on("mouseup.zoomRect", this.endRectSel.bind(this), true);

         this.zoom_rect = null;

         // disable tooltips in frame painter
         setPainterTooltipEnabled(this, false);

         evnt.stopPropagation();

         if (this.zoom_kind != 1)
            setTimeout(() => this.startLabelsMove(), 500);
      },

      startLabelsMove: function() {
         if (this.zoom_rect) return;

         let handle = this.zoom_kind == 2 ? this.x_handle : this.y_handle;

         if (!handle || (typeof handle.processLabelsMove != 'function') || !this.zoom_lastpos) return;

         if (handle.processLabelsMove('start', this.zoom_lastpos)) {
            this.zoom_labels = handle;
         }
      },

      moveRectSel: function(evnt) {

         if ((this.zoom_kind == 0) || (this.zoom_kind > 100)) return;

         evnt.preventDefault();
         let m = d3.pointer(evnt, this.getFrameSvg().node());

         if (this.zoom_labels)
            return this.zoom_labels.processLabelsMove('move', m);

         this.zoom_lastpos[0] = m[0];
         this.zoom_lastpos[1] = m[1];

         m[0] = Math.max(0, Math.min(this.getFrameWidth(), m[0]));
         m[1] = Math.max(0, Math.min(this.getFrameHeight(), m[1]));

         switch (this.zoom_kind) {
            case 1: this.zoom_curr[0] = m[0]; this.zoom_curr[1] = m[1]; break;
            case 2: this.zoom_curr[0] = m[0]; break;
            case 3: this.zoom_curr[1] = m[1]; break;
         }

         let x = Math.min(this.zoom_origin[0], this.zoom_curr[0]),
             y = Math.min(this.zoom_origin[1], this.zoom_curr[1]),
             w = Math.abs(this.zoom_curr[0] - this.zoom_origin[0]),
             h = Math.abs(this.zoom_curr[1] - this.zoom_origin[1]);

         if (!this.zoom_rect) {
            // ignore small changes, can be switching to labels move
            if ((this.zoom_kind != 1) && ((w < 2) || (h < 2))) return;

            this.zoom_rect = this.getFrameSvg()
                                 .append("rect")
                                 .attr("class", "zoom")
                                 .attr("pointer-events","none");
         }

         this.zoom_rect.attr("x", x).attr("y", y).attr("width", w).attr("height", h);
      },

      endRectSel: function(evnt) {
         if ((this.zoom_kind == 0) || (this.zoom_kind > 100)) return;

         evnt.preventDefault();

         d3.select(window).on("mousemove.zoomRect", null)
                          .on("mouseup.zoomRect", null);

         let m = d3.pointer(evnt, this.getFrameSvg().node()), kind = this.zoom_kind;

         if (this.zoom_labels) {
            this.zoom_labels.processLabelsMove('stop', m);
         } else {
            let changed = [true, true];
            m[0] = Math.max(0, Math.min(this.getFrameWidth(), m[0]));
            m[1] = Math.max(0, Math.min(this.getFrameHeight(), m[1]));

            switch (this.zoom_kind) {
               case 1: this.zoom_curr[0] = m[0]; this.zoom_curr[1] = m[1]; break;
               case 2: this.zoom_curr[0] = m[0]; changed[1] = false; break; // only X
               case 3: this.zoom_curr[1] = m[1]; changed[0] = false; break; // only Y
            }

            let xmin, xmax, ymin, ymax, isany = false,
                idx = this.swap_xy ? 1 : 0, idy = 1 - idx;

            if (changed[idx] && (Math.abs(this.zoom_curr[idx] - this.zoom_origin[idx]) > 10)) {
               xmin = Math.min(this.revertAxis("x", this.zoom_origin[idx]), this.revertAxis("x", this.zoom_curr[idx]));
               xmax = Math.max(this.revertAxis("x", this.zoom_origin[idx]), this.revertAxis("x", this.zoom_curr[idx]));
               isany = true;
            }

            if (changed[idy] && (Math.abs(this.zoom_curr[idy] - this.zoom_origin[idy]) > 10)) {
               ymin = Math.min(this.revertAxis("y", this.zoom_origin[idy]), this.revertAxis("y", this.zoom_curr[idy]));
               ymax = Math.max(this.revertAxis("y", this.zoom_origin[idy]), this.revertAxis("y", this.zoom_curr[idy]));
               isany = true;
            }

            if (isany) {
               this.zoomChangedInteractive("x", true);
               this.zoomChangedInteractive("y", true);
               this.zoom(xmin, xmax, ymin, ymax);
               kind = 0;
            }
         }

         let pnt =  (kind===1) ? { x: this.zoom_origin[0], y: this.zoom_origin[1] } : null;

         this.clearInteractiveElements();

         // if no zooming was done, select active object instead
         switch (kind) {
            case 1:
               this.processFrameClick(pnt);
               break;
            case 2: {
               let pp = this.getPadPainter();
               if (pp) pp.selectObjectPainter(this, null, "xaxis");
               break;
            }
            case 3: {
               let pp = this.getPadPainter();
               if (pp) pp.selectObjectPainter(this, null, "yaxis");
               break;
            }
         }

      },

      mouseDoubleClick: function(evnt) {
         evnt.preventDefault();
         let m = d3.pointer(evnt, this.getFrameSvg().node());
         this.clearInteractiveElements();

         let valid_x = (m[0] >= 0) && (m[0] <= this.getFrameWidth()),
             valid_y = (m[1] >= 0) && (m[1] <= this.getFrameHeight());

         if (valid_x && valid_y && this._dblclick_handler)
            if (this.processFrameClick({ x: m[0], y: m[1] }, true)) return;

         let kind = "xyz";
         if (!valid_x) kind = this.swap_xy ? "x" : "y"; else
         if (!valid_y) kind = this.swap_xy ? "y" : "x";
         if (this.unzoom(kind)) return;

         let pp = this.getPadPainter();
         let rect = this.getFrameRect();
         if (pp) pp.selectObjectPainter(pp, { x: m[0] + rect.x, y: m[1] + rect.y, dbl: true });
      },

      startTouchZoom: function(evnt) {
         // in case when zooming was started, block any other kind of events
         if (this.zoom_kind != 0) {
            evnt.preventDefault();
            evnt.stopPropagation();
            return;
         }

         let arr = d3.pointers(evnt, this.getFrameSvg().node());
         this.touch_cnt+=1;

         // normally double-touch will be handled
         // touch with single click used for context menu
         if (arr.length == 1) {
            // this is touch with single element

            let now = new Date(), diff = now.getTime() - this.last_touch.getTime();
            this.last_touch = now;

            if ((diff < 300) && this.zoom_curr
                && (Math.abs(this.zoom_curr[0] - arr[0][0]) < 30)
                && (Math.abs(this.zoom_curr[1] - arr[0][1]) < 30)) {

               evnt.preventDefault();
               evnt.stopPropagation();

               this.clearInteractiveElements();
               this.unzoom("xyz");

               this.last_touch = new Date(0);

               this.getFrameSvg().on("touchcancel", null)
                               .on("touchend", null, true);
            } else if (JSROOT.settings.ContextMenu) {
               this.zoom_curr = arr[0];
               this.getFrameSvg().on("touchcancel", this.endTouchSel.bind(this))
                               .on("touchend", this.endTouchSel.bind(this));
               evnt.preventDefault();
               evnt.stopPropagation();
            }
         }

         if ((arr.length != 2) || !JSROOT.settings.Zooming || !JSROOT.settings.ZoomTouch) return;

         evnt.preventDefault();
         evnt.stopPropagation();

         this.clearInteractiveElements();

         this.getFrameSvg().on("touchcancel", null)
                         .on("touchend", null);

         let pnt1 = arr[0], pnt2 = arr[1], w = this.getFrameWidth(), h = this.getFrameHeight();

         this.zoom_curr = [ Math.min(pnt1[0], pnt2[0]), Math.min(pnt1[1], pnt2[1]) ];
         this.zoom_origin = [ Math.max(pnt1[0], pnt2[0]), Math.max(pnt1[1], pnt2[1]) ];

         if ((this.zoom_curr[0] < 0) || (this.zoom_curr[0] > w)) {
            this.zoom_kind = 103; // only y
            this.zoom_curr[0] = 0;
            this.zoom_origin[0] = w;
         } else if ((this.zoom_origin[1] > h) || (this.zoom_origin[1] < 0)) {
            this.zoom_kind = 102; // only x
            this.zoom_curr[1] = 0;
            this.zoom_origin[1] = h;
         } else {
            this.zoom_kind = 101; // x and y
         }

         setPainterTooltipEnabled(this, false);

         this.zoom_rect = this.getFrameSvg().append("rect")
               .attr("class", "zoom")
               .attr("id", "zoomRect")
               .attr("x", this.zoom_curr[0])
               .attr("y", this.zoom_curr[1])
               .attr("width", this.zoom_origin[0] - this.zoom_curr[0])
               .attr("height", this.zoom_origin[1] - this.zoom_curr[1]);

         d3.select(window).on("touchmove.zoomRect", this.moveTouchSel.bind(this))
                          .on("touchcancel.zoomRect", this.endTouchSel.bind(this))
                          .on("touchend.zoomRect", this.endTouchSel.bind(this));
      },

      moveTouchSel: function(evnt) {
         if (this.zoom_kind < 100) return;

         evnt.preventDefault();

         let arr = d3.pointers(evnt, this.getFrameSvg().node());

         if (arr.length != 2)
            return this.clearInteractiveElements();

         let pnt1 = arr[0], pnt2 = arr[1];

         if (this.zoom_kind != 103) {
            this.zoom_curr[0] = Math.min(pnt1[0], pnt2[0]);
            this.zoom_origin[0] = Math.max(pnt1[0], pnt2[0]);
         }
         if (this.zoom_kind != 102) {
            this.zoom_curr[1] = Math.min(pnt1[1], pnt2[1]);
            this.zoom_origin[1] = Math.max(pnt1[1], pnt2[1]);
         }

         this.zoom_rect.attr("x", this.zoom_curr[0])
                        .attr("y", this.zoom_curr[1])
                        .attr("width", this.zoom_origin[0] - this.zoom_curr[0])
                        .attr("height", this.zoom_origin[1] - this.zoom_curr[1]);

         if ((this.zoom_origin[0] - this.zoom_curr[0] > 10)
              || (this.zoom_origin[1] - this.zoom_curr[1] > 10))
            setPainterTooltipEnabled(this, false);

         evnt.stopPropagation();
      },

      endTouchSel: function(evnt) {

         this.getFrameSvg().on("touchcancel", null)
                         .on("touchend", null);

         if (this.zoom_kind === 0) {
            // special case - single touch can ends up with context menu

            evnt.preventDefault();

            let now = new Date();

            let diff = now.getTime() - this.last_touch.getTime();

            if ((diff > 500) && (diff < 2000) && !this.isTooltipShown()) {
               this.showContextMenu('main', { clientX: this.zoom_curr[0], clientY: this.zoom_curr[1] });
               this.last_touch = new Date(0);
            } else {
               this.clearInteractiveElements();
            }
         }

         if (this.zoom_kind < 100) return;

         evnt.preventDefault();
         d3.select(window).on("touchmove.zoomRect", null)
                          .on("touchend.zoomRect", null)
                          .on("touchcancel.zoomRect", null);

         let xmin, xmax, ymin, ymax, isany = false,
             xid = this.swap_xy ? 1 : 0, yid = 1 - xid,
             changed = [true, true];
         if (this.zoom_kind === 102) changed[1] = false;
         if (this.zoom_kind === 103) changed[0] = false;

         if (changed[xid] && (Math.abs(this.zoom_curr[xid] - this.zoom_origin[xid]) > 10)) {
            xmin = Math.min(this.revertAxis("x", this.zoom_origin[xid]), this.revertAxis("x", this.zoom_curr[xid]));
            xmax = Math.max(this.revertAxis("x", this.zoom_origin[xid]), this.revertAxis("x", this.zoom_curr[xid]));
            isany = true;
         }

         if (changed[yid] && (Math.abs(this.zoom_curr[yid] - this.zoom_origin[yid]) > 10)) {
            ymin = Math.min(this.revertAxis("y", this.zoom_origin[yid]), this.revertAxis("y", this.zoom_curr[yid]));
            ymax = Math.max(this.revertAxis("y", this.zoom_origin[yid]), this.revertAxis("y", this.zoom_curr[yid]));
            isany = true;
         }

         this.clearInteractiveElements();
         this.last_touch = new Date(0);

         if (isany) {
            this.zoomChangedInteractive('x', true);
            this.zoomChangedInteractive('y', true);
            this.zoom(xmin, xmax, ymin, ymax);
         }

         evnt.stopPropagation();
      },

       /** @summary Analyze zooming with mouse wheel */
      analyzeMouseWheelEvent: function(event, item, dmin, test_ignore) {
         let handle = this[item.name + "_handle"];
         if (handle) return handle.analyzeWheelEvent(event, dmin, item, test_ignore);
         console.error('Fail to analyze zooming event for ', item.name);
      },

       /** @summary return true if default Y zooming should be enabled
         * @desc it is typically for 2-Dim histograms or
         * when histogram not draw, defined by other painters */
      isAllowedDefaultYZooming: function() {

         let pad_painter = this.getPadPainter();
         if (pad_painter && pad_painter.painters)
            for (let k = 0; k < pad_painter.painters.length; ++k) {
               let subpainter = pad_painter.painters[k];
               if (subpainter && (subpainter.wheel_zoomy !== undefined))
                  return subpainter.wheel_zoomy;
            }

         return false;
      },

      mouseWheel: function(evnt) {
         evnt.stopPropagation();

         evnt.preventDefault();
         this.clearInteractiveElements();

         let itemx = { name: "x", reverse: this.reverse_x, ignore: false },
             itemy = { name: "y", reverse: this.reverse_y, ignore: !this.isAllowedDefaultYZooming() },
             cur = d3.pointer(evnt, this.getFrameSvg().node()),
             w = this.getFrameWidth(), h = this.getFrameHeight();

         this.analyzeMouseWheelEvent(evnt, this.swap_xy ? itemy : itemx, cur[0] / w, (cur[1] >=0) && (cur[1] <= h));

         this.analyzeMouseWheelEvent(evnt, this.swap_xy ? itemx : itemy, 1 - cur[1] / h, (cur[0] >= 0) && (cur[0] <= w));

         this.zoom(itemx.min, itemx.max, itemy.min, itemy.max);

         if (itemx.changed) this.zoomChangedInteractive('x', true);
         if (itemy.changed) this.zoomChangedInteractive('y', true);
      },

      showContextMenu: function(kind, evnt, obj) {

         // ignore context menu when touches zooming is ongoing
         if (('zoom_kind' in this) && (this.zoom_kind > 100)) return;

         // this is for debug purposes only, when context menu is where, close is and show normal menu
         //if (!evnt && !kind && document.getElementById('root_ctx_menu')) {
         //   let elem = document.getElementById('root_ctx_menu');
         //   elem.parentNode.removeChild(elem);
         //   return;
         //}

         let menu_painter = this, exec_painter = null, frame_corner = false, fp = null; // object used to show context menu

         if (evnt.stopPropagation) {
            evnt.preventDefault();
            evnt.stopPropagation(); // disable main context menu

            if (kind == 'painter' && obj) {
               menu_painter = obj;
               kind = "";
            } else if (!kind) {
               let ms = d3.pointer(evnt, this.getFrameSvg().node()),
                   tch = d3.pointers(evnt, this.getFrameSvg().node()),
                   pp = this.getPadPainter(),
                   pnt = null, sel = null;

               fp = this;

               if (tch.length === 1) pnt = { x: tch[0][0], y: tch[0][1], touch: true }; else
               if (ms.length === 2) pnt = { x: ms[0], y: ms[1], touch: false };

               if ((pnt !== null) && (pp !== null)) {
                  pnt.painters = true; // assign painter for every tooltip
                  let hints = pp.processPadTooltipEvent(pnt), bestdist = 1000;
                  for (let n=0;n<hints.length;++n)
                     if (hints[n] && hints[n].menu) {
                        let dist = ('menu_dist' in hints[n]) ? hints[n].menu_dist : 7;
                        if (dist < bestdist) { sel = hints[n].painter; bestdist = dist; }
                     }
               }

               if (sel) menu_painter = sel; else kind = "frame";

               if (pnt) frame_corner = (pnt.x>0) && (pnt.x<20) && (pnt.y>0) && (pnt.y<20);

               fp.setLastEventPos(pnt);
            } else if (!this.v7_frame && ((kind=="x") || (kind=="y") || (kind=="z"))) {
               exec_painter = this.getMainPainter(); // histogram painter delivers items for axis menu
            }
         } else if (kind == 'painter' && obj) {
            // this is used in 3D context menu to show special painter
            menu_painter = obj;
            kind = "";
         }

         if (!exec_painter) exec_painter = menu_painter;

         if (!menu_painter || !menu_painter.fillContextMenu) return;

         this.clearInteractiveElements();

         jsrp.createMenu(evnt, menu_painter).then(menu => {
            let domenu = menu.painter.fillContextMenu(menu, kind, obj);

            // fill frame menu by default - or append frame elements when activated in the frame corner
            if (fp && (!domenu || (frame_corner && (kind!=="frame"))))
               domenu = fp.fillContextMenu(menu);

            if (domenu)
               exec_painter.fillObjectExecMenu(menu, kind).then(menu => {
                   // suppress any running zooming
                   setPainterTooltipEnabled(menu.painter, false);
                   menu.show().then(() => setPainterTooltipEnabled(menu.painter, true));
               });
         });
      },

     /** @summary Activate context menu handler via touch events
       * @private */
      startTouchMenu: function(kind, evnt) {
         // method to let activate context menu via touch handler

         let arr = d3.pointers(evnt, this.getFrameSvg().node());
         if (arr.length != 1) return;

         if (!kind || (kind=="")) kind = "main";
         let fld = "touch_" + kind;

         evnt.sourceEvent.preventDefault();
         evnt.sourceEvent.stopPropagation();

         this[fld] = { dt: new Date(), pos: arr[0] };

         let handler = this.endTouchMenu.bind(this, kind);

         this.getFrameSvg().on("touchcancel", handler)
                         .on("touchend", handler);
      },

      /** @summary Process end-touch event, which can cause content menu to appear
       * @private */
      endTouchMenu: function(kind, evnt) {
         let fld = "touch_" + kind;

         if (! (fld in this)) return;

         evnt.sourceEvent.preventDefault();
         evnt.sourceEvent.stopPropagation();

         let diff = new Date().getTime() - this[fld].dt.getTime();

         this.getFrameSvg().on("touchcancel", null)
                         .on("touchend", null);

         if (diff > 500) {
            let rect = this.getFrameSvg().node().getBoundingClientRect();
            this.showContextMenu(kind, { clientX: rect.left + this[fld].pos[0],
                                         clientY: rect.top + this[fld].pos[1] } );
         }

         delete this[fld];
      },

      clearInteractiveElements: function() {
         if (jsrp.closeMenu) jsrp.closeMenu();
         this.zoom_kind = 0;
         if (this.zoom_rect) { this.zoom_rect.remove(); delete this.zoom_rect; }
         delete this.zoom_curr;
         delete this.zoom_origin;
         delete this.zoom_lastpos;
         delete this.zoom_labels;


         // enable tooltip in frame painter
         setPainterTooltipEnabled(this, true);
      },

      /** Assign frame interactive methods */
      assign: function(painter) {
         JSROOT.extend(painter, this);

         /*
         painter.addBasicInteractivity = this.addBasicInteractivity;
         painter.addInteractivity = this.addInteractivity;
         painter.addKeysHandler = this.addKeysHandler;
         painter.processKeyPress = this.processKeyPress;
         painter.processFrameClick = this.processFrameClick;
         painter.startRectSel = this.startRectSel;
         painter.moveRectSel = this.moveRectSel;
         painter.endRectSel = this.endRectSel;
         painter.mouseDoubleClick = this.mouseDoubleClick;
         painter.startTouchZoom = this.startTouchZoom;
         painter.moveTouchSel = this.moveTouchSel;
         painter.endTouchSel = this.endTouchSel;
         painter.analyzeMouseWheelEvent = this.analyzeMouseWheelEvent;
         painter.isAllowedDefaultYZooming = this.isAllowedDefaultYZooming;
         painter.mouseWheel = this.mouseWheel;
         painter.showContextMenu = this.showContextMenu;
         painter.startTouchMenu = this.startTouchMenu;
         painter.endTouchMenu = this.endTouchMenu;
         painter.clearInteractiveElements = this.clearInteractiveElements;
         */
      }

   } // FrameInterative

   // some icons taken from http://uxrepo.com/
   let ToolbarIcons = {
      camera: { path: 'M 152.00,304.00c0.00,57.438, 46.562,104.00, 104.00,104.00s 104.00-46.562, 104.00-104.00s-46.562-104.00-104.00-104.00S 152.00,246.562, 152.00,304.00z M 480.00,128.00L 368.00,128.00 c-8.00-32.00-16.00-64.00-48.00-64.00L 192.00,64.00 c-32.00,0.00-40.00,32.00-48.00,64.00L 32.00,128.00 c-17.60,0.00-32.00,14.40-32.00,32.00l0.00,288.00 c0.00,17.60, 14.40,32.00, 32.00,32.00l 448.00,0.00 c 17.60,0.00, 32.00-14.40, 32.00-32.00L 512.00,160.00 C 512.00,142.40, 497.60,128.00, 480.00,128.00z M 256.00,446.00c-78.425,0.00-142.00-63.574-142.00-142.00c0.00-78.425, 63.575-142.00, 142.00-142.00c 78.426,0.00, 142.00,63.575, 142.00,142.00 C 398.00,382.426, 334.427,446.00, 256.00,446.00z M 480.00,224.00l-64.00,0.00 l0.00-32.00 l 64.00,0.00 L 480.00,224.00 z' },
      disk: { path: 'M384,0H128H32C14.336,0,0,14.336,0,32v448c0,17.656,14.336,32,32,32h448c17.656,0,32-14.344,32-32V96L416,0H384z M352,160   V32h32v128c0,17.664-14.344,32-32,32H160c-17.664,0-32-14.336-32-32V32h128v128H352z M96,288c0-17.656,14.336-32,32-32h256   c17.656,0,32,14.344,32,32v192H96V288z' },
      question: { path: 'M256,512c141.375,0,256-114.625,256-256S397.375,0,256,0S0,114.625,0,256S114.625,512,256,512z M256,64   c63.719,0,128,36.484,128,118.016c0,47.453-23.531,84.516-69.891,110.016C300.672,299.422,288,314.047,288,320   c0,17.656-14.344,32-32,32c-17.664,0-32-14.344-32-32c0-40.609,37.25-71.938,59.266-84.031   C315.625,218.109,320,198.656,320,182.016C320,135.008,279.906,128,256,128c-30.812,0-64,20.227-64,64.672   c0,17.664-14.336,32-32,32s-32-14.336-32-32C128,109.086,193.953,64,256,64z M256,449.406c-18.211,0-32.961-14.75-32.961-32.969   c0-18.188,14.75-32.953,32.961-32.953c18.219,0,32.969,14.766,32.969,32.953C288.969,434.656,274.219,449.406,256,449.406z' },
      undo: { path: 'M450.159,48.042c8.791,9.032,16.983,18.898,24.59,29.604c7.594,10.706,14.146,22.207,19.668,34.489  c5.509,12.296,9.82,25.269,12.92,38.938c3.113,13.669,4.663,27.834,4.663,42.499c0,14.256-1.511,28.863-4.532,43.822  c-3.009,14.952-7.997,30.217-14.953,45.795c-6.955,15.577-16.202,31.52-27.755,47.826s-25.88,32.9-42.942,49.807  c-5.51,5.444-11.787,11.67-18.834,18.651c-7.033,6.98-14.496,14.366-22.39,22.168c-7.88,7.802-15.955,15.825-24.187,24.069  c-8.258,8.231-16.333,16.203-24.252,23.888c-18.3,18.13-37.354,37.016-57.191,56.65l-56.84-57.445  c19.596-19.472,38.54-38.279,56.84-56.41c7.75-7.685,15.772-15.604,24.108-23.757s16.438-16.163,24.33-24.057  c7.894-7.893,15.356-15.33,22.402-22.312c7.034-6.98,13.312-13.193,18.821-18.651c22.351-22.402,39.165-44.648,50.471-66.738  c11.279-22.09,16.932-43.567,16.932-64.446c0-15.785-3.217-31.005-9.638-45.671c-6.422-14.665-16.229-28.504-29.437-41.529  c-3.282-3.282-7.358-6.395-12.217-9.325c-4.871-2.938-10.381-5.503-16.516-7.697c-6.121-2.201-12.815-3.992-20.058-5.373  c-7.242-1.374-14.9-2.064-23.002-2.064c-8.218,0-16.802,0.834-25.788,2.507c-8.961,1.674-18.053,4.429-27.222,8.271  c-9.189,3.842-18.456,8.869-27.808,15.089c-9.358,6.219-18.521,13.819-27.502,22.793l-59.92,60.271l93.797,94.058H0V40.91  l93.27,91.597l60.181-60.532c13.376-15.018,27.222-27.248,41.536-36.697c14.308-9.443,28.608-16.776,42.89-21.992  c14.288-5.223,28.505-8.74,42.623-10.557C294.645,0.905,308.189,0,321.162,0c13.429,0,26.389,1.185,38.84,3.562  c12.478,2.377,24.2,5.718,35.192,10.029c11.006,4.311,21.126,9.404,30.374,15.265C434.79,34.724,442.995,41.119,450.159,48.042z' },
      arrow_right: { path: 'M30.796,226.318h377.533L294.938,339.682c-11.899,11.906-11.899,31.184,0,43.084c11.887,11.899,31.19,11.893,43.077,0  l165.393-165.386c5.725-5.712,8.924-13.453,8.924-21.539c0-8.092-3.213-15.84-8.924-21.551L338.016,8.925  C332.065,2.975,324.278,0,316.478,0c-7.802,0-15.603,2.968-21.539,8.918c-11.899,11.906-11.899,31.184,0,43.084l113.391,113.384  H30.796c-16.822,0-30.463,13.645-30.463,30.463C0.333,212.674,13.974,226.318,30.796,226.318z' },
      arrow_up: { path: 'M295.505,629.446V135.957l148.193,148.206c15.555,15.559,40.753,15.559,56.308,0c15.555-15.538,15.546-40.767,0-56.304  L283.83,11.662C276.372,4.204,266.236,0,255.68,0c-10.568,0-20.705,4.204-28.172,11.662L11.333,227.859  c-7.777,7.777-11.666,17.965-11.666,28.158c0,10.192,3.88,20.385,11.657,28.158c15.563,15.555,40.762,15.555,56.317,0  l148.201-148.219v493.489c0,21.993,17.837,39.82,39.82,39.82C277.669,669.267,295.505,651.439,295.505,629.446z' },
      arrow_diag: { path: 'M279.875,511.994c-1.292,0-2.607-0.102-3.924-0.312c-10.944-1.771-19.333-10.676-20.457-21.71L233.97,278.348  L22.345,256.823c-11.029-1.119-19.928-9.51-21.698-20.461c-1.776-10.944,4.031-21.716,14.145-26.262L477.792,2.149  c9.282-4.163,20.167-2.165,27.355,5.024c7.201,7.189,9.199,18.086,5.024,27.356L302.22,497.527  C298.224,506.426,289.397,511.994,279.875,511.994z M118.277,217.332l140.534,14.294c11.567,1.178,20.718,10.335,21.878,21.896  l14.294,140.519l144.09-320.792L118.277,217.332z' },
      auto_zoom: { path: 'M505.441,242.47l-78.303-78.291c-9.18-9.177-24.048-9.171-33.216,0c-9.169,9.172-9.169,24.045,0.006,33.217l38.193,38.188  H280.088V80.194l38.188,38.199c4.587,4.584,10.596,6.881,16.605,6.881c6.003,0,12.018-2.297,16.605-6.875  c9.174-9.172,9.174-24.039,0.011-33.217L273.219,6.881C268.803,2.471,262.834,0,256.596,0c-6.229,0-12.202,2.471-16.605,6.881  l-78.296,78.302c-9.178,9.172-9.178,24.045,0,33.217c9.177,9.171,24.051,9.171,33.21,0l38.205-38.205v155.4H80.521l38.2-38.188  c9.177-9.171,9.177-24.039,0.005-33.216c-9.171-9.172-24.039-9.178-33.216,0L7.208,242.464c-4.404,4.403-6.881,10.381-6.881,16.611  c0,6.227,2.477,12.207,6.881,16.61l78.302,78.291c4.587,4.581,10.599,6.875,16.605,6.875c6.006,0,12.023-2.294,16.61-6.881  c9.172-9.174,9.172-24.036-0.005-33.211l-38.205-38.199h152.593v152.063l-38.199-38.211c-9.171-9.18-24.039-9.18-33.216-0.022  c-9.178,9.18-9.178,24.059-0.006,33.222l78.284,78.302c4.41,4.404,10.382,6.881,16.611,6.881c6.233,0,12.208-2.477,16.611-6.881  l78.302-78.296c9.181-9.18,9.181-24.048,0-33.205c-9.174-9.174-24.054-9.174-33.21,0l-38.199,38.188v-152.04h152.051l-38.205,38.199  c-9.18,9.175-9.18,24.037-0.005,33.211c4.587,4.587,10.596,6.881,16.604,6.881c6.01,0,12.024-2.294,16.605-6.875l78.303-78.285  c4.403-4.403,6.887-10.378,6.887-16.611C512.328,252.851,509.845,246.873,505.441,242.47z' },
      statbox: {
         path: 'M28.782,56.902H483.88c15.707,0,28.451-12.74,28.451-28.451C512.331,12.741,499.599,0,483.885,0H28.782   C13.074,0,0.331,12.741,0.331,28.451C0.331,44.162,13.074,56.902,28.782,56.902z' +
            'M483.885,136.845H28.782c-15.708,0-28.451,12.741-28.451,28.451c0,15.711,12.744,28.451,28.451,28.451H483.88   c15.707,0,28.451-12.74,28.451-28.451C512.331,149.586,499.599,136.845,483.885,136.845z' +
            'M483.885,273.275H28.782c-15.708,0-28.451,12.731-28.451,28.452c0,15.707,12.744,28.451,28.451,28.451H483.88   c15.707,0,28.451-12.744,28.451-28.451C512.337,286.007,499.599,273.275,483.885,273.275z' +
            'M256.065,409.704H30.492c-15.708,0-28.451,12.731-28.451,28.451c0,15.707,12.744,28.451,28.451,28.451h225.585   c15.707,0,28.451-12.744,28.451-28.451C284.516,422.436,271.785,409.704,256.065,409.704z'
      },
      circle: { path: "M256,256 m-150,0 a150,150 0 1,0 300,0 a150,150 0 1,0 -300,0" },
      three_circles: { path: "M256,85 m-70,0 a70,70 0 1,0 140,0 a70,70 0 1,0 -140,0  M256,255 m-70,0 a70,70 0 1,0 140,0 a70,70 0 1,0 -140,0  M256,425 m-70,0 a70,70 0 1,0 140,0 a70,70 0 1,0 -140,0 " },
      diamand: { path: "M256,0L384,256L256,511L128,256z" },
      rect: { path: "M80,80h352v352h-352z" },
      cross: { path: "M80,40l176,176l176,-176l40,40l-176,176l176,176l-40,40l-176,-176l-176,176l-40,-40l176,-176l-176,-176z" },
      vrgoggles: { size: "245.82 141.73", path: 'M175.56,111.37c-22.52,0-40.77-18.84-40.77-42.07S153,27.24,175.56,27.24s40.77,18.84,40.77,42.07S198.08,111.37,175.56,111.37ZM26.84,69.31c0-23.23,18.25-42.07,40.77-42.07s40.77,18.84,40.77,42.07-18.26,42.07-40.77,42.07S26.84,92.54,26.84,69.31ZM27.27,0C11.54,0,0,12.34,0,28.58V110.9c0,16.24,11.54,30.83,27.27,30.83H99.57c2.17,0,4.19-1.83,5.4-3.7L116.47,118a8,8,0,0,1,12.52-.18l11.51,20.34c1.2,1.86,3.22,3.61,5.39,3.61h72.29c15.74,0,27.63-14.6,27.63-30.83V28.58C245.82,12.34,233.93,0,218.19,0H27.27Z' },
      th2colorz: { recs: [{ x: 128, y: 486, w: 256, h: 26, f: 'rgb(38,62,168)' }, { y: 461, f: 'rgb(22,82,205)' }, { y: 435, f: 'rgb(16,100,220)' }, { y: 410, f: 'rgb(18,114,217)' }, { y: 384, f: 'rgb(20,129,214)' }, { y: 358, f: 'rgb(14,143,209)' }, { y: 333, f: 'rgb(9,157,204)' }, { y: 307, f: 'rgb(13,167,195)' }, { y: 282, f: 'rgb(30,175,179)' }, { y: 256, f: 'rgb(46,183,164)' }, { y: 230, f: 'rgb(82,186,146)' }, { y: 205, f: 'rgb(116,189,129)' }, { y: 179, f: 'rgb(149,190,113)' }, { y: 154, f: 'rgb(179,189,101)' }, { y: 128, f: 'rgb(209,187,89)' }, { y: 102, f: 'rgb(226,192,75)' }, { y: 77, f: 'rgb(244,198,59)' }, { y: 51, f: 'rgb(253,210,43)' }, { y: 26, f: 'rgb(251,230,29)' }, { y: 0, f: 'rgb(249,249,15)' }] },
      th2color: { recs: [{x:0,y:256,w:13,h:39,f:'rgb(38,62,168)'},{x:13,y:371,w:39,h:39},{y:294,h:39},{y:256,h:39},{y:218,h:39},{x:51,y:410,w:39,h:39},{y:371,h:39},{y:333,h:39},{y:294},{y:256,h:39},{y:218,h:39},{y:179,h:39},{y:141,h:39},{y:102,h:39},{y:64},{x:90,y:448,w:39,h:39},{y:410},{y:371,h:39},{y:333,h:39,f:'rgb(22,82,205)'},{y:294},{y:256,h:39,f:'rgb(16,100,220)'},{y:218,h:39},{y:179,h:39,f:'rgb(22,82,205)'},{y:141,h:39},{y:102,h:39,f:'rgb(38,62,168)'},{y:64},{y:0,h:27},{x:128,y:448,w:39,h:39},{y:410},{y:371,h:39},{y:333,h:39,f:'rgb(22,82,205)'},{y:294,f:'rgb(20,129,214)'},{y:256,h:39,f:'rgb(9,157,204)'},{y:218,h:39,f:'rgb(14,143,209)'},{y:179,h:39,f:'rgb(20,129,214)'},{y:141,h:39,f:'rgb(16,100,220)'},{y:102,h:39,f:'rgb(22,82,205)'},{y:64,f:'rgb(38,62,168)'},{y:26,h:39},{y:0,h:27},{x:166,y:486,h:14},{y:448,h:39},{y:410},{y:371,h:39,f:'rgb(22,82,205)'},{y:333,h:39,f:'rgb(20,129,214)'},{y:294,f:'rgb(82,186,146)'},{y:256,h:39,f:'rgb(179,189,101)'},{y:218,h:39,f:'rgb(116,189,129)'},{y:179,h:39,f:'rgb(82,186,146)'},{y:141,h:39,f:'rgb(14,143,209)'},{y:102,h:39,f:'rgb(16,100,220)'},{y:64,f:'rgb(38,62,168)'},{y:26,h:39},{x:205,y:486,w:39,h:14},{y:448,h:39},{y:410},{y:371,h:39,f:'rgb(16,100,220)'},{y:333,h:39,f:'rgb(9,157,204)'},{y:294,f:'rgb(149,190,113)'},{y:256,h:39,f:'rgb(244,198,59)'},{y:218,h:39},{y:179,h:39,f:'rgb(226,192,75)'},{y:141,h:39,f:'rgb(13,167,195)'},{y:102,h:39,f:'rgb(18,114,217)'},{y:64,f:'rgb(22,82,205)'},{y:26,h:39,f:'rgb(38,62,168)'},{x:243,y:448,w:39,h:39},{y:410},{y:371,h:39,f:'rgb(18,114,217)'},{y:333,h:39,f:'rgb(30,175,179)'},{y:294,f:'rgb(209,187,89)'},{y:256,h:39,f:'rgb(251,230,29)'},{y:218,h:39,f:'rgb(249,249,15)'},{y:179,h:39,f:'rgb(226,192,75)'},{y:141,h:39,f:'rgb(30,175,179)'},{y:102,h:39,f:'rgb(18,114,217)'},{y:64,f:'rgb(38,62,168)'},{y:26,h:39},{x:282,y:448,h:39},{y:410},{y:371,h:39,f:'rgb(18,114,217)'},{y:333,h:39,f:'rgb(14,143,209)'},{y:294,f:'rgb(149,190,113)'},{y:256,h:39,f:'rgb(226,192,75)'},{y:218,h:39,f:'rgb(244,198,59)'},{y:179,h:39,f:'rgb(149,190,113)'},{y:141,h:39,f:'rgb(9,157,204)'},{y:102,h:39,f:'rgb(18,114,217)'},{y:64,f:'rgb(38,62,168)'},{y:26,h:39},{x:320,y:448,w:39,h:39},{y:410},{y:371,h:39,f:'rgb(22,82,205)'},{y:333,h:39,f:'rgb(20,129,214)'},{y:294,f:'rgb(46,183,164)'},{y:256,h:39},{y:218,h:39,f:'rgb(82,186,146)'},{y:179,h:39,f:'rgb(9,157,204)'},{y:141,h:39,f:'rgb(20,129,214)'},{y:102,h:39,f:'rgb(16,100,220)'},{y:64,f:'rgb(38,62,168)'},{y:26,h:39},{x:358,y:448,h:39},{y:410},{y:371,h:39,f:'rgb(22,82,205)'},{y:333,h:39},{y:294,f:'rgb(16,100,220)'},{y:256,h:39,f:'rgb(20,129,214)'},{y:218,h:39,f:'rgb(14,143,209)'},{y:179,h:39,f:'rgb(18,114,217)'},{y:141,h:39,f:'rgb(22,82,205)'},{y:102,h:39,f:'rgb(38,62,168)'},{y:64},{y:26,h:39},{x:397,y:448,w:39,h:39},{y:371,h:39},{y:333,h:39},{y:294,f:'rgb(22,82,205)'},{y:256,h:39},{y:218,h:39},{y:179,h:39,f:'rgb(38,62,168)'},{y:141,h:39},{y:102,h:39},{y:64},{y:26,h:39},{x:435,y:410,h:39},{y:371,h:39},{y:333,h:39},{y:294},{y:256,h:39},{y:218,h:39},{y:179,h:39},{y:141,h:39},{y:102,h:39},{y:64},{x:474,y:256,h:39},{y:179,h:39}] },
      th2draw3d: {
         path: "M172.768,0H51.726C23.202,0,0.002,23.194,0.002,51.712v89.918c0,28.512,23.2,51.718,51.724,51.718h121.042   c28.518,0,51.724-23.2,51.724-51.718V51.712C224.486,23.194,201.286,0,172.768,0z M177.512,141.63c0,2.611-2.124,4.745-4.75,4.745   H51.726c-2.626,0-4.751-2.134-4.751-4.745V51.712c0-2.614,2.125-4.739,4.751-4.739h121.042c2.62,0,4.75,2.125,4.75,4.739 L177.512,141.63L177.512,141.63z "+
               "M460.293,0H339.237c-28.521,0-51.721,23.194-51.721,51.712v89.918c0,28.512,23.2,51.718,51.721,51.718h121.045   c28.521,0,51.721-23.2,51.721-51.718V51.712C512.002,23.194,488.802,0,460.293,0z M465.03,141.63c0,2.611-2.122,4.745-4.748,4.745   H339.237c-2.614,0-4.747-2.128-4.747-4.745V51.712c0-2.614,2.133-4.739,4.747-4.739h121.045c2.626,0,4.748,2.125,4.748,4.739 V141.63z "+
               "M172.768,256.149H51.726c-28.524,0-51.724,23.205-51.724,51.726v89.915c0,28.504,23.2,51.715,51.724,51.715h121.042   c28.518,0,51.724-23.199,51.724-51.715v-89.915C224.486,279.354,201.286,256.149,172.768,256.149z M177.512,397.784   c0,2.615-2.124,4.736-4.75,4.736H51.726c-2.626-0.006-4.751-2.121-4.751-4.736v-89.909c0-2.626,2.125-4.753,4.751-4.753h121.042 c2.62,0,4.75,2.116,4.75,4.753L177.512,397.784L177.512,397.784z "+
               "M460.293,256.149H339.237c-28.521,0-51.721,23.199-51.721,51.726v89.915c0,28.504,23.2,51.715,51.721,51.715h121.045   c28.521,0,51.721-23.199,51.721-51.715v-89.915C512.002,279.354,488.802,256.149,460.293,256.149z M465.03,397.784   c0,2.615-2.122,4.736-4.748,4.736H339.237c-2.614,0-4.747-2.121-4.747-4.736v-89.909c0-2.626,2.121-4.753,4.747-4.753h121.045 c2.615,0,4.748,2.116,4.748,4.753V397.784z"
      },

      createSVG: function(group, btn, size, title) {
         let svg = group.append("svg:svg")
                        .attr("class", "svg_toolbar_btn")
                        .attr("width", size + "px")
                        .attr("height", size + "px")
                        .attr("viewBox", "0 0 512 512")
                        .style("overflow", "hidden");

         if ('recs' in btn) {
            let rec = {};
            for (let n = 0; n < btn.recs.length; ++n) {
               JSROOT.extend(rec, btn.recs[n]);
               svg.append('rect').attr("x", rec.x).attr("y", rec.y)
                  .attr("width", rec.w).attr("height", rec.h)
                  .attr("fill", rec.f);
            }
         } else {
            svg.append('svg:path').attr('d', btn.path);
         }

         //  special rect to correctly get mouse events for whole button area
         svg.append("svg:rect").attr("x", 0).attr("y", 0).attr("width", 512).attr("height", 512)
            .style('opacity', 0).style('fill', "none").style("pointer-events", "visibleFill")
            .append("svg:title").text(title);

         return svg;
      }
   } // ToolbarIcons


   function getButtonSize(handler, fact) {
      return Math.round((fact || 1) * (handler.iscan || !handler.has_canvas ? 16 : 12));
   }

   function toggleButtonsVisibility(handler, action) {
      let group = handler.getLayerSvg("btns_layer"),
          btn = group.select("[name='Toggle']");

      if (btn.empty()) return;

      let state = btn.property('buttons_state');

      if (btn.property('timout_handler')) {
         if (action!=='timeout') clearTimeout(btn.property('timout_handler'));
         btn.property('timout_handler', null);
      }

      let is_visible = false;
      switch(action) {
         case 'enable': is_visible = true; break;
         case 'enterbtn': return; // do nothing, just cleanup timeout
         case 'timeout': is_visible = false; break;
         case 'toggle':
            state = !state;
            btn.property('buttons_state', state);
            is_visible = state;
            break;
         case 'disable':
         case 'leavebtn':
            if (!state) btn.property('timout_handler', setTimeout(() => toggleButtonsVisibility(handler, 'timeout'), 1200));
            return;
      }

      group.selectAll('svg').each(function() {
         if (this===btn.node()) return;
         d3.select(this).style('display', is_visible ? "" : "none");
      });
   }

   let PadButtonsHandler = {

      alignButtons:  function(btns, width, height) {
         let sz0 = getButtonSize(this, 1.25), nextx = (btns.property('nextx') || 0) + sz0, btns_x, btns_y;

         if (btns.property('vertical')) {
            btns_x = btns.property('leftside') ? 2 : (width - sz0);
            btns_y = height - nextx;
         } else {
            btns_x = btns.property('leftside') ? 2 : (width - nextx);
            btns_y = height - sz0;
         }

         btns.attr("transform","translate("+btns_x+","+btns_y+")");
      },

      findPadButton: function(keyname) {
         let group = this.getLayerSvg("btns_layer"), found_func = "";
         if (!group.empty())
            group.selectAll("svg").each(function() {
               if (d3.select(this).attr("key") === keyname)
                  found_func = d3.select(this).attr("name");
            });

         return found_func;
      },

      removePadButtons: function() {
         let group = this.getLayerSvg("btns_layer");
         if (!group.empty()) {
            group.selectAll("*").remove();
            group.property("nextx", null);
         }
      },

      showPadButtons: function() {
         let group = this.getLayerSvg("btns_layer");
         if (group.empty()) return;

         // clean all previous buttons
         group.selectAll("*").remove();

         let iscan = this.iscan || !this.has_canvas, ctrl,
             x = group.property('leftside') ? getButtonSize(this, 1.25) : 0, y = 0;

         if (this._fast_drawing) {
            ctrl = ToolbarIcons.createSVG(group, ToolbarIcons.circle, getButtonSize(this), "enlargePad");
            ctrl.attr("name", "Enlarge").attr("x", 0).attr("y", 0)
                .on("click", evnt => this.clickPadButton("enlargePad", evnt));
         } else {
            ctrl = ToolbarIcons.createSVG(group, ToolbarIcons.rect, getButtonSize(this), "Toggle tool buttons");

            ctrl.attr("name", "Toggle").attr("x", 0).attr("y", 0)
                .property("buttons_state", (JSROOT.settings.ToolBar!=='popup'))
                .on("click", () => toggleButtonsVisibility(this, 'toggle'))
                .on("mouseenter", () => toggleButtonsVisibility(this, 'enable'))
                .on("mouseleave", () => toggleButtonsVisibility(this, 'disable'));

            for (let k=0;k<this._buttons.length;++k) {
               let item = this._buttons[k];

               let btn = item.btn;
               if (typeof btn == 'string') btn = ToolbarIcons[btn];
               if (!btn) btn = ToolbarIcons.circle;

               let svg = ToolbarIcons.createSVG(group, btn, getButtonSize(this),
                           item.tooltip + (iscan ? "" : (" on pad " + this.this_pad_name)) + (item.keyname ? " (keyshortcut " + item.keyname + ")" : ""));

               if (group.property('vertical'))
                   svg.attr("x", y).attr("y", x);
               else
                  svg.attr("x", x).attr("y", y);

               svg.attr("name", item.funcname)
                  .style('display', (ctrl.property("buttons_state") ? '' : 'none'))
                  .on("mouseenter", () => toggleButtonsVisibility(this, 'enterbtn'))
                  .on("mouseleave", () => toggleButtonsVisibility(this, 'leavebtn'));

               if (item.keyname) svg.attr("key", item.keyname);

               svg.on("click", evnt => this.clickPadButton(item.funcname, evnt));

               x += getButtonSize(this, 1.25);
            }
         }

         group.property("nextx", x);

         this.alignButtons(group, this.getPadWidth(), this.getPadHeight());

         if (group.property('vertical'))
            ctrl.attr("y", x);
         else if (!group.property('leftside'))
            ctrl.attr("x", x);
      },

      assign: function(painter) {
         painter.alignButtons = this.alignButtons;
         painter.findPadButton = this.findPadButton;
         painter.removePadButtons = this.removePadButtons;
         painter.showPadButtons = this.showPadButtons;

      }
   } // PadButtonsHandler

   return {
      TooltipHandler: TooltipHandler,
      addDragHandler: addDragHandler,
      addMoveHandler: addMoveHandler,
      FrameInteractive: FrameInteractive,
      ToolbarIcons: ToolbarIcons,
      PadButtonsHandler: PadButtonsHandler
   };

})
