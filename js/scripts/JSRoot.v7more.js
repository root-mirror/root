/// @file JSRoot.v7more.js
/// JavaScript ROOT v7 graphics for different classes

JSROOT.define(['painter', 'v7gpad'], (jsrp) => {

   "use strict";

   function drawText() {
      let text      = this.getObject(),
          pp        = this.getPadPainter(),
          use_frame = false,
          p         = pp.getCoordinate(text.fPos),
          textFont  = this.v7EvalFont("text", { size: 12, color: "black", align: 22 });

      this.createG(use_frame);

      this.startTextDrawing(textFont, 'font');

      this.drawText({ x: p.x, y: p.y, text: text.fText, latex: 1 });

      return this.finishTextDrawing();
   }


   // =================================================================================

   function drawLine() {

       let line         = this.getObject(),
           pp           = this.getPadPainter(),
           p1           = pp.getCoordinate(line.fP1),
           p2           = pp.getCoordinate(line.fP2),
           line_width   = this.v7EvalAttr("line_width", 1),
           line_style   = this.v7EvalAttr("line_style", 1),
           line_color   = this.v7EvalColor("line_color", "black");

       this.createG();

       this.draw_g
           .append("svg:line")
           .attr("x1", p1.x)
           .attr("y1", p1.y)
           .attr("x2", p2.x)
           .attr("y2", p2.y)
           .style("stroke", line_color)
           .attr("stroke-width", line_width)
//        .attr("stroke-opacity", line_opacity)
           .style("stroke-dasharray", jsrp.root_line_styles[line_style]);
   }

   // =================================================================================

   function drawBox() {

       let box          = this.getObject(),
           pp           = this.getPadPainter(),
           p1           = pp.getCoordinate(box.fP1),
           p2           = pp.getCoordinate(box.fP2),
           line_width   = this.v7EvalAttr( "box_border_width", 1),
           line_style   = this.v7EvalAttr( "box_border_style", 1),
           line_color   = this.v7EvalColor( "box_border_color", "black"),
           fill_color   = this.v7EvalColor( "box_fill_color", "white"),
           fill_style   = this.v7EvalAttr( "box_fill_style", 1),
           round_width  = this.v7EvalAttr( "box_round_width", 0), // not yet exists
           round_height = this.v7EvalAttr( "box_round_height", 0); // not yet exists

    this.createG();

    if (fill_style == 0) fill_color = "none";

    this.draw_g
        .append("svg:rect")
        .attr("x", p1.x)
        .attr("width", p2.x-p1.x)
        .attr("y", p2.y)
        .attr("height", p1.y-p2.y)
        .attr("rx", round_width)
        .attr("ry", round_height)
        .style("stroke", line_color)
        .attr("stroke-width", line_width)
        .attr("fill", fill_color)
        .style("stroke-dasharray", jsrp.root_line_styles[line_style]);
   }

   // =================================================================================

   function drawMarker() {
       let marker       = this.getObject(),
           pp           = this.getPadPainter(),
           p            = pp.getCoordinate(marker.fP),
           marker_size  = this.v7EvalAttr( "marker_size", 1),
           marker_style = this.v7EvalAttr( "marker_style", 1),
           marker_color = this.v7EvalColor( "marker_color", "black"),
           att          = new JSROOT.TAttMarkerHandler({ style: marker_style, color: marker_color, size: marker_size }),
           path         = att.create(p.x, p.y);

       this.createG();

       if (path)
          this.draw_g.append("svg:path")
                     .attr("d", path)
                     .call(att.func);
   }

   // =================================================================================

   function drawLegendContent() {
      let legend     = this.getObject(),
          textFont  = this.v7EvalFont("legend_text", { size: 12, color: "black", align: 22 }),
          width      = this.pave_width,
          height     = this.pave_height,
          nlines     = legend.fEntries.length,
          pp         = this.getPadPainter();

      if (legend.fTitle) nlines++;

      if (!nlines || !pp) return Promise.resolve(this);

      let stepy = height / nlines, posy = 0, margin_x = 0.02 * width;

      textFont.setSize(height/(nlines * 1.2));
      this.startTextDrawing(textFont, 'font' );

      if (legend.fTitle) {
         this.drawText({ align: 22, latex: 1,
                         width: width - 2*margin_x, height: stepy, x: margin_x, y: posy, text: legend.fTitle });
         posy += stepy;
      }

      for (let i=0; i<legend.fEntries.length; ++i) {
         let objp = null, entry = legend.fEntries[i];

         this.drawText({ latex: 1, width: 0.75*width - 3*margin_x, height: stepy, x: 2*margin_x + width*0.25, y: posy, text: entry.fLabel });

         if (entry.fDrawableId != "custom") {
            objp = pp.findSnap(entry.fDrawableId, true);
         } else if (entry.fDrawable.fIO) {
            objp = new JSROOT.ObjectPainter(this.getDom(), entry.fDrawable.fIO);
            if (entry.fLine) objp.createv7AttLine();
            if (entry.fFill) objp.createv7AttFill();
            if (entry.fMarker) objp.createv7AttMarker();
         }

         if (objp && entry.fFill && objp.fillatt)
            this.draw_g
              .append("svg:rect")
              .attr("x", Math.round(margin_x))
              .attr("y", Math.round(posy + stepy*0.1))
              .attr("width", Math.round(width/4))
              .attr("height", Math.round(stepy*0.8))
              .call(objp.fillatt.func);

         if (objp && entry.fLine && objp.lineatt)
            this.draw_g
              .append("svg:line")
              .attr("x1", Math.round(margin_x))
              .attr("y1", Math.round(posy + stepy/2))
              .attr("x2", Math.round(margin_x + width/4))
              .attr("y2", Math.round(posy + stepy/2))
              .call(objp.lineatt.func);

         if (objp && entry.fMarker && objp.markeratt)
            this.draw_g.append("svg:path")
                .attr("d", objp.markeratt.create(margin_x + width/8, posy + stepy/2))
                .call(objp.markeratt.func);

         posy += stepy;
      }

      return this.finishTextDrawing();
   }

   function drawLegend(divid, legend, opt) {
      let painter = new JSROOT.v7.RPavePainter(divid, legend, opt, "legend");

      painter.drawContent = drawLegendContent;

     return jsrp.ensureRCanvas(painter, false).then(() => painter.drawPave());
   }

   // =================================================================================

   function drawPaveTextContent() {
      let pavetext  = this.getObject(),
          textFont  = this.v7EvalFont("pavetext_text", { size: 12, color: "black", align: 22 }),
          width     = this.pave_width,
          height    = this.pave_height,
          nlines    = pavetext.fText.length;

      if (!nlines) return;

      let stepy = height / nlines, posy = 0, margin_x = 0.02 * width;

      textFont.setSize(height/(nlines * 1.2))

      this.startTextDrawing(textFont, 'font');

      for (let i=0; i < pavetext.fText.length; ++i) {
         let line = pavetext.fText[i];

         this.drawText({ latex: 1, width: width - 2*margin_x, height: stepy, x: margin_x, y: posy, text: line });
         posy += stepy;
      }

      return this.finishTextDrawing();
   }

   function drawPaveText(divid, pave, opt) {
      let painter = new JSROOT.v7.RPavePainter(divid, pave, opt, "pavetext");

      painter.drawContent = drawPaveTextContent;

      return jsrp.ensureRCanvas(painter, false).then(() => painter.drawPave());
   }


   // ================================================================================

   JSROOT.v7.drawText     = drawText;
   JSROOT.v7.drawLine     = drawLine;
   JSROOT.v7.drawBox      = drawBox;
   JSROOT.v7.drawMarker   = drawMarker;
   JSROOT.v7.drawLegend   = drawLegend;
   JSROOT.v7.drawPaveText = drawPaveText;

   return JSROOT;

});
