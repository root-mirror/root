/// @file JSRoot.jq2d.js
/// Part of JavaScript ROOT, dependent from jQuery functionality

JSROOT.define(['d3', 'jquery', 'painter', 'hierarchy', 'jquery-ui', 'jqueryui-mousewheel', 'jqueryui-touch-punch'], (d3, $, jsrp) => {

   "use strict";

   JSROOT.loadScript('$$$style/jquery-ui');

   if (typeof jQuery === 'undefined') globalThis.jQuery = $;

   let BrowserLayout = JSROOT.BrowserLayout;

   /** @summary Set browser title text
     * @desc Title also used for dragging of the float browser */
   BrowserLayout.prototype.setBrowserTitle = function(title) {
      let main = d3.select("#" + this.gui_div + " .jsroot_browser");
      if (!main.empty())
         main.select(".jsroot_browser_title").text(title);
   }

   /** @summary Toggle browser kind */
   BrowserLayout.prototype.toggleBrowserKind = function(kind) {

      if (!this.gui_div) return;

      if (!kind) {
         if (!this.browser_kind) return;
         kind = (this.browser_kind === "float") ? "fix" : "float";
      }

      let main = d3.select("#"+this.gui_div+" .jsroot_browser"),
          jmain = $(main.node()),
          area = jmain.find(".jsroot_browser_area"),
          pthis = this;

      if (this.browser_kind === "float") {
          area.css('bottom', '0px')
              .css('top', '0px')
              .css('width','').css('height','')
              .toggleClass('jsroot_float_browser', false)
              .resizable("destroy")
              .draggable("destroy");
      } else if (this.browser_kind === "fix") {
         main.select(".jsroot_v_separator").remove();
         area.css('left', '0px');
         d3.select("#"+this.gui_div+"_drawing").style('left','0px'); // reset size
         main.select(".jsroot_h_separator").style('left','0px');
         d3.select("#"+this.gui_div+"_status").style('left','0px'); // reset left
         pthis.checkResize();
      }

      this.browser_kind = kind;
      this.browser_visible = true;

      if (kind==="float") {
         area.css('bottom', '40px')
           .toggleClass('jsroot_float_browser', true)
           .resizable({
              containment: "parent",
              minWidth: 100,
              resize: function(/* event, ui */) {
                 pthis.setButtonsPosition();
              },
              stop: function( event, ui ) {
                 let bottom = $(this).parent().innerHeight() - ui.position.top - ui.size.height;
                 if (bottom<7) $(this).css('height', "").css('bottom', 0);
              }
         })
         .draggable({
             containment: "parent",
             handle : $("#"+this.gui_div).find(".jsroot_browser_title"),
             snap: true,
             snapMode: "inner",
             snapTolerance: 10,
             drag: function(/* event, ui */) {
                pthis.setButtonsPosition();
             },
             stop: function(/* event, ui */) {
                let bottom = $(this).parent().innerHeight() - $(this).offset().top - $(this).outerHeight();
                if (bottom<7) $(this).css('height', "").css('bottom', 0);
             }
          });
         this.adjustBrowserSize();

     } else {

        area.css('left',0).css('top',0).css('bottom',0).css('height','');

        let vsepar =
           main.append('div')
               .classed("jsroot_separator", true).classed('jsroot_v_separator', true)
               .style('position', 'absolute').style('top',0).style('bottom',0);
        // creation of vertical separator
        $(vsepar.node()).draggable({
           axis: "x" , cursor: "ew-resize",
           containment: "parent",
           helper : function() { return $(this).clone().css('background-color','grey'); },
           drag: function(event,ui) {
              pthis.setButtonsPosition();
              pthis.adjustSeparators(ui.position.left, null);
           },
           stop: function(/* event,ui */) {
              pthis.checkResize();
           }
        });

        this.adjustSeparators(250, null, true, true);
     }

      this.setButtonsPosition();
   }

   /** @summary Set buttons position */
   BrowserLayout.prototype.setButtonsPosition = function() {
      if (!this.gui_div) return;

      let jmain = $("#"+this.gui_div+" .jsroot_browser"),
          btns = jmain.find(".jsroot_browser_btns"),
          top = 7, left = 7;

      if (!btns.length) return;

      if (this.browser_visible) {
         let area = jmain.find(".jsroot_browser_area"),
             off0 = jmain.offset(), off1 = area.offset();
         top = off1.top - off0.top + 7;
         left = off1.left - off0.left + area.innerWidth() - 27;
      }

      btns.css('left', left+'px').css('top', top+'px');
   }

   /** @summary Adjust browser size */
   BrowserLayout.prototype.adjustBrowserSize = function(onlycheckmax) {
      if (!this.gui_div || (this.browser_kind !== "float")) return;

      let jmain = $("#" + this.gui_div + " .jsroot_browser");
      if (!jmain.length) return;

      let area = jmain.find(".jsroot_browser_area"),
          cont = jmain.find(".jsroot_browser_hierarchy"),
          chld = cont.children(":first");

      if (onlycheckmax) {
         if (area.parent().innerHeight() - 10 < area.innerHeight())
            area.css('bottom', '0px').css('top','0px');
         return;
      }

      if (!chld.length) return;

      let h1 = cont.innerHeight(),
          h2 = chld.innerHeight();

      if ((h2!==undefined) && (h2<h1*0.7)) area.css('bottom', '');
   }

   /** @summary Toggle browser visibility */
   BrowserLayout.prototype.toggleBrowserVisisbility = function(fast_close) {
      if (!this.gui_div || (typeof this.browser_visible==='string')) return;

      let main = d3.select("#" + this.gui_div + " .jsroot_browser"),
          area = main.select('.jsroot_browser_area');

      if (area.empty()) return;

      let vsepar = main.select(".jsroot_v_separator"),
          drawing = d3.select("#" + this.gui_div + "_drawing"),
          tgt = area.property('last_left'),
          tgt_separ = area.property('last_vsepar'),
          tgt_drawing = area.property('last_drawing');

      if (!this.browser_visible) {
         if (fast_close) return;
         area.property('last_left', null).property('last_vsepar',null).property('last_drawing', null);
      } else {
         area.property('last_left', area.style('left'));
         if (!vsepar.empty()) {
            area.property('last_vsepar', vsepar.style('left'));
            area.property('last_drawing', drawing.style('left'));
         }
         tgt = (-$(area.node()).outerWidth(true)-10).toString() + "px";
         let mainw = $(main.node()).outerWidth(true);

         if (vsepar.empty() && ($(area.node()).offset().left > mainw/2)) tgt = (mainw+10) + "px";

         tgt_separ = "-10px";
         tgt_drawing = "0px";
      }

      let visible_at_the_end  = !this.browser_visible, _duration = fast_close ? 0 : 700;

      this.browser_visible = 'changing';

      area.transition().style('left', tgt).duration(_duration).on("end", () => {
         if (fast_close) return;
         this.browser_visible = visible_at_the_end;
         if (visible_at_the_end) this.setButtonsPosition();
      });

      if (!visible_at_the_end)
         main.select(".jsroot_browser_btns").transition().style('left', '7px').style('top', '7px').duration(_duration);

      if (!vsepar.empty()) {
         vsepar.transition().style('left', tgt_separ).duration(_duration);
         drawing.transition().style('left', tgt_drawing).duration(_duration).on("end", this.checkResize.bind(this));
      }

      if (this.status_layout && (this.browser_kind == 'fix')) {
         main.select(".jsroot_h_separator").transition().style('left', tgt_drawing).duration(_duration);
         main.select(".jsroot_status_area").transition().style('left', tgt_drawing).duration(_duration);
      }
   }

   /** @summary Toggle browser kind
     * @desc used together with browser buttons */
   BrowserLayout.prototype.toggleKind = function(browser_kind) {
      if (this.browser_visible!=='changing') {
         if (browser_kind === this.browser_kind) this.toggleBrowserVisisbility();
                                            else this.toggleBrowserKind(browser_kind);
      }
   }

   /** @summary Delete content */
   BrowserLayout.prototype.deleteContent = function() {
      let main = d3.select("#" + this.gui_div + " .jsroot_browser");
      if (main.empty()) return;

      this.createStatusLine(0, "delete");
      let vsepar = main.select(".jsroot_v_separator");
      if (!vsepar.empty())
         $(vsepar.node()).draggable('destroy');

      this.toggleBrowserVisisbility(true);

      main.selectAll("*").remove();
      delete this.browser_visible;
      delete this.browser_kind;

      this.checkResize();
   }

   /** @summary Creates status line */
   BrowserLayout.prototype.createStatusLine = function(height, mode) {

      let main = d3.select("#"+this.gui_div+" .jsroot_browser");
      if (main.empty())
         return Promise.resolve('');

      let id = this.gui_div + "_status",
          line = d3.select("#"+id),
          is_visible = !line.empty();

      if (mode==="toggle") { mode = !is_visible; } else
      if (mode==="delete") { mode = false; height = 0; delete this.status_layout; } else
      if (mode===undefined) { mode = true; this.status_layout = "app"; }

      if (is_visible) {
         if (mode === true)
            return Promise.resolve(id);

         let hsepar = main.select(".jsroot_h_separator");

         $(hsepar.node()).draggable("destroy");

         hsepar.remove();
         line.remove();

         if (this.status_layout !== "app")
            delete this.status_layout;

         if (this.status_handler && (jsrp.showStatus === this.status_handler)) {
            delete jsrp.showStatus;
            delete this.status_handler;
         }

         this.adjustSeparators(null, 0, true);
         return Promise.resolve("");
      }

      if (mode === false)
         return Promise.resolve("");

      let left_pos = d3.select("#" + this.gui_div + "_drawing").style('left');

      line = main.insert("div",".jsroot_browser_area").attr("id",id)
                 .classed("jsroot_status_area", true)
                 .style('position',"absolute").style('left',left_pos).style('height',"20px").style('bottom',0).style('right',0)
                 .style('margin',0).style('border',0);

      let hsepar = main.insert("div",".jsroot_browser_area")
                       .classed("jsroot_separator", true).classed("jsroot_h_separator", true)
                      .style('position','absolute').style('left',left_pos).style('right',0).style('bottom','20px').style('height','5px');

      let pthis = this;

      $(hsepar.node()).draggable({
         axis: "y" , cursor: "ns-resize", containment: "parent",
         helper: function() { return $(this).clone().css('background-color','grey'); },
         drag: function(event,ui) {
            pthis.adjustSeparators(null, -ui.position.top);
         },
         stop: function(/*event,ui*/) {
            pthis.checkResize();
         }
      });

      if (!height || (typeof height === 'string')) height = this.last_hsepar_height || 20;

      this.adjustSeparators(null, height, true);

      if (this.status_layout == "app")
         return Promise.resolve(id);

      this.status_layout = new JSROOT.GridDisplay(id, 'horizx4_1213');

      let frame_titles = ['object name','object title','mouse coordinates','object info'];
      for (let k=0;k<4;++k)
         d3.select(this.status_layout.getGridFrame(k)).attr('title', frame_titles[k]).style('overflow','hidden')
           .append("label").attr("class","jsroot_status_label");

      this.status_handler = this.showStatus.bind(this);

      jsrp.showStatus = this.status_handler;

      return Promise.resolve(id);
   }

   /** @summary Adjust separator positions */
   BrowserLayout.prototype.adjustSeparators = function(vsepar, hsepar, redraw, first_time) {

      if (!this.gui_div) return;

      let main = d3.select("#" + this.gui_div + " .jsroot_browser"), w = 5;

      if ((hsepar===null) && first_time && !main.select(".jsroot_h_separator").empty()) {
         // if separator set for the first time, check if status line present
         hsepar = main.select(".jsroot_h_separator").style('bottom');
         if ((typeof hsepar=='string') && (hsepar.indexOf('px')==hsepar.length-2))
            hsepar = hsepar.substr(0,hsepar.length-2);
         else
            hsepar = null;
      }

      if (hsepar!==null) {
         hsepar = parseInt(hsepar);
         let elem = main.select(".jsroot_h_separator"), hlimit = 0;

         if (!elem.empty()) {
            if (hsepar < 0) hsepar += ($(main.node()).outerHeight(true) - w);
            if (hsepar < 5) hsepar = 5;
            this.last_hsepar_height = hsepar;
            elem.style('bottom', hsepar+'px').style('height', w+'px');
            d3.select("#" + this.gui_div + "_status").style('height', hsepar+'px');
            hlimit = (hsepar+w) + 'px';
         }

         d3.select("#" + this.gui_div + "_drawing").style('bottom',hlimit);
      }

      if (vsepar!==null) {
         vsepar = parseInt(vsepar);
         if (vsepar < 50) vsepar = 50;
         main.select(".jsroot_browser_area").style('width',(vsepar-5)+'px');
         d3.select("#" + this.gui_div + "_drawing").style('left',(vsepar+w)+'px');
         main.select(".jsroot_h_separator").style('left', (vsepar+w)+'px');
         d3.select("#" + this.gui_div + "_status").style('left',(vsepar+w)+'px');
         main.select(".jsroot_v_separator").style('left',vsepar+'px').style('width',w+"px");
      }

      if (redraw) this.checkResize();
   }

   /** @summary Show status information inside special fields of browser layout
     * @private */
   BrowserLayout.prototype.showStatus = function(name, title, info, coordinates) {
      if (!this.status_layout) return;

      $(this.status_layout.getGridFrame(0)).children('label').text(name || "");
      $(this.status_layout.getGridFrame(1)).children('label').text(title || "");
      $(this.status_layout.getGridFrame(2)).children('label').text(coordinates || "");
      $(this.status_layout.getGridFrame(3)).children('label').text(info || "");

      if (!this.status_layout.first_check) {
         this.status_layout.first_check = true;
         let maxh = 0;
         for (let n=0;n<4;++n)
            maxh = Math.max(maxh, $(this.status_layout.getGridFrame(n)).children('label').outerHeight());
         if ((maxh>5) && ((maxh>this.last_hsepar_height) || (maxh<this.last_hsepar_height+5)))
            this.adjustSeparators(null, maxh, true);
      }
   }

   // =================================================================================================

   let HierarchyPainter = JSROOT.HierarchyPainter;

   /** @summary returns true if item is last in parent childs list
     * @private */
   HierarchyPainter.prototype.isLastSibling = function(hitem) {
      if (!hitem || !hitem._parent || !hitem._parent._childs) return false;
      let chlds = hitem._parent._childs, indx = chlds.indexOf(hitem);
      if (indx<0) return false;
      while (++indx < chlds.length)
         if (!('_hidden' in chlds[indx])) return false;
      return true;
   }

   /** @summary Create item html code
     * @private */
   HierarchyPainter.prototype.addItemHtml = function(hitem, d3prnt, arg) {
      if (!hitem || ('_hidden' in hitem)) return true;

      let isroot = (hitem === this.h),
          has_childs = ('_childs' in hitem),
          handle = jsrp.getDrawHandle(hitem._kind),
          img1 = "", img2 = "", can_click = false, break_list = false,
          d3cont, itemname = this.itemFullName(hitem);

      if (handle !== null) {
         if ('icon' in handle) img1 = handle.icon;
         if ('icon2' in handle) img2 = handle.icon2;
         if ((img1.length==0) && (typeof handle.icon_get == 'function'))
            img1 = handle.icon_get(hitem, this);
         if (('func' in handle) || ('execute' in handle) || ('aslink' in handle) ||
             (('expand' in handle) && (hitem._more !== false))) can_click = true;
      }

      if ('_icon' in hitem) img1 = hitem._icon;
      if ('_icon2' in hitem) img2 = hitem._icon2;
      if ((img1.length==0) && ('_online' in hitem))
         hitem._icon = img1 = "img_globe";
      if ((img1.length==0) && isroot)
         hitem._icon = img1 = "img_base";

      if (hitem._more || hitem._expand || hitem._player || hitem._can_draw)
         can_click = true;

      let can_menu = can_click;
      if (!can_menu && (typeof hitem._kind == 'string') && (hitem._kind.indexOf("ROOT.")==0))
         can_menu = can_click = true;

      if (img2.length==0) img2 = img1;
      if (img1.length==0) img1 = (has_childs || hitem._more) ? "img_folder" : "img_page";
      if (img2.length==0) img2 = (has_childs || hitem._more) ? "img_folderopen" : "img_page";

      if (arg === "update") {
         d3prnt.selectAll("*").remove();
         d3cont = d3prnt;
      } else {
         d3cont = d3prnt.append("div");
         if (arg && (arg >= (hitem._parent._show_limit || JSROOT.settings.HierarchyLimit))) break_list = true;
      }

      hitem._d3cont = d3cont.node(); // set for direct referencing
      d3cont.attr("item", itemname);

      // line with all html elements for this item (excluding childs)
      let d3line = d3cont.append("div").attr('class','h_line');

      // build indent
      let prnt = isroot ? null : hitem._parent;
      while (prnt && (prnt !== this.h)) {
         d3line.insert("div",":first-child")
               .attr("class", this.isLastSibling(prnt) ? "img_empty" : "img_line");
         prnt = prnt._parent;
      }

      let icon_class = "", plusminus = false;

      if (isroot) {
         // for root node no extra code
      } else
      if (has_childs && !break_list) {
         icon_class = hitem._isopen ? "img_minus" : "img_plus";
         plusminus = true;
      } else
      /*if (hitem._more) {
         icon_class = "img_plus"; // should be special plus ???
         plusminus = true;
      } else */ {
         icon_class = "img_join";
      }

      let h = this;

      if (icon_class.length > 0) {
         if (break_list || this.isLastSibling(hitem)) icon_class += "bottom";
         let d3icon = d3line.append("div").attr('class', icon_class);
         if (plusminus) d3icon.style('cursor','pointer')
                              .on("click", function(evnt) { h.tree_click(evnt, this, "plusminus"); });
      }

      // make node icons

      if (this.with_icons && !break_list) {
         let icon_name = hitem._isopen ? img2 : img1;

         let d3img;

         if (icon_name.indexOf("img_")==0)
            d3img = d3line.append("div")
                          .attr("class", icon_name)
                          .attr("title", hitem._kind);
         else
            d3img = d3line.append("img")
                          .attr("src", icon_name)
                          .attr("alt","")
                          .attr("title", hitem._kind)
                          .style('vertical-align','top').style('width','18px').style('height','18px');

         if (('_icon_click' in hitem) || (handle && ('icon_click' in handle)))
            d3img.on("click", function(evnt) { h.tree_click(evnt, this, "icon"); });
      }

      let d3a = d3line.append("a");
      if (can_click || has_childs || break_list)
         d3a.attr("class","h_item")
            .on("click", function(evnt) { h.tree_click(evnt, this); });

      if (break_list) {
         hitem._break_point = true; // indicate that list was broken here
         d3a.attr('title', 'there are ' + (hitem._parent._childs.length-arg) + ' more items')
            .text("...more...");
         return false;
      }

      if ('disp_kind' in h) {
         if (JSROOT.settings.DragAndDrop && can_click)
           this.enableDrag(d3a.node(), itemname);
         if (JSROOT.settings.ContextMenu && can_menu)
            d3a.on('contextmenu', function(evnt) { h.tree_contextmenu(evnt, this); });

         d3a.on("mouseover", function() { h.tree_mouseover(true, this); })
            .on("mouseleave", function() { h.tree_mouseover(false, this); });
      } else
      if (hitem._direct_context && JSROOT.settings.ContextMenu)
         d3a.on('contextmenu', function(evnt) { h.direct_contextmenu(evnt, this); });

      let element_name = hitem._name, element_title = "";

      if ('_realname' in hitem)
         element_name = hitem._realname;

      if ('_title' in hitem)
         element_title = hitem._title;

      if ('_fullname' in hitem)
         element_title += "  fullname: " + hitem._fullname;

      if (!element_title)
         element_title = element_name;

      d3a.attr('title', element_title)
         .text(element_name + ('_value' in hitem ? ":" : ""))
         .style('background', hitem._background ? hitem._background : null);

      if ('_value' in hitem) {
         let d3p = d3line.append("p");
         if ('_vclass' in hitem) d3p.attr('class', hitem._vclass);
         if (!hitem._isopen) d3p.html(hitem._value);
      }

      if (has_childs && (isroot || hitem._isopen)) {
         let d3chlds = d3cont.append("div").attr("class", "h_childs");
         for (let i = 0; i < hitem._childs.length; ++i) {
            let chld = hitem._childs[i];
            chld._parent = hitem;
            if (!this.addItemHtml(chld, d3chlds, i)) break; // if too many items, skip rest
         }
      }

      return true;
   }

   /** @summary Toggle open state of the item
     * @desc Used with "open all" / "close all" buttons in normal GUI
     * @param {boolean} isopen - if items should be expand or closed
     * @returns {boolean} true when any item was changed */
   HierarchyPainter.prototype.toggleOpenState = function(isopen, h) {
      let hitem = h || this.h;

      if (hitem._childs === undefined) {
         if (!isopen) return false;

         if (this.with_icons) {
            // in normal hierarchy check precisely if item can be expand
            if (!hitem._more && !hitem._expand && !this.canExpandItem(hitem)) return false;
         }

         this.expandItem(this.itemFullName(hitem));
         if (hitem._childs !== undefined) hitem._isopen = true;
         return hitem._isopen;
      }

      if ((hitem !== this.h) && isopen && !hitem._isopen) {
         // when there are childs and they are not see, simply show them
         hitem._isopen = true;
         return true;
      }

      let change_child = false;
      for (let i = 0; i < hitem._childs.length; ++i)
         if (this.toggleOpenState(isopen, hitem._childs[i]))
            change_child = true;

      if ((hitem !== this.h) && !isopen && hitem._isopen && !change_child) {
         // if none of the childs can be closed, than just close that item
         delete hitem._isopen;
         return true;
       }

      if (!h) this.refreshHtml();
      return false;
   }

   /** @summary Refresh HTML code of hierarchy painter
     * @returns {Promise} when done */
   HierarchyPainter.prototype.refreshHtml = function() {

      let d3elem = this.selectDom();
      if (d3elem.empty())
         return Promise.resolve(this);

      d3elem.html("")
            .style('overflow','hidden') // clear html - most simple way
            .style('display','flex')
            .style('flex-direction','column');

      let h = this, factcmds = [], status_item = null;
      this.forEachItem(item => {
         delete item._d3cont; // remove html container
         if (('_fastcmd' in item) && (item._kind == 'Command')) factcmds.push(item);
         if (('_status' in item) && !status_item) status_item = item;
      });

      if (!this.h || d3elem.empty())
         return Promise.resolve(this);

      if (factcmds.length) {
         let fastbtns = d3elem.append("div").attr("class","jsroot");
         for (let n=0;n<factcmds.length;++n) {
            let btn = fastbtns.append("button")
                       .text("")
                       .attr("class",'fast_command')
                       .attr("item", this.itemFullName(factcmds[n]))
                       .attr("title", factcmds[n]._title)
                       .on("click", function() { h.executeCommand(d3.select(this).attr("item"), this); } );

            if ('_icon' in factcmds[n])
               btn.append('img').attr("src", factcmds[n]._icon);
         }
      }

      let d3btns = d3elem.append("p").attr("class", "jsroot").style("margin-bottom","3px").style("margin-top",0);
      d3btns.append("a").attr("class", "h_button").text("open all")
            .attr("title","open all items in the browser").on("click", () => this.toggleOpenState(true));
      d3btns.append("text").text(" | ");
      d3btns.append("a").attr("class", "h_button").text("close all")
            .attr("title","close all items in the browser").on("click", () => this.toggleOpenState(false));

      if (typeof this.removeInspector == 'function') {
         d3btns.append("text").text(" | ");
         d3btns.append("a").attr("class", "h_button").text("remove")
               .attr("title","remove inspector").on("click", () => this.removeInspector());
      }

      if ('_online' in this.h) {
         d3btns.append("text").text(" | ");
         d3btns.append("a").attr("class", "h_button").text("reload")
               .attr("title","reload object list from the server").on("click", () => this.reload());
      }

      if ('disp_kind' in this) {
         d3btns.append("text").text(" | ");
         d3btns.append("a").attr("class", "h_button").text("clear")
               .attr("title","clear all drawn objects").on("click", () => this.clearHierarchy(false));
      }

      let maindiv =
         d3elem.append("div")
               .attr("class", "jsroot")
               .style('font-size', this.with_icons ? "12px" : "15px")
               .style("overflow","auto")
               .style("flex","1");

      if (this.background) // case of object inspector and streamer infos display
         maindiv.style("background-color", this.background)
                .style('margin', '2px').style('padding', '2px');

      this.addItemHtml(this.h, maindiv.append("div").attr("class","h_tree"));

      this.setTopPainter(); //assign hpainter as top painter

      if (status_item && !this.status_disabled && !JSROOT.decodeUrl().has('nostatus')) {
         let func = JSROOT.findFunction(status_item._status);
         if (typeof func == 'function')
            return this.createStatusLine().then(sdiv => {
               if (sdiv) func(sdiv, this.itemFullName(status_item));
            });
      }

      return Promise.resolve(this);
   }

   /** @summary Update item node
     * @private */
   HierarchyPainter.prototype.updateTreeNode = function(hitem, d3cont) {
      if ((d3cont===undefined) || d3cont.empty())  {
         d3cont = d3.select(hitem._d3cont ? hitem._d3cont : null);
         let name = this.itemFullName(hitem);
         if (d3cont.empty())
            d3cont = this.selectDom().select("[item='" + name + "']");
         if (d3cont.empty() && ('_cycle' in hitem))
            d3cont = this.selectDom().select("[item='" + name + ";" + hitem._cycle + "']");
         if (d3cont.empty()) return;
      }

      this.addItemHtml(hitem, d3cont, "update");

      if (this.brlayout) this.brlayout.adjustBrowserSize(true);
   }

   /** @summary Update item background
     * @private */
   HierarchyPainter.prototype.updateBackground = function(hitem, scroll_into_view) {

      if (!hitem || !hitem._d3cont) return;

      let d3cont = d3.select(hitem._d3cont);

      if (d3cont.empty()) return;

      let d3a = d3cont.select(".h_item");

      d3a.style('background', hitem._background ? hitem._background : null);

      if (scroll_into_view && hitem._background)
         d3a.node().scrollIntoView(false);
   }

   /** @summary Handler for click event of item in the hierarchy
     * @private */
   HierarchyPainter.prototype.tree_click = function(evnt, node, place) {
      if (!node) return;

      let d3cont = d3.select(node.parentNode.parentNode),
          itemname = d3cont.attr('item'),
          hitem = itemname ? this.findItem(itemname) : null;
      if (!hitem) return;

      if (hitem._break_point) {
         // special case of more item

         delete hitem._break_point;

         // update item itself
         this.addItemHtml(hitem, d3cont, "update");

         let prnt = hitem._parent, indx = prnt._childs.indexOf(hitem),
             d3chlds = d3.select(d3cont.node().parentNode);

         if (indx<0) return console.error('internal error');

         prnt._show_limit = (prnt._show_limit || JSROOT.settings.HierarchyLimit) * 2;

         for (let n=indx+1;n<prnt._childs.length;++n) {
            let chld = prnt._childs[n];
            chld._parent = prnt;
            if (!this.addItemHtml(chld, d3chlds, n)) break; // if too many items, skip rest
         }

         return;
      }

      let prnt = hitem, dflt = undefined;
      while (prnt) {
         if ((dflt = prnt._click_action) !== undefined) break;
         prnt = prnt._parent;
      }

      if (!place || (place=="")) place = "item";

      let sett = jsrp.getDrawSettings(hitem._kind), handle = sett.handle;

      if (place == "icon") {
         let func = null;
         if (typeof hitem._icon_click == 'function') func = hitem._icon_click; else
         if (handle && typeof handle.icon_click == 'function') func = handle.icon_click;
         if (func && func(hitem,this))
            this.updateTreeNode(hitem, d3cont);
         return;
      }

      // special feature - all items with '_expand' function are not drawn by click
      if ((place=="item") && ('_expand' in hitem) && !evnt.ctrlKey && !evnt.shiftKey) place = "plusminus";

      // special case - one should expand item
      if (((place == "plusminus") && !('_childs' in hitem) && hitem._more) ||
          ((place == "item") && (dflt === "expand"))) {
         return this.expandItem(itemname, d3cont);
      }

      if (place == "item") {

         if ('_player' in hitem)
            return this.player(itemname);

         if (handle && handle.aslink)
            return window.open(itemname + "/");

         if (handle && handle.execute)
            return this.executeCommand(itemname, node.parentNode);

         if (handle && handle.ignore_online && this.isOnlineItem(hitem)) return;

         let can_draw = hitem._can_draw,
             can_expand = hitem._more,
             dflt_expand = (this.default_by_click === "expand"),
             drawopt = "";

         if (evnt.shiftKey) {
            drawopt = (handle && handle.shift) ? handle.shift : "inspect";
            if ((drawopt==="inspect") && handle && handle.noinspect) drawopt = "";
         }
         if (handle && handle.ctrl && evnt.ctrlKey) drawopt = handle.ctrl;

         if (!drawopt) {
            for (let pitem = hitem._parent; !!pitem; pitem = pitem._parent) {
               if (pitem._painter) { can_draw = false; if (can_expand===undefined) can_expand = false; break; }
            }
         }

         if (hitem._childs) can_expand = false;

         if (can_draw === undefined) can_draw = sett.draw;
         if (can_expand === undefined) can_expand = sett.expand;

         if (can_draw && can_expand && !drawopt) {
            // if default action specified as expand, disable drawing
            if (dflt_expand || (handle && (handle.dflt === 'expand'))) can_draw = false; else
            if (this.isItemDisplayed(itemname)) can_draw = false; // if already displayed, try to expand
         }

         if (can_draw)
            return this.display(itemname, drawopt);

         if (can_expand || dflt_expand)
            return this.expandItem(itemname, d3cont);

         // cannot draw, but can inspect ROOT objects
         if ((typeof hitem._kind === "string") && (hitem._kind.indexOf("ROOT.")===0) && sett.inspect && (can_draw!==false))
            return this.display(itemname, "inspect");

         if (!hitem._childs || (hitem === this.h)) return;
      }

      if (hitem._isopen)
         delete hitem._isopen;
      else
         hitem._isopen = true;

      this.updateTreeNode(hitem, d3cont);
   }

   /** @summary Handler for mouse-over event
     * @private */
   HierarchyPainter.prototype.tree_mouseover = function(on, elem) {
      let itemname = d3.select(elem.parentNode.parentNode).attr('item');

      let hitem = this.findItem(itemname);
      if (!hitem) return;

      let painter, prnt = hitem;
      while (prnt && !painter) {
         painter = prnt._painter;
         prnt = prnt._parent;
      }

      if (painter && typeof painter.mouseOverHierarchy === 'function')
         painter.mouseOverHierarchy(on, itemname, hitem);
   }

   /** @summary alternative context menu, used in the object inspector
     * @private */
   HierarchyPainter.prototype.direct_contextmenu = function(evnt, elem) {
      evnt.preventDefault();
      let itemname = d3.select(elem.parentNode.parentNode).attr('item');
      let hitem = this.findItem(itemname);
      if (!hitem) return;

      if (typeof this.fill_context == 'function')
         jsrp.createMenu(evnt, this).then(menu => {
            this.fill_context(menu, hitem);
            if (menu.size() > 0) {
               menu.tree_node = elem.parentNode;
               menu.show();
            }
         });
   }

   /** @summary Handle context menu in the hieararchy
     * @private */
   HierarchyPainter.prototype.tree_contextmenu = function(evnt, elem) {
      // this is handling of context menu request for the normal objects browser

      evnt.preventDefault();

      let itemname = d3.select(elem.parentNode.parentNode).attr('item');

      let hitem = this.findItem(itemname);
      if (!hitem) return;

      let onlineprop = this.getOnlineProp(itemname),
          fileprop = this.getFileProp(itemname);

      function qualifyURL(url) {
         function escapeHTML(s) {
            return s.split('&').join('&amp;').split('<').join('&lt;').split('"').join('&quot;');
         }
         let el = document.createElement('div');
         el.innerHTML = '<a href="' + escapeHTML(url) + '">x</a>';
         return el.firstChild.href;
      }

      jsrp.createMenu(evnt, this).then(menu => {

         if ((itemname == "") && !('_jsonfile' in hitem)) {
            let files = [], addr = "", cnt = 0,
                separ = () => (cnt++ > 0) ? "&" : "?";

            this.forEachRootFile(item => files.push(item._file.fFullURL));

            if (!this.getTopOnlineItem())
               addr = JSROOT.source_dir + "index.htm";

            if (this.isMonitoring())
               addr += separ() + "monitoring=" + this.getMonitoringInterval();

            if (files.length==1)
               addr += separ() + "file=" + files[0];
            else
               if (files.length>1)
                  addr += separ() + "files=" + JSON.stringify(files);

            if (this.disp_kind)
               addr += separ() + "layout=" + this.disp_kind.replace(/ /g, "");

            let items = [];

            if (this.disp)
               this.disp.forEachPainter(p => {
                  if (p.getItemName())
                     items.push(p.getItemName());
               });

            if (items.length == 1) {
               addr += separ() + "item=" + items[0];
            } else if (items.length > 1) {
               addr += separ() + "items=" + JSON.stringify(items);
            }

            menu.add("Direct link", () => window.open(addr));
            menu.add("Only items", () => window.open(addr + "&nobrowser"));
         } else if (onlineprop) {
            this.fillOnlineMenu(menu, onlineprop, itemname);
         } else {
            let sett = jsrp.getDrawSettings(hitem._kind, 'nosame');

            // allow to draw item even if draw function is not defined
            if (hitem._can_draw) {
               if (!sett.opts) sett.opts = [""];
               if (sett.opts.indexOf("")<0) sett.opts.unshift("");
            }

            if (sett.opts)
               menu.addDrawMenu("Draw", sett.opts, arg => this.display(itemname, arg));

            if (fileprop && sett.opts && !fileprop.localfile) {
               let filepath = qualifyURL(fileprop.fileurl);
               if (filepath.indexOf(JSROOT.source_dir) == 0)
                  filepath = filepath.slice(JSROOT.source_dir.length);
               filepath = fileprop.kind + "=" + filepath;
               if (fileprop.itemname.length > 0) {
                  let name = fileprop.itemname;
                  if (name.search(/\+| |\,/)>=0) name = "\'" + name + "\'";
                  filepath += "&item=" + name;
               }

               menu.addDrawMenu("Draw in new tab", sett.opts, arg => {
                  window.open(JSROOT.source_dir + "index.htm?nobrowser&"+filepath +"&opt="+arg);
               });
            }

            if (sett.expand && !('_childs' in hitem) && (hitem._more || !('_more' in hitem)))
               menu.add("Expand", () => this.expandItem(itemname));

            if (hitem._kind === "ROOT.TStyle")
               menu.add("Apply", () => this.applyStyle(itemname));
         }

         if (typeof hitem._menu == 'function')
            hitem._menu(menu, hitem, this);

         if (menu.size() > 0) {
            menu.tree_node = elem.parentNode;
            if (menu.separ) menu.add("separator"); // add separator at the end
            menu.add("Close");
            menu.show();
         }

      }); // end menu creation

      return false;
   }

   /** @summary Creates configured JSROOT.MDIDisplay object
     * @returns {Promise} with created mdi object */
   HierarchyPainter.prototype.createDisplay = function() {

      if ('disp' in this) {
         if (this.disp.numDraw() > 0)
            return Promise.resolve(this.disp);
         this.disp.cleanup();
         delete this.disp;
      }

      // check that we can found frame where drawing should be done
      if (!document.getElementById(this.disp_frameid))
         return Promise.resolve(null);

      if (this.disp_kind == "tabs")
         this.disp = new TabsDisplay(this.disp_frameid);
      else if (this.disp_kind.indexOf("flex")==0)
         this.disp = new FlexibleDisplay(this.disp_frameid);
      else if (this.disp_kind.indexOf("coll")==0)
         this.disp = new CollapsibleDisplay(this.disp_frameid);
      else
         this.disp = new JSROOT.GridDisplay(this.disp_frameid, this.disp_kind);

      if (this.disp)
         this.disp.cleanupFrame = this.cleanupFrame.bind(this);

      return Promise.resolve(this.disp);
   }

   /** @summary Enable drag on the element
     * @private  */
   HierarchyPainter.prototype.enableDrag = function(element /*, itemname*/) {
      $(element).draggable({ revert: "invalid", appendTo: "body", helper: "clone" });
   }

   /** @summary Enable drop on the element
     * @private  */
   HierarchyPainter.prototype.enableDrop = function(frame, itemname) {
      let h = this;
      $(frame).droppable({
         hoverClass : "ui-state-active",
         accept: function(ui) {
            let dropname = ui.parent().parent().attr('item');
            if ((dropname == itemname) || !dropname) return false;

            let ditem = h.findItem(dropname);
            if (!ditem || (!('_kind' in ditem))) return false;

            return ditem._kind.indexOf("ROOT.")==0;
         },
         drop: function(event, ui) {
            let dropname = ui.draggable.parent().parent().attr('item');
            if (dropname) h.dropItem(dropname, $(this).attr("id"));
         }
      });
   }

   /** @summary Create browser elements */
   HierarchyPainter.prototype.createBrowser = function(browser_kind, update_html) {

      if (!this.gui_div || this.exclude_browser || !this.brlayout)
         return Promise.resolve(false);

      let main = d3.select("#" + this.gui_div + " .jsroot_browser"),
          jmain = $(main.node());

      // one requires top-level container
      if (main.empty())
         return Promise.resolve(false);

      if ((browser_kind==="float") && this.float_browser_disabled) browser_kind = "fix";

      if (!main.select('.jsroot_browser_area').empty()) {
         // this is case when browser created,
         // if update_html specified, hidden state will be toggled

         if (update_html) this.brlayout.toggleKind(browser_kind);

         return Promise.resolve(true);
      }

      let guiCode = "<p class='jsroot_browser_version'><a href='https://root.cern/js/'>JSROOT</a> version <span style='color:green'><b>" + JSROOT.version + "</b></span></p>";

      if (this.is_online) {
         guiCode +='<p> Hierarchy in <a href="h.json">json</a> and <a href="h.xml">xml</a> format</p>'
                 + '<div style="display:inline-block;">'
                 + '<label style="margin-right:5px; vertical-align:middle;">'
                 + '<input style="vertical-align:middle;" type="checkbox" name="monitoring" class="gui_monitoring"/>'
                 + 'Monitoring</label>';
      } else if (!this.no_select) {
         let myDiv = d3.select("#"+this.gui_div),
             files = myDiv.attr("files") || "../files/hsimple.root",
             path = JSROOT.decodeUrl().get("path") || myDiv.attr("path") || "",
             arrFiles = files.split(';');

         guiCode +=
            '<input type="text" value="" style="width:95%; margin:5px;border:2px;" class="gui_urlToLoad" title="input file name"/>'
            +'<div style="display:flex;flex-direction:row;padding-top:5px">'
            +'<select class="gui_selectFileName" style="flex:1;padding:2px;" title="select file name"'
            +'<option value="" selected="selected"></option>';
         for (let i in arrFiles)
            guiCode += '<option value = "' + path + arrFiles[i] + '">' + arrFiles[i] + '</option>';
         guiCode += '</select>'
            +'<input type="file" class="gui_localFile" accept=".root" style="display:none"/><output id="list" style="display:none"></output>'
            +'<input type="button" value="..." class="gui_fileBtn" style="min-width:3em;padding:3px;margin-left:5px;margin-right:5px;" title="select local file for reading"/><br/>'
            +'</div>'
            +'<p id="gui_fileCORS"><small><a href="https://github.com/root-project/jsroot/blob/master/docs/JSROOT.md#reading-root-files-from-other-servers">Read docu</a>'
            +' how to open files from other servers.</small></p>'
            +'<div style="display:flex;flex-direction:row">'
            +'<input style="padding:3px;margin-right:5px;"'
            +'       class="gui_ReadFileBtn" type="button" title="Read the Selected File" value="Load"/>'
            +'<input style="padding:3px;margin-right:5px;"'
            +'       class="gui_ResetUIBtn" type="button" title="Close all opened files and clear drawings" value="Reset"/>'
      } else if (this.no_select == "file") {
         guiCode += '<div style="display:flex;flex-direction:row">';
      }

      if (this.is_online || !this.no_select || this.no_select=="file")
         guiCode += '<select style="padding:2px;margin-right:5px;" title="layout kind" class="gui_layout"></select>'
                  + '</div>';

      guiCode += '<div id="' + this.gui_div+'_browser_hierarchy" class="jsroot_browser_hierarchy"></div>';

      this.brlayout.setBrowserContent(guiCode);

      this.brlayout.setBrowserTitle(this.is_online ? 'ROOT online server' : 'Read a ROOT file');

      let localfile_read_callback = null;

      if (!this.is_online && !this.no_select) {

         this.ReadSelectedFile = function() {
            let filename = main.select(".gui_urlToLoad").property('value').trim();
            if (!filename) return;

            if ((filename.toLowerCase().lastIndexOf(".json") == filename.length-5))
               this.openJsonFile(filename);
            else
               this.openRootFile(filename);
         }

         jmain.find(".gui_selectFileName").val("")
              .change(function() { jmain.find(".gui_urlToLoad").val($(this).val()); });
         jmain.find(".gui_fileBtn").button()
              .click(() => jmain.find(".gui_localFile").click());

         jmain.find(".gui_ReadFileBtn").button().click(() => this.ReadSelectedFile());

         jmain.find(".gui_ResetUIBtn").button().click(() => this.clearHierarchy(true));

         jmain.find(".gui_urlToLoad").keyup(e => {
            if (e.keyCode == 13) this.ReadSelectedFile();
         });

         jmain.find(".gui_localFile").change(evnt => {
            let files = evnt.target.files;

            for (let n=0;n<files.length;++n) {
               let f = files[n];
               main.select(".gui_urlToLoad").property('value', f.name);
               this.openRootFile(f).then(localfile_read_callback);
            }

            localfile_read_callback = null;
         });

         this.selectLocalFile = function() {
            return new Promise(resolveFunc => {
               localfile_read_callback = resolveFunc;
               $("#" + this.gui_div + " .jsroot_browser").find(".gui_localFile").click();
            });
         }
      }

      let jlayout = jmain.find(".gui_layout");
      if (jlayout.length) {
         let lst = ['simple', 'vert2', 'vert3', 'vert231', 'horiz2', 'horiz32', 'flex',
                     'grid 2x2', 'grid 1x3', 'grid 2x3', 'grid 3x3', 'grid 4x4', 'collapsible',  'tabs'];

         for (let k=0;k<lst.length;++k){
            let opt = document.createElement('option');
            opt.value = lst[k];
            opt.innerHTML = lst[k];
            jlayout.get(0).appendChild(opt);
         }

         let painter = this;
         jlayout.change(function() {
            painter.setDisplay($(this).val() || 'collapsible', painter.gui_div + "_drawing");
         });
      }

      this.setDom(this.gui_div + '_browser_hierarchy');

      if (update_html) {
         this.refreshHtml();
         this.initializeBrowser();
      }

      this.brlayout.toggleBrowserKind(browser_kind || "fix");

      return Promise.resolve(true);
   }

   /** @summary Initialize browser elements */
   HierarchyPainter.prototype.initializeBrowser = function() {

      let main = d3.select("#" + this.gui_div + " .jsroot_browser");
      if (main.empty() || !this.brlayout) return;
      let jmain = $(main.node());

      if (this.brlayout) this.brlayout.adjustBrowserSize();

      let selects = main.select(".gui_layout").node();

      if (selects) {
         let found = false;
         for (let i in selects.options) {
            let s = selects.options[i].text;
            if (typeof s !== 'string') continue;
            if ((s == this.getLayout()) || (s.replace(/ /g,"") == this.getLayout())) {
               selects.selectedIndex = i; found = true;
               break;
            }
         }
         if (!found) {
            let opt = document.createElement('option');
            opt.innerHTML = opt.value = this.getLayout();
            selects.appendChild(opt);
            selects.selectedIndex = selects.options.length-1;
         }
      }

      if (this.is_online) {
         if (this.h && this.h._toptitle)
            this.brlayout.setBrowserTitle(this.h._toptitle);
         let painter = this;
         jmain.find(".gui_monitoring")
           .prop('checked', this.isMonitoring())
           .click(function() {
               painter.enableMonitoring(this.checked);
               painter.updateItems();
            });
      } else if (!this.no_select) {
         let fname = "";
         this.forEachRootFile(item => { if (!fname) fname = item._fullurl; });
         jmain.find(".gui_urlToLoad").val(fname);
      }
   }

   /** @summary Enable monitoring mode */
   HierarchyPainter.prototype.enableMonitoring = function(on) {
      this.setMonitoring(undefined, on);

      let chkbox = d3.select("#" + this.gui_div + " .jsroot_browser .gui_monitoring");
      if (!chkbox.empty() && (chkbox.property('checked') !== on))
         chkbox.property('checked', on);
   }

   /** @summary Build main JSROOT GUI
     * @returns {Promise} when completed
     * @private  */
   JSROOT.buildGUI = function(gui_element, gui_kind) {
      let myDiv = (typeof gui_element == 'string') ? d3.select('#' + gui_element) : d3.select(gui_element);
      if (myDiv.empty()) return alert('no div for gui found');

      let online = false;
      if (gui_kind == "online") online = true;

      if (myDiv.attr("ignoreurl") === "true")
         JSROOT.settings.IgnoreUrlOptions = true;

      if (JSROOT.decodeUrl().has("nobrowser") || (myDiv.attr("nobrowser") && myDiv.attr("nobrowser")!=="false") || (gui_kind == "draw") || (gui_kind == "nobrowser"))
         return JSROOT.buildNobrowserGUI(gui_element, gui_kind);

      jsrp.readStyleFromURL();

      let hpainter = new JSROOT.HierarchyPainter('root', null);

      hpainter.is_online = online;

      return hpainter.startGUI(myDiv).then(() => {
         hpainter.initializeBrowser();
         return hpainter;
      });
   }

   // ==================================================

   class CollapsibleDisplay extends JSROOT.MDIDisplay {
      constructor(frameid) {
         super(frameid);
         this.cnt = 0; // use to count newly created frames
      }

      forEachFrame(userfunc,  only_visible) {
         let topid = this.frameid + '_collapsible';

         if (!document.getElementById(topid)) return;

         if (typeof userfunc != 'function') return;

         $('#' + topid + ' .collapsible_draw').each(function() {

            // check if only visible specified
            if (only_visible && $(this).is(":hidden")) return;

            userfunc($(this).get(0));
         });
      }

      getActiveFrame() {
         let found = super.getActiveFrame();
         if (found && !$(found).is(":hidden")) return found;

         found = null;
         this.forEachFrame(frame => { if (!found) found = frame; }, true);

         return found;
      }

      activateFrame(frame) {
         if ($(frame).is(":hidden")) {
            $(frame).prev().toggleClass("ui-accordion-header-active ui-state-active ui-state-default ui-corner-bottom")
                    .find("> .ui-icon").toggleClass("ui-icon-triangle-1-e ui-icon-triangle-1-s").end()
                    .next().toggleClass("ui-accordion-content-active").slideDown(0);
         }
         $(frame).prev()[0].scrollIntoView();
         // remember title
         this.active_frame_title = d3.select(frame).attr('frame_title');
      }

      createFrame(title) {

         this.beforeCreateFrame(title);

         let topid = this.frameid + '_collapsible';

         if (!document.getElementById(topid))
            $("#"+this.frameid).append('<div id="'+ topid  + '" class="jsroot ui-accordion ui-accordion-icons ui-widget ui-helper-reset" style="overflow:auto; overflow-y:scroll; height:100%; padding-left: 2px; padding-right: 2px"></div>');

         let mdi = this,
             hid = topid + "_sub" + this.cnt++,
             uid = hid + "h",
             entryInfo = "<h5 id=\"" + uid + "\">" +
                           "<span class='ui-icon ui-icon-triangle-1-e'></span>" +
                           "<a> " + title + "</a>&nbsp; " +
                           "<button type='button' class='jsroot_collaps_closebtn' style='float:right; width:1.4em' title='close canvas'/>" +
                           " </h5>\n" +
                           "<div class='collapsible_draw' id='" + hid + "'></div>\n";

         $("#" + topid).append(entryInfo);

         $('#' + uid)
               .addClass("ui-accordion-header ui-helper-reset ui-state-default ui-corner-top ui-corner-bottom")
               .hover(function() { $(this).toggleClass("ui-state-hover"); })
               .click( function() {
                        $(this).toggleClass("ui-accordion-header-active ui-state-active ui-state-default ui-corner-bottom")
                              .find("> .ui-icon").toggleClass("ui-icon-triangle-1-e ui-icon-triangle-1-s")
                              .end().next().toggleClass("ui-accordion-content-active").slideToggle(0);
                        let sub = $(this).next(), hide_drawing = sub.is(":hidden");
                        sub.css('display', hide_drawing ? 'none' : '');
                        if (!hide_drawing) JSROOT.resize(sub.get(0));
                     })
               .next()
               .addClass("ui-accordion-content ui-helper-reset ui-widget-content ui-corner-bottom")
               .hide();

         $('#' + uid).find(" .jsroot_collaps_closebtn")
              .button({ icons: { primary: "ui-icon-close" }, text: false })
              .click(function(){
                 mdi.cleanupFrame($(this).parent().next().attr('id'));
                 $(this).parent().next().remove(); // remove drawing
                 $(this).parent().remove();  // remove header
              });

         $('#' + uid)
               .toggleClass("ui-accordion-header-active ui-state-active ui-state-default ui-corner-bottom")
               .find("> .ui-icon").toggleClass("ui-icon-triangle-1-e ui-icon-triangle-1-s").end().next()
               .toggleClass("ui-accordion-content-active").slideToggle(0);

         return $("#" + hid).attr('frame_title', title).css('overflow','hidden')
                            .attr('can_resize','height') // inform JSROOT that it can resize height of the
                            .css('position','relative') // this required for correct positioning of 3D canvas in WebKit
                            .get(0);
       }

    } // class CollapsibleDisplay

   // ================================================

   class TabsDisplay extends JSROOT.MDIDisplay {

      constructor(frameid) {
         super(frameid);
         this.cnt = 0;
      }

      forEachFrame(userfunc, only_visible) {
         let topid = this.frameid + '_tabs';

         if (!document.getElementById(topid)) return;

         if (typeof userfunc != 'function') return;

         let cnt = -1;
         let active = $('#' + topid).tabs("option", "active");

         $('#' + topid + '> .tabs_draw').each(function() {
            cnt++;
            if (!only_visible || (cnt == active))
               userfunc($(this).get(0));
         });
      }

      getActiveFrame() {
         let found = null;
         this.forEachFrame(frame => { if (!found) found = frame; }, true);
         return found;
      }

      activateFrame(frame) {
         let cnt = 0, id = -1;
         this.forEachFrame(fr => {
            if ($(fr).attr('id') == $(frame).attr('id')) id = cnt;
            cnt++;
         });
         $('#' + this.frameid + "_tabs").tabs("option", "active", id);

         this.active_frame_title = d3.select(frame).attr('frame_title');
      }

      createFrame(title) {

         this.beforeCreateFrame(title);

         let mdi = this,
             topid = this.frameid + '_tabs',
             hid = topid + "_sub" + this.cnt++,
             li = '<li><a href="#' + hid + '">' + title
                    + '</a><span class="ui-icon ui-icon-close" style="float: left; margin: 0.4em 0.2em 0 0; cursor: pointer;" role="presentation">Remove Tab</span></li>',
            cont = '<div class="tabs_draw" id="' + hid + '"></div>';

         if (!document.getElementById(topid)) {
            $("#" + this.frameid).append('<div id="' + topid + '" class="jsroot">' + ' <ul>' + li + ' </ul>' + cont + '</div>');

            let tabs = $("#" + topid)
                          .css('overflow','hidden')
                          .tabs({
                             heightStyle : "fill",
                             activate : function (event,ui) {
                                $(ui.newPanel).css('overflow', 'hidden');
                                JSROOT.resize($(ui.newPanel).get(0));
                              }
                           });

            tabs.delegate("span.ui-icon-close", "click", function() {
               let panelId = $(this).closest("li").remove().attr("aria-controls");
               mdi.cleanupFrame(panelId);
               $("#" + panelId).remove();
               tabs.tabs("refresh");
               if ($('#' + topid + '> .tabs_draw').length == 0)
                  $("#" + topid).remove();

            });
         } else {
            $("#" + topid).find("> .ui-tabs-nav").append(li);
            $("#" + topid).append(cont);
            $("#" + topid).tabs("refresh");
            $("#" + topid).tabs("option", "active", -1);
         }
         $('#' + hid)
            .empty()
            .css('overflow', 'hidden')
            .attr('frame_title', title);

         return $('#' + hid).get(0);
      }

      checkMDIResize(frame_id, size) {
         $("#" + this.frameid + '_tabs').tabs("refresh");
         super.checkMDIResize(frame_id, size);
      }

   } // class TabsDisplay

   // ==================================================

   class FlexibleDisplay extends JSROOT.MDIDisplay {

      constructor(frameid) {
         super(frameid);
         this.cnt = 0; // use to count newly created frames
      }

      forEachFrame(userfunc,  only_visible) {
         let topid = this.frameid + '_flex';

         if (!document.getElementById(topid)) return;
         if (typeof userfunc != 'function') return;

         $('#' + topid + ' .flex_draw').each(function() {
            // check if only visible specified
            if (only_visible && $(this).is(":hidden")) return;

            userfunc($(this).get(0));
         });
      }

      getActiveFrame() {
         let found = super.getActiveFrame();
         if (found && !$(found).is(":hidden")) return found;

         found = null;
         this.forEachFrame(frame => { if (!found) found = frame; }, true);

         return found;
      }

      createFrame(title) {

         this.beforeCreateFrame(title);

         let topid = this.frameid + '_flex';

         if (!document.getElementById(topid))
            $("#" + this.frameid).append('<div id="'+ topid  + '" class="jsroot" style="overflow:none; height:100%; width:100%"></div>');

         let mdi = this,
             top = $("#" + topid),
             w = top.width(),
             h = top.height(),
             subid = topid + "_frame" + this.cnt;

         let entry ='<div id="' + subid + '" class="flex_frame" style="position:absolute">' +
                     '<div class="ui-widget-header flex_header">'+
                       '<p>'+title+'</p>' +
                       '<button type="button" style="float:right; width:1.4em"/>' +
                       '<button type="button" style="float:right; width:1.4em"/>' +
                       '<button type="button" style="float:right; width:1.4em"/>' +
                      '</div>' +
                     '<div id="' + subid + '_cont" class="flex_draw"></div>' +
                    '</div>';

         top.append(entry);

         function PopupWindow(div) {
            if (div === 'first') {
               div = null;
               $('#' + topid + ' .flex_frame').each(function() {
                  if (!$(this).is(":hidden") && ($(this).prop('state') != "minimal")) div = $(this);
               });
               if (!div) return;
            }

            div.appendTo(div.parent());

            if (div.prop('state') == "minimal") return;

            div = div.find(".flex_draw").get(0);
            jsrp.selectActivePad({ pp: jsrp.getElementCanvPainter(div), active: true });

            JSROOT.resize(div);
         }

         function ChangeWindowState(main, state) {
            let curr = main.prop('state');
            if (!curr) curr = "normal";
            main.prop('state', state);
            if (state==curr) return;

            if (curr == "normal") {
               main.prop('original_height', main.height());
               main.prop('original_width', main.width());
               main.prop('original_top', main.css('top'));
               main.prop('original_left', main.css('left'));
            }

            main.find(".jsroot_minbutton").find('.ui-icon')
                .toggleClass("ui-icon-triangle-1-s", state!="minimal")
                .toggleClass("ui-icon-triangle-2-n-s", state=="minimal");

            main.find(".jsroot_maxbutton").find('.ui-icon')
                .toggleClass("ui-icon-triangle-1-n", state!="maximal")
                .toggleClass("ui-icon-triangle-2-n-s", state=="maximal");

            switch (state) {
               case "minimal":
                  main.height(main.find('.flex_header').height()).width("auto");
                  main.find(".flex_draw").css("display","none");
                  main.find(".ui-resizable-handle").css("display","none");
                  break;
               case "maximal":
                  main.height("100%").width("100%").css('left','').css('top','');
                  main.find(".flex_draw").css("display","");
                  main.find(".ui-resizable-handle").css("display","none");
                  break;
               default:
                  main.find(".flex_draw").css("display","");
                  main.find(".ui-resizable-handle").css("display","");
                  main.height(main.prop('original_height'))
                      .width(main.prop('original_width'));
                  if (curr!="minimal")
                     main.css('left', main.prop('original_left'))
                         .css('top', main.prop('original_top'));
            }

            if (state !== "minimal")
               PopupWindow(main);
            else
               PopupWindow("first");
         }

         $("#" + subid)
            .css('left', parseInt(w * (this.cnt % 5)/10))
            .css('top', parseInt(h * (this.cnt % 5)/10))
            .width(Math.round(w * 0.58))
            .height(Math.round(h * 0.58))
            .resizable({
               helper: "jsroot-flex-resizable-helper",
               start: function(/* event, ui */) {
                  // bring element to front when start resizing
                  PopupWindow($(this));
               },
               stop: function(event, ui) {
                  let rect = { width : ui.size.width-1, height : ui.size.height - $(this).find(".flex_header").height()-1 };
                  JSROOT.resize($(this).find(".flex_draw").get(0), rect);
               }
             })
             .draggable({
               containment: "parent",
               start: function(event, ui) {
                  // bring element to front when start dragging
                  PopupWindow($(this));

                  let ddd = $(this).find(".flex_draw");

                  // block dragging when mouse below header
                  let elementMouseIsOver = document.elementFromPoint(event.clientX, event.clientY);
                  let isparent = false;
                  $(elementMouseIsOver).parents().map(function() { if ($(this).get(0) === ddd.get(0)) isparent = true; });
                  if (isparent) return false;
               }
            })
          .click(function() { PopupWindow($(this)); })
          .find('.flex_header')
            // .hover(function() { $(this).toggleClass("ui-state-hover"); })
            .click(function() {
               PopupWindow($(this).parent());
            })
            .dblclick(function() {
               let main = $(this).parent();
               if (main.prop('state') == "normal")
                  ChangeWindowState(main, "maximal");
               else
                  ChangeWindowState(main, "normal");
            })
           .find("button")
              .first()
              .attr('title','close canvas')
              .button({ icons: { primary: "ui-icon-close" }, text: false })
              .click(function() {
                 let main = $(this).parent().parent();
                 mdi.cleanupFrame(main.find(".flex_draw").get(0));
                 main.remove();
                 PopupWindow('first'); // set active as first window
              })
              .next()
              .attr('title','maximize canvas')
              .addClass('jsroot_maxbutton')
              .button({ icons: { primary: "ui-icon-triangle-1-n" }, text: false })
              .click(function() {
                 let main = $(this).parent().parent();
                 let maximize = $(this).find('.ui-icon').hasClass("ui-icon-triangle-1-n");
                 ChangeWindowState(main, maximize ? "maximal" : "normal");
              })
              .next()
              .attr('title','minimize canvas')
              .addClass('jsroot_minbutton')
              .button({ icons: { primary: "ui-icon-triangle-1-s" }, text: false })
              .click(function() {
                 let main = $(this).parent().parent();
                 let minimize = $(this).find('.ui-icon').hasClass("ui-icon-triangle-1-s");
                 ChangeWindowState(main, minimize ? "minimal" : "normal");
              });

         // set default z-index to avoid overlap of these special elements
         $("#" + subid).find(".ui-resizable-handle").css('z-index', '');

         this.cnt++;

         return $("#" + subid + "_cont").attr('frame_title', title).get(0);
      }

   } // class FlexibleDisplay

   // ================== new grid with flexible boundaries ========

   JSROOT.GridDisplay.prototype.createSeparator = function(handle, main, group) {
      let separ = $(main.append("div").node());

      separ.toggleClass('jsroot_separator', true)
           .toggleClass(handle.vertical ? 'jsroot_hline' : 'jsroot_vline', true)
           .prop('handle', handle)
           .attr('separator-id', group.id)
           .css('position','absolute')
           .css(handle.vertical ? 'top' : 'left', "calc(" + group.position+"% - 2px)")
           .css(handle.vertical ? 'width' : 'height', (handle.size || 100)+"%")
           .css(handle.vertical ? 'height' : 'width', '5px')
           .css('cursor', handle.vertical ? "ns-resize" : "ew-resize");

      separ.bind('changePosition', function(e, drag_ui) {
         let handle = $(this).prop('handle'),
             id = parseInt($(this).attr('separator-id')),
             pos = handle.groups[id].position;

         if (drag_ui === 'restore') {
            pos = handle.groups[id].position0;
         } else
         if (drag_ui && drag_ui.offset) {
            if (handle.vertical)
               pos = Math.round((drag_ui.offset.top+2-$(this).parent().offset().top)/$(this).parent().innerHeight()*100);
            else
               pos = Math.round((drag_ui.offset.left+2-$(this).parent().offset().left)/$(this).parent().innerWidth()*100);
         }

         let diff = handle.groups[id].position - pos;

         if (Math.abs(diff)<0.3) return; // if no significant change, do nothing

         // do not change if size too small
         if (Math.min(handle.groups[id-1].size-diff, handle.groups[id].size+diff) < 5) return;

         handle.groups[id-1].size -= diff;
         handle.groups[id].size += diff;
         handle.groups[id].position = pos;

         function SetGroupSize(prnt, grid) {
            let name = handle.vertical ? 'height' : 'width',
                size = handle.groups[grid].size+'%';
            prnt.children("[groupid='"+grid+"']").css(name, size)
                .children(".jsroot_separator").css(name, size);
         }

         $(this).css(handle.vertical ? 'top' : 'left', "calc("+pos+"% - 2px)");

         SetGroupSize($(this).parent(), id-1);
         SetGroupSize($(this).parent(), id);

         if (drag_ui === 'restore') {
            $(this).trigger('resizeGroup', id-1);
            $(this).trigger('resizeGroup', id);
         }
      });

      separ.bind('resizeGroup', function(e, grid) {
         let sel = $(this).parent().children("[groupid='"+grid+"']");
         if (!sel.hasClass('jsroot_newgrid')) sel = sel.find(".jsroot_newgrid");
         sel.each(function() { JSROOT.resize($(this).get(0)); });
      });

      separ.dblclick(function() {
         $(this).trigger('changePosition', 'restore');
      });

      separ.draggable({
         axis: handle.vertical ? "y" : "x",
         cursor: handle.vertical ? "ns-resize" : "ew-resize",
         containment: "parent",
         helper : function() { return $(this).clone().css('background-color','grey'); },
         start: function(/* event,ui */) {
            // remember start position
            let handle = $(this).prop('handle'),
                id = parseInt($(this).attr('separator-id'));
            handle.groups[id].startpos = handle.groups[id].position;
         },
         drag: function(event,ui) {
            $(this).trigger('changePosition', ui);
         },
         stop: function(/* event,ui */) {
            // verify if start position was changed
            let handle = $(this).prop('handle'),
               id = parseInt($(this).attr('separator-id'));
            if (Math.abs(handle.groups[id].startpos - handle.groups[id].position) < 0.5) return;

            $(this).trigger('resizeGroup', id-1);
            $(this).trigger('resizeGroup', id);
         }
      });
   }

   /** @summary Create painter to perform tree drawing on server side
     * @private */
   JSROOT.createTreePlayer = function(player) {

      player.draw_first = true;

      player.ConfigureOnline = function(itemname, url, askey, root_version, dflt_expr) {
         this.setItemName(itemname, "", this);
         this.url = url;
         this.root_version = root_version;
         this.askey = askey;
         this.dflt_expr = dflt_expr;
      }

      player.configureTree = function(tree) {
         this.local_tree = tree;
      }

      player.KeyUp = function(e) {
         if (e.keyCode == 13) this.PerformDraw();
      }

      player.ShowExtraButtons = function(args) {
         let main = $(this.selectDom().node());

          main.find(".treedraw_buttons")
             .append(" Cut: <input class='treedraw_cut ui-corner-all ui-widget' style='width:8em;margin-left:5px' title='cut expression'></input>"+
                     " Opt: <input class='treedraw_opt ui-corner-all ui-widget' style='width:5em;margin-left:5px' title='histogram draw options'></input>"+
                     " Num: <input class='treedraw_number' style='width:7em;margin-left:5px' title='number of entries to process (default all)'></input>" +
                     " First: <input class='treedraw_first' style='width:7em;margin-left:5px' title='first entry to process (default first)'></input>" +
                     " <button class='treedraw_clear' title='Clear drawing'>Clear</button>");

          let numentries = this.local_tree ? this.local_tree.fEntries : 0;

          main.find(".treedraw_cut").val(args && args.parse_cut ? args.parse_cut : "").keyup(this.keyup);
          main.find(".treedraw_opt").val(args && args.drawopt ? args.drawopt : "").keyup(this.keyup);
          main.find(".treedraw_number").val(args && args.numentries ? args.numentries : "").spinner({ numberFormat: "n", min: 0, page: 1000, max: numentries || 0 }).keyup(this.keyup);
          main.find(".treedraw_first").val(args && args.firstentry ? args.firstentry : "").spinner({ numberFormat: "n", min: 0, page: 1000, max: numentries || 0 }).keyup(this.keyup);
          main.find(".treedraw_clear").button().click(() => JSROOT.cleanup(this.drawid));
      }

      player.Show = function(args) {

         let main = $(this.selectDom().node());

         this.drawid = "jsroot_tree_player_" + JSROOT._.id_counter++ + "_draw";

         this.keyup = this.KeyUp.bind(this);

         let show_extra = args && (args.parse_cut || args.numentries || args.firstentry);

         main.html("<div class='treedraw_buttons' style='padding-left:0.5em'>" +
               "<button class='treedraw_exe' title='Execute draw expression'>Draw</button>" +
               " Expr:<input class='treedraw_varexp ui-corner-all ui-widget' style='width:12em;margin-left:5px' title='draw expression'></input> " +
               (show_extra ? "" : "<button class='treedraw_more'>More</button>") +
               "</div>" +
               "<hr/>" +
               "<div id='" + this.drawid + "' style='width:100%'></div>");

         // only when main html element created, one can painter
         // ObjectPainter allow such usage of methods from BasePainter
         this.setTopPainter();

         let p = this;

         if (this.local_tree)
            main.find('.treedraw_buttons').attr('title', "Tree draw player for: " + this.local_tree.fName);
         main.find('.treedraw_exe').button().click(() => p.PerformDraw());
         main.find('.treedraw_varexp')
              .val(args && args.parse_expr ? args.parse_expr : (this.dflt_expr || "px:py"))
              .keyup(this.keyup);

         if (show_extra) {
            this.ShowExtraButtons(args);
         } else {
            main.find('.treedraw_more').button().click(function() {
               $(this).remove();
               p.ShowExtraButtons();
            });
         }

         this.checkResize();
      }

      player.PerformLocalDraw = function() {
         if (!this.local_tree) return;

         let frame = $(this.selectDom().node()),
             args = { expr: frame.find('.treedraw_varexp').val() };

         if (frame.find('.treedraw_more').length==0) {
            args.cut = frame.find('.treedraw_cut').val();
            if (!args.cut) delete args.cut;

            args.drawopt = frame.find('.treedraw_opt').val();
            if (args.drawopt === "dump") { args.dump = true; args.drawopt = ""; }
            if (!args.drawopt) delete args.drawopt;

            args.numentries = parseInt(frame.find('.treedraw_number').val());
            if (isNaN(args.numentries)) delete args.numentries;

            args.firstentry = parseInt(frame.find('.treedraw_first').val());
            if (isNaN(args.firstentry)) delete args.firstentry;
         }

         if (args.drawopt) JSROOT.cleanup(this.drawid);

         let process_result = obj => JSROOT.redraw(this.drawid, obj);

         args.progress = process_result;

         this.local_tree.Draw(args).then(process_result);
      }

      player.PerformDraw = function() {

         if (this.local_tree) return this.PerformLocalDraw();

         let frame = $(this.selectDom().node()),
             url = this.url + '/exe.json.gz?compact=3&method=Draw',
             expr = frame.find('.treedraw_varexp').val(),
             hname = "h_tree_draw", option = "",
             pos = expr.indexOf(">>");

         if (pos<0) {
            expr += ">>" + hname;
         } else {
            hname = expr.substr(pos+2);
            if (hname[0]=='+') hname = hname.substr(1);
            let pos2 = hname.indexOf("(");
            if (pos2>0) hname = hname.substr(0, pos2);
         }

         if (frame.find('.treedraw_more').length==0) {
            let cut = frame.find('.treedraw_cut').val(),
                nentries = frame.find('.treedraw_number').val(),
                firstentry = frame.find('.treedraw_first').val();

            option = frame.find('.treedraw_opt').val();

            url += '&prototype="const char*,const char*,Option_t*,Long64_t,Long64_t"&varexp="' + expr + '"&selection="' + cut + '"';

            // provide all optional arguments - default value kMaxEntries not works properly in ROOT6
            if (nentries=="") nentries = (this.root_version >= 394499) ? "TTree::kMaxEntries": "1000000000"; // kMaxEntries available since ROOT 6.05/03
            if (firstentry=="") firstentry = "0";
            url += '&option="' + option + '"&nentries=' + nentries + '&firstentry=' + firstentry;
         } else {
            url += '&prototype="Option_t*"&opt="' + expr + '"';
         }
         url += '&_ret_object_=' + hname;

         let SubmitDrawRequest = () => {
            JSROOT.httpRequest(url, 'object').then(res => {
               JSROOT.cleanup(this.drawid);
               JSROOT.draw(this.drawid, res, option);
            });
         }

         if (this.askey) {
            // first let read tree from the file
            this.askey = false;
            JSROOT.httpRequest(this.url + "/root.json", 'text').then(SubmitDrawRequest);
         } else {
            SubmitDrawRequest();
         }
      }

      player.checkResize = function(/*arg*/) {
         let main = $(this.selectDom().node());

         $("#" + this.drawid).width(main.width());
         let h = main.height(),
             h0 = main.find(".treedraw_buttons").outerHeight(true),
             h1 = main.find("hr").outerHeight(true);

         $("#" + this.drawid).height(h - h0 - h1 - 2);

         JSROOT.resize(this.drawid);
      }

      return player;
   }

   /** @summary function used with THttpServer to assign player for the TTree object
     * @private */
   JSROOT.drawTreePlayer = function(hpainter, itemname, askey, asleaf) {

      let item = hpainter.findItem(itemname),
          top = hpainter.getTopOnlineItem(item),
          draw_expr = "", leaf_cnt = 0;
      if (!item || !top) return null;

      if (asleaf) {
         draw_expr = item._name;
         while (item && !item._ttree) item = item._parent;
         if (!item) return null;
         itemname = hpainter.itemFullName(item);
      }

      let url = hpainter.getOnlineItemUrl(itemname);
      if (!url) return null;

      let root_version = top._root_version ? parseInt(top._root_version) : 396545; // by default use version number 6-13-01

      let mdi = hpainter.getDisplay();
      if (!mdi) return null;

      let frame = mdi.findFrame(itemname, true);
      if (!frame) return null;

      let divid = d3.select(frame).attr('id'),
          player = new JSROOT.BasePainter(divid);

      if (item._childs && !asleaf)
         for (let n=0;n<item._childs.length;++n) {
            let leaf = item._childs[n];
            if (leaf && leaf._kind && (leaf._kind.indexOf("ROOT.TLeaf")==0) && (leaf_cnt<2)) {
               if (leaf_cnt++ > 0) draw_expr+=":";
               draw_expr+=leaf._name;
            }
         }

      JSROOT.createTreePlayer(player);
      player.ConfigureOnline(itemname, url, askey, root_version, draw_expr);
      player.Show();

      return player;
   }

   /** @summary function used with THttpServer when tree is not yet loaded
     * @private */
   JSROOT.drawTreePlayerKey = function(hpainter, itemname) {
      return JSROOT.drawTreePlayer(hpainter, itemname, true);
   }

   /** @summary function used with THttpServer when tree is not yet loaded
     * @private */
   JSROOT.drawLeafPlayer = function(hpainter, itemname) {
      return JSROOT.drawTreePlayer(hpainter, itemname, false, true);
   }

   return JSROOT;
});
