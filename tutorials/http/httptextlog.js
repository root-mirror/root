/// \file
/// \ingroup tutorial_http
/// JavaScript code for drawing TMsgList class from httptextlog.C macro
///
/// \macro_code
///
/// \author  Sergey Linev

JSROOT.require("painter").then(jsrp => {

   function MakeMsgListRequest(hitem, item) {
      // this function produces url for http request
      // here one provides id of last string received with previous request

      var arg = "&max=1000";
      if ('last-id' in item) arg+= "&id="+item['last-id'];
      return 'exe.json.gz?method=Select' + arg;
   }

   function AfterMsgListRequest(hitem, item, obj) {
      // after data received, one replaces typename for produced object

      if (item==null) return;

      if (obj==null) {
         delete item['last-id'];
         return;
      }
      // ignore all other classes
      if (obj._typename != 'TList') return;

      // change class name - it is only important for drawing
      obj._typename = "TMsgList";

      if (obj.arr.length>0) {
         item['last-id'] = obj.arr[0].fString;

         // add clear function for item
         if (!('clear' in item))
            item['clear'] = function() { delete this['last-id']; }
      }
   }


   function DrawMsgList(divid, lst, opt) {

      let painter = new JSROOT.BasePainter(divid);

      painter.Draw = function(lst) {
         if (!lst) return;

         let frame = this.selectDom();

         let main = frame.select("div");
         if (main.empty()) {
            main = frame.append("div")
                        .style('max-width','100%')
                        .style('max-height','100%')
                        .style('overflow','auto');
            // (re) set painter to first child element
            this.setTopPainter();
         }

         let old = main.selectAll("pre");
         let newsize = old.size() + lst.arr.length - 1;

         // in the browser keep maximum 1000 entries
         if (newsize > 1000)
            old.select(function(d,i) { return i < newsize - 1000 ? this : null; }).remove();

         for (let i = lst.arr.length - 1; i > 0; i--)
            main.append("pre").style('margin','2px').html(lst.arr[i].fString);
      }

      painter.redrawObject = function(obj) {
         this.Draw(obj);
         return true;
      }

      painter.Draw(lst);
      return Promise.resolve(painter);
   }

   // register draw function to JSROOT
   jsrp.addDrawFunc({
      name: "TMsgList",
      icon: "img_text",
      make_request: MakeMsgListRequest,
      after_request: AfterMsgListRequest,
      func: DrawMsgList,
      opt:"list"
   });

})
