/**
 * This is the JsMVA object, the JavaScript part of JPyInterface module. This is a bridge,
 * python code always produce an output, where one of this module's function will be called. These function calls
 * will redirect the date for specific visualization scripts.
 * It uses the following libraries:
 *                                - JsROOT
 *                                - d3js
 *                                - jquery
 * It uses the following submodules:
 *                                - NeuralNetwork
 *                                - DecisionTree
 *                                - NetworkDesigner
 * Author: Attila Bagoly <battila93@gmail.com>
 * Created: 5/14/16
 */

(function(factory){

    var JSROOT_source_dir = "https://root.cern.ch/js/notebook/scripts/";

    var url = "";
    if (requirejs.s.contexts.hasOwnProperty("_")) {
        url = requirejs.s.contexts._.config.paths["JsMVA"].replace("JsMVA.min","");
    }
    if ((console!==undefined) && (typeof console.log == 'function')) {
        if (url!=""){
            console.log("JsMVA source dir:" + url.substring(0, url.length-1));
        } else {
            console.log("JsMVA source dir can't be resolved, requireJS doesn't have context '_', this will be a problem!");
        }
    }

    require.config({
        paths: {
            'd3': JSROOT_source_dir+'d3.min',
            'JsRootCore': JSROOT_source_dir+'JSRoot.core.min',
            'nn': url+'NeuralNetwork.min',
            'dtree': url+'DecisionTree.min',
            'NetworkDesigner': url+'NetworkDesigner.min'
        }
    });

    define(['JsRootCore'], function(jsroot){
        return factory({}, jsroot);
    });

}(function(JsMVA, JSROOT){

    JsMVA.drawTH2 = function(divid, dat_json){
        var obj = JSROOT.parse(dat_json);
        JSROOT.draw(divid, obj, "colz;PAL50;text");
    };

    JsMVA.drawDNNMap = function(divid, dat_json){
        var obj = JSROOT.parse(dat_json);
        JSROOT.draw(divid, obj, "colz;PAL50");
    };

    JsMVA.draw = function(divid, dat_json){
        var obj = JSROOT.parse(dat_json);
        JSROOT.draw(divid, obj);
    };

    JsMVA.drawNeuralNetwork = function(divid, dat_json){
        var obj = JSON.parse(dat_json);
        require(['nn'], function(nn){
            nn.draw(divid, obj);
        });
    };

    JsMVA.drawDecisionTree = function(divid, dat_json){
        require(['dtree'], function(dtree){
            var obj = JSON.parse(dat_json);
            dtree.draw(divid, obj);
        });
    };

    var drawLabel = function(divid, obj){
        require(['d3'], function(d3){
            var csvg = d3.select("#"+divid+">.interactivePlot_Labels")[0][0];
            if (csvg!=null) return;
            var div = d3.select("#"+divid).style("position", "relative");
            var svg = div.append("svg").attr("class", "interactivePlot_Labels")
                .attr("width", "200px")
                .attr("height", "50px")
                .style({"position":"absolute", "top": "8px", "right": "8px"});
            var attr = {
                "pos": {"x": 150, "y": 0},
                "rect": {"width": 10, "height":10},
                "dy": 20,
                "padding": 10
            };
            canvas = {
                width:  160,
                height: 70
            };
            var container = svg.append("g").attr("id", "legend");
            container.selectAll("g")
                .data(obj.fGraphs.arr)
                .enter()
                .append("g")
                .each(function(d, i){
                    var g = d3.select(this);
                    g.append("rect")
                        .attr("x", canvas.width-attr.pos.x)
                        .attr("y", attr.pos.y+i*attr.dy)
                        .attr("width", attr.rect.width)
                        .attr("height", attr.rect.height)
                        .style("fill", function(d){return JSROOT.Painter.root_colors[d.fFillColor];});
                    g.append("text")
                        .attr("x", canvas.width-attr.pos.x+attr.rect.width+attr.padding)
                        .attr("y", attr.pos.y+i*attr.dy+attr.rect.height)
                        .text(function(d){return d.fTitle;})
                        .style("fill", function(d){return JSROOT.Painter.root_colors[d.fFillColor];});
                });
            div.append("svg").attr("width", "55px").attr("height", "20px")
                .style({"position":"absolute", "bottom": "15px", "right": "40px"})
                .append("text")
                .attr("x", "5px")
                .attr("y", "15px")
                .text(obj.fGraphs.arr[0].fTitle.indexOf("Error on training set")!=-1 ? "Epoch" : "#tree")
                .style({"font-size": "16px"});
        });
    };


    JsMVA.drawTrainingTestingErrors = function(divid, dat_json){
        var obj = JSROOT.parse(dat_json);
        JSROOT.draw(divid, obj);
        drawLabel(divid, obj);
    };

    JsMVA.updateTrainingTestingErrors = function(divid, dat_json){
        var obj = JSROOT.parse(dat_json);
        JSROOT.redraw(divid, obj);
        drawLabel(divid, obj);
    };

    JsMVA.NetworkDesigner = function(divid, dat_json){
      require(['NetworkDesigner'], function (nd) {
         nd.draw(divid);
      });
    };

    JsMVA.outputShowCorrelationMatrix = function(divid){
        require(['jquery', 'jquery-ui'], function($){
            var th2 = JSROOT.parse($("#"+divid).html());
            if (!$("#dialog_"+divid).length || $("#dialog_"+divid).length < 1) {
                $("#" + divid).parent().append("<div id='dialog_" + divid + "' title='" + th2.fTitle + "' style='width: 600px; height: 340px; z-index: 99;'></div>");
                JSROOT.draw("dialog_" + divid, th2, "colz;PAL50;text");
            }
            $("#dialog_" + divid).dialog({
                autoOpen: true,
                width: 600,
                show: {
                    effect: "blind",
                    duration: 1000
                },
                hide: {
                    effect: "explode",
                    duration: 500
                }
            });

        });
    };

    return JsMVA;
}));
