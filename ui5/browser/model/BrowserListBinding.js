sap.ui.define([
    "sap/ui/model/json/JSONListBinding"
], function(JSONListBinding) {
    "use strict";

    var hRootListBinding = JSONListBinding.extend("rootui5.browser.model.BrowserListBinding", {

        // called by the TreeTable to know the amount of entries
        getLength: function() {
           return this.getModel().getLength();
        },

        // function is called by the TreeTable when requesting the data to display
        getNodes: function(iStartIndex, iLength, iThreshold) {

           var args = { begin: iStartIndex, end: iStartIndex + iLength, threshold: iThreshold },
               nodes = this.getModel().buildFlatNodes(args);

           var aNodes = [];

           for (var i = args.begin; i < args.end; i++)
              aNodes.push(nodes && nodes[i] ? nodes[i] : null);

           return aNodes;
        },

        getContextByIndex: function(iIndex) {
           return this.getModel().getContext(this.getPath() + "/" + iIndex);
        },

        findNode: function() {
        },

        nodeHasChildren: function(oNode) {
           return oNode.type === "folder";
        },

        isExpanded: function(iIndex) {
            var elem = this.getModel().getElementByIndex(iIndex);
            var res = elem ? !!elem.expanded : false;

            return res;
        },

        expand: function(iIndex) {
        },

        collapse: function(iIndex) {
        },

        collapseToLevel: function(lvl) {
           // console.log('root.model.hListBinding#collapseToLevel', lvl);
        },

        expandToLevel: function(lvl) {
           // console.log('root.model.hListBinding#expandToLevel', lvl);
        },

        // called by the TreeTable when a node is expanded/collapsed
        toggleIndex: function(iIndex) {
            // console.log("root.model.hListBinding#toggleIndex(" + iIndex + ")");

            if (this.getModel().toggleNode(iIndex))
               this.checkUpdate(true);

            // QUESTION: why one should call checkUpdate?, should it be done automatically always?
        },

        getSelectedIndex: function() {
        },

        isIndexSelectable: function() {
        },

        attachSelectionChanged: function() {
           // dummy for compatibility with newest 1.70.0 version
        }

    });

    return hRootListBinding;

});