#!/bin/bash
# Scan libraries in ROOTSYS, and create a collaboration diagram for each. Put those in the
# doxygen html directory to be picked up by the class overview pages.

HTMLPATH=$DOXYGEN_OUTPUT_DIRECTORY/html
DOXYGEN_LDD=${DOXYGEN_LDD:=ldd}
dotFile=$(mktemp /tmp/libraries_XXXX.dot)

test -d "$HTMLPATH" || { echo "HTMLPATH '$HTMLPATH' not found."; exit 1; }
test -d "$ROOTSYS"  || { echo "ROOTSYS not set"; exit 1; }

for libname in ${ROOTSYS}/lib/lib[A-Z]*.so; do
   libsList=$(${DOXYGEN_LDD} ${libname} | grep -v ${libname})
   libname=${libname%.so}
   libname=${libname##*/}

   # Picture name containing the "coll graph"
   PICNAME="${HTMLPATH}/${libname}__coll__graph.svg"

   libsList=$(echo "$libsList" | sed -e "s/\.so.*$/\";/" | grep -v "\.dylib" | grep 'lib[A-Z]' |
	   sed -e "s/\(.*\)\(lib.*;\)$/   mainlib->\"\2/")

   # Generate the dot file describing the graph for libraries
   cat <<EOF > ${dotFile}
digraph G {
   rankdir=TB;
   node [shape=box, fontname=Arial];
   mainlib [label="${libname}"];
   ${libsList}
   mainlib [shape=box, fillcolor="#ABACBA", style=filled];
}
EOF

   # Generate the SVG image of the graph
   dot -Tsvg ${dotFile} -o $PICNAME
done
