#!/bin/sh

# Generates the HTML footer

echo '<!-- HTML footer for doxygen 1.8.14-->'
echo '<!-- start footer part -->'
echo '<!--BEGIN GENERATE_TREEVIEW-->'
echo '<div id="nav-path" class="navpath"><!-- id is needed for treeview function! -->'
echo '  <ul>'
echo '    $navpath'
echo '    <li class="footer">'
echo '    ROOT'$DOXYGEN_ROOT_VERSION' - Reference Guide Generated on $datetime (GVA Time) using Doxygen $doxygenversion &#160;&#160;'
echo '    <img class="footer" src="rootlogo_s.gif" alt="root"/></li>'
echo '  </ul>'
echo '</div>'
echo '<!--END GENERATE_TREEVIEW-->'
echo '<!--BEGIN !GENERATE_TREEVIEW-->'
echo '<hr class="footer"/><address class="footer"><small>'
echo 'ROOT'$DOXYGEN_ROOT_VERSION' - Reference Guide Generated on $datetime (GVA Time) using Doxygen $doxygenversion &#160;&#160;'
echo '<img class="footer" src="rootlogo_s.gif" alt="root"/>'
echo '</small></address>'
echo '<!--END !GENERATE_TREEVIEW-->'
echo '</body>'
echo '</html>'
