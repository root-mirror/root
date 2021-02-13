# Source this script to set up the ROOT build that this script is part of.
#
# This script is for the fish shell.
#
# Author: Axel Naumann, 2018-06-25

function update_path -d "Remove argv[2]argv[3] from argv[1] if argv[2], and prepend argv[4]"
   # Assert that we got enough arguments
   if test (count $argv) -ne 4
      echo "update_path: needs 4 arguments but have " (count $argv)
      return 1
   end

   set var $argv[1]

   set newpath $argv[4]
   for el in $$var
      if test "$argv[2]" = ""; or not test "$el" = "$argv[2]$argv[3]"
         set newpath $newpath $el
      end
   end

   set -xg $var $newpath
end

if set -q ROOTSYS
   set old_rootsys $ROOTSYS
end

set SOURCE (status -f)
# normalize path
set thisroot (dirname $SOURCE)
set -xg ROOTSYS (set oldpwd $PWD; cd $thisroot/.. > /dev/null;pwd;cd $oldpwd; set -e oldpwd)

if not set -q MANPATH
   # Grab the default man path before setting the path to avoid duplicates
   if which manpath > /dev/null ^ /dev/null
      set -xg MANPATH (manpath)
   else
      set -xg MANPATH (man -w 2> /dev/null)
   end
end

update_path PATH "$old_rootsys" "/bin" @bindir@
update_path LD_LIBRARY_PATH "$old_rootsys" "/lib" @libdir@
update_path DYLD_LIBRARY_PATH "$old_rootsys" "/lib" @libdir@
update_path PYTHONPATH "$old_rootsys" "/lib" @libdir@
update_path MANPATH "$old_rootsys" "/man" @mandir@
update_path CMAKE_PREFIX_PATH "$old_rootsys" "" $ROOTSYS
update_path JUPYTER_PATH "$old_rootsys" "/etc/notebook" ROOTSYS/etc/notebook
update_path JUPYTER_CONFIG_DIR "$old_rootsys" "/etc/notebook" ROOTSYS/etc/notebook

# Prevent Cppyy from checking the PCH (and avoid warning)
set -xg CLING_STANDARD_PCH none

functions -e update_path
set -e old_rootsys
set -e thisroot
set -e SOURCE
