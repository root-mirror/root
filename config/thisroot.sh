# Source this script to set up the ROOT build that this script is part of.
#
# Conveniently an alias like this can be defined in .bashrc:
#   alias thisroot=". bin/thisroot.sh"
#
# This script if for the bash like shells, see thisroot.csh for csh like shells.
#
# Author: Fons Rademakers, 18/8/2006

drop_from_path()
{
   # Assert that we got enough arguments
   if test $# -ne 2 ; then
      echo "drop_from_path: needs 2 arguments"
      return 1
   fi

   local p=$1
   local drop=$2

   newpath=`echo $p | sed -e "s;:${drop}:;:;g" \
                          -e "s;:${drop}\$;;g"   \
                          -e "s;^${drop}:;;g"   \
                          -e "s;^${drop}\$;;g"`
}

clean_environment()
{

   if [ -n "${old_rootsys}" ] ; then
      if [ -n "${PATH}" ]; then
         drop_from_path "$PATH" "${old_rootsys}/bin"
         PATH=$newpath
      fi
      if [ -n "${LD_LIBRARY_PATH}" ]; then
         drop_from_path "$LD_LIBRARY_PATH" "${old_rootsys}/lib"
         LD_LIBRARY_PATH=$newpath
      fi
      if [ -n "${DYLD_LIBRARY_PATH}" ]; then
         drop_from_path "$DYLD_LIBRARY_PATH" "${old_rootsys}/lib"
         DYLD_LIBRARY_PATH=$newpath
      fi
      if [ -n "${SHLIB_PATH}" ]; then
         drop_from_path "$SHLIB_PATH" "${old_rootsys}/lib"
         SHLIB_PATH=$newpath
      fi
      if [ -n "${LIBPATH}" ]; then
         drop_from_path "$LIBPATH" "${old_rootsys}/lib"
         LIBPATH=$newpath
      fi
      if [ -n "${PYTHONPATH}" ]; then
         drop_from_path "$PYTHONPATH" "${old_rootsys}/lib"
         PYTHONPATH=$newpath
      fi
      if [ -n "${MANPATH}" ]; then
         drop_from_path "$MANPATH" "${old_rootsys}/man"
         MANPATH=$newpath
      fi
      if [ -n "${CMAKE_PREFIX_PATH}" ]; then
         drop_from_path "$CMAKE_PREFIX_PATH" "${old_rootsys}"
         CMAKE_PREFIX_PATH=$newpath
      fi
      if [ -n "${JUPYTER_PATH}" ]; then
         drop_from_path "$JUPYTER_PATH" "${old_rootsys}/etc/notebook"
         JUPYTER_PATH=$newpath
      fi
      if [ -n "${JUPYTER_CONFIG_DIR}" ]; then
         drop_from_path "$JUPYTER_CONFIG_DIR" "${old_rootsys}/etc/notebook"
         JUPYTER_CONFIG_DIR=$newpath
      fi
   fi
   if [ -z "${MANPATH}" ]; then
      # Grab the default man path before setting the path to avoid duplicates
      if command -v manpath >/dev/null; then
         default_manpath=`manpath`
      elif command -v man >/dev/null; then
         default_manpath=`man -w 2> /dev/null`
      else
         default_manpath=""
      fi
   fi
}

set_environment()
{
   if [ -z "${PATH}" ]; then
      PATH=@bindir@; export PATH
   else
      PATH=@bindir@:$PATH; export PATH
   fi

   if [ -z "${LD_LIBRARY_PATH}" ]; then
      LD_LIBRARY_PATH=@libdir@
      export LD_LIBRARY_PATH       # Linux, ELF HP-UX
   else
      LD_LIBRARY_PATH=@libdir@:$LD_LIBRARY_PATH
      export LD_LIBRARY_PATH
   fi

   if [ -z "${DYLD_LIBRARY_PATH}" ]; then
      DYLD_LIBRARY_PATH=@libdir@
      export DYLD_LIBRARY_PATH       # Linux, ELF HP-UX
   else
      DYLD_LIBRARY_PATH=@libdir@:$DYLD_LIBRARY_PATH
      export DYLD_LIBRARY_PATH
   fi

   if [ -z "${SHLIB_PATH}" ]; then
      SHLIB_PATH=@libdir@
      export SHLIB_PATH       # Linux, ELF HP-UX
   else
      SHLIB_PATH=@libdir@:$SHLIB_PATH
      export SHLIB_PATH
   fi

   if [ -z "${LIBPATH}" ]; then
      LIBPATH=@libdir@
      export LIBPATH       # Linux, ELF HP-UX
   else
      LIBPATH=@libdir@:$LIBPATH
      export LIBPATH
   fi

   if [ -z "${PYTHONPATH}" ]; then
      PYTHONPATH=@libdir@
      export PYTHONPATH       # Linux, ELF HP-UX
   else
      PYTHONPATH=@libdir@:$PYTHONPATH
      export PYTHONPATH
   fi

   if [ -z "${MANPATH}" ]; then
      MANPATH=@mandir@:${default_manpath}; export MANPATH
   else
      MANPATH=@mandir@:$MANPATH; export MANPATH
   fi

   if [ -z "${CMAKE_PREFIX_PATH}" ]; then
      CMAKE_PREFIX_PATH=$ROOTSYS; export CMAKE_PREFIX_PATH       # Linux, ELF HP-UX
   else
      CMAKE_PREFIX_PATH=$ROOTSYS:$CMAKE_PREFIX_PATH; export CMAKE_PREFIX_PATH
   fi

   if [ -z "${JUPYTER_PATH}" ]; then
      JUPYTER_PATH=$ROOTSYS/etc/notebook; export JUPYTER_PATH       # Linux, ELF HP-UX
   else
      JUPYTER_PATH=$ROOTSYS/etc/notebook:$JUPYTER_PATH; export JUPYTER_PATH
   fi

   if [ -z "${JUPYTER_CONFIG_DIR}" ]; then
      JUPYTER_CONFIG_DIR=$ROOTSYS/etc/notebook; export JUPYTER_CONFIG_DIR # Linux, ELF HP-UX
   else
      JUPYTER_CONFIG_DIR=$ROOTSYS/etc/notebook:$JUPYTER_CONFIG_DIR; export JUPYTER_CONFIG_DIR
   fi
}


### main ###


if [ -n "${ROOTSYS}" ] ; then
   old_rootsys=${ROOTSYS}
fi


SOURCE=${BASH_ARGV[0]}
if [ "x$SOURCE" = "x" ]; then
   SOURCE=${(%):-%N} # for zsh
fi


if [ "x${SOURCE}" = "x" ]; then
   if [ -f bin/thisroot.sh ]; then
      ROOTSYS="$PWD"; export ROOTSYS
   elif [ -f ./thisroot.sh ]; then
      ROOTSYS=$(cd ..  > /dev/null; pwd); export ROOTSYS
   else
      echo ERROR: must "cd where/root/is" before calling ". bin/thisroot.sh" for this version of bash!
      ROOTSYS=; export ROOTSYS
      return 1
   fi
else
   # get param to "."
   thisroot=$(dirname ${SOURCE})
   ROOTSYS=$(cd ${thisroot}/.. > /dev/null;pwd); export ROOTSYS
fi


clean_environment
set_environment


# Prevent Cppyy from checking the PCH (and avoid warning)
export CLING_STANDARD_PCH=none

if [ "x`root-config --arch | grep -v win32gcc | grep -i win32`" != "x" ]; then
   ROOTSYS="`cygpath -w $ROOTSYS`"
fi

unset old_rootsys
unset thisroot
unset -f drop_from_path
unset -f clean_environment
unset -f set_environment
