#!/bin/bash
#
# PROOF setup script.
# This script distributes the ROOT binary on a cluster.
# Unpacks the binary, and start the necessary xrootd daemon.
#
# Author: Fons Rademakers

# On which platform do we run
ostype=`uname -s`

cwhich="type -path"
lsopt="--color=never"
ptimeout="-W"
if [ "$ostype" = "Darwin" ]; then
   lsopt=""
   ptimeout="-t"
fi

# Where to store values from a previous run
saverc="$HOME/.proofinstallrc"

#--- stored in $saverc
nodes=
domain=`hostname | awk -F '.' '{ print $2 "." $3 }'`
base=
lower=0
upper=0
width=1
rdomain="$domain"
rootbinary=
remoterootbinary=
autoremoterootbinary=
puser="$LOGNAME"
pmaster=
pport=1093
pdatadir="proof/data"
pworkers=0
pconfig=
#--- stored in $saverc

prompt_rootrepository="ftp://root.cern/root/"

prog=`basename $0`
progdir=`dirname $0`
runningdir=`pwd`
if echo $progdir | grep "^/" > /dev/null 2>& 1 || \
   echo $progdir | grep "^~" > /dev/null 2>& 1; then
   # absolute path
   fullpath=$progdir
else
   # relative path
   fullpath=$runningdir/$progdir
fi

# Helper programs
circle=$fullpath/circle.sh

# Set default argument values
batch="no"
ignore="no"
cluster=

# Process command line arguments
while [ $# -gt 0 ]; do
   case "$1" in
   -b|--batch)
      batch="yes"
      shift
      ;;
   -i|--ignore)
      ignore="yes"
      shift
      ;;
   -h|--help)
      echo "Usage: $prog [-b] [-i] [-h] [CLUSTER]"
      echo "Install and start PROOF on a cluster."
      echo ""
      echo "  -b, --batch    run in batch mode, requires a"
      echo "                 $saverc file, created"
      echo "                 by an interactive run"
      echo "                 mutual exclusive with the --ignore option"
      echo "  -i, --ignore   ignore an existing $saverc file"
      echo "                 mutual exclusive with the --batch option"
      echo "  -h, --help     display this help and exit"
      echo ""
      echo "CLUSTER, store and use the config for a specific cluster."
      exit 0
      ;;
   -*)
      echo "$prog: invalid option -- $1"
      echo "Try \`$prog --help' for more information."
      exit 1
      ;;
   *)
      cluster="$1"
      saverc="${saverc}-${cluster}"
      shift
      ;;
   esac
done

if [ "$ignore" = "yes" ] && [ "$batch" = "yes" ]; then
   echo "$prog: the options --ignore and --batch are mutual exclusive."
   exit 1
fi

if [ ! -r "$saverc" ] && [ "$batch" = "yes" ]; then
   echo "$prog: need a $saverc file with"
   echo "parameters to be able to run in batch mode."
   echo "First run interactively to create one, like"
   echo "\"$prog $cluster\"."
   exit 1
fi

if [ -r "$saverc" ] && [ "$ignore" = "no" ]; then
   source $saverc
fi


# Setup ssh-agent and load key using ssh-add subroutine.
### To be added.

if [ "$batch" = "no" ]; then

# Print welcome and instructions.
clear
echo ""
echo "              Welcome to the PROOF Installer"
echo ""
echo "To install PROOF on a cluster the following is needed:"
echo ""
echo "   - list of host names"
echo "   - ssh login capability on these machines"
echo "   - a binary ROOT distribution compatible with the cluster,"
echo "     or a download location (used by curl or wget)"
echo ""
echo "If any of the above is not true, please interrupt the installer"
echo "and re-run later."
echo ""
echo "If more than one answer is prompted, the one in upper case is the default."
echo ""
echo -n "Are you ready to continue [Y/n]: "
read answ &> /dev/null
if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
   echo "Exiting the install..."
   exit 0
fi

# Ask user for list of nodes to be used.
REPEAT_LOOP=1
while [ "$REPEAT_LOOP" = 1 ]; do
   clear
   echo ""
   echo "Give list of host names."
   echo ""
   echo "Please make your selection by entering an option:"
   echo ""
   echo "     1. Individual host names."
   echo "     2. A range of host names."
   echo "     h. Help."
   echo "     x. Exit."
   echo ""
   echo -n "Please type a selection: "
   read answ &> /dev/null
   case $answ in
   1)
      clear
      echo ""
      echo -n "Give list of host names [$nodes]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         nodes="$answ"
      fi
      if [ "x$nodes" = "x" ]; then
          echo -n "No hosts specified, try again [Y/n]: "
          read answ &> /dev/null
          if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
             echo "Exiting! Press <ENTER> to terminate install."
             read
             exit 0
          fi
          continue
      fi

      echo ""
      echo -n "Give domain name, if not already specified (- = no domain) [$domain]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         domain=`echo $answ | sed "s/^\.//"`
      fi
      if [ "$domain" != "-" ]; then
         n=
         for i in $nodes; do
            if [ "x$n" != "x" ]; then
               n="$n "
            fi
            n="$n$i.$domain"
         done
         nodes="$n"
      fi

      echo ""
      echo "The list of hosts is:"
      echo ""
      echo $nodes
      echo ""
      echo -n "Is this correct [Y/n]: "
      read answ &> /dev/null
      if [ "x$answ" = "x" ] || [ "$answ" = "Y" ] || [ "$answ" = "y" ]; then
         REPEAT_LOOP=0
      fi
      ;;
   2)
      nodes=
      clear
      echo ""
      echo -n "Give host base name [$base]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         base="$answ"
      fi
      if [ "x$base" = "x" ]; then
         echo -n "No base specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi

      echo ""
      echo -n "Give lower bound [$lower]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         lower="$answ"
      fi

      echo ""
      echo -n "Give upper bound [$upper]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         upper="$answ"
      fi
      if [ $lower -gt $upper ]; then
         echo -n "lower ($lower) > upper ($upper), try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
             echo "Exiting! Press <ENTER> to terminate install."
             read
             exit 0
          fi
          continue
      fi

      echo ""
      echo -n "Give width [$width]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         width="$answ"
      fi

      echo ""
      echo -n "Give domain name, if not already specified (- = no domain) [$rdomain]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         rdomain=`echo $answ | sed "s/^\.//"`
      fi

      form="%0${width}d"

      i=$lower
      while [ $i -le $upper ]; do
         ii=`printf $form $i`
         n="$base$ii"
         if [ "$rdomain" != "-" ]; then
            n="$n.$rdomain"
         fi
         if [ "x$nodes" != "x" ]; then
            nodes="$nodes "
         fi
         nodes="$nodes$n"
         i=$((i+1))
      done

      echo ""
      echo "The list of hosts is:"
      echo ""
      echo $nodes
      echo ""
      echo -n "Is this correct [Y/n]: "
      read answ &> /dev/null
      if [ "x$answ" = "x" ] || [ "$answ" = "Y" ] || [ "$answ" = "y" ]; then
         REPEAT_LOOP=0
      fi
      ;;
   h)
      echo "---------------------------------------------------------------------------"
      echo "Node selection options:"
      echo ""
      echo "-- Option 1 --"
      echo "Give a list of space separated host names. If they are in the same domain"
      echo "don't specify the domains. You will be asked later for the domain."
      echo "Example: thisnode1 thatnode2 ..."
      echo ""
      echo "-- Option 2 --"
      echo "Give the host base name, the begin and end of the range and the width."
      echo "Example:"
      echo "   node- 1 20 3"
      echo "Defines:"
      echo "   node-001 node-002 node-003 ... node-020"
      echo ""
      echo "-- Option h --"
      echo "Selecting option h displays this help message."
      echo " "
      echo "-- Option x --"
      echo "Selecting option x exits the PROOF installation."
      echo "---------------------------------------------------------------------------"
      echo -n "Press <Enter> to continue..."
      read answ &> /dev/null
      ;;
   x)
      echo "Exiting the install..."
      exit 0
      ;;
   *)
      echo "Invalid Choice. Please try again."
      REPEAT_LOOP=1
      ;;
   esac

   if [ "$REPEAT_LOOP" = 0 ]; then
      # Ask if we can ping the nodes
      echo ""
      echo -n "Check if hosts are alive using ping [Y/n]: "
      read answ &> /dev/null
      if [ "x$answ" = "x" ] || [ "$answ" = "Y" ] || [ "$answ" = "y" ]; then
         n=
         for i in $nodes; do
            ping -c 1 $ptimeout 1 $i > /dev/null 2>&1
            if [ $? -eq 0 ]; then
               if [ "x$n" != "x" ]; then
                  n="$n "
               fi
               n=$n$i
            else
              echo "Host $i not alive"
            fi
         done
         nodes=$n
      fi

      if [ "x$nodes" = "x" ]; then
          echo -n "No hosts alive, try again [Y/n]: "
          read answ &> /dev/null
          if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
             echo "Exiting! Press <ENTER> to terminate install."
             read
             exit 0
          fi
          REPEAT_LOOP=1
          continue
      fi

      echo ""
      echo "The list of available nodes is:"
      echo ""
      echo $nodes
      echo ""
      echo -n "Is this correct [Y/n]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
         REPEAT_LOOP=1
      fi
   fi
done

# Ask user which ROOT binary distribution to use.
# Either a local file, a remote file or a remote repository can be specified.
REPEAT_LOOP=1
while [ "$REPEAT_LOOP" = 1 ]; do
   bintype=0
   clear
   echo ""
   echo "Give ROOT binary distribution to be installed."
   echo ""
   echo "Please make your selection by entering an option:"
   echo ""
   echo "     1. Local ROOT binary distribution file."
   echo "     2. Remote ROOT binary distribution file."
   echo "     3. Remote ROOT binary distribution automatically"
   echo "        matching cluster node architecture."
   echo "     h. Help."
   echo "     x. Exit."
   echo ""
   echo -n "Please type a selection: "
   read answ &> /dev/null
   case $answ in
   1)
      clear
      echo ""
      echo -n "Give local ROOT binary distribution [$rootbinary]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         rootbinary=`eval echo $answ`
      fi
      if [ "x$rootbinary" = "x" ]; then
         echo -n "No local ROOT binary distribution file specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi
      REPEAT_LOOP=0
      bintype=1
      ;;
   2)
      clear
      echo ""
      echo -n "Give remote ROOT binary distribution [$remoterootbinary]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         remoterootbinary=`eval echo $answ`
      fi
      if [ "x$remoterootbinary" = "x" ]; then
         echo -n "No remote ROOT binary distribution file specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi
      REPEAT_LOOP=0
      bintype=2
      ;;
   3)
      clear
      echo ""
      echo -n "Give URL to ROOT binary repository [$prompt_rootrepository]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         autoremoterootbinary="$answ"
      else
         autoremoterootbinary="$prompt_rootrepository"
      fi
      if [ "x$autoremoterootbinary" = "x" ]; then
         echo -n "No URL to ROOT binary repository specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi
      REPEAT_LOOP=0
      bintype=3
      ;;
   h)
      echo "---------------------------------------------------------------------------"
      echo "ROOT binary distribution selection options:"
      echo ""
      echo "-- Option 1 --"
      echo "Give path to local ROOT binary distribution tar file."
      echo "This only makes sense if all the cluster nodes have the same architecture."
      echo "Example: ~/root_v5.17.04.Linux.slc4_amd64.gcc3.4.tar.gz"
      echo ""
      echo "-- Option 2 --"
      echo "Give URL to remote ROOT binary distribution tar file."
      echo "This only makes sense if all the cluster nodes have the same architecture."
      echo "Requires curl or wget on this machine."
      echo "Example: ftp://root.cern/root/root_v5.17.04.Linux.slc4_amd64.gcc3.4.tar.gz"
      echo ""
      echo "-- Option 3 --"
      echo "Give URL to remote ROOT binary distribution repository."
      echo "The remote nodes will try to download the binary version matching"
      echo "their architecture. Requires curl or wget on the remote nodes."
      echo "Example: ftp://root.cern/root/"
      echo ""
      echo "-- Option h --"
      echo "Selecting option h displays this help message."
      echo " "
      echo "-- Option x --"
      echo "Selecting option x exits the PROOF installation."
      echo "---------------------------------------------------------------------------"
      echo -n "Press <Enter> to continue..."
      read answ &> /dev/null
      ;;
   x)
      echo "Exiting the install..."
      exit 0
      ;;
   *)
      echo "Invalid Choice. Please try again."
      REPEAT_LOOP=1
      ;;
   esac

   if [ "$REPEAT_LOOP" = 0 ]; then
      echo ""
      case $bintype in
      1)
         echo "The local ROOT binary distribution to be used is:"
         echo ""
         echo $rootbinary
         ;;
      2)
         echo "The remote ROOT binary distribution to be used is:"
         echo ""
         echo $remoterootbinary
         ;;
      3)
         echo "The remote ROOT binary distribution repository to be used is:"
         echo ""
         echo $autoremoterootbinary
         ;;
      esac
      echo ""
      echo -n "Is this correct [Y/n]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
         REPEAT_LOOP=1
         continue
      fi

      case $bintype in
      1)
         if [ ! -r "$rootbinary" ]; then
            echo ""
            echo "File $rootbinary not found."
            echo -n "Try again [Y/n]: "
            read answ &> /dev/null
            if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
               echo "Exiting! Press <ENTER> to terminate install."
               read
               exit 0
            fi
            REPEAT_LOOP=1
         fi
         ;;
      2)
         if `$cwhich curl > /dev/null 2>&1`; then
            echo ""
            echo "Downloading $remoterootbinary..."
            echo ""
            rootbinary=`basename $remoterootbinary`.$$
            curl -o $rootbinary $remoterootbinary
            if [ "$?" -ne 0 ]; then
               echo ""
               echo "Remote file $remoterootbinary not found."
               echo -n "Try again [Y/n]: "
               read answ &> /dev/null
               if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
                  echo "Exiting! Press <ENTER> to terminate install."
                  read
                  exit 0
               fi
               REPEAT_LOOP=1
            fi
         elif `$cwhich wget > /dev/null 2>&1`; then
            echo ""
            echo "Downloading $remoterootbinary..."
            echo ""
            rootbinary=`basename $remoterootbinary`.$$
            wget -O $rootbinary $remoterootbinary
            if [ "$?" -ne 0 ]; then
               echo ""
               echo "Remote file $remoterootbinary not found."
               echo -n "Try again [Y/n]: "
               read answ &> /dev/null
               if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
                  echo "Exiting! Press <ENTER> to terminate install."
                  read
                  exit 0
               fi
               REPEAT_LOOP=1
            fi
         else
            echo ""
            echo "Neither curl nor wget available. Try getting the ROOT binary"
            echo "distribution file in some other way."
            echo ""
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
      esac
   fi
done

# Ask the user PROOF specific configuration parameters.
REPEAT_LOOP=1
while [ "$REPEAT_LOOP" = 1 ]; do
   conftype=0
   clear
   echo ""
   echo "Give PROOF specific configuration parameters."
   echo ""
   echo "Please make your selection by entering an option:"
   echo ""
   echo "     1. Basic PROOF configuration parameters."
   echo "     2. Advanced PROOF configuration parameters."
   echo "     h. Help."
   echo "     x. Exit."
   echo ""
   echo -n "Please type a selection: "
   read answ &> /dev/null
   case $answ in
   1)
      clear
      echo ""
      echo "Give basic PROOF specific configuration parameters."
      echo ""
      echo -n "Give name of user under which to run PROOF [$puser]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         puser="$answ"
      fi
      if [ "x$puser" = "x" ]; then
         echo -n "No user specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
             echo "Exiting! Press <ENTER> to terminate install."
             read
             exit 0
          fi
          continue
      fi

      echo ""
      nodearr=($nodes)
      pmaster=${nodearr[0]}
      echo -n "Give host on which to run the PROOF master [$pmaster]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         pmaster="$answ"
      fi
      if [ "x$pmaster" = "x" ]; then
         echo -n "No PROOF master host specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi

      echo ""
      echo -n "Give port on which the PROOF xrootd daemons listen [$pport]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         pport="$answ"
      fi
      if [ "x$pport" = "x" ]; then
         echo -n "No port specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi

      echo ""
      echo -n "Give the data file upload directory [$pdatadir]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         pdatadir=`eval echo $answ`
      fi
      if [ "x$pdatadir" = "x" ]; then
         echo -n "No data upload directory specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi

      echo ""
      echo -n "Give number of workers per host (0 = 1 worker per core) [$pworkers]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         pworkers="$answ"
      fi
      if [ "x$pworkers" = "x" ]; then
         echo -n "No workers per host specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi

      REPEAT_LOOP=0
      conftype=1
      ;;
   2)
      clear
      echo ""
      echo "Give advanced PROOF specific configuration parameters."
      echo ""
      echo -n "Give name of custom PROOF configuration file [$pconfig]: "
      read answ &> /dev/null
      if [ "x$answ" != "x" ]; then
         pconfig=`eval echo $answ`
      fi
      if [ "x$pconfig" = "x" ]; then
         echo -n "No config file specified, try again [Y/n]: "
         read answ &> /dev/null
         if [ "x$answ" != "x" ] && [ "$answ" != "Y" ] && [ "$answ" != "y" ]; then
            echo "Exiting! Press <ENTER> to terminate install."
            read
            exit 0
         fi
         continue
      fi

      REPEAT_LOOP=0
      conftype=2
      ;;
   h)
      echo "---------------------------------------------------------------------------"
      echo "PROOF specific configuration parameters:"
      echo ""
      echo "-- Option 1 --"
      echo "Give the following basic configuration parameters:"
      echo "   - the user name under which to run PROOF, in general this is the same"
      echo "     user name as used for ssh access."
      echo "   - host on which to run the PROOF master, you will be prompted with"
      echo "     the first host in the list of hosts you specified."
      echo "   - port on which the PROOF xrootd daemons will listen. Give some high"
      echo "     port number like 210210 or 9191. You will be prompted with a"
      echo "     port that is not being used on any of the cluster nodes."
      echo "   - data file upload directory, by default this is \"proof/data\","
      echo "     relative to your remote home directory."
      echo "   - number of workers per host, the default is as many workers"
      echo "     as there are cores (0 = 1 worker per core)."
      echo ""
      echo "-- Option 2 --"
      echo "Give a custom PROOF configuration file containing your own advanced"
      echo "settings, e.g. \"myxpd.conf\"."
      echo ""
      echo "-- Option h --"
      echo "Selecting option h displays this help message."
      echo " "
      echo "-- Option x --"
      echo "Selecting option x exits the PROOF installation."
      echo "---------------------------------------------------------------------------"
      echo -n "Press <Enter> to continue..."
      read answ &> /dev/null
      ;;
   x)
      echo "Exiting the install..."
      exit 0
      ;;
   *)
      echo "Invalid Choice. Please try again."
      REPEAT_LOOP=1
      ;;
   esac

   if [ "$REPEAT_LOOP" = 0 ]; then
      echo ""
      echo "The following PROOF configuration parameters have been specified:"
      echo ""
      case $conftype in
      1)
         echo "   - Name of user under which to run PROOF:        $puser"
         echo "   - Node on which to run the PROOF master:        $pmaster"
         echo "   - Port on which the PROOF daemons listen:       $pport"
         echo "   - Data file upload directory:                   $pdatadir"
         echo "   - Number of workers per host (0 = 1 per core):  $pworkers"
         ;;
      2)
         echo "   - Custom PROOF configuration file:  $pconfig"
         ;;
      esac
      echo ""
      echo -n "Is this correct [Y/n]: "
      read answ &> /dev/null
      if [ "x$answ" = "x" ] || [ "$answ" = "Y" ] || [ "$answ" = "y" ]; then
         REPEAT_LOOP=0
      fi
   fi
done

# End of info gathering


# Store settings in $saverc
echo "nodes=\"$nodes\"" > $saverc
echo "domain=\"-\""     >> $saverc
echo "base=\"$base\""   >> $saverc
echo "lower=\"$lower\"" >> $saverc
echo "upper=\"$upper\"" >> $saverc
echo "width=\"$width\"" >> $saverc
echo "rdomain=\"$rdomain\"" >> $saverc
echo "rootbinary=\"$rootbinary\"" >> $saverc
echo "remoterootbinary=\"$remoterootbinary\"" >> $saverc
echo "autoremoterootbinary=\"$autoremoterootbinary\"" >> $saverc
echo "puser=\"$puser\"" >> $saverc
echo "pmaster=\"$pmaster\"" >> $saverc
echo "pport=\"$pport\"" >> $saverc
echo "pdatadir=\"$pdatadir\"" >> $saverc
echo "pworkers=\"$pworkers\"" >> $saverc
echo "pconfig=\"$pconfig\"" >> $saverc

clear
echo ""
echo "Starting installation of PROOF on the cluster."
echo ""
echo "You might be asked to type your ssh password for user \"$puser\""
echo "several times. We try to use ssh-agent to store your ssh keys,"
echo "but if it does not work, you will have to type your password"
echo "up to three times per node. If you decide to quit and setup"
echo "automatic ssh access, re-running of the installer is simple as"
echo "all parameters have been saved. You can either run the installer"
echo "in batch (./proofinstall.sh -b) or you will be prompted with"
echo "your current choices."
echo ""
echo "Press <ENTER> to continue installation."
read

# end of "$batch" = "no"
fi

# Check if we can login to the nodes without passwd
# If we cannot copy public key to node


# Execute setup:
if [ -r "$rootbinary" ]; then
   echo $circle scp -u $puser -n \"$nodes\" "$rootbinary"
fi
echo $circle scp -u $puser -n \"$nodes\" remoteproofsetup.sh
echo $circle ssh -u $puser -n \"$nodes\" "remoteproofsetup.sh $pmaster $pport $autoremoterootbinary"

exit 0
