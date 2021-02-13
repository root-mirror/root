# Execute this script to set up the ROOT environment variables in Powershell.
#
# Author: Bertrand Bellenot, 17/11/2020

$scriptPath = split-path -parent $MyInvocation.MyCommand.Definition
$ROOTSYS = split-path -parent (get-item $scriptPath)
$env:PATH = $ROOTSYS + '\bin;' + $env:PATH
$env:CMAKE_PREFIX_PATH = $ROOTSYS + ';' + $env:CMAKE_PREFIX_PATH
$env:PYTHONPATH = $ROOTSYS + '\bin;' + $env:PYTHONPATH
$env:CLING_STANDARD_PCH = "none"
