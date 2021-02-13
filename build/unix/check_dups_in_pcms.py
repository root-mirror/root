#!/usr/bin/env python

import argparse
import glob
import os
import subprocess

try:
    ROOTSYS = os.environ["ROOTSYS"]
except KeyError:
    print("You need to set ROOTSYS to find rootcling")

parser = argparse.ArgumentParser(
    description="This script checks if the same headers are contained in multiple module files(pcms).\nDuplicates make sense only if they are textual, non-modular headers such as assert.h\nThe script is basic and reports the first two modules containing the duplicate headers."
)
parser.add_argument("--pcmfiles", nargs="+", help="List of PCM files", required=True)
args = parser.parse_args()
infiles = args.pcmfiles
pcmfiles = []

for globname in infiles:
    for pcmfile in glob.glob(globname):
        pcmfiles.append(pcmfile)

ROOTCLING_BINARY = os.path.join(ROOTSYS, "bin", "rootcling")

headerdict = {}

for pcmfile in pcmfiles:
    if not pcmfile.endswith(".pcm") or pcmfile.endswith("_rdict.pcm"):
        print("Ignoring ROOT pcm file ", pcmfile)
        continue
    rootcling_output = (
        subprocess.Popen(
            [ROOTCLING_BINARY, "bare-cling", "-module-file-info", pcmfile],
            stdout=subprocess.PIPE,
        )
        .communicate()[0]
        .decode("utf-8")
    )
    HEADERS_IN_PCM = []
    for line in rootcling_output.split("\n"):
        if "Input file:" in line:
            line = line[line.find("Input file:") + len("Input file:") :]
            if "[" and "]" in line:
                line = line[: line.rfind("[")]
            line = line.strip()
            HEADERS_IN_PCM.append(line)
    for header in HEADERS_IN_PCM:
        if header not in headerdict.keys():
            headerdict[header] = []
        headerdict[header].append(pcmfile)


headerpath = []
def shorten_path(filename):
    dirname = os.path.dirname(filename)
    for i in range(0, 3):
        dirname = os.path.dirname(dirname)
    # Do not shorten if the path has less than 3 components
    if dirname == '/' or dirname == '':
        return filename

    if dirname not in headerpath:
        headerpath.append(dirname)
    return "[{0}]{1}".format(headerpath.index(dirname), filename.replace(dirname, ''))

for header in sorted(headerdict, key=lambda header: len(headerdict[header]), reverse=True):
    pcmfiles = headerdict[header]
    if len(pcmfiles) > 1:
        short_header = shorten_path(header)
        short_pcmfiles = ','.join(shorten_path(x) for x in pcmfiles)
        print("Header {0} duplicated in {1} modules - {2}"
              .format(short_header, len(pcmfiles), short_pcmfiles))

print("Legend:")
for i in range(len(headerpath)):
    print("[{0}] - {1}".format(i, headerpath[i]))
