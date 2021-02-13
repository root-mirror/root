from __future__ import print_function

import os, sys, subprocess

target = sys.argv[1]
output = sys.argv[2]

def isokay(name):
 # filter standard symbols
    return name[0] != '_' and not name in {'memcpy', 'memmove', 'memset'}

popen = subprocess.Popen(['dumpbin', '/SYMBOLS', target+'.obj'],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

stdout, _ = popen.communicate()
stdout = stdout.decode('utf-8').strip()

outf = open(output+'.def', 'w')
outf.write('LIBRARY    %s.dll\nEXPORTS\n' % output)
for line in stdout.split('\r\n'):
    parts = line.split()
    if len(parts) < 8:
        continue
    if parts[7][0:4] in ['??_G', '??_E']:   # do not export deleting destructors
        continue
    if parts[4] == 'External':
        if isokay(parts[6]):
            outf.write('\t%s\tDATA\n' % parts[6])
    elif parts[4] == '()' and parts[5] == 'External':
        if isokay(parts[7]):
            outf.write('\t%s\n' % parts[7])

