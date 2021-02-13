import array, sys

N = 100000000 # 10^8


#- group: stl-vector ---------------------------------------------------------
looprange = range
if sys.hexversion < 0x3000000:
    looprange = xrange
global_vector = array.array('i', looprange(N))
