#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.012, 0.091, 0.178, 0.215, 0.234, 0.228, )
#consolidated.avx.run.log
yvals2 = (
0.357, 0.516, 0.347, 0.349, 0.321, 0.322, )
#consolidated.blocking_avx.run.log
yvals3 = (
0.000, 0.414, 0.376, 0.364, 0.355, 0.324, )
#consolidated.blocking.run.log
yvals4 = (
0.701, 0.358, 0.341, 0.365, 0.318, 0.387, )
#consolidated.inline2.run.log
yvals5 = (
0.005, 0.326, 0.356, 0.363, 0.345, 0.395, )
#consolidated.inline2x2.run.log
yvals6 = (
0.007, 0.485, 0.323, 0.358, 0.327, 0.408, )
#consolidated.inline2x4.run.log
yvals7 = (
0.513, 0.504, 0.327, 0.398, 0.414, 0.361, )
#consolidated.naive.run.log
yvals8 = (
0.003, 0.177, 0.279, 0.280, 0.272, 0.254, )
#consolidated.onestep_avx.run.log
yvals9 = (
0.005, 0.528, 0.583, 0.692, 0.403, 0.299, )
#consolidated.onestep.run.log
yvals10 = (
0.007, 0.429, 0.683, 0.749, 0.447, 0.313, )
#consolidated.store_grey.run.log
yvals11 = (
0.334, 0.372, 0.299, 0.349, 0.276, 0.217, )

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Cache miss rate',**yprops)
plt.xlabel('Total pixels per input image [pixels]')

plt.gca().tick_params(direction='in', length=3, color='k')

plt.plot(xvals, yvals1, 'ro-', linewidth=2) #matlab
plt.plot(xvals, yvals2, 'bo-.', linewidth=2) #inline2x4_AVX
plt.plot(xvals, yvals3, 'co-.', linewidth=2) #blocking_AVX
plt.plot(xvals, yvals4, 'co-', linewidth=2) #blocking
plt.plot(xvals, yvals5, 'co-', linewidth=2) #inline2
plt.plot(xvals, yvals6, 'co-', linewidth=2) #inline2x2
plt.plot(xvals, yvals7, 'bo-', linewidth=2) #inline2x4
plt.plot(xvals, yvals8, 'ko-', linewidth=2) #baseline
plt.plot(xvals, yvals9, 'mo-.', linewidth=2) #onestep_AVX
plt.plot(xvals, yvals10, 'mo-', linewidth=2) #onestep
plt.plot(xvals, yvals11, 'yo-', linewidth=2) #storegrey
plt.gca().set_axisbelow(True)

plt.ylim([0, 1.0]) 

plt.savefig('cache.png', dpi=300)

#plt.show()

