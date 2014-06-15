#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.138, 0.088, 0.172, 0.219, 0.205, 0.207, )
#consolidated.avx.run.log
yvals2 = (
0.278, 0.523, 0.356, 0.349, 0.309, 0.303, )
#consolidated.blocking_avx.run.log
yvals3 = (
0.290, 0.446, 0.380, 0.360, 0.340, 0.307, )
#consolidated.blocking.run.log
yvals4 = (
0.769, 0.521, 0.352, 0.349, 0.330, 0.340, )
#consolidated.inline2.run.log
yvals5 = (
0.369, 0.487, 0.377, 0.388, 0.362, 0.351, )
#consolidated.inline2x2.run.log
yvals6 = (
0.454, 0.489, 0.310, 0.327, 0.326, 0.323, )
#consolidated.inline2x4.run.log
yvals7 = (
0.649, 0.484, 0.336, 0.388, 0.383, 0.400, )
#consolidated.naive.run.log
yvals8 = (
0.105, 0.194, 0.269, 0.273, 0.264, 0.269, )
#consolidated.onestep_avx.run.log
yvals9 = (
0.428, 0.450, 0.582, 0.694, 0.404, 0.316, )
#consolidated.onestep.run.log
yvals10 = (
0.185, 0.526, 0.667, 0.746, 0.448, 0.318, )
#consolidated.store_grey.run.log
yvals11 = (
0.342, 0.332, 0.289, 0.311, 0.266, 0.223, )

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

