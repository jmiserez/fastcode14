#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.368, 0.354, 0.320, 0.278, 0.261, 0.249, )
#consolidated.avx.run.log
yvals2 = (
1.046, 0.947, 0.819, 0.656, 0.582, 0.509, )
#consolidated.blocking_avx.run.log
yvals3 = (
0.607, 0.632, 0.633, 0.505, 0.531, 0.491, )
#consolidated.blocking.run.log
yvals4 = (
0.510, 0.484, 0.436, 0.404, 0.378, 0.350, )
#consolidated.inline2.run.log
yvals5 = (
0.448, 0.410, 0.386, 0.363, 0.333, 0.304, )
#consolidated.inline2x2.run.log
yvals6 = (
0.426, 0.410, 0.374, 0.344, 0.320, 0.303, )
#consolidated.inline2x4.run.log
yvals7 = (
0.405, 0.386, 0.365, 0.323, 0.301, 0.289, )
#consolidated.naive.run.log
yvals8 = (
0.424, 0.418, 0.387, 0.346, 0.326, 0.305, )
#consolidated.onestep_avx.run.log
yvals9 = (
0.910, 0.713, 0.722, 0.625, 0.522, 0.486, )
#consolidated.onestep.run.log
yvals10 = (
0.442, 0.435, 0.427, 0.391, 0.328, 0.343, )
#consolidated.store_grey.run.log
yvals11 = (
0.474, 0.430, 0.404, 0.367, 0.346, 0.327, )

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Performance [flops/cycle]',**yprops)
plt.xlabel('Total pixels per input image [pixels]')

plt.gca().tick_params(direction='in', length=3, color='k')

plt.plot(xvals, yvals1, 'ro-', linewidth=2) #matlab
plt.plot(xvals, yvals2, 'bo-.', linewidth=2) #inline2x4_AVX
plt.plot(xvals, yvals3, 'co-.', linewidth=2) #blocking_AVX
plt.plot(xvals, yvals4, 'co-', linewidth=2) #blocking
#plt.plot(xvals, yvals5, 'co-', linewidth=2) #inline2
#plt.plot(xvals, yvals6, 'co-', linewidth=2) #inline2x2
#plt.plot(xvals, yvals7, 'bo-', linewidth=2) #inline2x4
plt.plot(xvals, yvals8, 'ko-', linewidth=2) #baseline
plt.plot(xvals, yvals9, 'mo-.', linewidth=2) #onestep_AVX
plt.plot(xvals, yvals10, 'mo-', linewidth=2) #onestep
#plt.plot(xvals, yvals11, 'yo-', linewidth=2) #storegrey
plt.gca().set_axisbelow(True)

plt.ylim([0, 1.2]) 

plt.savefig('perf.png', dpi=300)

#plt.show()

