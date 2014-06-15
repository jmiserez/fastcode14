#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.432, 0.431, 0.423, 0.412, 0.421, 0.415, )
#consolidated.avx.run.log
yvals2 = (
2.114, 2.091, 2.036, 2.096, 2.130, 2.117, )
#consolidated.blocking_avx.run.log
yvals3 = (
1.640, 1.620, 1.616, 1.625, 1.505, 1.594, )
#consolidated.blocking.run.log
yvals4 = (
0.672, 0.667, 0.657, 0.664, 0.667, 0.663, )
#consolidated.inline2.run.log
yvals5 = (
0.566, 0.562, 0.561, 0.557, 0.557, 0.555, )
#consolidated.inline2x2.run.log
yvals6 = (
0.551, 0.545, 0.532, 0.540, 0.540, 0.533, )
#consolidated.inline2x4.run.log
yvals7 = (
0.524, 0.510, 0.513, 0.514, 0.515, 0.513, )
#consolidated.naive.run.log
yvals8 = (
0.561, 0.510, 0.542, 0.543, 0.540, 0.545, )
#consolidated.onestep_avx.run.log
yvals9 = (
1.636, 1.622, 1.619, 1.613, 1.604, 1.588, )
#consolidated.onestep.run.log
yvals10 = (
0.655, 0.654, 0.642, 0.645, 0.645, 0.644, )
#consolidated.store_grey.run.log
yvals11 = (
0.617, 0.613, 0.602, 0.614, 0.615, 0.613, )

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
plt.plot(xvals, yvals5, 'co-', linewidth=2) #inline2
plt.plot(xvals, yvals6, 'co-', linewidth=2) #inline2x2
plt.plot(xvals, yvals7, 'bo-', linewidth=2) #inline2x4
plt.plot(xvals, yvals8, 'ko-', linewidth=2) #baseline
plt.plot(xvals, yvals9, 'mo-.', linewidth=2) #onestep_AVX
plt.plot(xvals, yvals10, 'mo-', linewidth=2) #onestep
plt.plot(xvals, yvals11, 'yo-', linewidth=2) #storegrey
plt.gca().set_axisbelow(True)

plt.ylim([0, 2.5]) 

plt.savefig('perf.png', dpi=300)

#plt.show()

