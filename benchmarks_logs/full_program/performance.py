#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.675, 0.660, 0.651, 0.646, 0.648, 0.645, )
#consolidated.avx.run.log
yvals2 = (
2.139, 2.169, 2.195, 2.245, 2.263, 2.165, )
#consolidated.blocking_avx.run.log
yvals3 = (
1.627, 1.622, 1.626, 1.629, 1.628, 1.617, )
#consolidated.blocking.run.log
yvals4 = (
0.662, 0.666, 0.667, 0.662, 0.666, 0.665, )
#consolidated.inline2.run.log
yvals5 = (
0.564, 0.563, 0.562, 0.564, 0.564, 0.562, )
#consolidated.inline2x2.run.log
yvals6 = (
0.547, 0.544, 0.542, 0.544, 0.544, 0.543, )
#consolidated.inline2x4.run.log
yvals7 = (
0.519, 0.510, 0.513, 0.514, 0.513, 0.512, )
#consolidated.naive.run.log
yvals8 = (
0.534, 0.537, 0.540, 0.538, 0.541, 0.540, )
#consolidated.onestep_avx.run.log
yvals9 = (
1.629, 1.620, 1.614, 1.609, 1.628, 1.628, )
#consolidated.onestep.run.log
yvals10 = (
0.658, 0.658, 0.654, 0.654, 0.643, 0.644, )
#consolidated.store_grey.run.log
yvals11 = (
0.620, 0.620, 0.618, 0.620, 0.611, 0.613, )

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

