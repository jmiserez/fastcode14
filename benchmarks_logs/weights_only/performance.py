#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.286, 0.285, 0.281, 0.279, 0.273, 0.276, )
#consolidated.avx.run.log
yvals2 = (
1.055, 1.362, 1.403, 1.430, 1.410, 1.404, )
#consolidated.blocking_avx.run.log
yvals3 = (
1.093, 1.077, 1.074, 1.075, 1.073, 1.076, )
#consolidated.blocking.run.log
yvals4 = (
0.446, 0.449, 0.449, 0.447, 0.449, 0.445, )
#consolidated.inline2.run.log
yvals5 = (
0.374, 0.372, 0.370, 0.371, 0.371, 0.366, )
#consolidated.inline2x2.run.log
yvals6 = (
0.354, 0.344, 0.347, 0.348, 0.348, 0.340, )
#consolidated.inline2x4.run.log
yvals7 = (
0.288, 0.294, 0.327, 0.316, 0.311, 0.319, )
#consolidated.naive.run.log
yvals8 = (
0.356, 0.354, 0.348, 0.348, 0.350, 0.351, )
#consolidated.onestep_avx.run.log
yvals9 = (
1.090, 1.075, 1.079, 1.071, 1.081, 1.080, )
#consolidated.onestep.run.log
yvals10 = (
0.441, 0.440, 0.439, 0.439, 0.439, 0.439, )
#consolidated.store_grey.run.log
yvals11 = (
0.389, 0.385, 0.390, 0.388, 0.390, 0.390, )

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

plt.ylim([0, 1.5]) 

plt.savefig('perf.png', dpi=300)

#plt.show()

