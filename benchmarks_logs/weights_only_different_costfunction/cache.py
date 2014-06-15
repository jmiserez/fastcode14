#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.004, 0.090, 0.177, 0.216, 0.209, 0.207, )
#consolidated.avx.run.log
yvals2 = (
0.001, 0.322, 0.347, 0.355, 0.320, 0.305, )
#consolidated.blocking_avx.run.log
yvals3 = (
0.001, 0.352, 0.377, 0.360, 0.378, 0.312, )
#consolidated.blocking.run.log
yvals4 = (
0.051, 0.365, 0.421, 0.372, 0.337, 0.346, )
#consolidated.inline2.run.log
yvals5 = (
0.004, 0.333, 0.333, 0.356, 0.334, 0.324, )
#consolidated.inline2x2.run.log
yvals6 = (
0.004, 0.315, 0.324, 0.337, 0.310, 0.327, )
#consolidated.inline2x4.run.log
yvals7 = (
0.005, 0.534, 0.332, 0.353, 0.317, 0.315, )
#consolidated.naive.run.log
yvals8 = (
0.002, 0.276, 0.268, 0.279, 0.276, 0.268, )
#consolidated.onestep_avx.run.log
yvals9 = (
0.011, 0.348, 0.584, 0.640, 0.373, 0.289, )
#consolidated.onestep.run.log
yvals10 = (
0.138, 0.635, 0.576, 0.607, 0.413, 0.333, )
#consolidated.store_grey.run.log
yvals11 = (
0.333, 0.304, 0.323, 0.321, 0.273, 0.220, )

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

