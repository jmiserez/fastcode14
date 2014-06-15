#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.088, 0.120, 0.179, 0.213, 0.210, 0.207, )
#consolidated.avx.run.log
yvals2 = (
0.290, 0.508, 0.343, 0.337, 0.309, 0.344, )
#consolidated.blocking_avx.run.log
yvals3 = (
0.262, 0.503, 0.372, 0.363, 0.332, 0.305, )
#consolidated.blocking.run.log
yvals4 = (
0.493, 0.481, 0.348, 0.430, 0.362, 0.360, )
#consolidated.inline2.run.log
yvals5 = (
0.257, 0.461, 0.372, 0.358, 0.349, 0.355, )
#consolidated.inline2x2.run.log
yvals6 = (
0.309, 0.489, 0.324, 0.337, 0.320, 0.315, )
#consolidated.inline2x4.run.log
yvals7 = (
0.207, 0.563, 0.324, 0.347, 0.359, 0.307, )
#consolidated.naive.run.log
yvals8 = (
0.106, 0.231, 0.239, 0.282, 0.270, 0.260, )
#consolidated.onestep_avx.run.log
yvals9 = (
0.166, 0.472, 0.585, 0.690, 0.402, 0.322, )
#consolidated.onestep.run.log
yvals10 = (
0.177, 0.523, 0.630, 0.664, 0.408, 0.347, )
#consolidated.store_grey.run.log
yvals11 = (
0.204, 0.320, 0.291, 0.309, 0.341, 0.251, )

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

