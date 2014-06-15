#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, )
#consolidated.avx.run.log
yvals2 = (
3.634, 3.277, 2.829, 2.265, 2.007, 1.757, )
#consolidated.blocking_avx.run.log
yvals3 = (
2.286, 2.380, 2.385, 1.903, 2.000, 1.849, )
#consolidated.blocking.run.log
yvals4 = (
1.789, 1.699, 1.530, 1.419, 1.326, 1.228, )
#consolidated.inline2.run.log
yvals5 = (
1.562, 1.431, 1.347, 1.266, 1.161, 1.060, )
#consolidated.inline2x2.run.log
yvals6 = (
1.466, 1.413, 1.287, 1.183, 1.101, 1.043, )
#consolidated.inline2x4.run.log
yvals7 = (
1.395, 1.329, 1.257, 1.115, 1.036, 0.998, )
#consolidated.naive.run.log
yvals8 = (
1.655, 1.630, 1.509, 1.352, 1.273, 1.192, )
#consolidated.onestep_avx.run.log
yvals9 = (
3.392, 2.661, 2.693, 2.331, 1.947, 1.814, )
#consolidated.onestep.run.log
yvals10 = (
1.554, 1.529, 1.502, 1.377, 1.154, 1.207, )
#consolidated.store_grey.run.log
yvals11 = (
1.763, 1.603, 1.507, 1.366, 1.288, 1.217, )

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Memory bandwidth [bytes/cycle]',**yprops)
plt.xlabel('Total pixels per input image [pixels]')

plt.gca().tick_params(direction='in', length=3, color='k')

#plt.plot(xvals, yvals1, 'ro-', linewidth=2) #matlab
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

plt.ylim([0, 4.0]) 

plt.savefig('bandwidth.png', dpi=300)

#plt.show()

