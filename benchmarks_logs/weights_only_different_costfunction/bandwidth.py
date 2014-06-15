#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, )
#consolidated.avx.run.log
yvals2 = (
0.924, 0.862, 0.814, 0.825, 0.831, 0.823, )
#consolidated.blocking_avx.run.log
yvals3 = (
1.471, 1.453, 1.449, 1.457, 1.350, 1.430, )
#consolidated.blocking.run.log
yvals4 = (
0.471, 0.468, 0.461, 0.466, 0.468, 0.465, )
#consolidated.inline2.run.log
yvals5 = (
0.333, 0.329, 0.327, 0.324, 0.324, 0.323, )
#consolidated.inline2x2.run.log
yvals6 = (
0.222, 0.214, 0.206, 0.208, 0.207, 0.204, )
#consolidated.inline2x4.run.log
yvals7 = (
0.215, 0.203, 0.201, 0.200, 0.200, 0.198, )
#consolidated.naive.run.log
yvals8 = (
0.438, 0.398, 0.423, 0.424, 0.421, 0.425, )
#consolidated.onestep_avx.run.log
yvals9 = (
1.409, 1.397, 1.394, 1.389, 1.381, 1.367, )
#consolidated.onestep.run.log
yvals10 = (
0.442, 0.442, 0.434, 0.435, 0.436, 0.435, )
#consolidated.store_grey.run.log
yvals11 = (
0.363, 0.361, 0.354, 0.361, 0.362, 0.361, )

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

#plt.ylim([0, 1.0]) 

plt.savefig('bandwidth.png', dpi=300)

#plt.show()

