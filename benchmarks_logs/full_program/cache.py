#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.111, 0.100, 0.174, 0.306, 0.389, 0.369, )
#consolidated.avx.run.log
yvals2 = (
0.138, 0.133, 0.163, 0.327, 0.570, 0.579, )
#consolidated.blocking_avx.run.log
yvals3 = (
0.252, 0.166, 0.183, 0.352, 0.567, 0.583, )
#consolidated.blocking.run.log
yvals4 = (
0.134, 0.139, 0.176, 0.332, 0.580, 0.582, )
#consolidated.inline2.run.log
yvals5 = (
0.137, 0.156, 0.190, 0.322, 0.576, 0.579, )
#consolidated.inline2x2.run.log
yvals6 = (
0.231, 0.127, 0.182, 0.325, 0.570, 0.578, )
#consolidated.inline2x4.run.log
yvals7 = (
0.184, 0.137, 0.169, 0.341, 0.578, 0.584, )
#consolidated.naive.run.log
yvals8 = (
0.216, 0.126, 0.167, 0.329, 0.568, 0.566, )
#consolidated.onestep_avx.run.log
yvals9 = (
0.164, 0.161, 0.169, 0.336, 0.586, 0.572, )
#consolidated.onestep.run.log
yvals10 = (
0.217, 0.148, 0.175, 0.351, 0.578, 0.568, )
#consolidated.store_grey.run.log
yvals11 = (
0.093, 0.159, 0.184, 0.333, 0.587, 0.581, )

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

