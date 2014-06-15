#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.817, 0.810, 0.793, 0.778, 0.792, 0.784, )
#consolidated.avx.run.log
yvals2 = (
2.908, 2.917, 2.965, 3.003, 3.061, 3.033, )
#consolidated.blocking_avx.run.log
yvals3 = (
2.173, 1.865, 2.151, 2.166, 2.141, 2.171, )
#consolidated.blocking.run.log
yvals4 = (
0.872, 0.883, 0.885, 0.887, 0.887, 0.885, )
#consolidated.inline2.run.log
yvals5 = (
0.756, 0.753, 0.752, 0.753, 0.755, 0.754, )
#consolidated.inline2x2.run.log
yvals6 = (
0.740, 0.738, 0.736, 0.738, 0.735, 0.735, )
#consolidated.inline2x4.run.log
yvals7 = (
0.700, 0.697, 0.695, 0.693, 0.694, 0.687, )
#consolidated.naive.run.log
yvals8 = (
0.762, 0.760, 0.752, 0.755, 0.754, 0.751, )
#consolidated.onestep_avx.run.log
yvals9 = (
2.159, 2.166, 2.167, 2.166, 2.164, 2.171, )
#consolidated.onestep.run.log
yvals10 = (
0.878, 0.878, 0.880, 0.880, 0.878, 0.878, )
#consolidated.store_grey.run.log
yvals11 = (
0.834, 0.847, 0.846, 0.847, 0.847, 0.845, )

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

plt.ylim([0, 3.5]) 

plt.savefig('perf.png', dpi=300)

#plt.show()

