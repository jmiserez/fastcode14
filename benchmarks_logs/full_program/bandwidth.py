#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, )
#consolidated.avx.run.log
yvals2 = (
7.431, 7.507, 7.584, 7.748, 7.805, 7.465, )
#consolidated.blocking_avx.run.log
yvals3 = (
6.127, 6.107, 6.123, 6.135, 6.131, 6.091, )
#consolidated.blocking.run.log
yvals4 = (
2.323, 2.338, 2.339, 2.322, 2.338, 2.334, )
#consolidated.inline2.run.log
yvals5 = (
1.969, 1.966, 1.961, 1.970, 1.971, 1.964, )
#consolidated.inline2x2.run.log
yvals6 = (
1.883, 1.874, 1.866, 1.875, 1.874, 1.871, )
#consolidated.inline2x4.run.log
yvals7 = (
1.787, 1.756, 1.767, 1.770, 1.767, 1.766, )
#consolidated.naive.run.log
yvals8 = (
2.081, 2.095, 2.108, 2.098, 2.113, 2.106, )
#consolidated.onestep_avx.run.log
yvals9 = (
6.075, 6.045, 6.023, 6.005, 6.076, 6.074, )
#consolidated.onestep.run.log
yvals10 = (
2.314, 2.316, 2.304, 2.302, 2.264, 2.268, )
#consolidated.store_grey.run.log
yvals11 = (
2.308, 2.309, 2.303, 2.308, 2.275, 2.285, )

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

