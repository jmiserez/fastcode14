#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, )
#consolidated.avx.run.log
yvals2 = (
7.442, 7.437, 7.546, 7.634, 7.777, 7.704, )
#consolidated.blocking_avx.run.log
yvals3 = (
6.119, 5.255, 6.064, 6.109, 6.037, 6.124, )
#consolidated.blocking.run.log
yvals4 = (
2.300, 2.330, 2.339, 2.343, 2.345, 2.337, )
#consolidated.inline2.run.log
yvals5 = (
1.966, 1.960, 1.959, 1.962, 1.968, 1.965, )
#consolidated.inline2x2.run.log
yvals6 = (
1.877, 1.871, 1.867, 1.872, 1.865, 1.864, )
#consolidated.inline2x4.run.log
yvals7 = (
1.776, 1.770, 1.764, 1.761, 1.763, 1.745, )
#consolidated.naive.run.log
yvals8 = (
2.175, 2.170, 2.150, 2.157, 2.156, 2.146, )
#consolidated.onestep_avx.run.log
yvals9 = (
6.021, 6.045, 6.053, 6.050, 6.046, 6.064, )
#consolidated.onestep.run.log
yvals10 = (
2.315, 2.318, 2.323, 2.326, 2.320, 2.320, )
#consolidated.store_grey.run.log
yvals11 = (
2.270, 2.307, 2.305, 2.309, 2.309, 2.304, )

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

