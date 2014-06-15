#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#consolidated.avx.run.log
yvals2 = (
2.287, 2.424, 2.501, 2.541, 2.562, 2.572, )
#consolidated.blocking_avx.run.log
yvals3 = (
1.115, 1.115, 1.115, 1.115, 1.115, 1.115, )
#consolidated.blocking.run.log
yvals4 = (
1.425, 1.425, 1.425, 1.425, 1.425, 1.425, )
#consolidated.inline2.run.log
yvals5 = (
1.702, 1.710, 1.714, 1.717, 1.718, 1.718, )
#consolidated.inline2x2.run.log
yvals6 = (
2.476, 2.543, 2.578, 2.596, 2.605, 2.610, )
#consolidated.inline2x4.run.log
yvals7 = (
2.435, 2.509, 2.549, 2.569, 2.580, 2.585, )
#consolidated.naive.run.log
yvals8 = (
1.281, 1.281, 1.281, 1.281, 1.281, 1.281, )
#consolidated.onestep_avx.run.log
yvals9 = (
1.161, 1.161, 1.161, 1.161, 1.161, 1.161, )
#consolidated.onestep.run.log
yvals10 = (
1.480, 1.480, 1.480, 1.480, 1.480, 1.480, )
#consolidated.store_grey.run.log
yvals11 = (
1.700, 1.700, 1.700, 1.700, 1.700, 1.700, )

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Operational intensity [flops/byte]',**yprops)
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

plt.ylim([0, 2.7]) 

plt.savefig('oi.png', dpi=300)

#plt.show()

