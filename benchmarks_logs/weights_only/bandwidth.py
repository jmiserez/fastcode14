#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, )
#consolidated.avx.run.log
yvals2 = (
0.717, 0.874, 0.874, 0.877, 0.857, 0.851, )
#consolidated.blocking_avx.run.log
yvals3 = (
1.477, 1.455, 1.452, 1.453, 1.450, 1.455, )
#consolidated.blocking.run.log
yvals4 = (
0.466, 0.470, 0.469, 0.468, 0.470, 0.465, )
#consolidated.inline2.run.log
yvals5 = (
0.333, 0.330, 0.328, 0.328, 0.328, 0.323, )
#consolidated.inline2x2.run.log
yvals6 = (
0.222, 0.211, 0.210, 0.209, 0.208, 0.203, )
#consolidated.inline2x4.run.log
yvals7 = (
0.184, 0.182, 0.200, 0.191, 0.188, 0.192, )
#consolidated.naive.run.log
yvals8 = (
0.438, 0.435, 0.429, 0.428, 0.430, 0.432, )
#consolidated.onestep_avx.run.log
yvals9 = (
1.414, 1.394, 1.399, 1.390, 1.403, 1.401, )
#consolidated.onestep.run.log
yvals10 = (
0.447, 0.445, 0.445, 0.445, 0.445, 0.445, )
#consolidated.store_grey.run.log
yvals11 = (
0.362, 0.359, 0.363, 0.361, 0.363, 0.363, )

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

