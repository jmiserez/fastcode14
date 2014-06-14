#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)



font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Memory bandwidth [accesses/cycle]',**yprops)
plt.xlabel('Total pixels per input image [pixels]')

plt.gca().tick_params(direction='in', length=3, color='k')

plt.plot(xvals, yvals2, 'go-', linewidth=2) #avx
plt.plot(xvals, yvals3, 'bo-', linewidth=2) #blocking
plt.plot(xvals, yvals4, 'bo-', linewidth=2) #inline2
plt.plot(xvals, yvals5, 'bo-', linewidth=2) #inline2x2
plt.plot(xvals, yvals6, 'bo-', linewidth=2) #inline2x4
plt.plot(xvals, yvals7, 'ko-.', linewidth=2) #baseline_opts
plt.plot(xvals, yvals8, 'ko-', linewidth=2) #baseline
plt.plot(xvals, yvals9, 'mo-', linewidth=2) #onestep
plt.plot(xvals, yvals10, 'co-', linewidth=2) #storegrey
plt.gca().set_axisbelow(True)

#plt.ylim([0, 1.0]) 

plt.savefig('memory.png', dpi=300)

#plt.show()

