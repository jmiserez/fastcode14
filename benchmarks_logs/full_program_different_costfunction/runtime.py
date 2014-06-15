#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
33995748, 135658003, 551770719, 2242189292, 8806033410, 35537104549, )
#consolidated.avx.run.log
yvals2 = (
6425606, 25604833, 100706255, 397713328, 1560615851, 6299990872, )
#consolidated.blocking_avx.run.log
yvals3 = (
8970440, 41895601, 145412222, 577744711, 2339044313, 9225255203, )
#consolidated.blocking.run.log
yvals4 = (
22722461, 89991086, 359095501, 1434764497, 5736518045, 23021665241, )
#consolidated.inline2.run.log
yvals5 = (
25566892, 102782294, 411846692, 1645138498, 6562646780, 26296095223, )
#consolidated.inline2x2.run.log
yvals6 = (
25227551, 101191623, 405696042, 1617817989, 6493660243, 25993853999, )
#consolidated.inline2x4.run.log
yvals7 = (
26723483, 107177982, 429865100, 1722356104, 6879083928, 27792822215, )
#consolidated.naive.run.log
yvals8 = (
24029447, 96624666, 390611621, 1558283676, 6238397293, 25071269437, )
#consolidated.onestep_avx.run.log
yvals9 = (
9029738, 36071506, 144292601, 577788148, 2313612407, 9227971678, )
#consolidated.onestep.run.log
yvals10 = (
22355574, 89532466, 357932665, 1431027166, 5739623601, 22959065709, )
#consolidated.store_grey.run.log
yvals11 = (
21869445, 86333737, 346119490, 1383367593, 5534847838, 22193508978, )

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Runtime [cycles]',**yprops)
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

#plt.ylim([0, 1.0]) 

plt.savefig('runtime.png', dpi=300)

#plt.show()

