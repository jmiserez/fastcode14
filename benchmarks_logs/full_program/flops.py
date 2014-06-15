#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
22857382, 90269122, 358781662, 1430576122, 5713255702, 22834990642, )
#consolidated.avx.run.log
yvals2 = (
13767260, 55017628, 219983132, 879773212, 3518789660, 14074567708, )
#consolidated.blocking_avx.run.log
yvals3 = (
14579644, 58453948, 234086332, 936886204, 3748626364, 14996668348, )
#consolidated.blocking.run.log
yvals4 = (
14907324, 59764668, 239329212, 957857724, 3832512444, 15332212668, )
#consolidated.inline2.run.log
yvals5 = (
14407356, 57716412, 231038652, 924499644, 3698687676, 14796127932, )
#consolidated.inline2x2.run.log
yvals6 = (
13751084, 54981100, 219881324, 879442540, 3517608044, 14070111340, )
#consolidated.inline2x4.run.log
yvals7 = (
13793612, 55062540, 220040588, 879757452, 3518234252, 14071360140, )
#consolidated.naive.run.log
yvals8 = (
13399996, 53735356, 215211964, 861388732, 3446636476, 13788708796, )
#consolidated.onestep_avx.run.log
yvals9 = (
14579644, 58453948, 234086332, 936886204, 3748626364, 14996668348, )
#consolidated.onestep.run.log
yvals10 = (
14710716, 58978236, 236183484, 945274812, 3782180796, 15130886076, )
#consolidated.store_grey.run.log
yvals11 = (
13334460, 53473212, 214163388, 857194428, 3429859260, 13721599932, )

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9}

plt.rc('font', **font)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=1, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Total flops',**yprops)
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

plt.savefig('flops.png', dpi=300)

#plt.show()

