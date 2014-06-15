#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
14598144, 58392576, 233570304, 934281216, 3737124864, 14948499456, )
#consolidated.avx.run.log
yvals2 = (
13805744, 55025392, 219732848, 878220400, 3511486064, 14043179632, )
#consolidated.blocking_avx.run.log
yvals3 = (
14614528, 58458112, 233832448, 935329792, 3741319168, 14965276672, )
#consolidated.blocking.run.log
yvals4 = (
14942208, 59768832, 239075328, 956301312, 3825205248, 15300820992, )
#consolidated.inline2.run.log
yvals5 = (
14442240, 57720576, 230784768, 922943232, 3691380480, 14764736256, )
#consolidated.inline2x2.run.log
yvals6 = (
13787168, 54986464, 219628640, 877887328, 3510302048, 14038720864, )
#consolidated.inline2x4.run.log
yvals7 = (
13832096, 55070304, 219790304, 878204640, 3510930656, 14039972064, )
#consolidated.naive.run.log
yvals8 = (
13434880, 53739520, 214958080, 859832320, 3439329280, 13757317120, )
#consolidated.onestep_avx.run.log
yvals9 = (
14614528, 58458112, 233832448, 935329792, 3741319168, 14965276672, )
#consolidated.onestep.run.log
yvals10 = (
14745600, 58982400, 235929600, 943718400, 3774873600, 15099494400, )
#consolidated.store_grey.run.log
yvals11 = (
13369344, 53477376, 213909504, 855638016, 3422552064, 13690208256, )

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

