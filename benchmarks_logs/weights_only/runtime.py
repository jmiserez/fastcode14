#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
33856540, 135781067, 551197888, 2221081868, 9082036826, 35956438065, )
#consolidated.avx.run.log
yvals2 = (
8421365, 25965728, 100536076, 394039070, 1598532038, 6415592819, )
#consolidated.blocking_avx.run.log
yvals3 = (
8872400, 36032107, 144435257, 577355968, 2314257838, 9226750507, )
#consolidated.blocking.run.log
yvals4 = (
22485689, 89303462, 357387950, 1435286264, 5715110983, 23097738121, )
#consolidated.inline2.run.log
yvals5 = (
25470575, 102404141, 410682732, 1638945074, 6553709779, 26586323029, )
#consolidated.inline2x2.run.log
yvals6 = (
25057666, 102673091, 406092494, 1617328561, 6468571068, 26470585879, )
#consolidated.inline2x4.run.log
yvals7 = (
30900928, 120495874, 431120791, 1785297991, 7238006034, 28221247684, )
#consolidated.naive.run.log
yvals8 = (
23924944, 96369115, 391316966, 1568430862, 6235772890, 24866060909, )
#consolidated.onestep_avx.run.log
yvals9 = (
8899407, 36103329, 143881218, 579479901, 2295926011, 9193878710, )
#consolidated.onestep.run.log
yvals10 = (
22284919, 89449029, 358049186, 1431914057, 5732177700, 22908192445, )
#consolidated.store_grey.run.log
yvals11 = (
21709762, 87730940, 346807892, 1394673534, 5547849567, 22203026562, )

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

