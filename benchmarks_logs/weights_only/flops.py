#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
9682944, 38731776, 154927104, 619708416, 2478833664, 9915334656, )
#consolidated.avx.run.log
yvals2 = (
8886944, 35360992, 141086048, 563644000, 2253191264, 9010011232, )
#consolidated.blocking_avx.run.log
yvals3 = (
9699328, 38797312, 155189248, 620756992, 2483027968, 9932111872, )
#consolidated.blocking.run.log
yvals4 = (
10027008, 40108032, 160432128, 641728512, 2566914048, 10267656192, )
#consolidated.inline2.run.log
yvals5 = (
9527040, 38059776, 152141568, 608370432, 2433089280, 9731571456, )
#consolidated.inline2x2.run.log
yvals6 = (
8870768, 35324464, 140984240, 563313328, 2252009648, 9005554864, )
#consolidated.inline2x4.run.log
yvals7 = (
8913296, 35405904, 141143504, 563628240, 2252635856, 9006803664, )
#consolidated.naive.run.log
yvals8 = (
8519680, 34078720, 136314880, 545259520, 2181038080, 8724152320, )
#consolidated.onestep_avx.run.log
yvals9 = (
9699328, 38797312, 155189248, 620756992, 2483027968, 9932111872, )
#consolidated.onestep.run.log
yvals10 = (
9830400, 39321600, 157286400, 629145600, 2516582400, 10066329600, )
#consolidated.store_grey.run.log
yvals11 = (
8454144, 33816576, 135266304, 541065216, 2164260864, 8657043456, )

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

