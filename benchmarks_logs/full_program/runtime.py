#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
33848822, 136822743, 551546688, 2213785661, 8810274154, 35383046474, )
#consolidated.avx.run.log
yvals2 = (
6435609, 25366088, 100202676, 391868363, 1555022309, 6501232400, )
#consolidated.blocking_avx.run.log
yvals3 = (
8958754, 36047121, 144002518, 575257375, 2303298567, 9276036310, )
#consolidated.blocking.run.log
yvals4 = (
22503741, 89694858, 359064645, 1447958253, 5752915795, 23052012078, )
#consolidated.inline2.run.log
yvals5 = (
25527694, 102465534, 411425169, 1638900626, 6552290938, 26310966746, )
#consolidated.inline2x2.run.log
yvals6 = (
25140915, 101048945, 405746277, 1615480467, 6464800579, 25899774482, )
#consolidated.inline2x4.run.log
yvals7 = (
26566941, 108015923, 429103543, 1712830432, 6863926069, 27464591391, )
#consolidated.naive.run.log
yvals8 = (
25111981, 100081698, 398301907, 1602410202, 6366153472, 25547769749, )
#consolidated.onestep_avx.run.log
yvals9 = (
8949083, 36073119, 145009488, 582111524, 2302113502, 9213706111, )
#consolidated.onestep.run.log
yvals10 = (
22358686, 89607588, 360870985, 1445963372, 5882436994, 23492170724, )
#consolidated.store_grey.run.log
yvals11 = (
21513512, 86250405, 346387581, 1383600849, 5616958266, 22372762719, )

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

