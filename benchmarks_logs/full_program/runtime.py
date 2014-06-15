#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
62038001, 255115597, 1122107295, 5143230843, 21850998719, 91588292062, )
#consolidated.avx.run.log
yvals2 = (
13157724, 58098999, 268590549, 1340625202, 6047651046, 27628312273, )
#consolidated.blocking_avx.run.log
yvals3 = (
24015747, 92503055, 369759572, 1854609318, 7062230601, 30558118926, )
#consolidated.blocking.run.log
yvals4 = (
29209766, 123378905, 548850308, 2369392042, 10142615893, 43829615924, )
#consolidated.inline2.run.log
yvals5 = (
32183732, 140798152, 598779510, 2549008760, 11123656746, 48732802247, )
#consolidated.inline2x2.run.log
yvals6 = (
32298636, 133975194, 588580648, 2560190614, 11004438487, 46451627922, )
#consolidated.inline2x4.run.log
yvals7 = (
34037969, 142765743, 603429482, 2719642725, 11707251833, 48616303524, )
#consolidated.naive.run.log
yvals8 = (
31592165, 128592603, 556617595, 2486420133, 10568446281, 45153116373, )
#consolidated.onestep_avx.run.log
yvals9 = (
16027908, 81932653, 324359022, 1499659231, 7185820509, 30847309401, )
#consolidated.onestep.run.log
yvals10 = (
33301161, 135721520, 553380896, 2416338959, 11536076115, 44122231725, )
#consolidated.store_grey.run.log
yvals11 = (
28153787, 124283911, 529554131, 2337210457, 9920658154, 42018381016, )

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

