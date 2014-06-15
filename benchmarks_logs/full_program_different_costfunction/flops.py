#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
27772624, 109929964, 437424904, 1745148964, 6971546944, 27868155484, )
#consolidated.avx.run.log
yvals2 = (
18686060, 74682028, 298629932, 1194349612, 4777084460, 19107736108, )
#consolidated.blocking_avx.run.log
yvals3 = (
19494844, 78114748, 312729532, 1251459004, 5006917564, 20029833148, )
#consolidated.blocking.run.log
yvals4 = (
19822524, 79425468, 317972412, 1272430524, 5090803644, 20365377468, )
#consolidated.inline2.run.log
yvals5 = (
19322556, 77377212, 309681852, 1239072444, 4956978876, 19829292732, )
#consolidated.inline2x2.run.log
yvals6 = (
18667484, 74643100, 298525724, 1194016540, 4775900444, 19103277340, )
#consolidated.inline2x4.run.log
yvals7 = (
18712412, 74726940, 298687388, 1194333852, 4776529052, 19104528540, )
#consolidated.naive.run.log
yvals8 = (
18315196, 73396156, 293855164, 1175961532, 4704927676, 18821873596, )
#consolidated.onestep_avx.run.log
yvals9 = (
19494844, 78114748, 312729532, 1251459004, 5006917564, 20029833148, )
#consolidated.onestep.run.log
yvals10 = (
19625916, 78639036, 314826684, 1259847612, 5040471996, 20164050876, )
#consolidated.store_grey.run.log
yvals11 = (
18249660, 73134012, 292806588, 1171767228, 4688150460, 18754764732, )

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

