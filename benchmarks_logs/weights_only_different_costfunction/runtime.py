#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128*128,256*256,512*512,1024*1024,2048*2048,4096*4096)

#01ref_matlab.run.log
yvals1 = (
33759731, 135570955, 552421919, 2269602461, 8885906130, 36015255316, )
#consolidated.avx.run.log
yvals2 = (
6530618, 26320746, 107908518, 419057938, 1648520783, 6632197319, )
#consolidated.blocking_avx.run.log
yvals3 = (
8912664, 36080342, 144703595, 575726267, 2486355131, 9388443278, )
#consolidated.blocking.run.log
yvals4 = (
22248602, 89602178, 363669177, 1439841159, 5739106180, 23081873546, )
#consolidated.inline2.run.log
yvals5 = (
25507247, 102644677, 411318480, 1657162263, 6629865282, 26606832934, )
#consolidated.inline2x2.run.log
yvals6 = (
25037628, 100949805, 412761160, 1625890323, 6496649209, 26344289782, )
#consolidated.inline2x4.run.log
yvals7 = (
26405665, 108007978, 428229955, 1709296069, 6816604481, 27383799330, )
#consolidated.naive.run.log
yvals8 = (
23938344, 105368702, 396348270, 1584037796, 6370643072, 25252242582, )
#consolidated.onestep_avx.run.log
yvals9 = (
8933014, 36032400, 144387730, 579849117, 2332319782, 9425765121, )
#consolidated.onestep.run.log
yvals10 = (
22518933, 90155507, 367458222, 1464133233, 5848577596, 23431624460, )
#consolidated.store_grey.run.log
yvals11 = (
21661717, 87248261, 355553508, 1393329455, 5565137182, 22331495194, )

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

