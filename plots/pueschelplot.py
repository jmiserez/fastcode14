#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (128,256,512,1024,2048,4096)

yvals1 = (
 0.262,
 0.241,
 0.219,
 0.190,
 0.177,
 0.165
)
#naive
yvals2 = (
 0.747,
 0.744,
 0.730,
 0.730,
 0.726,
 0.727

)
#store grey
yvals3 = (
 0.747,
 0.735,
 0.729,
 0.730,
 0.733,
 0.734

)
#onestep
yvals4 = (
 0.686,
 0.678,
 0.674,
 0.674,
 0.677,
 0.675
)
#blocking
yvals5 = (
 0.696,
 0.690,
 0.686,
 0.679,
 0.685,
 0.684
)
#inline2
yvals6 = (
 0.701,
 0.696,
 0.686,
 0.687,
 0.688,
 0.686
)
#inline2x2
yvals7 = (
 0.691,
 0.680,
 0.687,
 0.688,
 0.688,
 0.685
)
#inline2x4
yvals8 = (
 0.671,
 0.666,
 0.664,
 0.657,
 0.662,
 0.661
)
#avx
yvals9 = (
 2.220,
 2.147,
 2.295,
 2.383,
 2.395,
 2.411
)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#FFFFFF',alpha=0.1)
plt.grid(color='grey', alpha=0.5, linewidth=2, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Performance [flops/cycle]', **yprops)
plt.xlabel('Image edge dimension')

plt.gca().tick_params(direction='out', length=5, color='k')

plt.plot(xvals, yvals1, 'ro-', linewidth=2)
plt.plot(xvals, yvals2, 'ro-', linewidth=2)
plt.plot(xvals, yvals3, 'bo-', linewidth=2)
plt.plot(xvals, yvals4, 'bo-', linewidth=2)
plt.plot(xvals, yvals5, 'co-', linewidth=2)
plt.plot(xvals, yvals6, 'co-', linewidth=2)
plt.plot(xvals, yvals7, 'co-', linewidth=2)
plt.plot(xvals, yvals8, 'mo-', linewidth=2)
plt.plot(xvals, yvals9, 'go-', linewidth=2)
plt.gca().set_axisbelow(True)

plt.ylim([0, 1.0]) 

plt.savefig('test.png', dpi=300)

plt.show()

