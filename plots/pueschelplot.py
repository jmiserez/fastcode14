#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (10,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200)
yvals2 = (
  0.442
  ,0.446
  ,0.391
  ,0.394
  ,0.365
  ,0.380
  ,0.369
  ,0.360
  ,0.361
  ,0.361
  ,0.352
  ,0.358
  ,0.352
  ,0.347
  ,0.350
  ,0.348
  ,0.327
)
yvals3 = (
  0.590
  ,0.649
  ,0.605
  ,0.645
  ,0.623
  ,0.617
  ,0.611
  ,0.641
  ,0.643
  ,0.643
  ,0.632
  ,0.626
  ,0.643
  ,0.593
  ,0.630
  ,0.625
  ,0.624
)
yvals4 = (
  0.465
  ,0.387
  ,0.461
  ,0.387
  ,0.363
  ,0.332
  ,0.333
  ,0.329
  ,0.321
  ,0.319
  ,0.313
  ,0.309
  ,0.301
  ,0.298
  ,0.295
  ,0.285
  ,0.284
)
yvals5 = (
  0.482
  ,0.504
  ,0.527
  ,0.589
  ,0.575
  ,0.557
  ,0.559
  ,0.582
  ,0.583
  ,0.582
  ,0.576
  ,0.588
  ,0.576
  ,0.580
  ,0.588
  ,0.572
  ,0.567
)
yvals6 = (
  0.570
  ,0.580
  ,0.539
  ,0.440
  ,0.420
  ,0.403
  ,0.383
  ,0.380
  ,0.364
  ,0.358
  ,0.352
  ,0.343
  ,0.341
  ,0.335
  ,0.335
  ,0.331
  ,0.325
)

yvals7 = (
  0.778
  ,0.850
  ,0.834
  ,0.845
  ,0.836
  ,0.842
  ,0.841
  ,0.813
  ,0.828
  ,0.814
  ,0.793
  ,0.797
  ,0.805
  ,0.799
  ,0.759
  ,0.795
  ,0.793
)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#BBBBBB',alpha=0.1)
plt.grid(color='white', alpha=0.5, linewidth=2, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Performance [flops/cycle]', **yprops)
plt.xlabel('Image dim [pixels])')

plt.gca().tick_params(direction='out', length=0, color='k')

plt.plot(xvals, yvals2, 'ro-', linewidth=2)
plt.plot(xvals, yvals3, 'bo-', linewidth=2)
plt.plot(xvals, yvals4, 'go-', linewidth=2)
plt.plot(xvals, yvals5, 'bo-', linewidth=2)
plt.plot(xvals, yvals7, 'bo-', linewidth=2)
plt.gca().set_axisbelow(True)

plt.ylim([0, 1]) 

plt.savefig('test.png', dpi=300)

plt.show()

