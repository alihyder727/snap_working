#! /usr/bin/env python3
from netCDF4 import Dataset
from numpy import sum, sqrt
import sys

case1 = sys.argv[1]
case2 = sys.argv[2]

data1 = Dataset(case1, 'r')
theta1 = data1['theta'][-1,0,:,:]

data2 = Dataset(case2, 'r')
theta2 = data2['theta'][-1,0,:,:]

diff = sqrt(sum((theta2 - theta1)*(theta2 - theta1)))

if diff < 1.E-3:
  print('### 1D HJ RE test passed. ###')
else:
  raise ValueError('ERROR: 1D HJ RE test failed. L2-norm is %.2g' % diff)
