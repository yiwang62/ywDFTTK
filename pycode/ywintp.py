#!/global/project/projectdirs/m891/yiwang62/anaconda3/bin/python
#!/usr/bin/python -x

import sys
import numpy as np
from scipy import interpolate

def ywinterp(x):
  data = np.loadtxt(sys.stdin, comments="#", dtype=np.float)
  data[np.isnan(data)] = 0.0
  return float(interpolate.splev(at, interpolate.splrep(data[:,col0], data[:,col1])))

  #c_points = data[:,1] + _adjust

col0=0
col1=2
count = 1
while (count < len(sys.argv)):

  if (sys.argv[count] == "-col0"):
    count = count + 1
    if (count > len(sys.argv)):
      break
    col0 = int(sys.argv[count])
  elif (sys.argv[count] == "-col1"):
    count = count + 1
    if (count > len(sys.argv)):
      break
    col1 = int(sys.argv[count])
  elif (sys.argv[count] == "-at"):
    count = count + 1
    if (count > len(sys.argv)):
      break
    at = float(sys.argv[count])

  count = count + 1

y = ywinterp(at)/4
print('{:.2f},{:.4f}'.format(y, (np.log(y*(2*(2+1)-y)+1)-np.log(2))/26))
