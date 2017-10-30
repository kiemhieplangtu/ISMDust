import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import module            as md

## Common class ##
from restore             import restore # Read txt file, csv file

### MAIN ###
info = md.read_info_93src_csv(fname = '../doc/observed_src.csv',asarray=False)
src  = info['src']

print len(src)

for i in range(len(src)):
	print i, src[i]