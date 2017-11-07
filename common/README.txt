## Add Modules ##	
	- basicImport.py consists Basics Libs
	- Add Modules/Packages to myImport.py. 
	  Write Modules then add to myImport.py
	  Remeber to add: 
	      """
	      from basicImport import *
	      """
	  to new Modules

	- No need to change __init__.py

	- See ISMDust/references/someFolder as an example.


	- Every script needs these lines:

	"""
	import sys, os
	sys.path.insert(0, os.getenv("HOME")+'/ISMDust') # add folder of Class/Modules

	from common.myImport import *
	"""



## Add style ##
 # 1. Check dir: matplotlib.get_configdir() ; it's usually: '/home/vnguyen/.config/matplotlib'
 # 2. mkdir stylelib
 # 3. create files: stylelib/mystyle.mplstyle