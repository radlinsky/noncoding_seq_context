import pdb

try:
	print "hello"
	1/0
except:
	pdb.set_trace()
