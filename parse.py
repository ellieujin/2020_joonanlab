import os,sys

f = sys.argv[1]
with open(f) as fh:
	for l in fh:
		if l[0] != '#':
			info = l.split('\t')
			FILTER = info[6]
			gt_m = info[9].split(':')[0]
			gt_p = info[10].split(':')[0]
			gt_pm = info[11].split(':')[0]
			if FILTER == 'PASS' and gt_m == '0/0' and gt_p == '0/0' and gt_pm == '0/1':
				print (l)
