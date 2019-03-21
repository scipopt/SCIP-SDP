import os
import sys

def convert(filename):
	basename = os.path.basename(filename).split('.')[0]
	print(basename)
	os.system('cd /home/gally/gally_SCIPSDP/ \n ./bin/scipsdp -c "read '+ filename + '" -c "write prob /local/gally/cbffiles/' + basename + '.cbf" -c "quit"')
        os.system('gzip /local/gally/cbffiles/' + basename + '.cbf')
	

if __name__=="__main__":
	"""takes a text file as input with .dat-s files (full path) each line, for each of these SCIP-SDP is called to transform them to Cbf"""
	num_lines = sum(1 for line in open(sys.argv[1]))
	file=open(sys.argv[1])
	i=0				
	for line in file:
		convert(line.rstrip('\n'))
		i=i+1
		print str(i) + "/" + str(num_lines) + " finished"
	file.close()
