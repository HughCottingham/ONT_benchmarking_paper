import sys
import csv
import os
import glob

infile=sys.argv[1]
genome_size=int(sys.argv[2])
isolate_name=sys.argv[3]

with open(infile,'r') as input_file:
	csv_reader = csv.reader(input_file,delimiter='\t')
	IDs=[]
	depth_count=10
	cumsum=0
	print(f'Writing cumulative read lists for different subsets of {isolate_name}')
	for line in csv_reader:
		ID=line[0]
		length=int(line[1])
		IDs.append(ID)
		cumsum=cumsum+length
		depth=cumsum/genome_size
		with open(f'{isolate_name}_{depth_count}x_reads.txt','w') as outfile:
			if depth > 10:
				print(f'{depth_count}x read list complete.')
				depth_count=depth_count+10
				cumsum=0
				for ID in IDs:
					outfile.write(ID+'\n')
	print('Removing empty subset files (ie those that don\'t have enough reads)')
	subset_list = glob.glob(f'{isolate_name}*x*reads.txt')
	for subset in subset_list:
	    if os.path.getsize(subset) == 0:
	        os.remove(subset)
	print('All done!')
