#!/usr/bin/env python3
import os
import sys
lines_per_file = 9745
new_ped= None
current_directory = sys.argv[1]
final_directory = os.path.join(current_directory, r'files_ts')
if not os.path.exists(final_directory):
   os.makedirs(final_directory)
with open(os.path.join(current_directory,r'snp_dat_trans.tped')) as pedfile:
    for line_num, line in enumerate(pedfile):
        if line_num % lines_per_file == 0:
            if new_ped:
                new_ped.close()
            new_ped_filename = os.path.join(final_directory, 'fileR_{}.txt'.format(line_num))
            new_ped = open(new_ped_filename, "w")
        new_ped.write(line)
    if new_ped:
        new_ped.close()