#!/usr/bin/python
import optparse
import sys
import os
import subprocess
import multiprocessing

"""
Create distances from alignment files using Quicktree
"""
def run_quicktree(path_to_quicktree, path_to_alignments, filename, path_to_distances):
    dist_filename = filename[0:len(filename)-4]
    dist_filename += ".dist"
    alignment_path = os.path.join(path_to_alignments, filename)
    distance_path = os.path.join(path_to_distances, dist_filename)
    subprocess.call([path_to_quicktree, '-out',  'm', alignment_path], stdout=open(distance_path, 'w'))

def join_processes(my_processes):
    for proc in my_processes:
        proc.join()

def main():
    p = optparse.OptionParser(usage="usage: %prog [options] path_to_quicktree path_to_alignment_directory path_to_distance_directory begin end", 
                              version="%prog 0.0")
    p.add_option('--processes', '-p', help='Number of processes to use', default='6')
    options, args = p.parse_args()
    if len(sys.argv[1:]) == 0:
        print "%sA required argument is missing.  Use -h or --help for argument information."
        p.print_help()
        exit(1)
    num_processes = int(options.processes)
    path_to_program = args[0]
    mypath = args[1]
    path_to_distances = args[2]
    begin = int(args[3])
    end = int(args[4])
    onlyfiles = [ f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath,f)) ]
    onlyfiles.sort()
    #num_files = len(onlyfiles)
    count = begin
    my_processes = []
    while count < end:
        if len(my_processes) == num_processes:
            join_processes(my_processes)
            my_processes = []
        p1 = multiprocessing.Process(target=run_quicktree, args=(path_to_program, mypath, onlyfiles[count], path_to_distances))
        p1.start()
        my_processes.append(p1)
        count += 1
    join_processes(my_processes)
    print "Finished!"                                     

if __name__ == '__main__':
    main()