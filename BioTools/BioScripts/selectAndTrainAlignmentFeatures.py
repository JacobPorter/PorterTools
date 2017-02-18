#!/usr/bin/python
import optparse
import sys
import multiprocessing
import os
import subprocess

def main():
    p = optparse.OptionParser(usage="usage: %prog [options] path_to_MSVMpack features_file predict_file", 
                              version="%prog 0.0")
    p.add_option('--processes', '-p', help='Number of processes to use', default='7')
    options, args = p.parse_args()
    if len(sys.argv[1:]) == 0:
        print "%sA required argument is missing.  Use -h or --help for argument information."
        p.print_help()
        exit(1)    
    path_to_MSVMpack = args[0]
    features_file = args[1]
    predict_file = args[2]
    num_processes = int(options.processes)
    features_train = makeFeatureFiles(features_file)
    print "Finished creating training set"
    makeFeatureFiles(predict_file)
    print "Finished creating prediction set"
    count = 0
    my_processes = []
    num_files = len(features_train)
    while count < num_files:
        if len(my_processes) == num_processes:
            join_processes(my_processes)
            my_processes = []
        p1 = multiprocessing.Process(target=trainAndPredict, args=(path_to_MSVMpack, features_file, predict_file, features_train[count]))
        p1.start()
        my_processes.append(p1)
        count += 1
    join_processes(my_processes)

def join_processes(my_processes):
    for proc in my_processes:
        proc.join()

def trainAndPredict(path_to_MSVMpack, features_file, predict_file, feature):
    trainer = os.path.join(path_to_MSVMpack, "trainmsvm")
    predicter = os.path.join(path_to_MSVMpack, "predmsvm")
    feature_file_train = os.path.join(features_file + "." + feature)
    feature_file_predict = os.path.join(predict_file + "." + feature)
    output = open(feature_file_predict + ".outputs", 'w')
    subprocess.call([trainer, feature_file_train, feature_file_train+".model"])
    subprocess.call([predicter, feature_file_predict, feature_file_train+".model"], stdout=output)

def makeFeatureFiles(features_file):
    feature_file_desc = open(features_file)
    skipMe = True
    feature_dict = {}
    for line in feature_file_desc:
        if skipMe:
            header =  ["length", "seq_entropy", "first_20_entropy", 
                    "A_freq", "C_freq",  "G_freq", "T_freq", "N_freq",
                    "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                    "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", 
                    "length_X_seq_entropy", "mismatch_gap_percent" , "mapping_type"]
            
            for h in header:
                feature_dict[h] = []
            skipMe = False
        else:
            observation = line.strip("\n").split("\t")
            for i in xrange(len(observation)):
                feature_dict[header[i]].append(observation[i])
    for feature in feature_dict.keys():
        number_of_observations = len(feature_dict[feature])
        num_features = 1
        feature_write = open(features_file + "." + feature, 'w')
        feature_write.write(str(number_of_observations)+"\n")
        feature_write.write(str(num_features)+"\n") 
        for i in xrange(len(feature_dict[feature])):
            feature_write.write(feature_dict[feature][i]+"\t")
            feature_write.write(feature_dict["mapping_type"][i] + "\n")
    list1 = feature_dict.keys()
    list1.remove("mapping_type")
    return list1

if __name__ == '__main__':
    main()