#!/usr/bin/python
import optparse
import sys
from Utilities import SeqIterator
import multiprocessing
import math
import datetime

"""
@author: Jacob Porter
@since: 02/16/2016
@summary: AAn information based read trimmer
Written by Jacob Porter while at Virginia Tech
http://www.jacobporter.com

This program implements a read trimming strategy that maximizes information content.
It has expressions for a minimum threshold amount, length, correctness probability, and entropy.
For each read in the file, the linear time algorithm takes time O(n) where n is the read length.
This algorithm trims the read to a prefix of the original read.
The non-linear time algorithm searches for a substring that maximizes the score and returns that substring.
For each read in the file, the non-linear algorithm has time O(n^2).
This read trimmer is intended for Illumina style short DNA reads encoded as A, C, T, G in a fastq file format
with phred scores.  If the format is incorrect, consider converting to fastq format.
Multiprocessing is implemented with a producer / consumer design pattern using synchronized queues.

TODO: add a term or support for trimming out regions with homopolymer runs?
TODO: compare to other read trimmers, analyze on read data
TODO: Is the entropy only option redundant? Probably
TODO: Is balanced redundant? Probably
TODO: Should the scaling of the terms be changed because their contribution to the score will differ.
TODO: Should an exponential or logistic model or no model be used for biasing the non-linear version towards the start of the read.
"""

def calculate_read_trim_position_linear(read, m, alpha, beta, gamma, offset, balanced=False, entropyOnly=False):
    """
    This implements the linear version of the read trimmer.
    Read is a tuple.  It has the id string, the read string, and the phred string.
    alpha controls the influence of the length
    beta controls the influence of the correctness probability
    gamma controls the influence of entropy
    balanced uses the geometric mean
    entropyOnly finds the region of maximum entropy considering the length cutoff
    """
    num_terms = 3.0
    sequence = read[1]
    phred = read[2]
    total_length = len(read[1])
    max_position = 0
    max_score = float("-inf")
    initial_prob = 1
    frequencies = [0,0,0,0,0]
    for i in range(total_length):
        base = sequence[i].upper()
        prob_correct = 1 - 10**((ord(phred[i])-offset)/(-10.0))
        initial_prob *= prob_correct
        if base == 'A':
            frequencies[0] += 1
        elif base == 'C':
            frequencies[1] += 1
        elif base == 'T':
            frequencies[2] += 1
        elif base == 'G':
            frequencies[3] += 1
        else:
            frequencies[4] += 1
        entropy = 0
        for f in frequencies:
            if not (f == 0):
                my_freq = f / (i+1.0)
                entropy += my_freq*math.log(my_freq, 2)
        entropy = -entropy
        #if not entropyOnly:
        length_cutoff = 1 / (1+math.exp(m-(i+1.0)))
        if balanced and not entropyOnly:
            my_score = length_cutoff * math.pow( (i+1.0) * initial_prob * entropy, 1/3.0)
        elif not balanced and not entropyOnly:
            my_score = length_cutoff * math.pow((i+1.0), alpha/num_terms) * math.pow(initial_prob, beta/num_terms) * math.pow(entropy/2.0, gamma/num_terms)
        else:
            my_score = length_cutoff * entropy
        if my_score >= max_score:
            max_score = my_score
            max_position = i
    return (max_position, max_score)


def calculate_read_trim_position_dynamic(read, m, t, alpha, beta, gamma, offset, balanced=False, entropyOnly=False):
    """
    This implements the non-linear version of the read trimmer.
    It calls the linear version for each starting position.
    The position of the substring with the maximum score is returned.
    """
    sequence = read[1]
    phred = read[2]
    total_length = len(read[1])
    max_begin = 0
    max_end = 0
    max_score = float("-inf")
    for i in range(total_length):
        sliced_sequence = sequence[i:total_length]
        sliced_phred = phred[i:total_length]
        sliced_read = ("", sliced_sequence, sliced_phred)
        my_position, my_score = calculate_read_trim_position_linear(sliced_read, m, 
                                                                    alpha, beta, gamma, 
                                                                    offset, balanced=balanced, 
                                                                    entropyOnly=entropyOnly)
        if t != None and not entropyOnly:
            my_score = math.exp(-i/t) * my_score
        if my_score > max_score:
            max_begin = i
            max_end = my_position
            max_score = my_score
    return (max_begin, max_end+max_begin, max_score)


def trim_read(read, alpha, beta, gamma, m, t, offset, balanced, linear, entropyOnly, verbose):  
    """
    This function calls the functions that calculate the positions on the read to trim.
    It branches depending on whether the non-linear version is used.
    All other parameters are passed on.
    """       
    if linear:
        position, score = calculate_read_trim_position_linear(read, m, alpha, beta, gamma, offset,
                                                              balanced=balanced, entropyOnly=entropyOnly)
        position += 1
        if verbose:
            sys.stderr.write("%s  position: %s score: %s\n" % (read[0], str(position), str(score)))
            sys.stderr.flush()
        trimmed_read = (read[0], read[1][0:position], read[2][0:position])
    else:
        begin, end, score = calculate_read_trim_position_dynamic(read, m, t, alpha, beta, gamma, offset,
                                                                 balanced=balanced, entropyOnly=entropyOnly)
        end += 1
        if verbose:
            sys.stderr.write("%s  begin: %s end: %s score: %s\n" % (read[0], str(begin), str(end), str(score)))
            sys.stderr.flush()
        trimmed_read = (read[0], read[1][begin:end], read[2][begin:end])
    return trimmed_read
 


def trim_reads_single(reads, trimmed_output, alpha, beta, gamma, m, t, offset, balanced, linear, entropyOnly, verbose):
    """
    This is the main function that loops through the reads.
    The linear option is a boolean that selects for either the linear or the non-linear version.
    The floating point parameters r, s, m, and t control the tradeoffs in the scoring function.
    The balanced option is a boolean that uses the geometric mean of length, correctess probability, and entropy.  This option ignores s and r.
    The parameter t is only used in the non-linear version.  It controls the exponential model for weighting the starting position.
    reads is an iterator for the reads to trim and trimmed_output is an iterator for writing the reads.
    """
    for read in reads:
        trimmed_read = trim_read(read, alpha, beta, gamma, m, t, offset, balanced, linear, entropyOnly, verbose)
        trimmed_output.write(trimmed_read)
    

            
def trim_reads_multi(q_read, q_write, alpha, beta, gamma, m, t, offset, balanced, linear, entropyOnly, verbose):
    """
    A function for trimming reads with multiple processes.
    q_read is the queue for reading reads
    q_write is the queue for adding trimmed reads to.
    The other arguments are passed to the function that does all the work.
    A 'poison pill' is put onto the q_write queue so that the writer process knows when to quit.
    """
    #print str(os.getpid())
    while(True):
        read = q_read.get()
        if read is None:
            q_write.put(None)
            return
        else:
            trimmed_read = trim_read(read, alpha, beta, gamma, m, t, offset, balanced, linear, entropyOnly, verbose)
            q_write.put(trimmed_read)
            

def add_reads_to_queue_multi(q_read, input_file, num_processes):
    """
    This function puts the reads from the input_file into the Queue q_read.
    This function is called to add all the reads on a queue for other processes to consume.
    A 'poison pill' is added for all num_processes.  When other processes get this signal, they will stop.
    """
    reads = SeqIterator.SeqIterator(input_file, file_type='fastq')
    for read in reads:
        q_read.put(read)
    for i in range(num_processes):
        q_read.put(None)
    

def trim_reads(input_file, output_file, alpha, beta, gamma, m, t, offset=33, process_amount=1, balanced=False, linear=True, entropyOnly=False, verbose=False): 
    """
    This function is the entry point into the read trimmer.  It is the function that could be called
    by other programs to do read trimming.
    This function creates processes if multiprocessing is used.
    """
    trimmed_output = SeqIterator.SeqWriter(output_file, file_type='fastq')
    if process_amount == 1:
        try:
            reads = SeqIterator.SeqIterator(input_file, file_type='fastq')
        except IOError:
            sys.stderr.write("Something is wrong with the reads file.  Please check it.\n")
            return
        trim_reads_single(reads, trimmed_output, alpha, beta, gamma, m, t, offset, balanced, linear, entropyOnly, verbose)
        read_count = reads.records_processed()
    else:
        # Need to add records to a queue for the other processes.
        processes = []
        q_read = multiprocessing.Queue()
        q_write = multiprocessing.Queue()
        proc_read = multiprocessing.Process(target=add_reads_to_queue_multi, args=(q_read, input_file, process_amount))
        proc_read.start()
        # Start the processes.  Each one pulls a record from the queue and processes it.
        for i in range(process_amount):
            proc = multiprocessing.Process(target=trim_reads_multi, 
                                                   args=(q_read, q_write, alpha, beta, gamma, m, t, offset,  
                                                         balanced, linear, entropyOnly, verbose))
            
            processes.append(proc)
            proc.start()
        num_quit = 0
        read_count = 0
        # The calling process writes the reads to the output file.
        while(True):
            trimmed_read = q_write.get()
            if trimmed_read is None:
                num_quit += 1
                if num_quit == process_amount:
                    break
            else:
                read_count += 1
                trimmed_output.write(trimmed_read)
        proc_read.join()
        for proc in processes:
            proc.join()
    return read_count

def main():
    """
    This is the main file that defines the interface and executes the functions that do all
    of the read trimming.  Multiple processes are started and joined if the processes parameter
    is set to a value larger than one.
    """
    now = datetime.datetime.now()
    p = optparse.OptionParser()
    p.add_option('--reads', '-f', help='Fastq reads file location.', default=None)
#     p.add_option('--pair1', '-1', help='Fastq reads file for the first mate pair.', default=None)
#     p.add_option('--pair2', '-2', help='Fastq reads file for the second mate pair.', default=None)
    p.add_option('--processes', '-p', help='Number of processes to use.', default=1)
    p.add_option('--output','-o',help='The file name for the output file.', default=sys.stdout)
    p.add_option('--phred','-s',help='The tradeoff value for the phred score.', default=0.25)
    p.add_option('--entropy','-r',help='The tradeoff value for the entropy feature.', default=0.5)
    p.add_option('--min_length','-m',help='The minimum length target.', default=25)
    p.add_option('--start_weight', '-t', help='The weight for the exponential model for the starting position of the read for the non-linear option.  If set to 0, it will not be used.', default=0)
    #p.add_option('--balanced', '-b', help='Use the geometric mean of length, correctness probability, and entropy.', action="store_true", default=False)
    p.add_option('--non_linear', '-n', help='Use the non-linear approach to search for a substring that maximizes the score.', action="store_true", default=False)
    #p.add_option('--entropyOnly', '-e', help='Use the dynamic programming approach.', action="store_true", default=False)
    p.add_option('--ascii_offset', '-a', help='The offset value for the ascii encoding of phred scores.', default=33)
    p.add_option('--verbose', '-v', help='Run in verbose.  The score for each read is printed to stderr.', action="store_true", default=False)
    options, _ = p.parse_args()
    sys.stderr.write("Starting InfoTrim by Jacob Porter.\n")
    sys.stderr.flush()
    balanced = False #options.balanced
    linear = not options.non_linear
    entropyOnly = False #options.entropyOnly
    verbose = options.verbose
    try:
        r = float(options.entropy)
        s = float(options.phred)
        m = float(options.min_length)
        t = float(options.start_weight)
        if t == 0.0 or t == 0:
            t = None
        if not balanced and not entropyOnly:
            if r + s > 1 or r + s < 0:
                sys.stderr.write("The parameters r=%s (phred) and s=%s (entropy) need to add up to a value in [0,1].\n" % (str(s), str(r)))
                return
            #Is this formula correct?
            #sr_over_3 = s*r/3.0
            #alpha = (1-s)*(1-r) + sr_over_3
            #beta = s * (1-r) + sr_over_3
            #gamma = (1-s)*r + sr_over_3
            alpha = 1 - r - s
            beta = s
            gamma = r
        else:
            alpha = None
            beta = None
            gamma = None 
    except ValueError:
        sys.stderr.write("One of the parameters r, s, m, or t  was not set to a number.\n")
        sys.stderr.flush()
        return
    try:
        offset = int(options.ascii_offset)
    except ValueError:
        sys.stderr.write("The ascii offset must be an integer.\n")
        return
    input_file = options.reads
    if input_file == None:
        p.print_help()
        sys.exit(1)
    output_file = options.output
    output_file_string = "stdout"
    if not (output_file == sys.stdout):
        try:
            output_file_string = output_file
            output_file = open(output_file, 'w')
        except IOError:
            sys.stderr.write("Something is wrong with opening the output file.  Please check it.\n")
            return
    try:
        process_amount = int(options.processes)
    except ValueError:
        sys.stderr.write("The number of processes should be an integer.\n")
        return
    sys.stderr.write("Using file '" + str(input_file) + "'. Please wait.\n")
    sys.stderr.write("Current time: " + str(now) + " \n")
    sys.stderr.write("Writing output to file '" + output_file_string + "'.\n")
    sys.stderr.write("Using the following arguments.\n")
    sys.stderr.write("r=%s s=%s m=%s t=%s processes=%s ascii_offset=%s, balanced=%s linear=%s entropyOnly=%s verbose=%s\n"
                      % (str(r), str(s), str(m), str(t), str(process_amount), str(offset), 
                         str(balanced), str(linear), str(entropyOnly), str(verbose)))
    sys.stderr.flush()
    read_count = trim_reads(input_file, output_file, alpha, beta, gamma, m, t, offset=offset, process_amount=process_amount, balanced=balanced, linear=linear, entropyOnly=entropyOnly, verbose=verbose)
    sys.stderr.write("Finished trimming the reads.  There were " + str(read_count) + " reads processed.\n")
    now = datetime.datetime.now()
    sys.stderr.write("Current time: " + str(now) + " \n")
    sys.stderr.write("Have a good day!\n")
    sys.stderr.flush()

if __name__ == '__main__':
    main()