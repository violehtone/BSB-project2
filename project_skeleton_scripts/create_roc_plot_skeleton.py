#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, exports the coordinates and calculates the AUC.

import argparse
import numpy
import matplotlib
matplotlib.use('AGG')
import pylab

def parse_blast_results(filename):
    """
    Parse every protein pair's e-value out of a BLAST results file.
    :param filename: input file with BLAST results.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding e-value from BLAST (value).
            {(id, id) : e-value}
    """

    blast_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
    
            if len(arr) != 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])
                
            key = (arr[0], arr[1])
            if arr[2] == "NA":    # Substitute protein pairs whose e-value is
                value = 1e6       # not available with an e-value of 1 million
            else:
                value = float(arr[2])
            blast_results[key] = value

    return blast_results


def parse_benchmark_results(filename):
    """
    Parse every protein pair's classification out of the benchmark file.
    :param filename: input file with benchmark classifacations.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding call (value).
            {(id, id) : "similarity"}
    """

    benchmark_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
    
            if len(arr) < 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])

            # Benchmark classifications are symmetric, so add both possible keys:                
            key1 = (arr[0],arr[1])
            key2 = (arr[1],arr[0])
            value = arr[2]
            benchmark_results[key1] = value
            benchmark_results[key2] = value

    return benchmark_results


def integrate(x, y):
    """
    Calculate the Area Under the Curve (AUC) for a given list of coordinates
    :param x: a list of x-coordinates
    :param y: a list of y-coordinates
    :return: a float with the surface area under the curve described by x and y
    """
    auc = 0.
    last_x = x[0] # the i-1 x-coordinate in the list
    last_y = y[0] # the i-1 y-coordinate in the list
    for cur_x, cur_y in list(zip(x, y))[1:]:
        #########################
        ### START CODING HERE ###
        #########################
        auc += ((cur_x - last_x) * cur_y)
        #########################
        ###  END CODING HERE  ###
        #########################
        last_x = cur_x
        last_y = cur_y
    return auc


def roc_plot(blast_evalues, benchmark_dict, png_filename):
    """
    Draw the ROC plot for a given set of e-values and corresponding benchmark classifications.

    :param blast_evalues: the dictionary produced by parse_blast_results()
    :param benchmark_dict: the dictionary produced by parse_benchmark_results()
    """

    ### Create the lists of coordinates

    x = [0] # array of the ROC plot's x-coordinates: False Positive Rate = FP/(FP+TN)
    y = [0] # array of the ROC plot's y-coordinates: True  Positive Rate = TP/(TP+FN)
            # If the e-value is the same as the last one, only increase x or y of the last coordinate
    last_evalue = -1
    evalues = [(v, k) for k, v in blast_evalues.items()] # List of tuples consisting of (evalue, protein_pair)
    sorted_evalues = sorted(evalues)

    #initialize variables (i.e. fn = false negative)
    variable_map = {"FN":0.0, "TN":0.0, "FP":0.0, "TP":0.0}

    for evalue, protein_pair in sorted_evalues:
        # Iterate through the protein pairs, in order of ascending e-value
        # Determine whether it is 
        #    different -> actual negative, thus a false positive (x) 
        #    similar   -> actual positive, thus a true positive (y)
        # Increase the respective value and add a new coordinate for every unique e-value
        # If the e-value is the same as the last one, only increase x or y of the last coordinate
        # Ignore entries in the benchmark_dict classified as "ambiguous" and decide how to handle blast NA results
        try:
            ## Benchmark is either "different", "ambiguous", or "similar"
            benchmark = benchmark_dict[protein_pair]

            if (evalue != last_evalue):

                calculate_new_variable_values(evalue, benchmark, variable_map)
                false_pos_rate = FPR(variable_map["FP"], variable_map["TN"])
                true_pos_rate = TPR(variable_map["TP"], variable_map["FN"])

                x.append(false_pos_rate)
                y.append(true_pos_rate)

            ## if e-value is the same as last e-value
            else:
                calculate_new_variable_values(evalue, benchmark, variable_map)
                false_pos_rate = FPR(variable_map["FP"], variable_map["TN"])
                true_pos_rate = TPR(variable_map["TP"], variable_map["FN"])

                x[len(x)-1] = false_pos_rate
                y[len(y)-1] = true_pos_rate

        except KeyError:
            ## If benchmark was not found, continue to next e-value
            pass

        last_evalue = evalue

    # In order to get the rates for every coordinate we divide by the total number (last entry)
    print(x)
    print(y)
    
    x = numpy.array(x) / float(x[-1])
    y = numpy.array(y) / float(y[-1])

    ### Figure out the AUC
    auc = integrate(x, y)

    ### Draw the plot and write it to a file
    pylab.plot(x, y)
    pylab.plot([0,1],[0,1],'--k')
    pylab.xlabel('False Positive Rate')
    pylab.ylabel('True Positive Rate')
    pylab.title('AUC = %.3f' % auc)
    pylab.savefig(png_filename)

    ### Write coordinates to a file
    with open(png_filename.split('.')[0] + '_xy.tsv','w') as f:
        for a,b in zip(x,y):
            f.write(str(a) + '\t' + str(b) + '\n')


def FPR(fp, tn):
    if fp == tn == 0:
        return 0
    return fp / (fp+tn)


def TPR(tp, fn):
    if tp == fn == 0:
        return 0

    return tp / (tp + fn)


def calculate_new_variable_values(evalue, benchmark, variable_map):
    homolog = True

    ## Make 'NA's evaluate to non-homologs
    if evalue >= 1000000 or evalue == 'NA':
        homolog = False

    #FALSE NEGATIVE
    if (not homolog and benchmark == 'similar'):
        variable_map["FN"] += 1
    #TRUE NEGATIVE
    elif (not homolog and benchmark == 'different'):
        variable_map["TN"] += 1
    #FALSE POSITIVE
    elif (homolog and benchmark == 'different'):
        variable_map["FP"] += 1
    #TRUE POSITIVE
    elif (homolog and benchmark == 'similar'):
        variable_map["TP"] += 1


def main(blast_results_file, benchmark_results_file, png_file):
    # Parse the input files and retrieve every protein pair's e-value and benchmark classification.
    blast_evalues = parse_blast_results(blast_results_file)
    benchmark_results = parse_benchmark_results(benchmark_results_file)
    
    # Draw and save the ROC plot
    roc_plot(blast_evalues, benchmark_results, png_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw and save a ROC plot to a file")
    parser.add_argument("-iblast","--input_blast_results", help="tab-separated BLAST results file", required=True)
    parser.add_argument("-ibench","--input_benchmark_results", help="tab-separated benchmark classification file", required=True)
    parser.add_argument("-o", "--output_png", help="output png file", required=True)

    args = parser.parse_args()
    
    blast_file = args.input_blast_results
    benchmark_file = args.input_benchmark_results
    png_file = args.output_png

    main(blast_file,benchmark_file, png_file)
