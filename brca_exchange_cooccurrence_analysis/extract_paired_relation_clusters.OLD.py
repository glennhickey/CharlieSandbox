import argparse, sys, csv
from collections import defaultdict


def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input paired-relationship csv file.')
    parser.add_argument('-i', '--inPairCSV', type=str,
        help='Input paired-relation csv file.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options

def main(args):

    options = parse_args()
    
    clusters_1 = defaultdict(set)
    clusters_1_idx = 0
    
    # Cycle through each pair and assign a set cluster for each
    with open(options.inPairCSV, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            print("{}\t{}".format(row[0],row[1]))
            if "ID" in row[0]: continue
            clusters_1[clusters_1_idx].add(row[0])
            clusters_1[clusters_1_idx].add(row[1])
            clusters_1_idx += 1
    
    clusters_2 = { 0:clusters_1[0] }
    finished_merge = False
    while not finished_merge:
        merge = False
        print("Restart iteration")
        for idx1,set1 in clusters_1.iteritems():
            print("index_1: {}".format(idx1))
            for idx2,set2 in clusters_2.iteritems():
                print("index_2: {}".format(idx2))
                if bool(set1 & set2):
                    print("Merge clusters!")
                    clusters_2[idx2] |= set1
                else:
                    print("Add new cluster!")
                    clusters_2[idx2+1].append(set1)
        finished_merge = True
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

