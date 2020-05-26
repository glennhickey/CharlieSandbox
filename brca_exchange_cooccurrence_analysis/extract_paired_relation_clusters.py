import argparse, sys, csv, copy
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
    
    clusters = defaultdict(set)
    clusters_idx = 0
    
    # Cycle through each pair and assign a set cluster for each
    first_key = None
    first_value = None
    with open(options.inPairCSV, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            print("{}\t{}".format(row[0],row[1]))
            if "ID" in row[0]: continue
            if not first_key and not first_value:
                first_key = row[0]
                first_value = row[1]
            clusters[row[0]].add(row[1])
    
    clusters_2 = clusters.copy()
    for k1,v1 in clusters.items():
        for k2,v2 in clusters_2.items():
            if (k1 in v2) or (k2 in v1) or len(v1 & v2) >= 1:
                print('merging k1: {}, deleting k2: {}'.format(k1,k2))
                clusters[k1].update(v2)
                clusters[k1].add(k2)
                if clusters[k2]: del clusters[k2]
    print("Finished 1st pass")
    #import pdb; pdb.set_trace()
    with open(options.outReport, 'w') as cluster_output_file:
        for k,v in clusters.items():
            cluster_output_file.write("{}\t{}\n".format(k, '\t'.join(v)))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

