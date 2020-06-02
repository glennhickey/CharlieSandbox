import argparse, sys, csv, copy
from collections import defaultdict
import timeit

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
    with open(options.inPairCSV, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            print("{}\t{}".format(row[0],row[1]))
            if "ID" in row[0]: continue
            clusters[row[0]].add(row[1])
    
    clusters_2 = clusters.copy()
    key_del_list = list()
    for c, k1 in enumerate(sorted(clusters.keys())):
        v1 = clusters[k1]
        for k2 in sorted(clusters_2.keys()):
            v2 = clusters_2[k2]
            if (k1 in v2) or len(v1 & v2) >= 1 or (k1 == k2):
                clusters[k2].update(v1)
                clusters[k2].add(k1)
                if k1 != k2:
                    key_del_list.append(k1)
        print("line: {}/{}".format(c,'164351'))
    print("Finished 1st pass")
    clusters_final = clusters.copy()
    for key in key_del_list:
        if key in clusters_final:
            del clusters_final[key]
    print("Finished deletion pass")
    import pdb; pdb.set_trace()
    with open(options.outReport, 'w') as cluster_output_file:
        for k,v in clusters.items():
            cluster_output_file.write("{}\t{}\n".format(k, '\t'.join(v)))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

