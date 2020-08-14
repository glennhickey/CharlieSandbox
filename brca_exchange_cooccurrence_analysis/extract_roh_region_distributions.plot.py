import matplotlib
matplotlib.use('Agg')
import vcf, argparse, sys
import numpy as np
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import vcf, argparse, sys
import numpy as np
import pandas as pd
import math
import seaborn as sns
import matplotlib.pyplot as plt

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input bcftools roh tab-delimited file and output roh report and histogram.')
    parser.add_argument('-i', '--inROHdistA', type=str,
        help='Input 1st roh distribution filepath.')
    parser.add_argument('-j', '--inROHdistB', type=str,
        help='Input 2nd roh distribution filepath.')
    parser.add_argument('-l', '--regionLength', type=int,
        help='Input length of region used in calculating SROH.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output plot filename.')

    options = parser.parse_args()
    return options

def main(args):

    options = parse_args()
    
    roh_distribution_dict = defaultdict(list)
    with open(options.inROHdistA, 'r') as roh_file_a, open(options.inROHdistB, 'r') as roh_file_b:
        for line in roh_file_a:
            if 'sample_id' in line: continue
            parsed_line = line.strip().split('\t')
            roh_distribution_dict['SROH'].extend([float(parsed_line[2]),float(parsed_line[3]),float(parsed_line[4]),float(parsed_line[5]),float(parsed_line[6])])
            roh_distribution_dict['SROH_length'].extend(['all','100kb','1mb', '1500kb', '5mb'])
            roh_distribution_dict['group'].extend(['No']*5)
        for line in roh_file_b:
            if 'sample_id' in line: continue
            parsed_line = line.strip().split('\t')
            roh_distribution_dict['SROH'].extend([float(parsed_line[2]),float(parsed_line[3]),float(parsed_line[4]),float(parsed_line[5]),float(parsed_line[6])])
            roh_distribution_dict['SROH_length'].extend(['all','100kb','1mb', '1500kb', '5mb'])
            roh_distribution_dict['group'].extend(['Yes']*5)
    
    violin_df = pd.DataFrame(data=roh_distribution_dict) 
    sns.set(style="whitegrid", font_scale=1.5)
    
    fig, axes = plt.subplots(figsize=(10, 10)) 
    order=["all", "100kb", "1mb", "1500kb", "5mb"]
    sns.boxplot(
        x="SROH_length", y="SROH", hue="group", data=violin_df, 
        order=order,
        ax=axes
    )
    axes.set_xticklabels(["All", "100 (kb)", "1 (mb)", "1.5 (mb)", "5 (mb)"])
    axes.set_xlabel("Minimum ROH Length")
    axes.legend("")
    fig.savefig("roh_distribution_violin.{}.png".format(options.outReport)) 
    matplotlib.pyplot.close(fig)
     
if __name__ == "__main__":
    sys.exit(main(sys.argv))

