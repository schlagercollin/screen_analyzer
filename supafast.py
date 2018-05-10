import sys
import numpy as np
import csv
import pandas as pd
from collections import OrderedDict
import time
import math

#idfile = sys.argv[1]

def parse_lib(filename):
    """Creates a dictionary whose keys are the target sequence and whose values are the gene data"""
    print("Creating dict...")
    library_dict = {}
    with open(filename) as lib_file:
        reader = csv.DictReader(lib_file)  # read rows into a dictionary format
        for row in reader:
            library_dict[row["sgRNA Target Sequence"]] = 0
    print("Done. %d" % (len(library_dict)))
    return library_dict

def create_dataframe(library_file):
    print("Creating dataframe...")
    library_df = pd.read_csv(library_file) # read in library file as pandas data frame
    library_df.drop_duplicates(subset="sgRNA Target Sequence", keep="first", inplace=True)
    return library_df

def add_df_column(library_df, column_dict, column_name):
    library_df[column_name] = np.zeros(len(library_df), dtype="int") # create new column of zeros for Counts
    library_df.sort_values(by=["sgRNA Target Sequence"], inplace=True) # sort by sgRNA seq
    sorted_counts_dict = OrderedDict(sorted(column_dict.items(), key=lambda t: t[0])) # sort by sgRNA seq (align with df)
    counts_list = list(sorted_counts_dict.values()) # sort counts list
    assert (len(counts_list) == len(library_df))
    library_df[column_name] = counts_list # set column equal to list
    print("Done.")


def parse_qfast(qfast_file, library_file):
    print("Paring %r" % qfast_file)
    counts_dict = parse_lib(library_file)
    scan_counter = 0
    print("Iterating over file...")
    with open(qfast_file) as f:
        for line_number, sequence in enumerate(f):
            if line_number % (4) == 1:
                sequence = sequence[0:20]
                if sequence in counts_dict:
                    counts_dict[sequence] += 1
                scan_counter += 1
    print("Done.")
    return counts_dict

def compute_lfc(dataframe):
    log2 = lambda x: math.log2(x)
    sortedTotMatches = int(dataframe["Sorted Counts"].sum())
    unsortedTotMatches = int(dataframe["Unsorted Counts"].sum())
    normalized_sorted = (dataframe["Sorted Counts"]+1)/sortedTotMatches
    normalized_unsorted = (dataframe["Unsorted Counts"]+1)/unsortedTotMatches
    foldchange = normalized_sorted/normalized_unsorted
    dataframe["LFC"] = foldchange.apply(log2)
    return dataframe

def perform_analysis(sorted_fastq, unsorted_fastq, library_file, output_file):
    print("Beginning parsing...")
    library_df = create_dataframe(library_file)
    sorted_counts_dict = parse_qfast(sorted_file, library_file)
    add_df_column(library_df, sorted_counts_dict, "Sorted Counts")
    unsorted_counts_dict = parse_qfast(unsorted_file, library_file)
    add_df_column(library_df, unsorted_counts_dict, "Unsorted Counts")
    compute_lfc(library_df)
    library_df.sort_values(by=["LFC"], ascending=False, inplace=True)
    library_df.to_csv(output_file)
    return library_df

if __name__ == "__main__":
    start = time.time()
    sorted_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/fastq/100mM-rep1-Bot5_S2_L001_R1_001.fastq"
    unsorted_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/fastq/100mM-rep1-Unsorted_S1_L001_R1_001.fastq"
    library_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/library/Mouse_kinome_list_brie_updated.csv"
    output_file = "output.csv"

    result = perform_analysis(sorted_file, unsorted_file, library_file, output_file)

    print("All done!")
    end = time.time()
    duration = (end - start) # duration is second
    print("Took %r seconds." % duration)
