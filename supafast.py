import sys
import numpy as np
import csv
import os
import pandas as pd
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as compTool
from collections import OrderedDict
import time
import math
import threading
import subprocess

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
    print("Parsing %r" % qfast_file)
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
    print("Computing LFC...")
    log2 = lambda x: math.log2(x)
    sortedTotMatches = int(dataframe["Sorted Counts"].sum())
    unsortedTotMatches = int(dataframe["Unsorted Counts"].sum())
    normalized_sorted = (dataframe["Sorted Counts"]+1)/sortedTotMatches
    normalized_unsorted = (dataframe["Unsorted Counts"]+1)/unsortedTotMatches
    foldchange = normalized_sorted/normalized_unsorted
    dataframe["LFC"] = foldchange.apply(log2)
    print("Done.")
    return dataframe

def fischer_function(row, sortedTot = None, unsortedTot = None):
    oddsratio, pValue = stats.fisher_exact([[row["Unsorted Counts"], unsortedTot], [row["Sorted Counts"], sortedTot]]);
    return pValue

def compute_fischer(dataframe):
    print("Computing Fischer...")
    sortedTotMatches = int(dataframe["Sorted Counts"].sum())
    unsortedTotMatches = int(dataframe["Unsorted Counts"].sum())
    dataframe["Fischer P-Values"] = dataframe.apply(fischer_function, axis=1, sortedTot=sortedTotMatches, unsortedTot=unsortedTotMatches)
    adjEnrichPVals = compTool.multipletests(list(dataframe["Fischer P-Values"]), method='fdr_bh', is_sorted=False)[1]
    assert (len(adjEnrichPVals) == len(dataframe))
    dataframe["FDR-Corrected P-Values"] = adjEnrichPVals
    print("Done.")
    return dataframe

def perform_analysis(sorted_fastq, unsorted_fastq, library_file):
    print("Beginning parsing...")
    library_df = create_dataframe(library_file)
    sorted_counts_dict = parse_qfast(sorted_fastq, library_file)
    add_df_column(library_df, sorted_counts_dict, "Sorted Counts")
    unsorted_counts_dict = parse_qfast(unsorted_fastq, library_file)
    add_df_column(library_df, unsorted_counts_dict, "Unsorted Counts")
    return library_df

def collapse_to_gene_level(dataframe):
    print("Collapsing to Gene Level...")
    aggregation_rule = {c : "sum" if (c == "Sorted Counts" or c == "Unsorted Counts") else "first" for c in dataframe.columns}
    aggregation_rule.pop("Target Gene Symbol")
    print(aggregation_rule)
    gene_level_df = dataframe.groupby("Target Gene Symbol", as_index=False).agg(aggregation_rule)
    gene_level_df.to_csv("gene_level_output.csv")
    print("Done.")
    return gene_level_df

def create_mageck_input_file(dataframe):
    """ Mageck input format is tab-delimited file """
    """ sgRNA   \tGene  \tControl1 (unsorted)  \tSorted1   \n"""
    mageck_data_frame = dataframe[["sgRNA Target Sequence", "Target Gene Symbol", "Unsorted Counts", "Sorted Counts"]]
    return mageck_data_frame



class parseThread(threading.Thread):
    def __init__(self, sorted_file, unsorted_file, output_prefix, guides_file, output_dir, control_file="Brie_Kinome_controls.txt"):
        self.sorted_file = sorted_file
        self.unsorted_file = unsorted_file
        self.output_prefix = output_prefix
        self.output_dir = output_dir
        self.guides_file = guides_file
        self.control_file = control_file
        self.count_status = 0
        self.output = 0
        self.status = "Thread Initialized"
        print("Initializing thread information...")
        print(self.output_prefix)
        super().__init__()

    def run(self):
        self.status = "Parsing files..."
        self.parse_files()

        self.status = "Collapsing to Gene Level..."
        self.genes_result = collapse_to_gene_level(self.guides_result)

        self.status = "Computing LFC..."
        compute_lfc(self.guides_result)
        compute_lfc(self.genes_result)

        self.status = "Computing fischer for guides..."
        #compute_fischer(self.guides_result)
        self.status = "Computing fischer for genes..."
        #compute_fischer(self.genes_result)

        # Create mageck input file
        self.status = "Creating Mageck input file..."
        mageck_data_frame = create_mageck_input_file(self.guides_result.sort_values(by=["sgRNA Target Sequence"]))
        self.mageck_path = self.output_prefix+"_mageck_input_file.txt"
        mageck_data_frame.to_csv(self.mageck_path, sep="\t", index=False)

        # Perform mageck Test
        self.perform_mageck_test()

        # self.guides_result.sort_values(by=["FDR-Corrected P-Values"], ascending=True, inplace=True)
        self.guides_result.sort_values(by=["LFC"], ascending=False, inplace=True)
        guide_path = self.output_prefix+"_guide_enrichment_calculation.csv"
        print(guide_path)
        self.guides_result.to_csv(guide_path)

        self.genes_result.sort_values(by=["LFC"], ascending=False, inplace=True)
        genes_path = self.output_prefix+"_gene_enrichment_calculation.csv"
        print(genes_path)
        self.genes_result.to_csv(genes_path)


        self.status = "Analysis Complete"
        self.output = self.guides_result, self.genes_result
        return self.guides_result, self.genes_result

    def parse_files(self):
        self.guides_result = create_dataframe(self.guides_file)
        self.status = "Parsing sorted file..."
        sorted_counts_dict = parse_qfast(self.sorted_file, self.guides_file)
        add_df_column(self.guides_result, sorted_counts_dict, "Sorted Counts")

        self.status = "Parsing unsorted file..."
        unsorted_counts_dict = parse_qfast(self.unsorted_file, self.guides_file)
        add_df_column(self.guides_result, unsorted_counts_dict, "Unsorted Counts")

    def perform_mageck_test(self):
        command = "mageck test -k " + self.mageck_path + " -t 1 -c 0 -n " + str(self.output_dir)+"_mageck" + " --sort-criteria pos --gene-lfc-method alphamean --control-sgrna " + str(self.control_file)
        self.status = "Running mageck..."
        #os.chdir(self.output_dir)
        os.system(command)
        # p = subprocess.Popen(command, shell=True, cwd=self.output_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # for line in p.stdout.readlines():
        #    self.status = str(line)
        #    print("MAGECK: ", str(line))
        # p.wait()

def main():
    start = time.time()
    sorted_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/fastq/100mM-rep1-Bot5_S2_L001_R1_001.fastq"
    unsorted_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/fastq/100mM-rep1-Unsorted_S1_L001_R1_001.fastq"
    library_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/library/Mouse_kinome_list_brie_updated.csv"
    output_file = "output.csv"

    guides_result = perform_analysis(sorted_file, unsorted_file, library_file, output_file)
    genes_result = collapse_to_gene_level(guides_result)

    compute_lfc(guides_result)
    compute_lfc(genes_result)
    compute_fischer(guides_result)
    compute_fischer(genes_result)

    guides_result.sort_values(by=["FDR-Corrected P-Values"], ascending=True, inplace=True)
    guides_result.to_csv(output_file)
    genes_result.sort_values(by=["FDR-Corrected P-Values"], ascending=True, inplace=True)
    genes_result.to_csv("gene_level_"+output_file)


    print("All done!")
    end = time.time()
    duration = (end - start) # duration is second
    print("Took %r seconds." % duration)
    return guides_result, genes_result, duration

if __name__ == "__main__":
    guides, genes, time_elapsed = main()
