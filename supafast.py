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
from pdb import set_trace as bp

#idfile = sys.argv[1]

def complete_merge(df1, df2):
    column_set_1 = set(df1.columns)
    column_set_2 = set(df2.columns)
    intersection = list(column_set_1.intersection(column_set_2))
    return pd.merge(df1, df2, on=intersection)

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
    sortedTotMatches = int(dataframe["Top Sorted Counts"].sum())
    unsortedTotMatches = int(dataframe["Unsorted Counts"].sum())
    normalized_sorted = (dataframe["Top Sorted Counts"]+1)/sortedTotMatches
    normalized_unsorted = (dataframe["Unsorted Counts"]+1)/unsortedTotMatches
    foldchange = normalized_sorted/normalized_unsorted
    dataframe["LFC"] = foldchange.apply(log2)
    print("Done.")
    return dataframe

def fischer_function(row, sortedTot = None, unsortedTot = None):
    oddsratio, pValue = stats.fisher_exact([[row["Unsorted Counts"], unsortedTot], [row["Top Sorted Counts"], sortedTot]]);
    return pValue

def compute_fischer(dataframe):
    print("Computing Fischer...")
    sortedTotMatches = int(dataframe["Top Sorted Counts"].sum())
    unsortedTotMatches = int(dataframe["Unsorted Counts"].sum())
    dataframe["Fischer P-Values"] = dataframe.apply(fischer_function, axis=1, sortedTot=sortedTotMatches, unsortedTot=unsortedTotMatches)
    adjEnrichPVals = compTool.multipletests(list(dataframe["Fischer P-Values"]), method='fdr_bh', is_sorted=False)[1]
    assert (len(adjEnrichPVals) == len(dataframe))
    dataframe["FDR-Corrected P-Values"] = adjEnrichPVals
    dataframe["-log(FDR-Corrected P-Values)"] = dataframe["FDR-Corrected P-Values"].apply(lambda x: -1*math.log(x,10))
    print("Done.")
    return dataframe

def perform_analysis(sorted_fastq, unsorted_fastq, library_file):
    print("Beginning parsing...")
    library_df = create_dataframe(library_file)
    sorted_counts_dict = parse_qfast(sorted_fastq, library_file)
    add_df_column(library_df, sorted_counts_dict, "Top Sorted Counts")
    unsorted_counts_dict = parse_qfast(unsorted_fastq, library_file)
    add_df_column(library_df, unsorted_counts_dict, "Unsorted Counts")
    return library_df

def collapse_to_gene_level(dataframe):
    """Sums up top sorted counts, bot sorted counts, and unsorted counts based on Target Gene Symbol
    All other stats should be the same across a given Target Gene Symbol, so the 'first' value is preserved
    as opposed to a sum (see aggregation_rule)"""

    print("Collapsing to Gene Level...")
    aggregation_rule = {c : "sum" if (c == "Top Sorted Counts" or c == "Bot Sorted Counts" or c == "Unsorted Counts") else "first" for c in dataframe.columns}
    aggregation_rule.pop("Target Gene Symbol")
    print(aggregation_rule)
    dataframe.sort_values(by=["sgRNA Target Sequence"], inplace=True)
    gene_level_df = dataframe.groupby("Target Gene Symbol", as_index=False).agg(aggregation_rule)
    gene_level_df.to_csv("gene_level_output.csv")
    print("Done.")
    return gene_level_df

def create_mageck_input_file(dataframe):
    """ Mageck input format is tab-delimited file """
    """ sgRNA   \tGene  \tControl1 (unsorted)  \tSorted1   \n"""
    mageck_data_frame = dataframe[["sgRNA Target Sequence", "Target Gene Symbol", "Unsorted Counts", "Top Sorted Counts"]]
    return mageck_data_frame

class parseThread(threading.Thread):
    def __init__(self, top_sorted_file, bot_sorted_file, unsorted_file, output_prefix, guides_file, output_dir, control_file="Brie_Kinome_controls.txt"):
        self.config_analysis = {"Mageck": True, "Fischer": False, "Ratio": False} # default values; set in controller.py
        self.top_sorted_file = top_sorted_file
        self.bot_sorted_file = bot_sorted_file
        self.unsorted_file = unsorted_file
        self.output_prefix = output_prefix
        self.output_dir = output_dir
        self.output_dir = os.path.join(output_dir, output_prefix)
        self.guides_file = guides_file
        self.control_file = control_file
        self.count_status = 0
        self.output = 0
        self.status = "Thread Initialized"
        print("Initializing thread information...")
        print(self.output_prefix)
        super().__init__()

    def run(self, lfc=False, fischer=False, mageck=False):
        self.status = "Parsing files..."
        self.parse_files()

        self.status = "Collapsing to Gene Level..."
        self.genes_result = collapse_to_gene_level(self.guides_result)

        if self.config_analysis["Mageck"] == True:

            self.mageck_guides_result = self.guides_result.copy()
            self.mageck_genes_result = self.genes_result.copy()

            # Create mageck input file
            self.status = "Creating Mageck input file..."
            mageck_data_frame = create_mageck_input_file(self.mageck_guides_result.sort_values(by=["sgRNA Target Sequence"]))
            self.mageck_path = self.output_prefix+"_mageck_input_file.txt"
            mageck_data_frame.to_csv(self.mageck_path, sep="\t", index=False)

            # Perform mageck Test
            self.perform_mageck_test()
            # Ensure that the mageck tests are sorted the same as the genes_result
            # Is this supposed to be "genes"
            # Debug something here regarding index needing to be target gene symbol
            self.mageck_genes_result.set_index('Target Gene Symbol', drop=False, inplace=True)
            self.mageck_result_df.set_index('id', drop=False, inplace=True)
            #self.genes_result.drop("Non-Targeting-Control", inplace=True) #temporary!
            self.mageck_genes_result.sort_values(by=["Target Gene Symbol"], ascending=False, inplace=True)
            self.mageck_result_df.sort_values(by=["id"], ascending=False, inplace=True)
            self.mageck_genes_result["pos|lfc"] = self.mageck_result_df["pos|lfc"]
            self.mageck_genes_result["pos|p-value"] = self.mageck_result_df["pos|p-value"]
            self.mageck_genes_result["-log(pos|p-value)"] = self.mageck_result_df["-log(pos|p-value)"]

            # Sort by mageck p-value and then write out mageck output file
            self.mageck_genes_result.sort_values(by=["-log(pos|p-value)"], ascending=False, inplace=True)
            genes_path = self.output_prefix+"_mageck_gene.csv"
            print(genes_path)
            self.mageck_genes_result.to_csv(genes_path)

            # Append new analysis to master file
            self.genes_result = complete_merge(self.genes_result, self.mageck_genes_result)



        if self.config_analysis["Ratio"] == True:
            print(self.guides_result.columns)
            for_ratio_test = self.guides_result[['sgRNA Target Sequence', 'Target Gene Symbol', 'Bot Sorted Counts', 'Top Sorted Counts']].copy()
            self.ratio_test(for_ratio_test)



        self.genes_result.sort_values(by=["Target Gene Symbol"], inplace=True)
        genes_path = self.output_prefix+"_gene_enrichment_calculation.csv"
        print(genes_path)
        self.genes_result.to_csv(genes_path)
        self.guides_result.sort_values(by=["sgRNA Target Sequence"], ascending=False, inplace=True)
        guide_path = self.output_prefix+"_guide_enrichment_calculation.csv"
        self.guides_result.to_csv(guide_path)


        self.status = "Analysis Complete"
        self.output = self.guides_result, self.genes_result
        return self.guides_result, self.genes_result

    def parse_files(self):
        self.guides_result = create_dataframe(self.guides_file)

        if self.top_sorted_file:
            self.status = "Parsing top sorted file..."
            top_sorted_counts_dict = parse_qfast(self.top_sorted_file, self.guides_file)
            add_df_column(self.guides_result, top_sorted_counts_dict, "Top Sorted Counts")

        if self.bot_sorted_file:
            self.status = "Parsing top sorted file..."
            bot_sorted_counts_dict = parse_qfast(self.bot_sorted_file, self.guides_file)
            add_df_column(self.guides_result, bot_sorted_counts_dict, "Bot Sorted Counts")

        if self.unsorted_file:
            self.status = "Parsing unsorted file..."
            unsorted_counts_dict = parse_qfast(self.unsorted_file, self.guides_file)
            add_df_column(self.guides_result, unsorted_counts_dict, "Unsorted Counts")

    def perform_mageck_test(self):
        command = "mageck test -k " + self.mageck_path + " -t 1 -c 0 -n " + str(self.output_dir)+"_mageck" + " --sort-criteria pos --gene-lfc-method alphamean --control-sgrna " + str(self.control_file)
        self.status = "Running mageck..."
        #os.chdir(self.output_dir)
        os.system(command)
        time.sleep(5); #wait 5 seconds for mageck to be done?
        mageck_prefix = self.output_prefix+"_mageck.gene_summary.txt"
        self.mageck_result_df = pd.read_csv(mageck_prefix, sep="\t")
        self.mageck_result_df["-log(pos|p-value)"] = self.mageck_result_df["pos|p-value"].apply(lambda x: -1*math.log(x, 10))
        # self.mageck_result_df["pos|lfc"] = self.mageck_result_df["pos|lfc"]
        # p = subprocess.Popen(command, shell=True, cwd=self.output_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # for line in p.stdout.readlines():
        #    self.status = str(line)
        #    print("MAGECK: ", str(line))
        # p.wait()

    def ratio_test(self, df):
        #df = df.rename(columns={'Sorted Counts': 'bot.Rep1.Reads', 'Unsorted Counts': 'top.Rep1.Reads'})

        # Add 1 to reads
        df['Bot Sorted Counts+1'] = df['Bot Sorted Counts']+1
        df['Top Sorted Counts+1'] = df['Top Sorted Counts']+1

        forQuant = df[['Bot Sorted Counts+1','Top Sorted Counts+1']].copy()
        quantNorm = self.quantileNormalize(forQuant)
        # rank_mean = forQuant.stack().groupby(forQuant.rank(method='first').stack().astype(int)).mean()
        # quantNorm = forQuant.rank(method='min').stack().astype(int).map(rank_mean).unstack()


        df['bot.Rep1.Quant'] = quantNorm.loc[:,'Top Sorted Counts+1']
        df['top.Rep1.Quant'] = quantNorm.loc[:,'Bot Sorted Counts+1']

        ratios = df[['sgRNA Target Sequence', 'Target Gene Symbol', 'bot.Rep1.Quant']].copy()
        ratios.rename(columns = {'bot.Rep1.Quant':'Bot5.Geom.Mean'}, inplace=True)
        ratios['Top5.Geom.Mean'] = df['top.Rep1.Quant']
        ratios['ratio'] = ratios['Bot5.Geom.Mean'] / ratios['Top5.Geom.Mean']
        ratios['log_ratio'] = ratios['ratio'].apply(np.log2)

        ratios['mean_abundance'] = forQuant.apply(self.geomMean, axis=1)
        ratios['log_MA'] = ratios['mean_abundance'].apply(np.log2)

        ratios.sort_values(by=['log_MA'], ascending=True, inplace=True)
        ratios['ZScore'] = np.zeros(len(ratios['log_MA']))

        bins=3

        guidesPerBin = math.ceil(len(ratios)/bins)
        ratios.reset_index(inplace=True)

        for i in range(bins):
            start = (i)*guidesPerBin
            end = start+guidesPerBin
            if end > len(ratios):
                end = len(ratios)
            guides = ratios[start:end]
            guides['ZScore'] = stats.zscore(guides.loc[:,'log_ratio'])
            # ratios['ZScore'][start:end] = guides.loc[:,'ZScore']
            ratios.loc[start:end,'ZScore'] = guides.loc[:,'ZScore']

        print("Collapsing to Gene Level...")
        ratios.sort_values(by=['sgRNA Target Sequence'], inplace=True)
        aggregation_rule = {c : "mean" if (c == "ZScore" or c == "log_ratio") else "first" for c in ratios.columns}
        aggregation_rule.pop("Target Gene Symbol")
        gene_level_df = ratios.groupby("Target Gene Symbol", as_index=False).agg(aggregation_rule)
        gene_level_df.sort_values(by=['ZScore'], ascending=False, inplace=True)
        gene_level_df.drop(columns=["index"], inplace=True)
        gene_level_df.reset_index(drop=True, inplace=True)
        ratios_path = self.output_prefix+"_ratio_test.csv"
        gene_level_df.to_csv(ratios_path)

        self.genes_result = complete_merge(self.genes_result, gene_level_df)

        return ratios

    def geomMean(self, x):
        product = np.prod(x)
        numElem = len(x)
        result = product**(1/numElem)
        return result

    def quantileNormalize(self, df_input):
        df = df_input.copy()
        #compute rank
        dic = {}
        for col in df:
            dic.update({col : sorted(df[col])})
        sorted_df = pd.DataFrame(dic)
        rank = sorted_df.mean(axis = 1).tolist()
        #sort
        for col in df:
            t = np.searchsorted(np.sort(df[col]), df[col])
            df[col] = [rank[i] for i in t]
        return df

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
