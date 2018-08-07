import sys
import numpy as np
import csv
import os
import copy
import pandas as pd
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as compTool
from collections import OrderedDict
import time
import math
import threading
import subprocess
from pdb import set_trace as debugger
import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')

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

def collapse_to_gene_level(data):
    """Sums up top sorted counts, bot sorted counts, and unsorted counts based on Target Gene Symbol
    All other stats should be the same across a given Target Gene Symbol, so the 'first' value is preserved
    as opposed to a sum (see aggregation_rule). sgRNA Target Sequences are all included for
    a given gene, and they are separated by commas.

    Do you need the information alongside the frequencies like this? It takes a significant amount of time to
    aggregate the descriptor data, and you might be able to just reference a separate dataframe to get the
    annotation information.
    """

    gene_level = pd.pivot_table(data, index="Target Gene Symbol", columns=["replicate","condition"], values="frequency", aggfunc="sum")
    info_columns = list(set(data.columns.tolist()) - set(["condition","frequency","replicate"]))
    get_info = data[info_columns]
    joiner = lambda x: ", ".join(x)
    aggregation_rule = {c : joiner if (c.startswith("sgRNA Target")) else "first" for c in get_info.columns}

    get_info = get_info.groupby("Target Gene Symbol").agg(aggregation_rule)
    get_info.head()
    gene_level[get_info.columns] = get_info[get_info.columns]
    gene_level.head()

    return gene_level

def create_mageck_input_files(raw_data, conditions):
    """ Mageck input format is tab-delimited file """
    """ sgRNA   \tGene  \tControl1 \tControl2 \tRep1 \tRep2   \n"""
    """ Note: separate inputs for separate sorts [top versus bottom] (for now)"""

    pivot_table = pd.pivot_table(raw_data, index=["sgRNA Target Sequence", "Target Gene Symbol"], \
                                columns=["condition","replicate"], values="frequency", aggfunc="first", \
                                fill_value=0)
    mageck_input = pivot_table[conditions]
    return mageck_input

# class parseThread(threading.Thread):
class analysisThread():
    def __init__(self, input_files, output_prefix, output_dir):
        self.config_analysis = {
                                "Mageck": {"Top_Sorted": True,
                                            "Bottom_Sorted": True},
                                "Fischer": False,
                                 "Ratio": False
                                }

        self.output = copy.deepcopy(self.config_analysis)

        self.control_file = input_files["Control"]
        self.library_file = input_files["Library"]

        self.output_prefix = output_prefix
        self.output_dir = os.path.join(output_dir, output_prefix)
        self.input_files = input_files
        self.count_status = 0
        self.status = "Thread Initialized"
        print("Initializing thread information...")
        # super().__init__()

    def run(self):

        logging.info("Begun run sequence.")

        # Parse the input files to generate conglomerate dataframes
        integrated_output_guide_level, integrated_output_gene_level = self.getCountData(self.input_files)

        # Save count data before performing any sort of analyses
        guide_path = self.output_dir+"_guide_counts.csv"
        genes_path = self.output_dir+"_gene_counts.csv"
        integrated_output_guide_level.to_csv(guide_path)
        integrated_output_gene_level.to_csv(genes_path)

        if self.config_analysis["Mageck"]:
            for condition, perform_condition_analysis in self.config_analysis["Mageck"].items():
                if perform_condition_analysis == True:
                     self.output["Mageck"][condition] = self.perform_mageck_analysis(integrated_output_guide_level, condition)

        if self.config_analysis["Ratio"] == True:
            for_ratio_test = self.guides_result[['sgRNA Target Sequence', 'Target Gene Symbol', 'Bot Sorted Counts', 'Top Sorted Counts']].copy()
            ratio_result_df = self.ratio_test(for_ratio_test)

        return "Complete"

    def getCountData(self, input_files):
        logging.info("Parsing files...")
        self.status = "Parsing files..."
        # Create master guides_result dataframe
        integrated_output_guide_level = self.parse_input_files(input_files)

        # Collapse to gene level
        logging.info("Collapsing to gene level...")
        self.status = "Collapsing to Gene Level..."
        integrated_output_gene_level = collapse_to_gene_level(integrated_output_guide_level)

        return integrated_output_guide_level, integrated_output_gene_level

    def parse_input_files(self, input_files):

        # Extract fastq files from input_files dictionary
        fastq_files = input_files["Fastq"]

        # Create main dataframe with library information
        library_annotations = create_dataframe(self.library_file)

        integrated_output_df = pd.DataFrame([])


        # Iterate over fastq files for parsing
        os.mkdir(self.output_prefix)
        self.output_path = os.path.join(self.output_prefix, self.output_prefix)

        for replicate_name, replicate_value in fastq_files.items():
            logging.info("Parsing replicate %s" % replicate_name)

            for fastq_description, fastq_file in replicate_value.items():
                logging.info("Parsing file %s" % fastq_description)

                # Parse the fastq file and convert it to a dataframe
                parsed_output_dict = parse_qfast(fastq_file, self.library_file)
                parsed_output_df = pd.DataFrame.from_dict(parsed_output_dict, orient="index")
                parsed_output_df = parsed_output_df.reset_index(drop=False)

                parsed_output_df = parsed_output_df.rename(columns={"index":"sgRNA Target Sequence",0:"frequency"})

                replicate_column = pd.Series(np.zeros(len(parsed_output_df)))
                replicate_column[:] = replicate_name
                parsed_output_df["replicate"] = replicate_column
                condition_column = pd.Series(np.zeros(len(parsed_output_df)))
                condition_column[:] = fastq_description
                parsed_output_df["condition"] = condition_column
                annotated_parsed_output_df = pd.merge(parsed_output_df, library_annotations.copy(), how="left", on=["sgRNA Target Sequence"])

                # Append the dataframe
                integrated_output_df = integrated_output_df.append(annotated_parsed_output_df).reset_index(drop=True)

        return integrated_output_df

    def perform_mageck_analysis(self, input_dataframe, condition):
        # Create mageck input file
        logging.info("Beginning mageck analysis.")
        self.status = "Creating Mageck input file..."
        mageck_df = create_mageck_input_files(input_dataframe, ["Unsorted Population", "Top Sorted Population"])
        mageck_path = self.output_path+"_"+condition+"_mageck_input_file.txt"

        output_column_values = mageck_df.columns.get_level_values(1)+":"+mageck_df.columns.get_level_values(0)
        mageck_df.to_csv(mageck_path, sep="\t", index=True, header=output_column_values)

        # Perform mageck test, passing the file path location of the input file
        mageck_result_df = self.run_mageck_command(mageck_path, condition)

        # Sort by mageck p-value and then write out mageck output file
        mageck_result_df = mageck_result_df.sort_values(by=["-log(pos|p-value)"], ascending=False)

        mageck_output_path = self.output_path+"_"+condition+"_mageck_gene_results.csv"
        mageck_result_df.to_csv(mageck_output_path)

        return mageck_result_df

    def run_mageck_command(self, mageck_path, condition):
        command = "mageck test -k " + mageck_path + " -t 2,3 -c 0,1 -n " + \
                str(self.output_path)+"_mageck_" + condition + \
                " --sort-criteria pos --gene-lfc-method alphamean --control-sgrna "\
                 + str(self.control_file)

        self.status = "Running mageck..."

        os.system(command)

        mageck_prefix = self.output_path+"_mageck_"+condition+".gene_summary.txt"

        mageck_result_df = pd.read_csv(mageck_prefix, sep="\t")
        mageck_result_df["-log(pos|p-value)"] = mageck_result_df["pos|p-value"].apply(lambda x: -1*math.log(x, 10))

        return mageck_result_df

    def ratio_test(self, df):

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

        # Add more below

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

def testing():
    raw_data = pd.read_csv("output.csv")
    return me.perform_mageck_analysis(raw_data)

if __name__ == "__main__":

    unsorted_rep1_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/Fastq Files/Replicate Files/fastq/Genome-Pos-3T3-Unsorted_S1_L001_R1_001.fastq"
    top_sorted_rep1_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/Fastq Files/Replicate Files/fastq/Genome-Pos-3T3-SortedTop5_S3_L001_R1_001.fastq"
    bot_sorted_rep1_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/Fastq Files/Replicate Files/fastq/Genome-Pos-3T3-SortedBot10_S2_L001_R1_001.fastq"

    unsorted_rep2_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/Fastq Files/Replicate Files/fastq/Genome-Pos2-3T3-Unsorted_S1_L001_R1_001.fastq"
    top_sorted_rep2_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/Fastq Files/Replicate Files/fastq/Genome-Pos2-3T3-Top5_S3_L001_R1_001.fastq"
    bot_sorted_rep2_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/Fastq Files/Replicate Files/fastq/Genome-Pos2-3T3-Bot10_S2_L001_R1_001.fastq"

    library_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/screen_analyzer/tmp/data/library/Brie_crispr_library_with_controls_for_analysis_updated.csv"
    control_file = "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/screen_analyzer/tmp/Brie-library-controls.txt"

    input_files = {"Fastq":
                        {"Rep1":
                            {
                            "Unsorted Population" : unsorted_rep1_file,
                            "Top Sorted Population" : top_sorted_rep1_file,
                            "Bottom Sorted Population" : bot_sorted_rep1_file
                            },
                        "Rep2":
                            {
                            "Unsorted Population" : unsorted_rep2_file,
                            "Top Sorted Population" : top_sorted_rep2_file,
                            "Bottom Sorted Population" : bot_sorted_rep2_file
                            }
                        },
                    "Library": library_file,
                    "Control": control_file
    }

    output_file = "output.csv"

    me = analysisThread(input_files, "test1", "/Users/collinschlager/Documents/Rohatgi_Lab/Informatics/screen_analyzer")
