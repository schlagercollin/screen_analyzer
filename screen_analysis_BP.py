#Script to analyze sequencing results post-selection of cells. This script takes in sequencing reads from the sorted and unsorted (control) populations, counts the number of reads per guide, and outputs some statistics about the number of matched guides and the distribution to a few output files. It also outputs a read count per gene to another csv file. It also outputs a files that contain the perfectly-matched and non-matched sequences for further analysis. If a fastq file is provided for both sorted and unsorted populations, a statistical analysis is conducted to look for genes with an enrichment of guides/reads.

from Bio import SeqIO
from collections import OrderedDict
import numpy as np
import sys, math, os, csv
import xlwt
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as compTool
import argparse
import threading
import copy
import subprocess
#import annotate_functions

#this method creates the dictionaries that map guide sequence to number to reads and guide sequence to gene name
def createDictionaries(input_file):
    countDict = {}; #dictionary to keep track of the read counts for each guide
    guideGeneDict = {}; #dictionary to keep track of which gene each guide belongs to
    geneCountDict = {}; #dictionary to keep track of the read counts for a gene
    geneToGuideDict = {}; #dictionary to keep track of the sgRNAs for a gene
    geneToDescriptionDict = {}; #dictionary to keep track of the description and summary for a gene
    # open library of guide sequences (.csv file) and create dictionary of read counts for each guide where count is initialized to 0
    #try:
    with open(input_file, mode='rU') as infile:
        reader = csv.DictReader(infile);
        #reader.next();
        for row in reader:
            guideSeq = row["sgRNA Target Sequence"];
            description = row["Description"];
            summary = row["Summary"];
            countDict[guideSeq] = 0 #initialize sequence read count in countDict
            gene = row["Target Gene Symbol"];
            guideGeneDict[guideSeq] = gene; #create guide->gene mapping in guideGeneDict
            geneCountDict[gene] = 0; #initialize read count for each gene
            #create gene->list of sgRNAs in geneToGuideDict
            if geneToGuideDict.get(gene) == None:
                geneToGuideDict[gene] = [guideSeq];
            else:
                (geneToGuideDict[gene]).append(guideSeq);
            if geneToDescriptionDict.get(gene) == None:
                geneToDescriptionDict[gene] = {}
                geneToDescriptionDict[gene]["Description"] = description
                geneToDescriptionDict[gene]["Summary"] = summary
            else:
                pass
    #except Exception as e:
    #    print('could not open', input_file);
    #    print(e)
    return countDict, guideGeneDict, geneCountDict, geneToGuideDict, geneToDescriptionDict;


GUIDE_START = 0 #start index of guide sequence
GUIDE_END = 20 #end index of guide sequence
BEFORE_GUIDE = "CACCG" #identifies sequence after guide

#print reads counts per guide and per gene for each population
def geneReadCounts(outputFileName, geneCountDict, geneToGuideDict, countDict, perfect_matches):
    #create ordered dictionary with genes and respective counts based on read counts and output as a csv file
    dict_sorted = OrderedDict(sorted(geneCountDict.items(), key=lambda t: t[1]));
    #open Excel workbook
    workbook = xlwt.Workbook();
    sheet = workbook.add_sheet("Reads Per Gene");
    #create headers
    sheet.write(0,0,"Gene Name");
    sheet.write(0,1,"Read Count");
    sheet.write(0,2,"Fraction of Total Reads");
    #input data
    runningSum = 0.0;
    rowNum = 1;
    for gene in dict_sorted:
        count = dict_sorted[gene];
        runningSum += count;
        fraction = runningSum/perfect_matches;
        sheet.write(rowNum, 0, gene);
        sheet.write(rowNum, 1, count);
        sheet.write(rowNum, 2, fraction);
        rowNum += 1;
    #create second sheet to hold data on preservation of guides for each gene
    sheet = workbook.add_sheet("Guides Per Gene");
    #create headers
    sheet.write(0,0,"Gene Name");
    sheet.write(0,1,"Number of Guides");
    sheet.write(0,2,"Number of Guides Preserved");
    sheet.write(0,3,"Fraction of Guides Preserved");
    rowNum = 1;
    geneGuidesPreserved = 0; #track the number of genes for which all guides are preserved in library
    for gene in geneToGuideDict:
        guides = geneToGuideDict[gene]; #get list of guides
        guidesPreserved = 0; # holds count of guides preserved for this gene
        for guide in guides: #count the number of guides preserved
            if countDict[guide] > 0:
                guidesPreserved += 1;
        if guidesPreserved == len(guides): #increment total count of genes with all guides preserved if all guides are preserved for current gene
            geneGuidesPreserved += 1;
        fractionPreserved = float(guidesPreserved)/len(guides)*100;
        sheet.write(rowNum, 0, gene);
        sheet.write(rowNum, 1, len(guides));
        sheet.write(rowNum, 2, guidesPreserved);
        sheet.write(rowNum, 3, fractionPreserved);
        rowNum += 1;
    #calculate fraction of genes for which all guides were preserved
    numGenes = len(geneToGuideDict.keys());
    sheet.write(0,5,"Percentage of Genes with all Guides:");
    sheet.write(0,6,(float(geneGuidesPreserved)/numGenes*100));
    #save workbook
    workbook.save(outputFileName+"_gene_count.xls");
    print("Percentage of Genes with all Guides:", (float(geneGuidesPreserved)/numGenes))

#class to hold information about gene
class Gene:
    def __init__(self, name, unsortedReads, unsortedTotal, sortedReads, sortedTotal, fold_change, summary, description):
        self.name = name; #keeps track of gene name
        self.description = description;
        self.summary = summary;
        self.unsortedReads = unsortedReads;
        self.unsortedTotal = unsortedTotal;
        #sorted stats
        self.sortedReads = sortedReads;
        self.sortedTotal = sortedTotal;
        #p-values
        self.enrichPVal = 0;
        self.adjEnrichPVal = 0;
        #fold change
        self.foldChange = fold_change;

#class to hold information about each guide
class Guide:
    def __init__(self, sequence, gene, unsortedReads, unsortedTotal, sortedReads, sortedTotal):
        self.sequence = sequence; #keeps track of the guide sequence
        self.gene = gene; #keeps track of the gene the guide sequence targets
        self.unsortedReads = unsortedReads;
        self.unsortedTotal = unsortedTotal;
        #sorted stats
        self.sortedReads = sortedReads;
        self.sortedTotal = sortedTotal;
        #p-values
        self.enrichPVal = 0;
        self.adjEnrichPVal = 0;
        #store fold change
        logFoldChange = 0;

def printAnalysesInputs(guideStatsDict, output_file):
    guidesList = list(guideStatsDict.values());
    guidesList.sort(key = lambda guide: guide.sequence); #sort on sequence so replicates can be easily combined later
    #setup MaGeCK file
    mageckoutputName = str(output_file+"_mageck_input.txt");
    mageckFile = open(mageckoutputName, 'w');
    mageckFile.write("sgRNA\tGene\tControl1\tSorted1\n");
    count = 1; #need to track guide count for RIGER file
    for guide in guidesList:
        mageckFile.write(guide.sequence+"\t"+guide.gene+"\t"+str(guide.unsortedReads)+"\t"+str(guide.sortedReads)+"\n");
        count += 1; #update count

    mageckFile.close();
    print("Printed MAGeCK input files");

#function to be called to analyze a single replicate for a screen
def run_screen_analysis(unsorted_fastq, sorted_fastq, guides_file, output_file):
    #necessary dictionaries to keep track of guide/gene info
    countDict, guideGeneDict, geneCountDict, geneToGuideDict, geneToDescriptionDict = createDictionaries(guides_file);
    print("created required dictionaries");

    #count reads per guides and per gene for unsorted population if the fastq file is provided
    unsortedGuideCountDict = dict(countDict);
    unsortedGeneCountDict = dict(geneCountDict);
    unsortedTotMatches = 0;
    print("Counting reads for unsorted population");
    unsortedOutputName = str(output_file+"_unsorted");
    perfect_matches = count_spacers(unsorted_fastq, unsortedOutputName, unsortedGuideCountDict, guideGeneDict, unsortedGeneCountDict);
    unsortedTotMatches = perfect_matches;
    #print gene related info for unsorted population
    geneReadCounts(unsortedOutputName, unsortedGeneCountDict, geneToGuideDict, unsortedGuideCountDict, perfect_matches);
    print("Done counting reads for unsorted population");

    #count reads per guides and per gene for sorted population if the fastq file is provided
    sortedGuideCountDict = dict(countDict);
    sortedGeneCountDict = dict(geneCountDict);
    sortedTotMatches = 0;
    print("Counting reads for sorted population");
    sortedOutputName = str(output_file+"_sorted");
    perfect_matches = count_spacers(sorted_fastq, sortedOutputName, sortedGuideCountDict, guideGeneDict, sortedGeneCountDict);
    sortedTotMatches = perfect_matches;
    #print gene related info for sorted population
    geneReadCounts(sortedOutputName, sortedGeneCountDict, geneToGuideDict, sortedGuideCountDict, perfect_matches);
    print("Done counting reads for sorted population");

    #perform enrichment analysis for genes if both unsorted and sorted files are provided
    print("Performing enrichment analysis");
    #look for genes enriched for reads
    geneStatsDict = calcGeneEnrich(unsortedGeneCountDict, unsortedTotMatches, sortedGeneCountDict, sortedTotMatches);
    printGeneEnrichStats(geneStatsDict, output_file);
    #look for guides enriched for reads
    guideStatsDict = calcGuideEnrich(unsortedGuideCountDict, unsortedTotMatches, sortedGuideCountDict, sortedTotMatches, guideGeneDict);
    printGuideEnrichStats(guideStatsDict, output_file);
    #print MaGeCK input files
    printAnalysesInputs(guideStatsDict, output_file);
    #perform mageck analysis
    # mageckFile = str(output_file+"_mageck_input.txt");
    # os.system("cp /Volumes/G-DRIVE\ slim/CRISPR\ Screens/Brie-library-controls.txt ."); #need to copy in control guides file
    # command = "mageck test -k " + mageckFile + " -t 1 -c 0 -n " + str(output_file+"_mageck" + " --sort-criteria pos --gene-lfc-method alphamean --control-sgrna Brie-library-controls.txt");
    # os.system(command);
    # os.system("rm Brie-library-controls.txt"); #rm control guides file
    # #annotate files
    # print "Annotating Guides file and MaGeCK results";
    # geneFuncsDict = annotate_functions.readFunctions();
    # #annoate guide enrichment file
    # guidesEnrichFile = output_file+"_guide_enrichment_calculation.csv";
    # annotate_functions.annotateGuidesFile(guidesEnrichFile, geneFuncsDict);
    # #annotate mageck results file
    # mageckResults = output_file+"_mageck.gene_summary.txt";
    # annotate_functions.annotateMageckFile(mageckResults, geneFuncsDict);

class parseThread(threading.Thread):
    def __init__(self, sorted, unsorted, output, guides, control):
        self.sorted = sorted
        self.unsorted = unsorted
        self.output = output
        self.guides = guides
        self.control = control
        self.count_status = 0
        self.status = "Thread Initialized"
        print("Initializing thread information...")
        super().__init__()

    def run(self):
        self.status = "Beginning analysis..."
        print("Running thread...")
        #necessary dictionaries to keep track of guide/gene info
        countDict, guideGeneDict, geneCountDict, geneToGuideDict, geneToDescriptionDict = createDictionaries(self.guides);
        print("created required dictionaries");
        self.status = "Created dictionaries..."
        #count reads per guides and per gene for unsorted population if the fastq file is provided
        #print("countDict = ", countDict)
        unsortedGuideCountDict = copy.deepcopy(countDict);
        unsortedGeneCountDict = copy.deepcopy(geneCountDict);
        unsortedTotMatches = 0;
        if self.unsorted != None:
            print("Counting reads for unsorted population...");
            self.status = "Counting reads for unsorted population"
            unsortedOutputName = str(self.output+"_unsorted");
            perfect_matches = self.count_spacers(self.unsorted, unsortedOutputName, unsortedGuideCountDict, guideGeneDict, unsortedGeneCountDict);
            unsortedTotMatches = perfect_matches;
            #print gene related info for unsorted population
            geneReadCounts(unsortedOutputName, unsortedGeneCountDict, geneToGuideDict, unsortedGuideCountDict, perfect_matches);
            print("Done counting reads for unsorted population");
            self.status = "Task completed."
        #count reads per guides and per gene for sorted population if the fastq file is provided
        sortedGuideCountDict = copy.deepcopy(countDict);
        sortedGeneCountDict = copy.deepcopy(geneCountDict);
        sortedTotMatches = 0;
        if self.sorted != None:
            print("Counting reads for sorted population");
            self.status = "Counting reads for sorted population"
            sortedOutputName = str(self.output+"_sorted");
            perfect_matches = self.count_spacers(self.sorted, sortedOutputName, sortedGuideCountDict, guideGeneDict, sortedGeneCountDict);
            sortedTotMatches = perfect_matches;
            #print gene related info for sorted population
            geneReadCounts(sortedOutputName, sortedGeneCountDict, geneToGuideDict, sortedGuideCountDict, perfect_matches);
            print("Done counting reads for sorted population");
            self.status = "Task completed."

        #perform enrichment analysis for genes if both unsorted and sorted files are provided
        if self.sorted != None and self.unsorted != None:
            print("Performing enrichment analysis");
            self.status = "Performing enrichment analysis on genes..."
            #look for genes enriched for reads
            geneStatsDict = self.calcGeneEnrich(unsortedGeneCountDict, unsortedTotMatches, sortedGeneCountDict, sortedTotMatches, geneToDescriptionDict);
            self.printGeneEnrichStats(geneStatsDict, self.output);
            #look for guides enriched for reads
            self.status = "Performing enrichment analysis on guides..."
            guideStatsDict = self.calcGuideEnrich(unsortedGuideCountDict, unsortedTotMatches, sortedGuideCountDict, sortedTotMatches, guideGeneDict);
            self.printGuideEnrichStats(guideStatsDict, self.output);

        self.printAnalysesInputs(guideStatsDict, self.output);
        mageckFile = str(self.output+"_mageck_input.txt");
        command = "mageck test -k " + mageckFile + " -t 1 -c 0 -n " + str(self.output)+"_mageck" + " --sort-criteria pos --gene-lfc-method alphamean --control-sgrna " + str(self.control)
        print("Running command!")
        #command = command.split(" ")
        #print(command)
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            self.status = line
            print(line)
        p.wait()
        #runcommand(command)

        print("DONE!")
        self.status = "Analysis Complete"
        self.myResult = self.output

    #matches the reads to the guides and prints statistics regarding the library. Returns the number of perfectly matched reads
    def count_spacers(self, fastq_file, output_file, countDict, guideGeneDict, geneCountDict):
        """
        creates a dictionary with guide counts from fastq_file, writes to output_file
        fastq_file: forward read fastq file
        output_file: csv file to write guide dictionary to
        dictionary: guide sequence as key, guide count as entry
        """
        print("CountDict = ", countDict)
        self.count_status = 0 #total number of reads processed
        perfect_matches = 0 # guides with perfect match to library
        non_perfect_matches = 0 #number of guides without a perfect match to the library
        perfMatchList = []; #array to hold sequences that perfectly match to library
        nonPerfList = []; #array to hold sequences that do not perfectly match to the library

        # open fastq file
        try:
            handle = open(fastq_file, "rU")
        except:
            print("could not find fastq file")
            return

        # process reads in fastq file
        readiter = SeqIO.parse(handle, "fastq")
        for record in readiter: #contains the seq and Qscore etc.
            sys.stdout.write("\r%d" % self.count_status)
            sys.stdout.flush()
            self.count_status += 1
            read_sequence = str.upper(str(record.seq))
            guide = read_sequence[GUIDE_START:GUIDE_END] #get sgRNA sequence from read
            if guide in countDict: #only count reads that match a guide sequence exactly
                countDict[guide] += 1; #add one to count for guide
                perfect_matches += 1;
                gene = guideGeneDict[guide]; #find gene associated with guide
                geneCountDict[gene] += 1; #add one to count for gene
                perfMatchList.append(read_sequence); #add read to perfect-match list
            else:
                non_perfect_matches += 1;
                nonPerfList.append(read_sequence);  #add read to non-perfect match list
        numGuidesOverHund = 0;
        # create ordered dictionary with guides and respective counts based on read counts and output as a csv file
        dict_sorted = OrderedDict(sorted(countDict.items(), key=lambda t: t[1]));
        outputName = str(output_file+"_guide_count.csv")
        with open(outputName, 'w') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',');
            mywriter.writerow(["Gene Name", "sgRNA", "Read Count", "Fraction of Total Reads"]);
            runningSum = 0.0;
            for guide in dict_sorted:
                count = dict_sorted[guide];
                if count >= 100:
                    numGuidesOverHund += 1;
                gene = guideGeneDict[guide];
                runningSum += count;
                fraction = runningSum/perfect_matches;
                mywriter.writerow([gene, guide, count, fraction]);
        #print perfect-matches to an output file
        outputName = str(output_file+"_perfect_matches.csv");
        with open(outputName, 'w') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',');
            for seq in perfMatchList:
                mywriter.writerow([seq]);
        #print non-perfect matches to an output file
        outputName = str(output_file+"_not_matched.csv");
        with open(outputName, 'w') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',');
            for seq in nonPerfList:
                mywriter.writerow([seq]);
        print("Printed perfectly matched and not-matched sequences");
        print(countDict)

        # percentage of guides that matched perfectly
        percent_matched = round(perfect_matches/float(perfect_matches + non_perfect_matches) * 100, 1)
        # percentage of undetected guides with no read counts
        guides_with_reads = np.count_nonzero(list(countDict.values()))
        guides_no_reads = len(list(countDict.values())) - guides_with_reads
        percent_no_reads = round(guides_no_reads/float(len(list(countDict.values()))) * 100, 1);
        #percent of genes with no read counts
        numGenes = len(list(geneCountDict.values()));
        genes_with_reads = np.count_nonzero(list(geneCountDict.values()));
        genes_no_reads = numGenes - genes_with_reads;
        percent_genes_no_reads = round(genes_no_reads/float(numGenes) * 100, 1);
        # skew ratio of top 10% to bottom 10% of guide counts
        top_10 = np.percentile(list(countDict.values()), 90)
        bottom_10 = np.percentile(list(countDict.values()), 10)
        if top_10 != 0 and bottom_10 != 0:
            skew_ratio = top_10/bottom_10
        else:
            skew_ratio = 'Not enough perfect matches to determine skew ratio'
        #percent of guides with more than 100 reads
        percent_guides_greater_hund = numGuidesOverHund / float(len(countDict.keys())) *100;

        # write analysis statistics to statistics.txt
        outputName = output_file+"_statistics.txt";
        with open(outputName, 'w') as infile:
            infile.write('Number of perfect read matches: ' + str(perfect_matches) + '\n')
            infile.write('Number of nonperfect read matches: ' + str(non_perfect_matches) + '\n')
            #infile.write('Number of reads where key was not found: ' + str(key_not_found) + '\n')
            infile.write('Number of reads processed: ' + str(self.count_status) + '\n')
            infile.write('Percentage of reads that matched perfectly: ' + str(percent_matched) + '\n')
            infile.write('Percentage of undetected guides: ' + str(percent_no_reads) + '\n')
            infile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n\n')
            infile.write('Number of guides with more than 100 reads: '+str(numGuidesOverHund)+'\n');
            infile.write('Percent of guides with more than 100 reads: '+str(percent_guides_greater_hund)+'\n');
            infile.write('Percent of genes with no reads: '+str(percent_genes_no_reads)+'\n');

        handle.close()
        outputName = str(output_file+"_statistics.csv")
        with open(outputName, 'w') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',');
            mywriter.writerow(["Number of perfect read matches", perfect_matches]);
            mywriter.writerow(["Number of nonperfect read matches", non_perfect_matches]);
            mywriter.writerow(["Number of reads processed", self.count_status]);
            mywriter.writerow(["Percentage of reads that matched perfectly", percent_matched]);
            mywriter.writerow(["Percentage of undetected guides", percent_no_reads]);
            mywriter.writerow(["Skew ratio of top 10% to bottom 10%", skew_ratio]);
            mywriter.writerow(["Number of guides with more than 100 reads", numGuidesOverHund]);
            mywriter.writerow(["Percent of guides with more than 100 reads", percent_guides_greater_hund]);
            mywriter.writerow(["Percent of genes with no reads", percent_genes_no_reads]);
            print("Printed Gene Enrichment stats")
        return perfect_matches;

    #extract gene names and read counts from unsorted and sorted populations, and perform Fisher's exact test to determine enrichment p-value
    def calcGeneEnrich(self, unsortedGeneCountDict, unsortedTotMatches, sortedGeneCountDict, sortedTotMatches, descriptionDict):
        geneStatsDict = {}; #new dictionary to hold Gene objects
        enrichPVals = [];
        total_keys = len(unsortedGeneCountDict.keys())
        keys = unsortedGeneCountDict.keys();
        keys_total = len(keys)
        print(keys_total)
        for counter, key in enumerate(keys):
            self.count_status = counter
            #make new gene object with the correct counts
            unsorted_count = unsortedGeneCountDict[key]
            sorted_count = sortedGeneCountDict[key]
            fold_change = math.log2(((sorted_count+1)/sortedTotMatches)/((unsorted_count+1)/unsortedTotMatches))
            gene = Gene(key, unsortedGeneCountDict[key], unsortedTotMatches, sortedGeneCountDict[key], sortedTotMatches, fold_change, descriptionDict[key]["Summary"], descriptionDict[key]["Description"]);
            geneStatsDict[key] = gene;
            oddsratio, pValue = stats.fisher_exact([[gene.unsortedReads, gene.unsortedTotal], [gene.sortedReads, gene.sortedTotal]]);
            gene.enrichPVal = pValue;
            enrichPVals.append(pValue);
        #get FDR-corrected pvalues
        adjEnrichPVals = compTool.multipletests(enrichPVals, method='fdr_bh', is_sorted=False)[1];
        for i, key in enumerate(keys):
            gene = geneStatsDict[key];
            gene.adjEnrichPVal = adjEnrichPVals[i];

        return geneStatsDict;

    #prints the enrichment p-values for each gene to a .csv file. The file is ranked from lowest FDR-corrected p-value to highest
    def printGeneEnrichStats(self, geneStatsDict, output_file):
        genesList = list(geneStatsDict.values());
        genesList.sort(key = lambda gene: gene.adjEnrichPVal);
        outputName = str(output_file+"_gene_enrichment_calculation.csv");
        with open(outputName, 'w') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',');
            mywriter.writerow(["Gene Name", "FDR-corrected p-value", "p-value", "Reads in Unsorted Population", "Total Reads in Unsorted Population", "Reads in Sorted Population", "Total Reads in Sorted Population","Fold Change", "Summary","Description"]);
            for gene in genesList:
                mywriter.writerow([gene.name, gene.adjEnrichPVal, gene.enrichPVal, gene.unsortedReads, gene.unsortedTotal, gene.sortedReads, gene.sortedTotal, gene.foldChange, gene.summary, gene.description]);
            print("Printed Gene Enrichment stats")

    #extract guides sequences and read counts from unsorted and sorted populations, and perform Fisher's exact test to determine enrichment p-value
    def calcGuideEnrich(self, unsortedGuideCountDict, unsortedTotMatches, sortedGuideCountDict, sortedTotMatches, guideGeneDict):
        print("Performing guide enrichment...")
        guideStatsDict = {}; #new dictionary to hold guide objects
        enrichPVals = [];
        keys = unsortedGuideCountDict.keys();
        keys_total = len(keys)
        for counter, key in enumerate(keys):
            #make new guide object with the correct counts
            self.count_status = counter
            guide = Guide(key, guideGeneDict[key], unsortedGuideCountDict[key], unsortedTotMatches, sortedGuideCountDict[key], sortedTotMatches);
            guideStatsDict[key] = guide;
            oddsratio, pValue = stats.fisher_exact([[guide.unsortedReads, guide.unsortedTotal], [guide.sortedReads, guide.sortedTotal]]);
            guide.enrichPVal = pValue;
            enrichPVals.append(pValue);
        #get FDR-corrected pvalues
        adjEnrichPVals = compTool.multipletests(enrichPVals, method='fdr_bh', is_sorted=False)[1];
        for i, key in enumerate(keys):
            guide = guideStatsDict[key];
            guide.adjEnrichPVal = adjEnrichPVals[i];

        return guideStatsDict;

    #prints the enrichment p-values for each guide to a .csv file. The file is ranked from lowest FDR-corrected p-value to highest
    def printGuideEnrichStats(self, guideStatsDict, output_file):
        guidesList = list(guideStatsDict.values());
        guidesList.sort(key = lambda guide: guide.adjEnrichPVal);
        outputName = str(output_file+"_guide_enrichment_calculation.csv");
        with open(outputName, 'w') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',');
            mywriter.writerow(["Guide Sequence", "Gene", "FDR-corrected p-value", "p-value", "Reads in Unsorted Population", "Total Reads in Unsorted Population", "Reads in Sorted Population", "Total Reads in Sorted Population","","Reads in Unsorted+1","Final Total Unsorted Reads","Reads in Sorted + 1", "Final Total Sorted Reads", "Fold-change", "Log2 (Fold-change)"]);
            plusUnsortedReads = ((guidesList[0]).unsortedTotal)+len(guidesList); #get the total number of reads in the unsorted population and add 1 additional read per guide
            plusSortedReads = ((guidesList[0]).sortedTotal)+len(guidesList);
            for guide in guidesList:
                #calculate normalized fold change of each guide
                normSorted = (guide.sortedReads+1.0)/(plusSortedReads);
                normUnsorted = (guide.unsortedReads+1.0)/(plusUnsortedReads);
                foldChange = normSorted/normUnsorted;
                logFoldChange = math.log(foldChange, 2); #get log2(fold-change)
                guide.logFoldChange = logFoldChange;
                mywriter.writerow([guide.sequence, guide.gene, guide.adjEnrichPVal, guide.enrichPVal, guide.unsortedReads, guide.unsortedTotal, guide.sortedReads, guide.sortedTotal,"",(guide.unsortedReads+1), plusUnsortedReads, (guide.sortedReads+1), plusSortedReads, foldChange, logFoldChange]);
            print("Printed Guide Enrichment stats");

    def printAnalysesInputs(self, guideStatsDict, output_file):
        guidesList = list(guideStatsDict.values());
        guidesList.sort(key = lambda guide: guide.sequence); #sort on sequence so replicates can be easily combined later
        #setup MaGeCK file
        mageckoutputName = str(output_file+"_mageck_input.txt");
        mageckFile = open(mageckoutputName, 'w');
        mageckFile.write("sgRNA\tGene\tControl1\tSorted1\n");
        count = 1; #need to track guide count for RIGER file
        for guide in guidesList:
            mageckFile.write(guide.sequence+"\t"+guide.gene+"\t"+str(guide.unsortedReads)+"\t"+str(guide.sortedReads)+"\n");
            count += 1; #update count

        mageckFile.close();
        print("Printed MAGeCK input files");


def runcommand(command):
    pro = os.system(command)

#main method that ties processes together
def main(argv):
    parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution');
    parser.add_argument('-s', '--sorted', type=str, dest='sorted_fastq', help="fastq file for the sorted population", default=None);
    parser.add_argument('-u', '--unsorted', type=str, dest='unsorted_fastq', help='fastq file for unsorted population', default=None);
    parser.add_argument('-o', '--output', type=str, dest='output_file',
                        help='output file name', default='library');
    parser.add_argument('-g', '--guides', type=str, dest='guides_file',
                        help='input file name', required=True);
    args = parser.parse_args();
    #check to make sure there is at least one fastq file provided
    if args.sorted_fastq == None and args.unsorted_fastq == None:
        print("Must provide at least one fastq file");
        sys.exit();

    #necessary dictionaries to keep track of guide/gene info
    countDict, guideGeneDict, geneCountDict, geneToGuideDict, geneToDescriptionDict = createDictionaries(args.guides_file);
    print("created required dictionaries");
    #count reads per guides and per gene for unsorted population if the fastq file is provided
    unsortedGuideCountDict = dict(countDict);
    unsortedGeneCountDict = dict(geneCountDict);
    unsortedTotMatches = 0;
    if args.unsorted_fastq != None:
        print("Counting reads for unsorted population");
        unsortedOutputName = str(args.output_file+"_unsorted");
        perfect_matches = count_spacers(args.unsorted_fastq, unsortedOutputName, unsortedGuideCountDict, guideGeneDict, unsortedGeneCountDict);
        unsortedTotMatches = perfect_matches;
        #print gene related info for unsorted population
        geneReadCounts(unsortedOutputName, unsortedGeneCountDict, geneToGuideDict, unsortedGuideCountDict, perfect_matches);
        print("Done counting reads for unsorted population");
    #count reads per guides and per gene for sorted population if the fastq file is provided
    sortedGuideCountDict = dict(countDict);
    sortedGeneCountDict = dict(geneCountDict);
    sortedTotMatches = 0;
    if args.sorted_fastq != None:
        print("Counting reads for sorted population");
        sortedOutputName = str(args.output_file+"_sorted");
        perfect_matches = count_spacers(args.sorted_fastq, sortedOutputName, sortedGuideCountDict, guideGeneDict, sortedGeneCountDict);
        sortedTotMatches = perfect_matches;
        #print gene related info for sorted population
        geneReadCounts(sortedOutputName, sortedGeneCountDict, geneToGuideDict, sortedGuideCountDict, perfect_matches);
        print("Done counting reads for sorted population");

    #perform enrichment analysis for genes if both unsorted and sorted files are provided
    if args.sorted_fastq != None and args.unsorted_fastq != None:
        print("Performing enrichment analysis");
        #look for genes enriched for reads
        geneStatsDict = calcGeneEnrich(unsortedGeneCountDict, unsortedTotMatches, sortedGeneCountDict, sortedTotMatches, geneToDescriptionDict);
        printGeneEnrichStats(geneStatsDict, args.output_file);
        #look for guides enriched for reads
        guideStatsDict = calcGuideEnrich(unsortedGuideCountDict, unsortedTotMatches, sortedGuideCountDict, sortedTotMatches, guideGeneDict);
        printGuideEnrichStats(guideStatsDict, args.output_file);
        # #print RIGER and MaGeCK input files
        # printAnalysesInputs(guideStatsDict, args.output_file);
        # #perform mageck analysis
        # mageckFile = str(args.output_file+"_mageck_input.txt");
        # os.system("cp /Volumes/G-DRIVE\ slim/CRISPR\ Screens/Brie-library-controls.txt ."); #need to copy in control guides file
        # command = "mageck test -k " + mageckFile + " -t 1 -c 0 -n " + str(args.output_file+"_mageck"
        #             + " --sort-criteria pos --gene-lfc-method alphamean --control-sgrna Brie-library-controls.txt");
        # os.system(command);
        # os.system("rm Brie-library-controls.txt"); #rm control guides file
        # #annotate files
        # print "Annotating Guides file and MaGeCK results";
        # geneFuncsDict = annotate_functions.readFunctions();
        # #annoate guide enrichment file
        # guidesEnrichFile = args.output_file+"_guide_enrichment_calculation.csv";
        # annotate_functions.annotateGuidesFile(guidesEnrichFile, geneFuncsDict);
        # #annotate mageck results file
        # mageckResults = args.output_file+"_mageck.gene_summary.txt";
        # annotate_functions.annotateMageckFile(mageckResults, geneFuncsDict);





#run main method
#main(sys.argv[1:]);
