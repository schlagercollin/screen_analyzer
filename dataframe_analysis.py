import pandas as pd
import numpy as np
import sys

def parse_qfast(qfast_file, library_file, output_file):
    library = pd.read_csv(library_file) # read in library file as pandas data frame
    library.set_index("sgRNA Target Sequence", inplace=True) # set search index as guide sequences
    library["Counts"] = np.zeros(len(library), dtype="int") # create new column of zeros for Counts

    with open(qfast_file) as f: #open fastq file
        scan_counter = 0
        for line_number, sequence in enumerate(f): #iterate over fastq file
            if line_number % (4) == 1: #only care about the sequence (every 4 lines)
                sys.stdout.write("\r%d" % scan_counter) # console status info
                sys.stdout.flush() # console status info
                sequence = sequence[0:20]
                if sequence in library.index:
                    library.at[sequence,"Counts"] += 1
                scan_counter += 1

    print("Saving to Output File")
    library.to_csv(output_file)
    print("Saved!")

if __name__ == "__main__":
    qfast_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_data/100mM-rep1-Unsorted_S1_L001_R1_001.fastq"
    library_file = "/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/library/Brie_crispr_library_with_controls_for_analysis_updated.csv"
    output_file = "output.csv"
    print("Beginning parsing...")
    parse_qfast(qfast_file, library_file, output_file)
    print("All done!")
