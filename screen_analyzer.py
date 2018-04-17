"""
Performs analysis of qfast file given a library file
"""

import csv, sys, json, threading, math, os, copy, time, datetime
from collections import OrderedDict


class parseThread(threading.Thread):
    def __init__(self, qfast_file, library_file, frequency_output_file, unmatched_output_file, stats_output_file):
        self.qfast_file = qfast_file
        self.library_file = library_file
        self.frequency_output_file = frequency_output_file
        self.unmatched_output_file = unmatched_output_file
        self.stats_output_file = stats_output_file
        self.counter = 0
        super().__init__()

    def run(self):
        library_dict = parse_lib(self.library_file)
        #print(library_dict)
        streamlined_result_dict = OrderedDict()
        unmatched = {}
        number_matched = 0
        number_unmatched = 0
        stats = {'time': str(datetime.date.today()),
                    'total_reads': 0,
                    'number_aligned': 0,
                    'percent_aligned': 0,
                    'number_unaligned': 0,
                    'percent_guides_zero_reads': 0,
                    'hundred_plus_reads': 0,
                    'percent_hundred_plus_reads': 0,
                    'percent_genes_zero_reads': "To be determined",
                    'output_file': 0,
                    'output_unmatched_file': 0
                }

        with open(self.qfast_file) as f:
            for line_number, sequence in enumerate(f):
                if line_number % (4) == 1:
                    sys.stdout.write("\r%d" % self.counter)
                    sys.stdout.flush()
                    sequence = sequence[0:20]
                    if sequence in library_dict:
                        try:
                            library_dict[sequence]["Frequency"] += 1
                        except KeyError:
                            library_dict[sequence]["Frequency"] = 1
                        try:
                            new_key = library_dict[sequence]["Target Gene Symbol"]
                            streamlined_result_dict[new_key]["Frequency"] += 1
                        except KeyError:
                            streamlined_result_dict[new_key] = {}
                            streamlined_result_dict[new_key]["Frequency"] = 1
                            try:
                                streamlined_result_dict[new_key]["Description"] = library_dict[sequence]["Description"]
                            except KeyError:
                                streamlined_result_dict[new_key]["Description"] = None
                            try:
                                streamlined_result_dict[new_key]["Summary"] = library_dict[sequence]["Summary"]
                            except KeyError:
                                streamlined_result_dict[new_key]["Summary"] = None
                        number_matched += 1
                    else:
                        try:
                            unmatched[sequence]["Frequency"] += 1
                        except KeyError:
                            unmatched[sequence] = {}
                            unmatched[sequence]["Unmatched Sequence"] = sequence
                            unmatched[sequence]["Frequency"] = 1
                        number_unmatched += 1
                    self.counter += 1

            guides_zero_reads = 0
            greater_than_100_reads = 0
            for sequence in library_dict:
                try:
                    if library_dict[sequence]["Frequency"] == 0:
                        guides_zero_reads += 1
                    elif library_dict[sequence]["Frequency"] >= 100:
                        greater_than_100_reads += 1
                except KeyError:
                    library_dict[sequence]["Frequency"] = 0
                    guides_zero_reads += 1


        stats['total_reads'] = self.counter
        stats['number_aligned'] = number_matched
        stats['percent_aligned'] = (number_matched/self.counter)*100
        stats['number_unaligned'] = number_unmatched
        stats['percent_guides_zero_reads'] = (guides_zero_reads / len(library_dict))*100
        stats['hundred_plus_reads'] = greater_than_100_reads
        stats['percent_hundred_plus_reads'] = (greater_than_100_reads / number_matched)*100
        stats['output_file'] = self.frequency_output_file
        stats['output_unmatched_file'] = self.unmatched_output_file

        write_to_file(library_dict, self.frequency_output_file)
        write_to_file(unmatched, self.unmatched_output_file)

        sorted_results = OrderedDict(sorted(list(streamlined_result_dict.items()), key=lambda x: (streamlined_result_dict[x[0]]['Frequency']), reverse=True))

        with open(self.stats_output_file, "w") as outfile:
            json.dump([sorted_results, unmatched, stats], outfile)

        self.myResult = json.dumps([sorted_results,unmatched,stats])
        self.counter = "Complete"
    def status(self):
        return self.counter

def get_files():
    dict1 = load_from_file(os.path.join('tmp/data/output', os.listdir('tmp/data/output')[2]), output="file")
    dict2 = load_from_file(os.path.join('tmp/data/output', os.listdir('tmp/data/output')[-1]), output="file")
    return dict1, dict2

def compare(first, second, output_file_json, output="json"):
    first,blank,blank = load_from_file(first, output="file")
    second,blank,blank = load_from_file(second, output="file")
    compare_dict = {}
    for gene in first:
        if gene in second:
            value1 = first[gene]["Frequency"]
            value2 = second[gene]["Frequency"]
            ratio = value2 / value1 # sorted / unsorted
            compare_dict[gene] = copy.deepcopy(first[gene])
            compare_dict[gene]["logRatio"] = math.log(ratio, 2)
    #sorted_results = compare_dict
    sorted_results = OrderedDict(sorted(list(compare_dict.items()), key=lambda x: (compare_dict[x[0]]["logRatio"]), reverse=True))
    with open(output_file_json, 'w') as output_file:
        output_file.write(json.dumps([sorted_results]))
    if output=="json":
        return json.dumps([sorted_results])
    else:
        return sorted_results


def parse_lib(filename):
    """Creates a dictionary whose keys are the target sequence and whose values are the gene data"""
    library_dict = {}
    with open(filename) as lib_file:
        reader = csv.DictReader(lib_file)  # read rows into a dictionary format
        for row in reader:
            library_dict[row["sgRNA Target Sequence"]] = row
    return library_dict

def parse_qfast(qfast_file, library_file, frequency_output_file, unmatched_output_file, stats_output_file):
    """Cross references fastq_file with library_file and outputs a dictionary
    containing various information (sequence, frequency, description, summary) about
    a given matched gene."""
    library_dict = parse_lib(library_file)
    unmatched = {}
    number_matched = 0
    number_unmatched = 0
    scan_counter = 0
    stats = {'time': time.time(),
                'total_reads': 0,
                'number_aligned': 0,
                'percent_aligned': 0,
                'number_unmatched': 0,
                'percent_guides_zero_reads': 0,
                'hundred_plus_reads': 0,
                'percent_hundred_plus_reads': 0,
                'percent_genes_zero_reads': 0,
                'output_file': 0,
                'output_unmatched_file': 0
            }

    with open(qfast_file) as f:
        for line_number, sequence in enumerate(f):
            if line_number % (4) == 1:
                sys.stdout.write("\r%d" % scan_counter)
                sys.stdout.flush()
                sequence = sequence[0:20]
                if sequence in library_dict:
                    try:
                        library_dict[sequence]["Frequency"] += 1
                    except KeyError:
                        library_dict[sequence]["Frequency"] = 1
                    number_matched += 1
                else:
                    try:
                        unmatched[sequence]["Frequency"] += 1
                    except KeyError:
                        unmatched[sequence] = {}
                        unmatched[sequence]["Unmatched Sequence"] = sequence
                        unmatched[sequence]["Frequency"] = 1
                    number_unmatched += 1
                scan_counter += 1

        guides_zero_reads = 0
        greater_than_100_reads = 0
        for sequence in library_dict:
            try:
                if library_dict[sequence]["Frequency"] == 0:
                    guides_zero_reads += 1
                elif library_dict[sequence]["Frequency"] >= 100:
                    greater_than_100_reads += 1
            except KeyError:
                library_dict[sequence]["Frequency"] = 0
                guides_zero_reads += 1


    stats['total_reads'] = scan_counter
    stats['number_aligned'] = number_matched
    stats['percent_aligned'] = (number_matched/scan_counter)*100
    stats['number_unaligned'] = number_unmatched
    stats['percent_guides_zero_reads'] = (guides_zero_reads / len(library_dict))*100
    stats['hundred_plus_reads'] = greater_than_100_reads
    stats['percent_hundred_plus_reads'] = (greater_than_100_reads / len(library_dict))*100
    stats['output_file'] = frequency_output_file
    stats['output_unmatched_file'] = unmatched_output_file

    write_to_file(library_dict, frequency_output_file)
    write_to_file(unmatched, unmatched_output_file)

    with open(stats_output_file, "w") as outfile:
        json.dump([sorted_results, unmatched, stats], outfile)

    return library_dict, stats


def write_to_file(library_dict, file_name):
    with open(file_name, "w") as output_file:
        try:
            fieldnames = list((next(iter(library_dict.values()))).keys()) #get fieldnames from libary_dict
        except:
            fieldnames = list((library_dict.keys()))
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        for item in library_dict:
            writer.writerow(library_dict[item])
    print("Done writing output!")

def load_from_csv(file_name):
    library_dict = {}
    unmatched = {}
    streamlined_result_dict = OrderedDict()
    with open(file_name, "r") as lib_file:
        reader = csv.DictReader(lib_file)  # read rows into a dictionary format
        for row in reader:
            if int(row["Frequency"]) > 0:
                #print(row["Target Gene Symbol"], row["Frequency"])
                library_dict[row["sgRNA Target Sequence"]] = row
                symbol_key = row["Target Gene Symbol"]
                try:
                    streamlined_result_dict[symbol_key]["Frequency"] += int(row["Frequency"])
                except KeyError:
                    streamlined_result_dict[symbol_key] = {}
                    streamlined_result_dict[symbol_key]["Frequency"] = int(row["Frequency"])
                    streamlined_result_dict[symbol_key]["Description"] = row["Description"]
                    streamlined_result_dict[symbol_key]["Summary"] = row["Summary"]
    stats = {'time': time.time(),
                'total_reads': 0,
                'number_aligned': 0,
                'percent_aligned': 0,
                'number_unmatched': 0,
                'percent_guides_zero_reads': 0,
                'hundred_plus_reads': 0,
                'percent_hundred_plus_reads': 0,
                'percent_genes_zero_reads': 0,
                'output_file': 0,
                'output_unmatched_file': 0
    }
    sorted_results = OrderedDict(sorted(list(streamlined_result_dict.items()), key=lambda x: (streamlined_result_dict[x[0]]['Frequency']), reverse=True))
    print("Ready to return.")
    return json.dumps([sorted_results, unmatched, stats])

def load_from_file(file_name, output="json"):
    """Loads json file produced by 'parse_qfast' and returns the same information as
    'parse_qfast' as if you analyzed the files for the first time."""
    print(file_name[-4:])
    if file_name[-4:] == ".csv":
        return load_from_csv(file_name)
    else:
        with open(file_name, "r") as data_file:
            loaded_dictionaries = json.loads(data_file.read())
        print(loaded_dictionaries)
        sorted_results = loaded_dictionaries[0]
        unmatched = loaded_dictionaries[1]
        stats = loaded_dictionaries[2]
        if output=="json":
            return json.dumps([sorted_results, unmatched, stats])
        else:
            return sorted_results, unmatched, stats


if __name__ == "__main__":
    dictionary, stats = parse_qfast("/Users/collinschlager/Documents/Rohatgi_Lab/screen_data/100mM-rep1-Unsorted_S1_L001_R1_001.fastq","/Users/collinschlager/Documents/Rohatgi_Lab/screen_data/Mouse_Kinome_list_Brie.csv","/Users/collinschlager/Documents/Rohatgi_Lab/screen_data/output_frequency.csv", "/Users/collinschlager/Documents/Rohatgi_Lab/screen_data/output_unmatched.csv", "/Users/collinschlager/Documents/Rohatgi_Lab/screen_data/output_stats.json")
    #library_dict = "/Users/collinschlager/PycharmProjects/screen_analyzer/tmp/data/library/Mouse_Kinome_list_Brie.csv"
    #result = parse_qfast("/Users/collinschlager/PycharmProjects/screen_analyzer/tmp/data/fastq/100mM-rep1-Unsorted_S1_L001_R1_001.fastq", library_dict, "output.json")
    #print(result)
