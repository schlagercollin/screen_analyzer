"""
Performs analysis of qfast file given a library file
"""

import csv, sys, json, threading
from collections import OrderedDict


class parseThread(threading.Thread):
    def __init__(self, qfast_file, library_file, output_file_json):
        self.qfast_file = qfast_file
        self.library_file = library_file
        self.output_file_json = output_file_json
        self.counter = 0
        super().__init__()

    def run(self):
        library_dict = parse_lib(self.library_file)
        self.counter = 0
        results = {}
        unmatched = {}
        number_matched = 0
        number_unmatched = 0
        sorted_results = OrderedDict()
        stats = {}
        with open(self.qfast_file) as f:
            for line_number, sequence in enumerate(f):
                if line_number % (4) == 1:
                    sys.stdout.write("\r%d" % self.counter)
                    sys.stdout.flush()
                    sequence = sequence[0:20]
                    if sequence in library_dict:
                        try:
                            results[library_dict[sequence]["Target Gene Symbol"]]["frequency"] += 1
                        except KeyError:
                            results[library_dict[sequence]["Target Gene Symbol"]] = {}
                            results[library_dict[sequence]["Target Gene Symbol"]]["frequency"] = 1
                            try:
                                results[library_dict[sequence]["Target Gene Symbol"]]["description"] = library_dict[sequence]["Description"]
                            except KeyError:
                                results[library_dict[sequence]["Target Gene Symbol"]]["description"] = None
                            try:
                                results[library_dict[sequence]["Target Gene Symbol"]]["summary"] = library_dict[sequence]["Summary"]
                            except KeyError:
                                results[library_dict[sequence]["Target Gene Symbol"]]["summary"] = None
                        number_matched += 1
                    else:
                        try:
                            unmatched[sequence] += 1
                        except KeyError:
                            unmatched[sequence] = 1
                        number_unmatched += 1
                    self.counter += 1
        sorted_results = OrderedDict(sorted(list(results.items()), key=lambda x: (results[x[0]]['frequency']), reverse=True))
        stats["number_matched"] = number_matched
        stats["gene_count"] = len(sorted_results)
        stats["number_unmatched"] = number_unmatched
        with open(self.output_file_json, 'w') as output:
            output.write(json.dumps([sorted_results, unmatched, stats]))
        self.counter = "Complete"
        self.myResult = json.dumps([sorted_results, stats])
        return self.myResult
    def status(self):
        return self.counter

        

def parse_lib(filename):
    """Creates a dictionary whose keys are the target sequence and whose values are the gene data"""
    library_dict = {}
    with open(filename) as lib_file:
        reader = csv.DictReader(lib_file)  # read rows into a dictionary format
        for row in reader:
            library_dict[row["sgRNA Target Sequence"]] = row
    return library_dict

def parse_qfast(qfast_file, library_file, output_file_json):
    """Cross references fastq_file with library_file and outputs a dictionary
    containing various information (sequence, frequency, description, summary) about
    a given matched gene."""
    library_dict = parse_lib(library_file)
    counter = 0
    results = {}
    unmatched = {}
    number_matched = 0
    number_unmatched = 0
    sorted_results = OrderedDict()
    stats = {}
    with open(qfast_file) as f:
        for line_number, sequence in enumerate(f):
            if line_number % (4) == 1:
                sys.stdout.write("\r%d" % counter)
                sys.stdout.flush()
                sequence = sequence[0:20]
                if sequence in library_dict:
                    try:
                        results[library_dict[sequence]["Target Gene Symbol"]]["frequency"] += 1
                    except KeyError:
                        results[library_dict[sequence]["Target Gene Symbol"]] = {}
                        results[library_dict[sequence]["Target Gene Symbol"]]["frequency"] = 1
                        try:
                            results[library_dict[sequence]["Target Gene Symbol"]]["description"] = library_dict[sequence]["Description"]
                        except KeyError:
                            results[library_dict[sequence]["Target Gene Symbol"]]["description"] = None
                        try:
                            results[library_dict[sequence]["Target Gene Symbol"]]["summary"] = library_dict[sequence]["Summary"]
                        except KeyError:
                            results[library_dict[sequence]["Target Gene Symbol"]]["summary"] = None
                    number_matched += 1
                else:
                    try:
                        unmatched[sequence] += 1
                    except KeyError:
                        unmatched[sequence] = 1
                    number_unmatched += 1
                counter += 1
    sorted_results = OrderedDict(sorted(list(results.items()), key=lambda x: (results[x[0]]['frequency']), reverse=True))
    stats["number_matched"] = number_matched
    stats["gene_count"] = len(sorted_results)
    stats["number_unmatched"] = number_unmatched
    with open(output_file_json, 'w') as output:
        output.write(json.dumps([sorted_results, unmatched, stats]))
    return json.dumps([sorted_results, stats])


def load_from_file(file_name):
    """Loads json file produced by 'parse_qfast' and returns the same information as
    'parse_qfast' as if you analyzed the files for the first time."""
    with open(file_name, "r") as data_file:
        loaded_dictionaries = json.loads(data_file.read())
    matched = loaded_dictionaries[0]
    unmatched = loaded_dictionaries[1]
    stats = loaded_dictionaries[2]
    return json.dumps([matched, stats])


if __name__ == "__main__":
    library_dict = "/Users/collinschlager/PycharmProjects/screen_analyzer/tmp/data/library/Mouse_Kinome_list_Brie.csv"
    result = parse_qfast("/Users/collinschlager/PycharmProjects/screen_analyzer/tmp/data/fastq/100mM-rep1-Unsorted_S1_L001_R1_001.fastq", library_dict, "output.json")
    print(result)
