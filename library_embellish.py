"""Takes a library file and polls NCBI to get some more information
Uses four cores to split up the API requests.
For the library file give, this process takes 886.94 seconds (14 minutes)

TODO: batch queries to speed this up.
"""

import csv
import time
import urllib.request
import itertools
from xml.dom import minidom
from multiprocessing import Pool


API_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%d&retmode=xml'

class Gene:
    def __init__(self, gene_id):
        url = API_URL % gene_id
        dom = minidom.parse(urllib.request.urlopen(url))
        try:
            self.summary = dom.getElementsByTagName('Entrezgene_summary')[0].firstChild.nodeValue
        except:
            self.summary = None
        self.id = gene_id
        self.description = dom.getElementsByTagName('Gene-ref_desc')[0].firstChild.nodeValue


def parse_lib(filename):
    library_dict = {}
    with open(filename) as lib_file:
        reader = csv.DictReader(lib_file)  # read rows into a dictionary format
        for row in reader:
            library_dict[row["sgRNA Target Sequence"]] = row
    return library_dict

def add_column(reader , output_file):
    with open(input_file,'r') as csvinput:
        with open(output_file, 'w') as csvoutput:
            reader = csv.DictReader(csvinput)

            new_file = []
            counter = 0

            for row in reader:
                counter += 1
                gene = Gene(int(row["Target Gene ID"]))
                row["Summary"] = gene.summary
                row["Description"] = gene.description
                new_file.append(row)
                print(counter)

            fieldnames = reader.fieldnames
            fieldnames.append("Summary")
            fieldnames.append("Description")

            writer = csv.DictWriter(csvoutput, fieldnames=fieldnames, lineterminator='\n')

            writer.writeheader()
            writer.writerows(new_file)

def process_part(reader):
    new_file = []
    counter = 0

    for row in reader:
        counter += 1
        gene = Gene(int(row["Target Gene ID"]))
        row["Summary"] = gene.summary
        row["Description"] = gene.description
        new_file.append(row)
        print(counter)
    return new_file


def process(rows):
    result = []
    try:
        gene = Gene(int(rows[1]))
        rows.append(gene.summary)
        rows.append(gene.description)
    except:
        pass
    result.append(rows)
    return result
    # for row in rows:
    #     try:
    #         gene = Gene(int(row[1]))
    #         row.append(gene.summary)
    #         row.append(gene.description)
    #     except ValueError:
    #         pass
    #     result.append(row)
    #
    # return result

def gen_chunks(reader, chunksize=100):
    """
    Chunk generator. Take a CSV `reader` and yield
    `chunksize` sized slices.
    """
    chunk = []
    for index, line in enumerate(reader):
        if (index % chunksize == 0 and index > 0):
            yield chunk
            del chunk[:]
        chunk.append(line)
    yield chunk

def embellish(FILE):
    print("Running...")
    start = time.time()
    N = 4
    #FILE = "tmp/data/library/Mouse_Kinome_list_Brie.csv"
    #FILE = "test.csv"
    filename = FILE.split('.')[0]
    OUTPUT_FILE = filename+"_updated"+".csv"
    results = []
    output = []
    pool = Pool(N)
    with open(FILE) as source_file:
        reader = csv.reader(source_file)
        fieldnames = reader.__next__() #pass over header
        #chunks = gen_chunks(reader, chunksize=int(file_length/N)+1)
        result = pool.map(process, reader)
        print("Complete!")
    fieldnames.append("Summary")
    fieldnames.append("Description")
    pool.close()
    pool.join()
    flat_list = [item for sublist in result for item in sublist]
    output_rows = [fieldnames] + flat_list
    with open(OUTPUT_FILE, 'w') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(output_rows)
    end = time.time()
    duration = end - start
    print("Time Elapsed: ", duration)
    return OUTPUT_FILE
