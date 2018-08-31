from Bio import Entrez
from Bio import SeqIO
import json
import pandas as pd
Entrez.email = "schlager@stanford.edu" # email for identification purposes

class Protein_Entry:
    def __init__(self, accession):
        self.accession = accession

    def annotate(self, gene, gene_synonyms, db_xref, description):
        self.gene = gene
        self.gene_synonyms = gene_synonyms
        self.description = description
        self.extractXREF(db_xref)


    def extractXREF(self, db_xref):
        if db_xref == None:
            self.CCDS = None
            self.GeneID = None
            self.MGI = None
        else:
            for xref in db_xref:
                if xref.startswith("CCDS:"):
                    self.CCDS = xref.split(":")[-1]
                elif xref.startswith("GeneID:"):
                    self.GeneID = xref.split(":")[-1]
                elif xref.startswith("MGI:"):
                    self.MGI = ":".join(xref.split(":")[-2:])

def get_NCBI_data(accession_number_list):
    """Direct interface with NCBI Entrez database.
        Returned more than mygene did."""
    number_of_accession_numbers = len(accession_number_list)
    number_of_unique_accession_numbers = len(list(set(accession_number_list)))
    print("Fetching data for %d entries (%d unique accession numbers)..." % (number_of_accession_numbers, number_of_unique_accession_numbers))
    acc_list = list(map(lambda x: x+"[accn]", accession_number_list))
    query = (" OR ").join(acc_list)

    with Entrez.esearch(db="protein", retmax=3000, term=query, idtype="acc", usehistory="y") as handle:
        record = Entrez.read(handle)
        errors_list = record["ErrorList"]['PhraseNotFound']
        print("%d entries returned information." % len(record["IdList"]))
        print("%d entries returned no records." % len(errors_list))
        # print(errors_list)

    with Entrez.efetch(db="protein", query_key=record['QueryKey'], webenv=record['WebEnv'], rettype="gp") as handle:
        with open("output.gp","w") as output_file:
            output_file.write(handle.read())

        records = SeqIO.parse("output.gp","gb")

        entries = []

        for record in records:
            new_entry = Protein_Entry(record.id)

            description = record.description
            feature = [feature for feature in record.features if feature.type == "CDS"][0]

            try:
                gene = feature.qualifiers['gene'][0]
            except KeyError:
                gene = None
            try:
                gene_synonyms = feature.qualifiers['gene_synonym'][0].split("; ")
            except KeyError:
                gene_synonyms = None
            try:
                db_xref = feature.qualifiers['db_xref']
            except KeyError:
                db_xref = None

            new_entry.annotate(gene, gene_synonyms, db_xref, description)
            entries.append(new_entry)

        return pd.DataFrame([vars(entry) for entry in entries], columns=["accession","gene","GeneID","description","MGI","CCDS","gene_synonyms"]), errors_list
