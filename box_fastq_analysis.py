from __future__ import print_function, unicode_literals
import os
import io
import csv
from boxsdk import Client
from boxsdk.exception import BoxAPIException
from boxsdk.object.collaboration import CollaborationRole
from auth import authenticate
import pandas as pd

APP_ROOT_FOLDER_ID = "53019807615"

class BoxAccess():
    def __init__(self, root_folder_id):
        oauth, _, _ = authenticate()
        self.client = Client(oauth)
        self.root = self.client.folder(folder_id=root_folder_id)
        self.curdir = self.root

    def search(self, query, limit=100, offset=0, **kwargs):
        search_result = self.client.search(query, ancestor_folders=[self.root], limit=limit, offset=offset, **kwargs)
        return search_result

    def ls(self, limit=100, offset=0):
        return self.curdir.get_items(limit=limit, offset=offset)

    def cd(self, folder_name):
        if folder_name == "":
            self.curdir = self.root
        else:
            try:
                my_folder_id = self.create_id_mapping()[folder_name]
                self.curdir = self.client.folder(folder_id=my_folder_id)
            except KeyError:
                print("%s was not found as a directory." % folder_name)

    def create_id_mapping(self, limit=100, offset=0):
        folder_contents = self.curdir.get_items(limit=limit, offset=offset)
        mapping = {}
        for folder in folder_contents:
            if folder.type == "folder":
                mapping[folder.name] = folder.id
        return mapping

    def upload_file_from_stream(self, stream, filename):
        """
        Upload contents of stream to filename in self.curdir
        """
        self.curdir.upload_stream(stream, filename)

    def stream_by_lines(self, fileObj):
        """
        Custom stream file. Input file object.
        Output an iterable that streams out the file contents.
        """
        url = fileObj.get_url('content')
        box_response = fileObj._session.get(url, expect_json_response=False, stream=True)
        streamObj = box_response.network_response._request_response.iter_lines() #get Response object
        return streamObj

def parse_lib(library_generator):
    """Creates a dictionary whose keys are the target sequence and whose values are the gene data"""
    print("Creating dict...")
    library_dict = {}
    reader = csv.DictReader(library_generator)  # read rows into a dictionary format
    for row in reader:
        library_dict[row["sgRNA Target Sequence"].encode()] = 0
    print("Done. %d" % (len(library_dict)))
    return library_dict

def parse_fastq_from_box(fastq_generator, library_generator):
    counts_dict = parse_lib(library_generator)
    print("Parsing file...")
    scan_counter = 0
    for index, line in enumerate(fastq_generator):
        if index % (4) == 1: #sequence lines are every four
            sequence = line[0:20]
            if sequence in counts_dict:
                counts_dict[sequence] += 1
    print("Done")
    return counts_dict

def string_stream(stream):
    """Converts bytes generator to string generator"""
    for item in stream:
        yield item.decode("utf-8")


if __name__ == "__main__":
    me = BoxAccess(APP_ROOT_FOLDER_ID)
    me.cd("Library")
    library_file = me.ls()[0]
    library_stream = string_stream(me.stream_by_lines(library_file))
    me.cd("")
    me.cd("Fastq")
    fastq_file = me.ls()[0]
    fastq_stream = me.stream_by_lines(fastq_file)
    counts_dict = parse_fastq_from_box(fastq_stream, library_stream)
    me.cd("")
    with open('dict.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in counts_dict.items():
           writer.writerow([key, value])
    with open("dict.csv", "rb") as csv_file:
        me.upload_file_from_stream(csv_file, "counts.csv")
    print("All done.")
