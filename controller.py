#!/Users/collinschlager/anaconda/bin/python3
"""
Flask controller for parse fastq webapp.
"""

import os, json, threading, sys, csv, shutil
from flask import Flask, render_template, request, jsonify, url_for, send_from_directory
from werkzeug.utils import secure_filename
#import screen_analysis_BP
import supafast
import library_embellish
import time
import logging
import pandas as pd
#logging.basicConfig()
#logging.getLogger().setLevel(logging.DEBUG)

curdir = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(curdir, 'tmp/data')
print(UPLOAD_FOLDER)
ALLOWED_EXTENSIONS = ['csv', 'fastq']

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    """check if uploaded file is allowed"""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def upload_file(myfile, file_type):
    """upload file to corresponding location"""
    if allowed_file(myfile.filename):
        filename = secure_filename(myfile.filename)
        if file_type=="fastq":
            myfile.save(os.path.join(app.config['UPLOAD_FOLDER'], "fastq", filename))
        elif file_type=="library":
            myfile.save(os.path.join(app.config['UPLOAD_FOLDER'], "library", filename))
        elif file_type=="result":
            myfile.save(os.path.join(app.config['UPLOAD_FOLDER'], "output", filename))
        return filename
    else:
        return False

def check_data_files():
    """get data files and return them"""
    data_files = {}

    data_files['fastq'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'fastq')) if not f.startswith('.')]
    data_files['library'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'library')) if not f.startswith('.')]
    data_files['output'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'output')) if not f.startswith('.')]
    return data_files


@app.route('/analysis/submit', methods=['POST'])
def analysis_submit():
    """on analysis click: fetches stashed files or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':

        sorted_fastq = request.values['sorted']
        if sorted_fastq == "Upload your own":
            sorted_fastq = request.files['sorted']
            sorted_fastq = upload_file(sorted_fastq, "fastq")

        unsorted_fastq = request.values['unsorted']
        if unsorted_fastq == "Upload your own":
            unsorted_fastq = request.files['unsorted']
            unsorted_fastq = upload_file(unsorted_fastq, "fastq")

        guides = request.values['guides']
        if 'guides' in request.files and request.files['guides'].filename != '':
            guides = request.files['guides']
            print("File: ",guides)
            guides = upload_file(guides, "library")

        output_file_name = request.values['output']

        output = analyze_data(sorted_fastq, unsorted_fastq, output_file_name, guides)
        return jsonify(result=output)

@app.route('/analysis/load', methods=['POST'])
def analysis_load():
    """on load click: fetches stashed result file or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        result_file = request.values['result_file']
        print(result_file)
        if result_file == "Upload your own":
            result_file = request.files['result_file']
            result_file = upload_file(result_file, "result")

        result_file_prefix = os.path.join(UPLOAD_FOLDER, "output", result_file, result_file)
        gene_enrichment, guide_enrichment = load_from_dir(result_file_prefix)
        print("made it here")
        return gene_enrichment

@app.route('/downloads/<path:dir_name>', methods=['GET', 'POST'])
def download(dir_name):
    dir_name_path = os.path.join(UPLOAD_FOLDER, "output", dir_name)
    #print("Making Archive...")
    #shutil.make_archive(dir_name_path, 'zip', dir_name_path)
    #print("Done.")
    file_name = dir_name+".zip"
    directory_path = os.path.join(UPLOAD_FOLDER, "output")
    print(directory_path, file_name)
    shutil.make_archive(dir_name_path, 'zip', dir_name_path)
    return send_from_directory(directory=directory_path, filename=file_name)

def load_from_dir(result_file_prefix):
    # function to load relevant portions from a directory
    #sorted_stats_file = result_file_prefix + "_sorted_statistics.csv"
    #unsorted_stats_file = result_file_prefix + "_unsorted_statistics.csv"
    gene_enrichment_file = result_file_prefix + "_gene_enrichment_calculation.csv"
    guide_enrichment_file = result_file_prefix + "_guide_enrichment_calculation.csv"
    gene_enrichment_df = pd.read_csv(gene_enrichment_file)
    guide_enrichment_df = pd.read_csv(guide_enrichment_file)
    gene_enrichment = str(list(gene_enrichment_df.T.to_dict().values()))
    guide_enrichment = str(list(guide_enrichment_df.T.to_dict().values()))
    return gene_enrichment, guide_enrichment
    #mageck_file = result_file_prefix + "_mageck.sgrna_summary.txt"
    #gene_enrichment = load_csv_as_list(gene_enrichment_file, skip_header=True)
    #sorted_stats = load_csv_as_list(sorted_stats_file, skip_header=False)
    #unsorted_stats = load_csv_as_list(unsorted_stats_file, skip_header=False)
    #mageck_results = load_csv_as_list(mageck_file, skip_header=True, delimiter="\t")
    #return gene_enrichment, sorted_stats, unsorted_stats, mageck_results

def load_csv_as_list(file_name, skip_header=False, delimiter=','):
    with open(file_name, "r") as csv_file:
        reader = csv.reader(csv_file, delimiter=delimiter)
        if skip_header == True:
            next(reader)
        return list(reader)

@app.route('/analysis/status', methods=['POST'])
def analysis_status():
    if request.method == 'POST':
        result = None
        count_status, status = check_status()
        if status == "Analysis Complete":
            print("Analysis Confirmed Complete!")
        return jsonify(count_status=count_status, status=status)

@app.route('/analysis/status', methods=['POST'])
def analysis_get_result():
    if request.method == 'POST':
        output = get_result()
        return jsonify(result=output)

global myThread
myThread = None

def analyze_data(sorted, unsorted, output, guides, control="Brie_Kinome_controls.txt"):
    global myThread
    """wrapper for parse_qfast function. handles some path information"""
    print(UPLOAD_FOLDER, sorted, unsorted, output, guides)

    os.mkdir(os.path.join(UPLOAD_FOLDER, 'output', output)) # Make output directory
    output = os.path.join(UPLOAD_FOLDER, 'output', output, output) # Set output path to new dir

    guides = os.path.join(UPLOAD_FOLDER, 'library', guides) # Get input file paths
    sorted = os.path.join(UPLOAD_FOLDER,'fastq', sorted)
    unsorted = os.path.join(UPLOAD_FOLDER,'fastq', unsorted)
    control = os.path.join(UPLOAD_FOLDER, control)
    print("About to start the thread...")

    myThread = supafast.parseThread(sorted, unsorted, output, guides, control_file=None)
    myThread.start()

    return True

def check_status():
    # Get count status and status of analysis thread
    global myThread
    count_status = myThread.count_status
    status = myThread.status
    return count_status, status

def get_result():
    global myThread
    myResult = myThread.output
    return myResult

@app.route('/compare')
def compare():
    data_files = check_data_files()
    return render_template("compare.html", data_files=data_files)

@app.route('/compare/submit', methods=['POST'])
def compare_submit():
    """on analysis click: fetches stashed files or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        unsorted_pop = request.values['unsorted_pop']
        if unsorted_pop == "Upload your own":
            unsorted_pop = request.files['unsorted_pop']
            unsorted_pop = upload_file(unsorted_pop, "fastq")
        sorted_pop = request.values['sorted_pop']
        if sorted_pop == "Upload your own":
            sorted_pop = request.files['sorted_pop']
            sorted_pop = upload_file(sorted_pop, "fastq")
        output_file = request.values['output_file_name']
        if output_file == '':
            output_file = "comparison"+ str(int(time.time())) + ".json"
        output = compare_data(unsorted_pop, sorted_pop, output_file)
        return jsonify(result=output)

def compare_data(unsorted_pop, sorted_pop, output_file):
    """wrapper for compare function. handles some path information"""
    output_file = os.path.join(UPLOAD_FOLDER, "comparisons", output_file)
    unsorted_pop = os.path.join(UPLOAD_FOLDER, 'output', unsorted_pop)
    sorted_pop = os.path.join(UPLOAD_FOLDER,'output', sorted_pop)
    output = screen_analyzer.compare(unsorted_pop, sorted_pop, output_file)
    return output

@app.route('/')
@app.route('/index')
@app.route('/analysis')
def index():
    """main index route. just updates current data files and returns them for drop down menu purposes"""
    data_files = check_data_files()
    return render_template("index.html", data_files=data_files)

class embellishThread(threading.Thread):
    def __init__(self, FILE):
        self.file = FILE
        super().__init__()
    def run(self):
        print("Embellish started.")
        library_embellish.embellish(self.file)

@app.route('/embellish/submit', methods=['POST'])
def embellish_load():
    """on load click: fetches stashed result file or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        result_file = request.values['library']
        if result_file == "Upload your own":
            result_file = request.files['library']
            result_file = upload_file(result_file, "library")
        result_file = os.path.join(UPLOAD_FOLDER, "library", result_file)
        #output = library_embellish.embellish(result_file)
        myThread = embellishThread(result_file)
        myThread.start()
        output = "Embellishing process started. Check back for embellished file."
        return jsonify(result=output)


@app.route('/embellish')
def embellish():
    """main route for embellish"""
    data_files = check_data_files()
    return render_template("embellish.html", data_files=data_files)




if __name__ == '__main__':
    app.run(debug=True)
