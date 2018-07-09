#!/usr/local/bin/python3
"""
Flask controller for parse fastq webapp.
"""

import os, json, threading, sys, csv, shutil, subprocess
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
ALLOWED_EXTENSIONS = ['csv', 'fastq', 'zip']

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

def is_valid_data_file(filename, include_dirs=False):
    if not filename.startswith('.'):
        try:
            if filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS:
                if include_dirs == True:
                    return False
                else:
                    return True
        except IndexError:
            if include_dirs == True:
                return True
            else:
                return False
    else:
        return False

def check_data_files():
    """get data files and return them"""
    data_files = {}

    data_files['fastq'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'fastq')) if is_valid_data_file(f)]
    data_files['library'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'library')) if is_valid_data_file(f)]
    data_files['output'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'output')) if is_valid_data_file(f, include_dirs=True)]
    return data_files


@app.route('/analysis/submit', methods=['POST'])
def analysis_submit():
    """on analysis click: fetches stashed files or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':

        top_sorted_pop = request.values['top_sorted']
        if top_sorted_pop == "Upload your own":
            top_sorted_pop = request.files['top_sorted']
            top_sorted_pop = upload_file(top_sorted_pop, "fastq")

        bot_sorted_pop = request.values['bot_sorted']
        if bot_sorted_pop == "Upload your own":
            bot_sorted_pop = request.files['bot_sorted']
            bot_sorted_pop = upload_file(bot_sorted_pop, "fastq")

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

        if "Mageck" in request.values:
            print("Mageck enabled!")
            mageck = True
        else:
            mageck = False

        if "Fischer" in request.values:
            print("Fischer enabled!")
            fischer = True
        else:
            fischer = False

        if "Ratio" in request.values:
            print("Ratio enabled.")
            ratio = True
        else:
            ratio = False

        config = {"Mageck": mageck, "Fischer": fischer, "Ratio": ratio}
        print(config)

        output = analyze_data(top_sorted_pop, bot_sorted_pop, unsorted_fastq, output_file_name, guides, config)
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

        result_file_dir = os.path.join(UPLOAD_FOLDER, "output", result_file)
        result_file_prefix = os.path.join(UPLOAD_FOLDER, "output", result_file, result_file)

        try:
            contents = os.listdir(result_file_dir)
            date = time.strftime("%D %H:%M", time.localtime(int(os.path.getctime(result_file_dir))))
        except:
            contents = None
            date = "Unable to fetch date information."
        gene_enrichment, guide_enrichment = load_from_dir(result_file_prefix)
        print("made it here")
        return jsonify(data=gene_enrichment, contents=contents, date=date)

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
        status, output_dir = check_status()
        output_dir = str(output_dir)
        if status == "Analysis Complete":
            print("Analysis Confirmed Complete!")
        return jsonify(status=status, output_dir=output_dir)

@app.route('/analysis/status', methods=['POST'])
def analysis_get_result():
    if request.method == 'POST':
        output = get_result()
        return jsonify(result=output)

global myThread
myThread = None

def analyze_data(top_sorted, bot_sorted, unsorted, output, guides, config, control="Brie_Kinome_controls.txt"):
    global myThread
    """wrapper for parse_qfast function. handles some path information"""
    print(UPLOAD_FOLDER, sorted, unsorted, output, guides)

    os.mkdir(os.path.join(UPLOAD_FOLDER, 'output', output)) # Make output directory
    output_dir = os.path.join(UPLOAD_FOLDER, 'output', output) # Set output path to new dir
    output = os.path.join(UPLOAD_FOLDER, 'output', output, output) # Set output path to new dir

    guides = os.path.join(UPLOAD_FOLDER, 'library', guides) # Get input file paths
    top_sorted = os.path.join(UPLOAD_FOLDER,'fastq', top_sorted)
    bot_sorted = os.path.join(UPLOAD_FOLDER, 'fastq', bot_sorted)
    unsorted = os.path.join(UPLOAD_FOLDER,'fastq', unsorted)
    control = os.path.join(UPLOAD_FOLDER, control)
    print("About to start the thread...")

    myThread = supafast.parseThread(top_sorted, bot_sorted, unsorted, output, guides, output_dir, control_file="Brie_Kinome_controls.txt")
    myThread.config_analysis = config
    myThread.start()

    return True

def check_status():
    # Get count status and status of analysis thread
    global myThread
    status = myThread.status
    output_dir = myThread.output_prefix
    return status, output_dir

def get_result():
    global myThread
    myResult = myThread.output
    return myResult

@app.route('/')
@app.route('/index')
@app.route('/analysis')
def index():
    """main index route. just updates current data files and returns them for drop down menu purposes"""
    data_files = check_data_files()
    return render_template("index.html", data_files=data_files)

@app.route('/load')
def load_route():
    data_files = check_data_files()
    return render_template("analysis_load.html", data_files=data_files)

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

@app.route('/about')
def about():
    """main route for about"""
    return render_template("about.html")




if __name__ == '__main__':
    app.run(debug=True)
