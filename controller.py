#!/Users/collinschlager/anaconda/bin/python3
"""
Flask controller for parse fastq webapp.
"""

import os, json, threading, sys, csv
from flask import Flask, render_template, request, jsonify, url_for
from werkzeug.utils import secure_filename
import screen_analyzer
import library_embellish
import time

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
        fastq = request.values['fastq']
        if fastq == "Upload your own":
            fastq = request.files['fastq']
            fastq = upload_file(fastq, "fastq")
        library = request.values['library']
        print("Request: ", request.files)
        output_file_name = request.values['file_name']
        if 'library' in request.files and request.files['library'].filename != '':
            library = request.files['library']
            print("File: ",library)
            library = upload_file(library, "library")
        output = analyze_data(fastq, library, output_file_name)
        return jsonify(result=output)

@app.route('/analysis/load', methods=['POST'])
def analysis_load():
    """on load click: fetches stashed result file or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        result_file = request.values['result_file']
        if result_file == "Upload your own":
            result_file = request.files['result_file']
            result_file = upload_file(result_file, "result")
        result_file = os.path.join(UPLOAD_FOLDER, "output", result_file)
        output = load_from_file(result_file)
        return jsonify(result=output)

def load_from_file(result_file):
    print(result_file)
    with open(result_file, "r") as csv_file:
        reader = csv.reader(csv_file)
        next(reader)
        unsorted_list = list(reader)
    print(unsorted_list)
    #unsorted_list.sort(key=lambda x: float(x[7]), reverse=True)
    return unsorted_list

@app.route('/analysis/status', methods=['POST'])
def analysis_status():
    if request.method == 'POST':
        result = None
        output = check_status()
        print(output)
        if output == "Complete":
            result = get_result()
        return jsonify(myStatus=output, result=result)

@app.route('/analysis/status', methods=['POST'])
def analysis_get_result():
    if request.method == 'POST':
        output = get_result()
        return jsonify(result=output)

global myThread
myThread = None

def analyze_data(fastq, library, file_name):
    global myThread
    """wrapper for parse_qfast function. handles some path information"""
    print(UPLOAD_FOLDER, fastq, library)
    matched_file = os.path.join(UPLOAD_FOLDER, "output", file_name+"_matched.csv")
    unmatched_file = os.path.join(UPLOAD_FOLDER, "output", file_name+"_unmatched.csv")
    stats_file = os.path.join(UPLOAD_FOLDER, "output", file_name+".json")
    library = os.path.join(UPLOAD_FOLDER, 'library', library)
    fastq = os.path.join(UPLOAD_FOLDER,'fastq', fastq)
    myThread = screen_analyzer.parseThread(fastq, library, matched_file, unmatched_file, stats_file)
    myThread.start()
    return True

def check_status():
    global myThread
    myStatus = myThread.status()
    return myThread.status()

def get_result():
    global myThread
    myResult = myThread.myResult
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
