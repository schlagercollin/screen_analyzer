"""
Flask controller for parse fastq webapp.
"""

import os, json, threading, sys
from flask import Flask, render_template, request, jsonify, url_for
from werkzeug.utils import secure_filename
from screen_analyzer import *

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

def upload_file(myfile):
    """upload file to corresponding location"""
    if allowed_file(myfile.filename):
        filename = secure_filename(myfile.filename)
        if type=="fastq":
            myfile.save(os.path.join(app.config['UPLOAD_FOLDER'], "fastq", filename))
        elif type=="library":
            myfile.save(os.path.join(app.config['UPLOAD_FOLDER'], "library", filename))
        return filename
    else:
        return False

def check_data_files():
    """get data files and return them"""
    data_files = {}
    data_files['fastq'] = os.listdir(os.path.join(UPLOAD_FOLDER, 'fastq'))
    data_files['library'] = os.listdir(os.path.join(UPLOAD_FOLDER, 'library'))
    data_files['output'] = os.listdir(os.path.join(UPLOAD_FOLDER, 'output'))
    return data_files


@app.route('/analysis/submit', methods=['POST'])
def analysis_submit():
    """on analysis click: fetches stashed files or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        fastq = request.values['fastq']
        if fastq == "Upload your own":
            fastq = request.files['fastq']
            fastq = upload_file(fastq)
        library = request.values['library']
        print("Request: ", request.files)
        if 'library' in request.files and request.files['library'].filename != '':
            library = request.files['library']
            print("File: ",library)
            library = upload_file(library)
        output = analyze_data(fastq, library)
        return jsonify(result=output)

@app.route('/analysis/load', methods=['POST'])
def analysis_load():
    """on load click: fetches stashed result file or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        result_file = request.values['result_file']
        if result_file == "Upload your own":
            result_file = request.files['result_file']
            result_file = upload_file(result_file)
        result_file = os.path.join(UPLOAD_FOLDER, "output", result_file)
        output = load_from_file(result_file)
        return jsonify(result=output)

def analyze_data(fastq, library):
    """wrapper for parse_qfast function. handles some path information"""
    print(UPLOAD_FOLDER, fastq, library)
    output_file = os.path.join(UPLOAD_FOLDER, "output", fastq+library+".json")
    library = os.path.join(UPLOAD_FOLDER, 'library', library)
    fastq = os.path.join(UPLOAD_FOLDER,'fastq', fastq)
    return parse_qfast(fastq, library, output_file)


@app.route('/')
@app.route('/index')
@app.route('/analysis')
def index():
    """main index route. just updates current data files and returns them for drop down menu purposes"""
    data_files = check_data_files()
    return render_template("index.html", data_files=data_files)




if __name__ == '__main__':
    print(sys.version)
    app.run(debug=True)

