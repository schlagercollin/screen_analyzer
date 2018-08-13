#!/usr/local/bin/python3
"""
Flask controller for parse fastq webapp.
"""

import os, json, threading, sys, csv, shutil, subprocess
from flask import Flask, render_template, request, jsonify, url_for, send_from_directory
from werkzeug.utils import secure_filename
# import supafast
import collin_screen_analysis as screen_analysis
import library_embellish
import time
import logging
import pandas as pd
import copy
import yaml
import shutil
#
# werkzeug_log = logging.getLogger('werkzeug')
# werkzeug_log.setLevel(logging.ERROR)

my_logger = logging.getLogger(__name__)
my_logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
my_logger.addHandler(ch)

curdir = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(curdir, 'tmp/data')
ANALYSIS_FOLDER = os.path.join(curdir, 'tmp/Screen_Analyses')
#print(UPLOAD_FOLDER)
ALLOWED_EXTENSIONS = ['csv', 'fastq', 'zip', 'txt']

application = Flask(__name__)
application.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def debug(message):
    return my_logger.debug(message)

def allowed_file(filename):
    """check if uploaded file is allowed"""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def upload_file(myfile, file_type):
    """upload file to corresponding location"""
    if allowed_file(myfile.filename):
        filename = secure_filename(myfile.filename)
        if file_type=="fastq":
            myfile.save(os.path.join(application.config['UPLOAD_FOLDER'], "Fastq", filename))
        elif file_type=="library":
            myfile.save(os.path.join(application.config['UPLOAD_FOLDER'], "Library", filename))
        elif file_type=="control":
            myfile.save(os.path.join(application.config['UPLOAD_FOLDER'], "Control", filename))
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

    data_files['fastq'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'Fastq')) if is_valid_data_file(f)]
    data_files['library'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'Library')) if is_valid_data_file(f)]
    data_files['output'] = [f for f in os.listdir(ANALYSIS_FOLDER) if is_valid_data_file(f, include_dirs=True)]
    data_files['control'] = [f for f in os.listdir(os.path.join(UPLOAD_FOLDER, 'Control')) if is_valid_data_file(f)]
    return data_files


@application.route('/analysis/submit', methods=['POST'])
def analysis_submit():
    """on analysis click: fetches stashed files or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':

        raw_values, analysis_input_config = get_analysis_files(request)

        # Run analyze_data wrapper function
        output = analyze_data(analysis_input_config)
        return jsonify(result=output)

def get_analysis_files(request_object):
    """ Extract analysis files from frontend request object"""

    # Begin to create input files nested dictionary object
    analysis_input_config = dict()

    # Store how many replicates so we can form our analysis_input_config object
    number_of_replicates = int(request_object.values["replicate_number"])
    # For each replicate, there will be a dictionary with the corresponding
    # fastq files stored in a overarching dictionary called 'Fastq' (below)
    input_files = dict()
    input_files["Fastq"] = dict()
    for replicate_number in range(number_of_replicates):
        rep = str(replicate_number+1) #so its not zero-indexed
        replicate_files = dict()
        replicate_files["Bottom Sorted Population"] = retrieve_request_info(request_object, "bot_sorted_"+rep)
        replicate_files["Top Sorted Population"] = retrieve_request_info(request_object, "top_sorted_"+rep)
        replicate_files["Unsorted Population"] = retrieve_request_info(request_object, "unsorted_"+rep)
        input_files["Fastq"]["Rep"+rep] = replicate_files

    # Get True or False for whether each analysis type was checked
    analyses_queued = dict()
    for analysis_type in ["Mageck", "Ratio", "Fischer"]:
        analyses_queued[analysis_type] = bool(retrieve_request_info(request_object, analysis_type))
    if analyses_queued["Mageck"]:
        mageck_analysis_types = dict()
        mageck_analysis_types["Top_Sorted"] = \
                            bool(replicate_files["Top Sorted Population"])
        mageck_analysis_types["Bottom_Sorted"] = \
                            bool(replicate_files["Bottom Sorted Population"])
        analyses_queued["Mageck"] = mageck_analysis_types

    analysis_input_config["Analyses Queued"] = analyses_queued


    # Extracting remaining file data from request object passed from frontend
    input_files["Library"] = retrieve_request_info(request_object, "library")
    input_files["Control"] = retrieve_request_info(request_object, "control")
    analysis_input_config["Input Files"] = input_files

    # Extracting meta data from request object
    metadata = dict()
    metadata["Notes"] = retrieve_request_info(request_object, "notes")
    metadata["Timestamp"] = retrieve_request_info(request_object, "timestamp")
    metadata["Replicates"] = number_of_replicates
    analysis_input_config["Metadata"] = metadata

    folder_name = retrieve_request_info(request_object, "analysis_name")

    # Create analysis folder
    analysis_folder = create_analysis_folders(folder_name, analysis_input_config["Analyses Queued"])
    metadata["Analysis Path"] = analysis_folder
    metadata["Analysis Name"] = os.path.basename(analysis_folder)

    return request.values, analysis_input_config

def save_config_file(config_file, analysis_folder, filename):
    analysis_config_path = os.path.join(analysis_folder, filename)
    with open(analysis_config_path, "w") as outfile:
        # json.dump(analysis_input_config, outfile, sort_keys=True, indent=4, separators=(',',':'))
        yaml.dump(config_file, outfile, default_flow_style=False)

def retrieve_request_info(request_object, request_item):
    try:
        item_value = request_object.values[request_item]
    except KeyError:
        debug("%s item not found in request." % request_item)
        return None

    if item_value == "Upload your own":
        return request_object.files[request_item]
    elif item_value == "--" or item_value == "":
        debug("%s item not found in request." % request_item)
        return None
    else:
        return request_object.values[request_item]

def create_analysis_folders(analysis_name, analysis_types):
    if os.path.exists(os.path.join(ANALYSIS_FOLDER, analysis_name)):
        debug("Analysis folder %s already exists." % analysis_name)
        incremental_value = 1
        new_analysis_name = analysis_name+"_"+str(incremental_value)
        while os.path.exists(os.path.join(ANALYSIS_FOLDER, new_analysis_name)):
            incremental_value += 1
            new_analysis_name = analysis_name+"_"+str(incremental_value)
        debug("Creating '%s' instead." % new_analysis_name)

        dir_path = os.path.join(ANALYSIS_FOLDER, new_analysis_name)
        os.mkdir(dir_path)

    else:
        debug("Analysis folder %s created." % analysis_name)

        dir_path = os.path.join(ANALYSIS_FOLDER, analysis_name)
        os.mkdir(dir_path)

    for analysis, checked in analysis_types.items():
        if checked:
            os.mkdir(os.path.join(dir_path, analysis))

    return dir_path

def delete_analysis_folder(analysis_name):
    """Function to delete an analysis folder and all of its contents"""
    filepath = os.path.join(ANALYSIS_FOLDER, analysis_name)
    if os.path.isdir(filepath):

        shutil.rmtree(os.path.join(ANALYSIS_FOLDER, analysis_name))
        return "Deleted"
    else:
        return "Folder does not exist"

def list_analysis_folders():
    """Function to return analysis folder names"""
    return os.listdir(ANALYSIS_FOLDER)

@application.route('/admin/analysis/remove/<path:dir_name>')
def delete_analysis_folder_route(dir_name):
    result = delete_analysis_folder(dir_name)
    return jsonify(action_result=result, folders=list_analysis_folders())

@application.route('/admin/analysis/removeall')
def delete_all_analysis_folders_route():
    for folder in list_analysis_folders():
        delete_analysis_folder(folder)
    return jsonify(folders=list_analysis_folders())

@application.route('/admin/analysis/list')
def list_analysis_folder_route():
    folders = list_analysis_folders()
    return jsonify(folders=folders)


class Analysis:
    def __init__(self, analysis_name):
        self.path = os.path.join(ANALYSIS_FOLDER, analysis_name)
        self.name = analysis_name
        self.path_and_prefix = os.path.join(self.path, analysis_name)

        self.get_config_info()

    def get_config_info(self):
        config_file_path = os.path.join(self.path, "analysis_config.yaml")
        with open(config_file_path, "r") as config_file:
            self.config_info = yaml.load(config_file)

        self.analyses_performed = self.config_info["Analyses Queued"]
        self.config_name = self.config_info["Metadata"]["Analysis Name"]
        self.timestamp = self.config_info["Metadata"]["Timestamp"]
        self.replicates = self.config_info["Metadata"]["Replicates"]
        self.notes = self.config_info["Metadata"]["Notes"]

        return self.config_info

    def load_specific_data(self, analysis_type, level, columns, sort_by, ascending, datapoints=1000):
        df = self.load_specific_df(analysis_type, level)
        if df is not None:
            print(sort_by, ascending)
            print(type(ascending))
            df = df.sort_values(by=sort_by, ascending=ascending)
            df = df.head(n=int(datapoints))
            if columns == "all":
                return df.tolist()
            else:
                return df[columns].tolist()
        else:
            return None

    def load_specific_df(self, analysis_type, level):
        if analysis_type == "Counts":
            df = self.get_counts_result(level=level)
            df = df.drop(index="Non-Targeting-Control") # Drop Non-Targeting Control
        elif analysis_type == "Mageck Top":
            df = self.get_mageck_result(level, "Top Sorted")
        elif analysis_type == "Mageck Bottom":
            df = self.get_mageck_result(level, "Bottom Sorted")
        elif analysis_type == "Ratio":
            df = self.get_ratio_result(level)
        else:
            raise "Analysis type %s not recognized in load_specific_df." % analysis_type
        return df

    def get_counts_result(self, level="Gene"):
        if level == "Gene":
            counts_file = self.path_and_prefix + "_gene_counts.csv"
        else:
            counts_file = self.path_and_prefix + "_guide_counts.csv"
        counts_df = pd.read_csv(counts_file, header=[0,1], index_col=[0])
        counts_df.columns = flattenHierarchicalCol(counts_df.columns)
        counts_df = counts_df.fillna("--")
        return counts_df

    def get_ratio_result(self, level):
        if level == "Gene":
            ratio_file = os.path.join(self.path, "Ratio")
            ratio_file = os.path.join(ratio_file, self.name+"_ratio_genes.csv")
        else:
            ratio_file = os.path.join(self.path, "Ratio")
            ratio_file = os.path.join(ratio_file, self.name+"_ratio_guides.csv")
        try:
            ratio_df = pd.read_csv(ratio_file)
        except:
            ratio_df = None
        return ratio_df

    def get_mageck_result(self, level, condition):
        # Build correct mageck file path based on level and condition
        mageck_analysis_directory = os.path.join(self.path, "Mageck")
        mageck_path_and_prefix = os.path.join(mageck_analysis_directory, self.name)
        mageck_file_path = mageck_path_and_prefix
        if condition == "Bottom Sorted":
            mageck_file_path += "_Bottom_Sorted_mageck"
        elif condition == "Top Sorted":
            mageck_file_path += "_Top_Sorted_mageck"
        if level == "Gene":
            mageck_file_path += "_gene_results.csv"
        elif level == "Guide":
            mageck_file_path += "_guide_results.csv"

        try:
            mageck_result_df = pd.read_csv(mageck_file_path)
            mageck_result_df = mageck_result_df.fillna("--")
        except:
            mageck_result_df = None
        return mageck_result_df

    def fetch_combined_data(self, datapoints=1000):
        print("Fetching combined data...")
        gene_counts_df = self.get_counts_result()
        combined_df = gene_counts_df
        mageck_statistics_df = self.get_mageck_result("Gene", "Top Sorted")
        if mageck_statistics_df is not None:
            mageck_statistics_df = mageck_statistics_df.set_index("id", drop=True)
            combined_df = combined_df.join(mageck_statistics_df)
        ratio_statistics_df = self.get_ratio_result("Gene")
        if ratio_statistics_df is not None:
            ratio_statistics_df = ratio_statistics_df.set_index("Target Gene Symbol", drop=True)
            combined_df = combined_df.join(ratio_statistics_df)

        combined_df = combined_df.fillna("--")

        initial_col_order = [
            'Target Gene Symbol'
        ]

        combined_df = self.truncate_df(combined_df, initial_col_order, datapoints)

        json_data, columns = self.return_json(combined_df)

        return json_data, columns

    def truncate_df(self, dataframe, initial_column_order, datapoints):
        remaining_cols = list(set(dataframe.columns) - set(initial_column_order))
        columns = initial_column_order + remaining_cols

        dataframe = dataframe[ columns ]
        dataframe = dataframe.head(n=int(datapoints))

        return dataframe


    def fetch_general_data(self, datapoints=1000):
        gene_counts_df = self.get_counts_result()
        # Defines the order of the first couple of columns
        my_columns = [
            'Target Gene Symbol',
            'Target Gene ID',
            'Description',
            'Summary',
            'sgRNA Target Sequence',
            'Rep1: Top Sorted Population',
            'Rep1: Unsorted Population',
        ]
        my_columns += list(set(gene_counts_df.columns) - set(my_columns))
        gene_counts_df = gene_counts_df[ my_columns ]
        gene_counts_df = gene_counts_df.head(n=int(datapoints))
        json_data, columns = return_json(gene_counts_df)
        return json_data, columns

    def return_json(self, dataframe):
        columns = dataframe.columns.tolist()
        json_data = eval(dataframe.to_json(orient="records"))
        return json_data, columns



@application.route('/analysis/load_specific', methods=['POST'])
def analysis_load_specific():
    if request.method == 'POST':

        analysis_name = request.values['analysis_name']
        analysis_type = request.values['analysis_type']
        analysis_level = request.values['analysis_level']
        analysis_columns = request.values['analysis_columns']
        analysis_sort_by = request.values['analysis_sort_by']
        analysis_ascending = request.values['analysis_ascending']
        analysis_datapoints = request.values['analysis_datapoints']

        print(analysis_ascending, type(analysis_ascending))

        if analysis_ascending == "false":
            analysis_ascending = False
        else:
            analysis_ascending = True

        print(analysis_ascending, type(analysis_ascending))

        analysis = Analysis(analysis_name)

        data = analysis.load_specific_data(analysis_type, analysis_level, \
                                            analysis_columns, analysis_sort_by, \
                                            analysis_ascending, \
                                            datapoints=analysis_datapoints)

        return jsonify(data=data)


@application.route('/analysis/load', methods=['POST'])
def analysis_load():
    """on load click: fetches stashed result file or uploads and uses the new one
    returns json object for display (see index.html)"""
    if request.method == 'POST':
        # Fetch form value: analysis name to load
        analysis_name = request.values['analysis_name']
        try:
            datapoints = request.values['datapoints']
        except:
            datapoints = 1000
        analysis = Analysis(analysis_name)
        data, columns = analysis.fetch_combined_data(datapoints=datapoints)

        print("Passing data to frontend")
        return jsonify(analysis_info=analysis.config_info, data=data, columns=columns, name=analysis_name)

def load_general_info(output_path, output_name):
    """function to load count data from a result directory"""
    config_file_path = os.path.join(output_path, "analysis_config.yaml")
    with open(config_file_path, "r") as config_file:
        analysis_config = yaml.load(config_file)
    return analysis_config, gene_counts_df, guide_counts_df

def load_analysis_data(output_path, output_name, analysis_type, number_of_datapoints=1000):
    """function to load slice of analysis dataframes for quick transfer to frontend"""
    result_file_prefix = os.path.join(output_path, output_name)
    gene_counts_file = result_file_prefix + "_gene_counts.csv"
    guide_counts_file = result_file_prefix + "_guide_counts.csv"

    gene_counts_df = pd.read_csv(gene_counts_file, header=[0,1], index_col=[0])
    guide_counts_df = pd.read_csv(guide_counts_file)

def get_extracted_data(dataframe, number=1000):
    """function to sort dataframe and return top <number> of hits for
    easy transfer to frontend"""

    analyses = {'mageck': '-log(pos|p-value)',
                'ratio': 'ZScore',
                'fischer':'-log(FDR-Corrected P-Values)',
                'top_sorted': 'Top Sorted Counts',
                'bot_sorted': 'Bot Sorted Counts',
                'unsorted': 'Unsorted Counts'}

    print("Retrieving %d datapoints", number)
    outputs = copy.deepcopy(analyses)

    for name, sort_value in analyses.items():
        if sort_value in dataframe.columns:
            extracted = extract_data(dataframe, sort_value, number=number)
            outputs[name] = extracted
        else:
            outputs[name] = False

    print(outputs)

    return outputs

def extract_data(dataframe, analysis_type, number=1000):
    """Sorts given dataframe by analysis_type. Slices top 1000 in format passable
    to frontend javascript"""
    df = dataframe
    df.columns = load_analysis_panelrarchicalCol(df.columns)
    #assert (analysis_type in df.columns), "Analysis type not found in df columns"
    df.sort_values(by=analysis_type, ascending=False, inplace=True)
    df = df.head(n=number)
    return str(list(df.T.to_dict().values()))

def load_csv_as_list(file_name, skip_header=False, delimiter=','):
    with open(file_name, "r") as csv_file:
        reader = csv.reader(csv_file, delimiter=delimiter)
        if skip_header == True:
            next(reader)
        return list(reader)

def flattenHierarchicalCol(columns, sep = ': '):
    """Flattens multi-index columns by joining tuple with separator
    e.g. ("Rep1", "Unsorted Population") --> "Rep1_Unsorted Population"""
    flattenedCols = []
    for column_name in columns.values:
        if isinstance(column_name, tuple):
            if column_name[0].startswith("Rep"):
                flattenedCols.append(sep.join(column_name).rstrip(sep))
            else:
                flattenedCols.append(column_name[0])
        else:
            flattenedCols.append(column_name)
    return flattenedCols

def prep_for_javascript(df):
    """Prepares dataframe in a format that is enterable into dataframe-js"""
    df.columns = flattenHierarchicalCol(df.columns, sep = '_')
    return df

@application.route('/downloads/<path:dir_name>', methods=['GET', 'POST'])
def download(dir_name):
    result_folder_path = os.path.join(ANALYSIS_FOLDER, dir_name)
    if os.path.isdir(result_folder_path):
        file_name = dir_name+".zip"
        shutil.make_archive(result_folder_path, 'zip', result_folder_path)
        deleteTimer = threading.Timer(10.0, deleteFile, [result_folder_path+".zip"])
        deleteTimer.start()
        return send_from_directory(directory=ANALYSIS_FOLDER, filename=file_name)
    else:
        error_msg = "Analysis folder %s not found." % dir_name
        return jsonify(Error=error_msg)

def deleteFile(filename):
    if os.path.exists(filename):
        print("DEL %r" % filename)
        os.remove(filename)
    return "Deleted"

@application.route('/analysis/status', methods=['POST'])
def analysis_status():
    if request.method == 'POST':
        result = None
        status, name = check_status()
        if status == "Analysis Complete":
            print("Analysis Confirmed Complete!")
        return jsonify(status=status, name=name)

@application.route('/analysis/status', methods=['POST'])
def analysis_get_result():
    if request.method == 'POST':
        output = get_result()
        return jsonify(result=output)

global myThread
myThread = None

def convert_to_path(d):
    """Goes through input file dictionary and converts filenames to filepaths
    Returns a new dictionary object with the same structure as the input"""
    new_dict = copy.deepcopy(d)
    for key, value in d.items():
        if key == "Fastq":
            for rep, conditions in d["Fastq"].items():
                for condition, filename in conditions.items():
                    if filename:
                        new_dict["Fastq"][rep][condition] = os.path.join(UPLOAD_FOLDER, "Fastq", filename)
        else:
            if value:
                new_dict[key] = os.path.join(UPLOAD_FOLDER, key, value)
    return new_dict

def analyze_data(analysis_input):
    global myThread
    """wrapper for parse_qfast function. handles some path information"""

    # Convert input file names to path name
    input_file_paths = convert_to_path(analysis_input["Input Files"])
    analysis_input["Input File Paths"] = input_file_paths

    analysis_name = analysis_input["Metadata"]["Analysis Name"]
    output_dir = analysis_input["Metadata"]["Analysis Path"]
    output_prefix = analysis_name

    save_config_file(analysis_input, output_dir, "analysis_config.yaml")

    myThread = screen_analysis.analysisThread(analysis_input, output_prefix, output_dir)
    # myThread.run()
    myThread.start()

    return "Thread Started"

def check_status():
    # Get count status and status of analysis thread
    global myThread
    status = myThread.status
    name = myThread.analysis_name
    if status == "Analysis Complete":
        print("Joining thread...")
        myThread.join()
        print("Joined!")
    return status, name

def get_result():
    global myThread
    myResult = myThread.output
    return myResult

@application.route('/')
@application.route('/index')
@application.route('/analysis')
def index():
    """main index route. just updates current data files and returns them for drop down menu purposes"""
    data_files = check_data_files()
    return render_template("index.html", data_files=data_files)

@application.route('/test')
def test():
    return render_template("test.html")

@application.route('/load')
@application.route('/load/<analysis_name>')
def load_route(analysis_name=None):
    data_files = check_data_files()
    return render_template("analysis_load.html", data_files=data_files, analysis_name=analysis_name)

class embellishThread(threading.Thread):
    def __init__(self, FILE):
        self.file = FILE
        super().__init__()
    def run(self):
        print("Embellish started.")
        library_embellish.embellish(self.file)

@application.route('/embellish/submit', methods=['POST'])
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


@application.route('/embellish')
def embellish():
    """main route for embellish"""
    data_files = check_data_files()
    return render_template("embellish.html", data_files=data_files)

@application.route('/about')
def about():
    """main route for about"""
    return render_template("about.html")




if __name__ == '__main__':
    application.run(debug=True)
