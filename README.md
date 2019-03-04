# screen_analyzer
A simple program to analyze genetic screen data.
Author: Collin Schlager

### Dependencies:
```python3, python3-flask```

### Install Instructions

##### Cloning the Source Code from the Repository

First, clone the github repository into a directory of your choosing. 

Here, I show how you can make a directory in your Documents folder called "Screens"

```cd```
```cd Documents```
```mkdir Screens```

With the directory created, you can now "change directory" into the new folder.

```cd Screens```

Next, you are going to "clone" the github repository (this downloads the source code of the webapp into a new folder)

```git clone https://github.com/schlagercollin/screen_analyzer.git```

Listing the files, you should now see a new folder titled "screen_analyzer."

```ls```

Let's checkout the stuff that downloaded by changing directory into the new folder.

```cd screen_analyzer```




### Instructions:
Run with ```python3 controller.py```
Then go to web browser and enter ```localhost:5000```.

Enter fastq file and library file (.csv) and click analyze. (Note, while you can upload both, it is faster to place the files in /tmp/data/fastq and /tmp/data/library, respectively, to skip upload times.
The program will create an output (.json) file in ```/tmp/data/output``` whose name is a concatenation of the fastq and library filenames.

The result of the analyis will be displayed below. To load a previous analysis, select the file and click load.
=======
A Flask-based web app for CRISPR/Cas9 screen analysis.
Author: Collin Schlager

### Dependencies:
```python3, python3-flask, python3-numpy, python3-pandas, python3-pyaml, python3-scipy```

### Instructions:
- Fork or download source files
- Change directory into source file root folder: ```cd <insert-path-where-you-put-the-code>/screen_analyzer```
- Run ```python3 application.py```
- Then go to web browser and enter ```localhost:5000```.

Enter an analysis name, your screen .fastq files, your library file, and a control file. Then, select the anlayses you wish to perform and hit run.

The anlaysis will take some time, as the .fastq files need to be parsed and then cross-referenced with your library file. Next, the counts are normalized and put through the chosen statistical pipelines.

The result will redirect you to the ```analysis/load``` url, where plots and a table representing the data will be rendered. You may download a .zip file of the result directory on this page.

Clicking on a gene will display more information about it, including its symbol, frequency in the screen run, name, and a summary (if available). In order to obtain this additional information, please see the section below.

### library_embellisher.py

In order to display the additional gene information such as name and summary, the library file used in the screen analysis must be "embellished" using __library_embellisher.py__. This program polls NCBI's gene database for the additional information and updates the library file accordingly. This process can take some time due to the API polling. For the given library file, this process takes ~14 minutes. However, once your library file is "embellished," there is no need to undergo this process again. 

Embellish your library file with the following python3 commands:

```
from library_embellish import embellish

embellish(YOUR_FILE_PATH_HERE)
```

A new updated (.csv) file will be created with the additional information. Use this file for future analysis.
