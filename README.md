# screen_analyzer
A simple program to analyze genetic screen data.
Author: Collin Schlager

### Dependencies:
```python3, python3-flask```

### Instructions:
Run with ```python3 controller.py```
Then go to web browser and enter ```localhost:5000```.

<img src="https://user-images.githubusercontent.com/23715298/34709947-3a9910e2-f4ce-11e7-8adf-be76f3ac27ab.png" width=600px>

Enter fastq file and library file (.csv) and click analyze. (Note, while you can upload both, it is faster to place the files in /tmp/data/fastq and /tmp/data/library, respectively, to skip upload times.
The program will create an output (.json) file in ```/tmp/data/output``` whose name is a concatenation of the fastq and library filenames.

The result of the analyis will be displayed below. To load a previous analysis, select the file and click load.

<img src="https://user-images.githubusercontent.com/23715298/34710029-9348e8c0-f4ce-11e7-80e1-8a744e13fd3f.png" width=600px>

Clicking on a gene will display more information about it, including its symbol, frequency in the screen run, name, and a summary (if available). In order to obtain this additional information, please see the section below.

<img src="https://user-images.githubusercontent.com/23715298/34710046-a60c30f2-f4ce-11e7-8674-35ef8a4a0535.png" width=600px>

### library_embellisher.py
In order to display the additional gene information such as name and summary, the library file used in the screen analysis must be "embellished" using __library_embellisher.py__. This program polls NCBI's gene database for the additional information and updates the library file accordingly. This process can take some time due to the API polling. For the given library file, this process takes ~14 minutes. However, once your library file is "embellished," there is no need to undergo this process again. 

Embellish your library file with the following python3 commands:

```
from library_embellish import embellish

embellish(YOUR_FILE_PATH_HERE)
```

A new updated (.csv) file will be created with the additional information. Use this file for future analysis.
