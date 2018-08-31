# screen_analyzer
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

** note: this section has not been updated recently. While is still works as designed, a new and improved library embellisher is on its way!

In order to display the additional gene information such as name and summary, the library file used in the screen analysis must be "embellished" using __library_embellisher.py__. This program polls NCBI's gene database for the additional information and updates the library file accordingly. This process can take some time due to the API polling. For the given library file, this process takes ~14 minutes. However, once your library file is "embellished," there is no need to undergo this process again. 

Embellish your library file with the following python3 commands:

```
from library_embellish import embellish

embellish(YOUR_FILE_PATH_HERE)
```

A new updated (.csv) file will be created with the additional information. Use this file for future analysis.
