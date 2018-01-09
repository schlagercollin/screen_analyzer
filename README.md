# screen_analyzer
A simple program to analyze genetic screen data.
Author: Collin Schlager

### Dependencies:
```python3, python3-flask```

### Instructions:
Run with ```python3 controller.py```
Then go to web browser and enter ```localhost:5000```.

Enter fastq file and library file (.csv) and click analyze. Wait a few seconds.
The program will create an output (.json) file in ```/tmp/data/output``` whose name is a concatenation of the fastq and library filenames.

The result of the analyis will be displayed below. To load a previous analysis, select the file and click load.
