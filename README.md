# screen_analyzer
A simple program to analyze genetic screen data.
Author: Collin Schlager

## Dependencies:
```python3, pip, virtualenv, and Mageck```

First, make sure you have pip and virtualenv. These will help download the other packages that the web app depends on.

```python3 -m pip install --user --upgrade pip```

```python3 -m pip --version```

```python3 -m pip install --user virtualenv```

You also need Mageck. Download the compressed source files here: [Sourceforge link](https://sourceforge.net/projects/mageck/files/latest/download)

Assuming it downloaded to your Downloads folder, go into there and uncompress the source file.

```tar xvzf mageck-0.5.8.tar.gz```
```cd mageck-0.5.8```

Install with the following command.

```python3 setup.py install```

You should now be able to run Mageck. Check that you can by entering

```mageck -v```

>If something like "command not found" is returned instead, you will need to add the command to your PATH variable. Let me know if you run into this problem. Essentially, you just need to tell your computer where to find the program named "mageck."

This should display the version number (0.5.8).

Okay, now that we have the pre-requisite software, we can download the webapp and get it running!

## Install Instructions

### Cloning the Source Code from the Repository

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

Now, we will create a virtual environment to keep all the software dependencies in one place. Note, if you haven't installed virtualenv or pip (described above), you should do so now.

```python3 -m virtualenv env```

This should create a folder called "env."

To enter the virtual environment just created, set your source using the following command.

```source env/bin/activate```

Now that we are in the virtual environment, use pip (package manager) to install the dependencies under requirements.txt. This is done with the following command

```pip install -r requirements.txt```

After this installs the necessary packages, you should be able to run webapp! To start-up the webapp local server, run the application.py file.

```./application.py```

Then, you should be able to access the webapp by going to your favorite web-browser, and going to `localhost:5000` in the address bar.



### library_embellisher.py

In order to display the additional gene information such as name and summary, the library file used in the screen analysis must be "embellished" using __library_embellisher.py__. This program polls NCBI's gene database for the additional information and updates the library file accordingly. This process can take some time due to the API polling. For the given library file, this process takes ~14 minutes. However, once your library file is "embellished," there is no need to undergo this process again. 

Embellish your library file with the following python3 commands:

```
from library_embellish import embellish

embellish(YOUR_FILE_PATH_HERE)
```

A new updated (.csv) file will be created with the additional information. Use this file for future analysis.
