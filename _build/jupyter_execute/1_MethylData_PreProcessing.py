#!/usr/bin/env python
# coding: utf-8

# # Data Download and Initial Processing

# ## Overview

# 
# This is a comprehensive chapter on the pipeline utilized to download and process raw methylation array data from _Infinium MethylationEPIC_ and _Infinium HumanMethylation450_:
# 
# 1. Download raw datasets from GEO or GDC using ```methylprep```
# 2. Pre-process raw data using ```methylprep```, 
#     which is a python library that offers R's [SeSAMe](https://www.bioconductor.org/packages/release/bioc/html/sesame.html) as the backend analysis suite.
# 2. Process data exclusion criteria using ```methylcheck```
# 3. Run Illumina's Quality Control
# 4. Prepare data for analysis
# 
# ```{note}
# This pipeline has been adapted from python's [methylsuite](https://pypi.org/project/methylsuite/).
# ```

# ## Download and Pre-Process Raw Data from GEO

# 
# 1. Install Python
# 
#     - Go to python.org and download the appropriate version of python
#     - During install, click on the box to enable it on PATH
# 
# 2. Open a terminal (Windows: git bash, Linux/Mac: terminal)
# 
#     - Move to the directory/folder where your code files will be stored
# 
# 
# 3. Create a virtual environment with python 3.7 using either pip or conda:
# 
# ```bash
# python37 -m venv <venv_name>
# ```
# or, if python37 location cannot be found:
# 
# ```bash
#  ~/AppData/Local/Programs/Python/Python37/python.exe -m venv <venv_name>
# ```
# 
# 4. Then activate the environment:
#     - If using Mac/Linux:
#     ```bash
#     source activate <venv_name>
#     ```
#     - If using Windows (git bash):
#     ```bash
#     source <venv_name>/Scripts/activate
#     ```
# 
# 5. Install methylsuite and dependencies:
# 
# ```bash
# pip install statsmodels
# pip install methylsuite
# ```
# 
# 6. Download dataset directly from GEO:
# 
# ```bash
# python -m methylprep beta_bake -i <GEO series> -d <directory> 
# ```
# 
# 5. Process dataset:
# 
# ```bash
# python -m methylprep -v process --array_type <epic or 450k> -d <directory> --all --no_sample_sheet
# ```
# 
# ```{note}
# If files are larger than RAM, specify batches
# Full ```methylprep``` API can be found [here](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/docs/cli.html)

# 
# ## Download and Pre-Process Raw Data from Genomic Data Commons (GDC)

# 
# - Go to [GDC](https://portal.gdc.cancer.gov/)
# - Click on _Repository_
# - In _Files_ on the left side of the screen, find _Experimental Strategy_
# - Select _Methylation Array_, which will prompt you to see all cases and files containing methylation array data available
# - For our project, we will select the following filter/search strategy, which you can copy and paste into _advanced search_ (or search manually):
# ```
#     cases.project.program.name in ["TARGET"] and cases.project.project_id 
#     in ["TARGET-AML"] and files.data_format in ["idat"] and 
#     files.experimental_strategy in ["Methylation Array"] and 
#     files.platform in ["illumina human methylation 450"]
# ```
# ```{note}
# Most (if not all?) methylation data available through GDC is publicly available, so no access to controlled data is necessary!
# ```
# - Once you have a selection of your samples, add them all to _cart_, which you can do by clicking on the _cart_ symbol in the header of the first column
# - Then go to the top right corner of the screen and click on _cart_
# - Download the _Manisfest_ file (under the _Download_ tab), as well as the _Sample Sheet_, and _Clinical_. Do not download the _Cart_ - files are too large.
# - To download the data, you will need the GDC Data Transfer Tool (command-line) that you can find [here](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
# - Download the software, unzip it, and place the .exe in a directory of easy access
# - In your terminal, go to the directory containing the .exe file and type the following:
# ```
# .\gdc-client download -m <path to the manisfest file you downloaded>
# ```
# - When the .idat files are downloaded, open the sample sheet file you downloaded previously. It should be a .tsv file readable in excell.
# - Save the sample sheet file as _samplesheet.csv_ in the directory where the data are located
# - Congrats! You are ready to process the files using ```methylprep``` (see chapter above for how to install methylprep and dependencies).
# - Like before, we will need to open methylprep on the terminal:
# ```
# source activate methylsuite
# ```
# - And run the command to process the files:
# ```
# python -m methylprep -v process --array_type 450k -d <path to data and samplesheet.csv> --all
# ```
# ```{note}
# If files are larger than your RAM memory, you will need to run in batches (check methylprep's documentation).
# Also, be mindful of \ and / differences if working in a windows machine through git bash. Some of methylprep code may break on Windows. Use a linux, Mac, or WSL if possible.
# ```
# - Fantastic, files are now processed and ready for ```methylcheck```!
