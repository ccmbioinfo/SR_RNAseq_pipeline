#!/bin/bash

#this script is used to set up the required python packages for the clinical pipeline
#make sure you have a vitural environemnt (python_venv) set up in the project directory
#if not, set up by running command in the project directory: python -m venv ./python_venv

module load python/3.12.4

if [ -d "./python_venv/bin" ]; then
  ### Take action if $DIR exists ###
  echo "venv exist, python ready"
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "Creating python vitrual environment"
  python -m venv ./python_venv

fi

source python_venv/bin/activate

pip install -r requirements.txt
