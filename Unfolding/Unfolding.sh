#!/bin/bash

# Run Python script
python run_notebooks.py

# Run first ROOT macro
root -l -b -q Unfolding.cpp

# Run second ROOT macro
root -l -b -q PlotUnfoldingCorrections.cpp