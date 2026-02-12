#!/bin/bash

# This script generates the ROOT dictionary for the MyStruct class and compiles it into a shared library.

rootcling -f My_DictOutput.cxx -s My_libMyDict.so My_LinkDef.h
g++ -fPIC -shared My_DictOutput.cxx $(root-config --cflags --libs) -o My_libMyDict.so