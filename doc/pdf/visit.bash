#!/bin/bash

# For serial FronTier VTK output, this script generates a .visit file.
# It must be run in the same directory as the FronTier output. Its one argument
# is the run name for the collection of files to use.

# Adjust this line for 2d/3d files
if [ -z $1 ]; then
        echo Usage: $(basename $0) run_name
        exit -1
fi
INTFC=""
FILES=$(ls | grep $1-vtk.ts)
for i in ${FILES}; do
        if [ -z ${INTFC} ]; then
                INTFC=$(ls $i | grep intfc.vtk)
        elif [ ! ${INTFC} = $(ls $i | grep intfc.vtk) ]; then
                echo INCONSISTENT VTK FILES
                exit -1
        fi
done
echo ${FILES} | sed -e "s/ /\/${INTFC}\n/g" -e "s/$/\/${INTFC}/" > $1.visit
echo Successfully created $1.visit
