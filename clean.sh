#!/bin/bash

if (( $# != 1 )); then
    echo "Illegal number of arguments. There should be 1 argument: case name."
    exit 1
fi

caseName=$1
rm -f ./pdf.txt *.out
rm -fr ./outputs/$1/cache ./outputs/$1/DOFSelection ./outputs/$1/PDFs

read -p "Delete also derivatives? (y/n)" -n 1 -r
echo    
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -fr ./outputs/$1/derivatives
fi
