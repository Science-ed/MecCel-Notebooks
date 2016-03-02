#!/bin/bash
echo "Converting Jupyter python scripts to python..."
for file in *.py.html
do
    script=$(echo $file |sed -e "s/\.html//")
    echo "Converting $file to $script..."
    sed -e "s/get_/#get_/" $file > $script
    rm -r $file
done
