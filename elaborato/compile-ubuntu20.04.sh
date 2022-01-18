#!/bin/bash

#compile documents from md to pdf or html with pandoc

#katex base url
url='https://cdn.jsdelivr.net/npm/katex@0.15.1/dist/'
disable_float='disable_float.tex'
#flags
while getopts t:o:n: flags

do
    case "${flags}" in
    n) filename=${OPTARG};;
    t) format=${OPTARG};;
    o) outfile=${OPTARG};;
    esac
done

#echo "$format";

if [ "$format" = "html" ]
then
pandoc $filename -f markdown+tex_math_dollars --filter pandoc-fignos --citeproc --katex=$url -t html5 -s -o ${outfile}".html"
elif [ "$format" = "pdf" ]
then
pandoc $filename -f markdown+tex_math_dollars -H $disable_float --filter pandoc-fignos --citeproc -t pdf  --pdf-engine=pdflatex  -s -o ${outfile}".pdf";
else
echo "wrong format"
fi

