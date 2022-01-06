#!/bin/bash

#compile documents from md to pdf or html

#katex base url
url='https://cdn.jsdelivr.net/npm/katex@0.15.1/dist/'
#flags
while getopts t:o:n: flags

do
    case "${flags}" in
    n) filename=${OPTARG};;
    t) format=${OPTARG};;
    o) outfile=${OPTARG};;
    esac
done

echo "$format";

if [ "$format" = "html" ]
then
pandoc $filename -f markdown+tex_math_dollars --katex=$url -t html5 -s -o ${outfile}".html"
elif [ "$format" = "pdf" ]
then
pandoc $filename -f markdown+tex_math_dollars -t pdf  --pdf-engine=xelatex  -s -o ${outfile}".pdf";
else
echo "wrong format"
fi

