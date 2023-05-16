#! /bin/bash
cat samplesheetRADsExample.csv | while read LINE
do

i=$(awk '{print $6}' <<< $LINE)
b=$(awk '{print $4}' <<< $LINE)
n=$(awk '{print $7}' <<< $LINE)
a=$(awk '{print $1}' <<< $LINE)
a2=$(awk '{print $5}' <<< $LINE)
a3=$(awk '{print $3}' <<< $LINE)

rename "s/$i\..{4}\.$b/$a3\.$a\.$n\.$a2/" *.fasta

done
