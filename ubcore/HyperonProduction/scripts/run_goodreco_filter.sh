#!/bin/bash

#list of files to run filter over
input=/pnfs/uboone/scratch/users/cthorpe/hyperon/v08_00_00_46/all_fhc/make_events_numi_reco/files.list

#create output directory

#read input file

while IFS= read -r line
do
  echo "$line"

#strip off the "root" from the filename
common=${line%".root"}

#create output file name
outfilename="${common}_filtered.root"

echo Output file: $outfilename

lar -c run_HyperonGoodRecoFilter.fcl -s ${line} -o ${outfilename}

done < "$input"


