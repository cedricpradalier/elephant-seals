#!/bin/bash

for f in ../2*.labels
do
	dname=`basename $f .labels`
	if test -d output.$dname
	then
		echo "Skipping $dname"
	else
		echo "Processing $dname"
		./bin/segment_labels $f ${f/labels/csv}
		mv output output.$dname
	fi
done

