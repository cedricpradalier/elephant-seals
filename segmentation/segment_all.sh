#!/bin/bash

for f in ../15734_inter.labels
do
	dname=`basename $f .labels`
	echo $dname
	rm -rf output output.$dname
	./bin/segment_labels $f ${f/labels/csv}
	mv output output.$dname
done

