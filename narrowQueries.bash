#!/bin/bash
# enter Path of the file to read from
myCSV='./query.csv'
# enter Path to all the sequence zip files
sequencesPath='/media/evan/PROJECT/Data'
while read line
do
  echo $line
done < $myCSV
# unfinished, maybe delete
