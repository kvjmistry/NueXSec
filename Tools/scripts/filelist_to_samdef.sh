#!/bin/bash

# This script will take a fileist and convert it to a samdef
# first argument is the sam definition name
# second argument is the filelist

samdef=$1
filelist=$2

files=`awk '{print $1}' $filelist | paste -s -d, -`
files=$files

command1="samweb create-definition "
command2=$samdef
command3=" \"file_name $files\""
command=$command1$command2$command3

echo $command

eval "$command"
