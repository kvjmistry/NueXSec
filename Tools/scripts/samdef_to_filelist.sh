# This script will take a sam definition and convert it into a file list. 
# The name of the output list will be the samdef name + _file.list
# USAGE: source samdef_to_file.list <the samdef>

echo "SAMDEF: $1"
samdef=$1

for files in `samweb list-files defname: ${samdef}`
do
  # Locate the file
  newfile=`samweb locate-file "${files}"`

  # get rid of enstore in the name
  modnewfile=${newfile#"enstore:"}

  # get rid of anything after the bracket
  modnewfile=${modnewfile%(*}

  # create the full path
  modnewfile="${modnewfile}/${files}"
  echo $modnewfile

  # Put into the new file
  echo ${modnewfile} >> ${samdef}_file.list

done

echo "${samdef}_file.list CREATED"
