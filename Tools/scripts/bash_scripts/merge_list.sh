# script to split a filelist of many root files, hadd them together and then merge those files to one file
# USAGE: source merge_list.sh <input list of files> <output root file name>

split -dl 500 $1 filelist_to_merge

for file in $(ls | grep filelist_to_merge); do
	echo "Number of lines to merge: $(cat $file | wc -l)"
	echo "hadd "${file}.root" @$file"
	hadd "${file}.root" @$file
	echo "rm $file"
	rm $file
	echo
done

echo "now merging files into one..."

echo "hadd $2 filelist_to_merge*.root"
hadd $2 filelist_to_merge*.root

# Cleanup
echo "rm filelist_to_merge*.root"

