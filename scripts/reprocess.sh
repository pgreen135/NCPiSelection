#!/bin/bash

# Number of expected command-line arguments
num_expected=2

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./reprocess.sh OUTPUT_DIRECTORY NTUPLE_LIST_FILE"
  #exit 1
fi

output_dir=$1
ntuple_list_file=$2

if [ ! -f "$ntuple_list_file" ]; then
  echo "Ntuple list file \"${ntuple_list_file}\" not found"
  #exit 1
fi

if [ ! -d "${output_dir}" ]; then
  echo "Output directory \"${output_dir}\" not found"
  #exit 2
fi

# Loop over each line of the ntuple list file
while read line; do
  # Select lines that do not begin with a '#' character and contain at least
  # one non-whitespace character. These are assumed to be input file names
  if [[ ! $line = \#* ]] && [[ $line = *[^[:space:]]* ]]; then
    # Process the next input ntuple file
    input_file_name=$( echo $line | awk '{print $1}' )
    type_enum=$( echo $line | awk '{print $2}' )
    run_enum=$( echo $line | awk '{print $3}' )
    output_file_name="${output_dir}/xsec-$(basename ${input_file_name})"
    echo "PROCESSING ${input_file_name} --> ${output_file_name}"
    time ../build/ProcessNTuple ${input_file_name} ${output_file_name} ${type_enum} ${run_enum}
  fi
done < "${ntuple_list_file}"