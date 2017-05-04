#!/usr/bin/env bash
# Loop through all f90 files in a directory and extract test data from them.
# Compile and run the test afterwards.

# Call this script from the JOREK directory, with a directory as argument
# If the second argument is -x or --xml output results to result.xml

outfile=`sed 's|/|_|g' <<< $1`
if [ -z "$outfile" ]; then
  exit 1
fi
rm -f "$outfile" "$outfile.{f90,mods,tests}"
has_setup=0
has_teardown=0
xml=""
junit=0
# A bit ugly this
if [ "$#" -eq 2 ]; then
  if [ "$2" == "-x" ] || [ "$2" == "--xml" ]; then
    xml="_xml"
  fi
  if [ "$2" == "-j" ] || [ "$2" == "--junit" ]; then
    junit=1
  fi
fi


# For each file in the requested directory
for file in $1/*.f90; do
  echo "Scanning $file"
  # Look for any lines named module *
  grep -o '^ *module [a-Z_0-9]*' $file | sed 's/module/use/g' >> $outfile.mods
  # Look for setup subroutine (global, can be only one)
  if grep -q 'subroutine \bsetup\b' $file; then
    has_setup=1
  fi
  # Look for any subroutines called setup_* or *_setup
  grep -o 'subroutine setup_[^ ]*\|subroutine [^ ]*_setup' $file >> $outfile.tests
  # Look for any subroutines called test_* or *_test
  # might contain duplicates due to end subroutine. Remove those later
  grep -o 'subroutine test_[^ ]*' $file >> $outfile.tests
  # Look for any subroutines called teardown_* or *_teardo
  grep -o 'subroutine teardown_[^ ]*\|subroutine [^ ]*_teardown' $file >> $outfile.tests

  # Look for teardown subroutine (global, can be only one)
  if grep -q 'subroutine \bteardown\b' $file; then
    has_teardown=1
  fi
done

# Write out the test driver
echo "program $outfile" > $outfile.f90
echo "use fruit" >> $outfile.f90 # use fruit_mpi if you want mpi support
sort $outfile.mods | uniq >> $outfile.f90
echo "implicit none" >> $outfile.f90
echo "call init_fruit$xml" >> $outfile.f90
if [ $has_setup -eq 1 ]; then
  echo "call setup" >> $outfile.f90
fi
uniq $outfile.tests | sed -e 's/subroutine \([^ ]*\)/call run_test_case(\1,"\1")/g' >> $outfile.f90
echo "call fruit_summary$xml" >> $outfile.f90
echo "call fruit_finalize" >> $outfile.f90
if [ $has_teardown -eq 1 ]; then
  echo "call teardown" >> $outfile.f90
fi
echo "end program $outfile" >> $outfile.f90

rm "$outfile.tests"
rm "$outfile.mods"

# build the test with a new invokation of make to get the dependencies right
make $outfile
if [ $? -eq 0 ]; then
  if [ "$junit" -eq 1 ]; then
    util/fruit2junit.sh $outfile
  else
    ./$outfile
  fi
fi
