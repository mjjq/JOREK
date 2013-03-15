if [ $# -ne 2 ] || [ ! -d "$1" ] ; then
  printf "\nUsage: $0 <Source Directory> <Target Directory>\n\n"
  printf "Copy data from <Source directory> to <Target Directory>\n"
  printf "in order to restart the Non regression test in the\n"
  printf "directory <Target Directory>.\n\n"
  exit 0
fi

mkdir $2
cp -r $1/jorek_rst* $1/jorek_equil* $2/
cp -r $1/??[0-9][0-9][0-9] $1/out_* $2/
cd $2
if [ -f jorek_rst3.rst ]; then
   cp jorek_rst3.rst jorek_restart.rst
   rm out_loop[4-9]
else if [ -f jorek_rst2.rst ]; then
   cp jorek_rst2.rst jorek_restart.rst
   rm out_loop[3-9]
else if [ -f jorek_rst1.rst ]; then
   cp jorek_rst1.rst jorek_restart.rst
   rm out_loop[2-9]
fi
fi
fi