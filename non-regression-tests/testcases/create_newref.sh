#!/bin/bash

startdir=`readlink -f $(dirname $0)`

TARBALL="\$FORGEUSER@scm.gforge.inria.fr:/home/groups/aster/nrt/jorek_rst_files.tgz"
FILE=`basename $TARBALL`
DISTANTDIR=`dirname $TARBALL`

printf "\n  Creating the archive $FILE with the files:\n\n"
LIST="
ballooning_xpoint_303/jorek_restart.rst 
tearing_limiter_199/jorek_restart.rst 
tearing_limiter_303/jorek_restart.rst
h5_tearing_limiter_199/jorek00089_export.rst
h5_tearing_limiter_199/jorek00090_export.h5
"
tar cvzf new_$FILE $LIST
returncode=$?

if [ $returncode -eq 0 ]; then
mv new_$FILE $FILE
cat<<EOF

  To export this archive and put it as reference archive,
  please copy it to the shared directory yourself (change the login):

   rsync --progress -av "${FILE}" "${DISTANTDIR}"
EOF
fi