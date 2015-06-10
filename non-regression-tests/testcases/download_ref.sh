#!/bin/bash

startdir=`readlink -f $(dirname $0)`

TARBALL="latu@scm.gforge.inria.fr:/home/groups/aster/nrt/jorek_rst_files.tgz"
FILE=`basename $TARBALL`
rsync --progress -av --delete-after "${TARBALL}" ${startdir}/
returncode=$?
if [ $returncode -ne 0 ]; then
cat <<EOF

###################################################################################
  Failed to download from INRIA gforge.
  Copy jorek_rst_files.tgz yourself into testcases directory
  and launch this script "$0" again to decompress the archive.
###################################################################################
EOF
else
  (cd $startdir; tar xvzf  $FILE)
fi