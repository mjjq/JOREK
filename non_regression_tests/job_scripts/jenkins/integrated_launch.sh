startdir=`readlink -f $(dirname $0)`
nonregdir=`readlink -f ${startdir}/../..` 
cd $nonregdir || exit 1
export JOREK_HOST=jenkins
./launch_all.sh 