startdir=`readlink -f $(dirname $0)`
nonregdir=`readlink -f ${startdir}/../..` 
cd $nonregdir || exit 1
rm -f testcases/*/*tgz 2>/dev/null
cp ../Make.inc/Makefile.jenkins_linux_gnu.inc ../Makefile.inc || exit 1
./get_all_data.sh || exit 1
export compilethreads=2
export JOREK_HOST=jenkins
./compile_all.sh || exit 1