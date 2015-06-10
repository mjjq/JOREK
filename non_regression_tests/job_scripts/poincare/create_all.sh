for i in $(ls ../../testcases/ | grep tear); do 
   ./copy.sh tearing_limiter_199.job $i.job; 
done
