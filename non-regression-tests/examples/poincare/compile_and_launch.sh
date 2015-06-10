. env.sh
LIST="xpoint_equil xpoint_grid polar_grid tearing_limiter_199 tearing_limiter_303 ballooning_xpoint_303"
for case in ${LIST}; do
  echo $case
  ../run_test.sh -r $case
  llsubmit ${case}.cmd
done
