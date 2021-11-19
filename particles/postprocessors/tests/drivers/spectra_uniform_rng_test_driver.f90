! This program is the driver of the spectra_uniform_rng_test tests
program spectra_uniform_rng_test
use fruit
use mod_spectra_uniform_rng_test, only: run_fruit_spectra_uniform_rng_test
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the spectra uniform rng test basket
  call run_fruit_spectra_uniform_rng_test

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program spectra_uniform_rng_test
