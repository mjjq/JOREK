! This program is the driver of the spectra_deterministic_test tests
program spectra_deterministic_test
use fruit
use mod_spectra_deterministic_test, only: run_fruit_spectra_deterministic_test
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the spectra deterministic test basket
  call run_fruit_spectra_deterministic_test

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program spectra_deterministic_test
