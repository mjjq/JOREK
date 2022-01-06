! This program is the driver of the vertices tests
program vertices_test
use fruit
use mod_vertices_test, only: run_fruit_vertices
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the vertices test basket
  call run_fruit_vertices

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program vertices_test
