program test

  USE netcdfmod
  
  implicit none
  
  character(len = 100) :: file_output
  
  print*, 'Enter file name:'
  read*, file_output
  
  call read_dim (file_output)
  
end program
