#!/bin/bash
printf "\n\n===========================\n"
printf "===========================\n"
printf "===========================\n"
printf "Run all pFUnit tests\n"
printf " --------------------------\n"


printf "\n\n--------------------------\n"
printf "test_activation\n\n"
/pFUnitTests/test_activation

printf "\n\n --------------------------\n"
printf "test_shared_thermo\n\n"
pFUnitTests/test_shared_thermo
