#!/bin/bash

# compile the program
make distclean &> /dev/null
make &> /dev/null

#executing tha main using as parameters the file parametersL2.pot 
#which impose to use the L2 norm as the stopping critetion
./main -p parametersL2.pot
#output file: resultL2.dot

#executing tha main using as parameters the file parametersH1.pot 
#which impose to use the H1 norm as the stopping critetion
./main -p parametersH1.pot
#output file: resultH1.dot

exit
