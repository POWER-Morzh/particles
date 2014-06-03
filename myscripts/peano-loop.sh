#!/bin/bash

 

for((f = 1; f<=2; f=$f*2))
do
	for((m = 7; m <=10; m=$m+1))
	do
		(cd ..; make clean; make MANTISSA=$m MUL_FACTOR=$f)
		(cd ../output; ../bin/peano-Particles-F$f-M$m pit video 250 1000000 0.05 0.1)
	done
done

