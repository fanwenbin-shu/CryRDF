

make:
	ifort -c module.f90 -mkl
	ifort -c main.f90
	ifort -o CryRDF main.o module.o -mkl
	time ./CryRDF input

clean:
	rm -f main.o *.mod *.out CryRDF
