gp_lck: gp_lck.o grid.o init.o rhs.o time.o
	gfortran *.o -O3 -march=native -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -L${HOME}/json-fortran-8.3.0/lib -ljsonfortran -I${HOME}/json-fortran-8.3.0/lib -o gp_lck
grid.mod: grid.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -L${HOME}/json-fortran-8.3.0/lib -ljsonfortran -I${HOME}/json-fortran-8.3.0/lib
init.mod: init.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -L${HOME}/json-fortran-8.3.0/lib -ljsonfortran -I${HOME}/json-fortran-8.3.0/lib
rhs.mod: rhs.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -L${HOME}/json-fortran-8.3.0/lib -ljsonfortran -I${HOME}/json-fortran-8.3.0/lib
time.mod: time.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -L${HOME}/json-fortran-8.3.0/lib -ljsonfortran -I${HOME}/json-fortran-8.3.0/lib
gp_lck.o: gp_lck.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -L${HOME}/json-fortran-8.3.0/lib -ljsonfortran -I${HOME}/json-fortran-8.3.0/lib

clean:
	rm -r ./*.o ./*.mod
