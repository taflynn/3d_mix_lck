gp_lck: gp_lck.o grid.o init.o rhs.o time.o
	gfortran *.o -O3 -march=native -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp -o gp_lck
grid.mod: grid.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp
init.mod: init.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp
init.mod: rhs.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp
init.mod: time.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp
gp_lck.o: gp_lck.f03
	gfortran -c FFTW3.f03 grid.f03 init.f03 rhs.f03 time.f03 gp_lck.f03 -L/usr/lib/x86_64-linux-gnu/ -I/usr/include -lfftw3_threads -lfftw3 -lm -fopenmp

clean:
	rm -r ./*.o ./*.mod
