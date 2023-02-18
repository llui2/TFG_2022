gfortran -c spin_model.f
gfortran -c metropolis.f
chmod +x spin_model.o
chmod +x metropolis.o
gfortran spin_model.o metropolis.o
./a.out

gfortran -c performance.f
chmod +x performance.o
gfortran spin_model.o performance.o
./a.out