prefix=/Data/gfi/users/local
gfortran example.f90 -I$prefix/include -L$prefix/lib -Wl,-R$prefix/lib -ldynfor -o example.x
