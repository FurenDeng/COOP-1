TabZ1: ellipse_collapse.f90 make_zvir1_table.f90
	gfortran -cpp --free-line-length-none -fmax-identifier-length=63 -O3 -fimplicit-none -o $@ ellipse_collapse.f90 make_zvir1_table.f90
.PHONY : clean
clean:
	rm -f *.o *.*~ *.mod \#* *.a makefile~ TabZ1 *.dat
