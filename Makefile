SHELL=/bin/bash

default: build

build:
	ar -r lib/libcoop.a typedef/*.o utils/*.o background/*.o firstorder/*.o nonlinear/*.o forecast/*.o mapio/*.o lib/*.o

clean:
	rm -f *.*~ Makefile~ include/*.*~ data/*.*~ data/\#*\# include/\#*\#
