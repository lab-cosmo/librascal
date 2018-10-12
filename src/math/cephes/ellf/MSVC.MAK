#
# Microsoft Visual C++ makefile for ellf.exe
#
# Open a command line DOS window.
# Set directory paths by running the Microsoft batch file vcvars32.bat
# located at \program files\microsoft visual studio\vc98\bin\vcvars32.bat
# Then build the program with nmake.exe by the command
#   nmake -f msvc.mak
# To run the examples, type
#   ellf < ellf.que > temp.tmp
# and compare the results in temp.tmp with ellf.ans


ellf.exe: ellf.obj ellik.obj ellpj.obj ellpk.obj cmplx.obj const.obj polevl.obj mtherr.obj
	link ellf ellik ellpj ellpk cmplx polevl const mtherr

ellf.obj: ellf.c
	cl /c ellf.c

ellik.obj: ellik.c
	cl /c ellik.c

ellpj.obj: ellpj.c
	cl /c ellpj.c

ellpk.obj: ellpk.c mconf.h
	cl /c ellpk.c

cmplx.obj: cmplx.c mconf.h
	cl /c cmplx.c

polevl.obj: polevl.c mconf.h
	cl /c polevl.c

mtherr.obj: mtherr.c mconf.h
	cl /c mtherr.c

const.obj: const.c mconf.h
	cl /c const.c
