#----------------------------------------------------------------------------

NEXE="exe2"
echo $NEXE

#----------------------------------------------------------------------------

OPT="-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -gen-interfaces -warn interfaces \
     -debug all -implicitnone -assume realloc_lhs -fstack-protector -assume protect_parens  "
OPT="-O3 -traceback"
echo $OPT

# ulimit -s unlimited

#----------------------------------------------------------------------------

MKL="-mkl=parallel"
#echo $MKL

#----------------------------------------------------------------------------

rm -f $NEXE
rm -f nohup.out

ifort $OPT -o $NEXE  src/Tools_mod.f90 \
                    src/Fem2D_mod.f90 \
                    src/Models_mod.f90 \
                    src/Stress_mod.f90 \
                    src/Equations_mod.f90 \
                    src/Newton_mod.f90 \
                    src/PostProcess_mod.f90 \
                    src/Regression_mod.f90 \
                    src/Continuation_mod.f90 \
                    Fem2D_prg.f90 	\
			   $MKL -qopenmp
 
#----------------------------------------------------------------------------
#rm -f *.out
rm -f OUTFE.OUT
rm -f *.mod
rm -f *.DAT
rm -f *.dat
rm -f *.plt
rm -f *.PLT
rm -f fort.*
rm -f *.SOL
rm -f *.FRS
rm -f *.STL
rm -f *genmod.f90
#----------------------------------------------------------------------------
