if (( $# != 1 ))
then
  echo "Usage: ./import.sh gromacsdir"
  exit 1
fi

GRO="$1"

cp "$GRO"/src/gromacs/linearalgebra/gmx_lapack/lapack_copyright COPYRIGHT

cp "$GRO"/src/gromacs/linearalgebra/gmx_lapack/lapack_limits.h .


sed 's|"gromacs/utility/real.h"|"real.h"|' "$GRO"/src/gromacs/linearalgebra/gmx_lapack.h |
  sed 's|F77_FUNC|PLUMED_BLAS_F77_FUNC|' |
  grep -v config.h |
  grep -v gromacs/utility/basedefinitions.h |
  awk '{
         if($1=="#ifdef" && $2=="__cplusplus"){
           inside=1;
           a++;
           if(a==1){
             print "#ifndef __PLUMED_LAPACK_RETURNS_FLOAT"
             print "#define __PLUMED_LAPACK_RETURNS_FLOAT float"
             print "#endif"
             print "#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)"
             print "#include \"def_internal.h\""
             print "namespace PLMD{"
             print "namespace lapack{"
             print "#else"
             print "#include \"def_external.h\""
             print "extern \"C\"{"
             print "#endif"
           }
           if(a==2){
             print "}"
             print "#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)"
             print "}"
             print "#endif"
           }
         }
         if(!inside){
           if(NF==1 && $1=="float") print "__PLUMED_LAPACK_RETURNS_FLOAT"
           else print
         }
         if(inside && $1=="#endif") inside=0;
       }' > lapack.h

grep PLUMED_BLAS_F77_FUNC lapack.h  | sed 's/(/ /' | sed 's/,/ /' | sed 's/)/ /' |
  awk '{print "/** \\ingroup internal-lapack */"; print "#define plumed_lapack_"$2" PLMD::lapack::PLUMED_BLAS_F77_FUNC("$2","$3")"}' > def_internal.h
grep PLUMED_BLAS_F77_FUNC lapack.h  | sed 's/(/ /' | sed 's/,/ /' | sed 's/)/ /' |
  awk '{print "#define plumed_lapack_"$2" PLUMED_BLAS_F77_FUNC("$2","$3")"}' > def_external.h

cat << EOF > real.h
#ifndef PLUMED_lapack_real_h
#define PLUMED_lapack_real_h
    /*! \brief Double precision accuracy */
#define PLUMED_GMX_DOUBLE_EPS   2.2204460492503131e-16

    /*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define PLUMED_GMX_DOUBLE_MAX   1.7976931348623157e+308

    /*! \brief Minimum double precision value */
#define PLUMED_GMX_DOUBLE_MIN   2.2250738585072014e-308

    /*! \brief Single precision accuracy */
#define PLUMED_GMX_FLOAT_EPS    1.19209290e-07F

    /*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define PLUMED_GMX_FLOAT_MAX    3.40282346E+38F

    /*! \brief Minimum single precision value */
#define PLUMED_GMX_FLOAT_MIN    1.175494351E-38F

#if defined(F77_NO_UNDERSCORE)
#define PLUMED_BLAS_F77_FUNC(lower,upper) lower
#else
#define PLUMED_BLAS_F77_FUNC(lower,upper) lower ## _
#endif

#endif
EOF

{
echo "#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)"
for file in "$GRO"/src/gromacs/linearalgebra/gmx_lapack/*.cpp
do
  cat "$file" |
  sed 's|"gromacs/utility/real.h"|"real.h"|' |
  sed 's|"../gmx_lapack.h"|"lapack.h"|' |
  sed 's|"../gmx_blas.h"|"blas/blas.h"|' |
  sed 's| gmx_unused | |' |
  awk '{
    if(match($0,"F77_FUNC") && !done){
      print "#include \"blas/blas.h\""
      print "namespace PLMD{"
      print "namespace lapack{"
      print "using namespace blas;"
      done=1;
    }
    if(NR>1)print save;
    save=$0;
  }END{
    print save
    print "}"
    print "}"
  }' |
 sed 's|F77_FUNC|PLUMED_BLAS_F77_FUNC|' |
 sed 's|GMX_|PLUMED_GMX_|g' 
done
echo "#endif"
# removing static seems to be still necessary
# even though the problem was reported upstream https://gerrit.gromacs.org/#/c/4651/
} | sed 's/static //' > lapack.cpp

cd ../
./header.sh lapack


