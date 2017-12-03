if (( $# != 1 ))
then
  echo "Usage: ./import.sh gromacsdir"
  exit 1
fi

GRO="$1"

cp "$GRO"/src/gmxlib/gmx_lapack/lapack_copyright COPYRIGHT

cp "$GRO"/src/gmxlib/gmx_lapack/lapack_limits.h .


sed 's|"types/simple.h"|"simple.h"|' "$GRO"/include/gmx_lapack.h |
  sed 's|F77_FUNC|PLUMED_BLAS_F77_FUNC|' |
  grep -v visibility.h |
  awk '{
         if($1=="#ifdef" && $2=="__cplusplus"){
           inside=1;
           a++;
           if(a==1){
             print "#include \"simple.h\""
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

cat << EOF > simple.h
#ifndef PLUMED_lapack_simple_h
#define PLUMED_lapack_simple_h
#define _simple_h
    /*! \brief Double precision accuracy */
#define PLUMED_GMX_DOUBLE_EPS   1.11022302E-16

    /*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define PLUMED_GMX_DOUBLE_MAX   1.79769312E+308

    /*! \brief Minimum double precision value */
#define PLUMED_GMX_DOUBLE_MIN   2.22507386E-308

    /*! \brief Single precision accuracy */
#define PLUMED_GMX_FLOAT_EPS    5.96046448E-08

    /*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define PLUMED_GMX_FLOAT_MAX    3.40282346E+38

    /*! \brief Minimum single precision value */
#define PLUMED_GMX_FLOAT_MIN    1.17549435E-38

#if defined(F77_NO_UNDERSCORE)
#define PLUMED_BLAS_F77_FUNC(lower,upper) lower
#else
#define PLUMED_BLAS_F77_FUNC(lower,upper) lower ## _
#endif

#define GMX_LIBGMX_EXPORT

#endif
EOF

{
echo "#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)"
for file in "$GRO"/src/gmxlib/gmx_lapack/*.c
do
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
  }' "$file" |
 sed 's|F77_FUNC|PLUMED_BLAS_F77_FUNC|' |
 sed 's|GMX_|PLUMED_GMX_|g' |
 sed 's|gmx_lapack|lapack|g' |
 sed 's|gmx_blas|blas/blas|g' |
 sed 's|<types/simple.h>|"simple.h"|'
done
echo "#endif"
} | sed 's/static //' > lapack.cpp

cd ../
./header.sh lapack


