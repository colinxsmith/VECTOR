PATH=$PATH:/usr/bin
mkdir lapack
cp ../safeqp/lapack-release/SRC/ieeeck.f lapack
cp ../safeqp/lapack-release/SRC/ilaenv.f lapack
cp ../safeqp/lapack-release/SRC/dsptrf.f lapack
cp ../safeqp/lapack-release/SRC/dsptrs.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/drotg.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/idamax.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dswap.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dspr.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dscal.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dger.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dgemv.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dgemm.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/daxpy.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/dcopy.f lapack
cp ../safeqp/lapack-release/BLAS/SRC/ddot.f lapack
cp $(grep -i BUNch ../safeqp/lapack-release/SRC/d*.f | sed "s/:.*//") lapack
cp $(grep -i dsytrf ../safeqp/lapack-release/SRC/d*.f | sed "s/:.*//") lapack
cd lapack
for i in *.f
do
f2c $i
done
rm -rf *.f