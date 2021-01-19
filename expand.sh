PATH=$PATH:/usr/bin
cp ../safeqp/lapack-release/SRC/dsptrf.f .
cp ../safeqp/lapack-release/SRC/dsptrs.f .
cp ../safeqp/lapack-release/BLAS/SRC/idamax.f .
cp ../safeqp/lapack-release/BLAS/SRC/dswap.f .
cp ../safeqp/lapack-release/BLAS/SRC/dspr.f .
cp ../safeqp/lapack-release/BLAS/SRC/dscal.f .
cp ../safeqp/lapack-release/BLAS/SRC/dger.f .
cp ../safeqp/lapack-release/BLAS/SRC/dgemv.f .
cp ../safeqp/lapack-release/BLAS/SRC/dgemm.f .
for i in *.f
do
f2c $i
done
rm -rf *.f