#include <ldefns.h>
extern "C" short Optimise_internalCVP(dimen n,long nfac,char** names,vector w_opt,dimen m,
                                    vector A,vector L,vector U,vector alpha,
                                    vector benchmark,vector Q,double gamma,vector initial,
                                    double delta,vector buy,vector sell,double kappa,long basket,
                                    long trades,int revise,int costs,double min_holding,double min_trade,
                                    int m_LS,int Fully_Invested,double Rmin,double Rmax,
                                    int m_Round,vector min_lot,vector size_lot,int* shake,dimen ncomp,vector Composites,
									double LSValue,dimen npiece,vector hpiece,vector pgrad);
int main()
{
    dimen n = 10;
    dimen m = 1;
    double x[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double c[] = {-1, -2, -3, -4, -5, -6, -17, -8, -9, -10};
    double A[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double L[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    double U[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double hess[] = {1,
                     0, 1,
                     0, 0, 1,
                     0, 0, 0, 1,
                     0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    short back = 2000;

    back = Optimise_internalCVP(n, -1, 0, x, m, A, L, U, c, 0, hess, 1, 0, -1, 0, 0, -1, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    for(int i=0;i<n;++i){
        printf("%d %f\n",i+1,x[i]);
    }
    printf("%f\n",-ddotvec(n,c,x));
}