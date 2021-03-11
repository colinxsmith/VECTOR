
//Make 64 bit version on Windows in safeqp64 with
//cl -D__SYSNT__ -EHsc -I . safeqp.lib ..\VECTOR\LP.cpp
#include <valarray>
#include <ldefns.h>
extern "C" short Optimise_internalCVP(dimen n, long nfac, char **names, vector w_opt, dimen m,
                                      vector A, vector L, vector U, vector alpha,
                                      vector benchmark, vector Q, double gamma, vector initial,
                                      double delta, vector buy, vector sell, double kappa, long basket,
                                      long trades, int revise, int costs, double min_holding, double min_trade,
                                      int m_LS, int Fully_Invested, double Rmin, double Rmax,
                                      int m_Round, vector min_lot, vector size_lot, int *shake, dimen ncomp, vector Composites,
                                      double LSValue, dimen npiece, vector hpiece, vector pgrad);
int main()
{
    dimen n = 10;
    dimen m = 2;
    double x[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double c[] = {-1, -2, -3, -4, -5, -6, -17, -8, -9, -10};
    double A[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  0, 0, 1, 1, 1, 0, 0, 0, 0, 0};
    double L[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.1};
    double U[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5};
    dmx_transpose(n, m, A, A);
    double hess[] = {0.07622384475840693,
                     -0.0016365991417207626,
                     0.08604333317714313,
                     -0.011316860823247121,
                     -0.030434772808639987,
                     0.11211083919582854,
                     -0.005442735249764186,
                     0.004464673502705685,
                     0.02884916845203872,
                     0.07030671358796253,
                     0.011709823440838929,
                     0.02086299454025381,
                     -0.015170796975208789,
                     -0.0005067292359139386,
                     0.059822876920856805,
                     0.030315371151201392,
                     0.0034082127351833247,
                     -0.022933127136383458,
                     0.0016861298409444059,
                     0.012600705278736524,
                     0.08316673826220947,
                     0.0115607066361883,
                     -0.0034806816621613668,
                     -0.0010519351552628065,
                     -0.0028576346296797506,
                     0.008573189337515441,
                     -0.025115408356197216,
                     0.08293672602298946,
                     0.014699977532219966,
                     -0.019735088940840084,
                     0.05551477101220931,
                     -0.019246450916202917,
                     -0.01412419584221336,
                     -0.0025011076826836065,
                     -0.013239288762594226,
                     0.0907984789770602,
                     -0.013526556333058826,
                     0.0098303768770443,
                     -0.042483694444829995,
                     -0.035685449490480525,
                     0.015581944899869915,
                     -0.0008483921768689395,
                     -0.006279985059017501,
                     -0.00375529364246191,
                     0.10266907898364636,
                     0.03999356589115988,
                     0.045353250040677695,
                     -0.0019274282309346968,
                     0.010712719440722496,
                     0.024816489058402724,
                     0.03416835576827826,
                     -0.001184721225989449,
                     -0.011735959703553567,
                     -0.04135384837511391,
                     0.12585651441173307};
    dscalvec(n * (n + 1) / 2, 1e3, hess);
    short back = 2000;
    double gamma = 0.5;
    back = Optimise_internalCVP(n, -1, 0, x, m, A, L, U, c, 0, hess, gamma, 0, -1, 0, 0, -1, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    for (int i = 0; i < n; ++i)
    {
        printf("%d %f\n", i + 1, x[i]);
    }
    std::valarray<double> implied(n);
    dsmxmulv(n, hess, x, &implied[0]);
    if (gamma == 1)
        printf("back=%d U=%f\n", back, -ddotvec(n, c, x));
    else
        printf("back=%d U=%f\n", back, -ddotvec(n, c, x) * gamma / (1 - gamma) + 0.5 * ddotvec(n, x, &implied[0]));
}