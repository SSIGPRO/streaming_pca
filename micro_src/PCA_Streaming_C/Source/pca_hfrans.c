#include "pca_streaming_c.h"

#define MU 1.0

void pca_hfrans_minor(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_n_buff,    /* n sized vector buffer */
                FLOAT_T * p_m_buff)    /* m sized vector buffer */
{
    uint16_t n = pSrcU->numRows, m = pSrcU->numCols;
    FLOAT_T * px = pSrcx;
    FLOAT_T * py = p_m_buff;
    FLOAT_T * pp = p_n_buff;
    FLOAT_T * pu = p_n_buff;
    FLOAT_T * pv = p_m_buff;

    ARM_MATRIX_INSTANCE U, V;
    ARM_MAT_INIT(&U, 1, n, pu);
    ARM_MAT_INIT(&V, 1, m, pv);

    FLOAT_T xnorm2, ynorm2, beta, delta, rho, tau, unorm;

    ext_vec_mat_mult(px, pSrcU, py); // y = x@U
    ARM_MAT_VEC_MULT(pSrcU, py, pp); // p = U@y
    ARM_DOT_PROD(px, px, n, &xnorm2); // xnorm2 = px.T@px
    ARM_DOT_PROD(py, py, m, &ynorm2); // ynorm2 = py.T@py
    beta = MU/(2*xnorm2);
    delta = 4*beta*(1.0-beta*xnorm2)*ynorm2;
    ARM_SQRT(1.0-delta, &rho); // rho = sqrt(1-delta)
    tau = (1.0/rho-1.0)/ynorm2;
    ARM_SCALE(px, 2*(1.0+tau*ynorm2), px, n); // x = 2*(1+tau*ynorm2)*x
    ARM_SCALE(pp, tau/beta, pp, n); // p = tau/beta*p
    ARM_SUB(px, pp, pu, n); // u = x - p
    ARM_DOT_PROD(pu, pu, n, &unorm); // unorm^2 = u.T@u
    ARM_SQRT(unorm, &unorm); // unorm = sqrt(unorm^2)
    ARM_SCALE(pu, 1.0/unorm, pu, n); // u = u/unorm
    ARM_MAT_MULT(&U, pSrcU, &V); // v = u@U
    ARM_SCALE(pv, 2, pv, m); // v = 2v
    ext_sub_outer_prod(pSrcU, pu, pv, pSrcU); // U = U - 2u@v.T
}
