#include "pca_streaming_c.h"

#define BETA 0.99

void pca_past(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
    ARM_MATRIX_INSTANCE * pSrcPT,      /* m x m matrix P  */
                FLOAT_T * p_n_buff,    /* n sized vector buffer */
                FLOAT_T * p_m_buff)    /* m sized vector buffer */
{
    uint16_t m = pSrcU->numCols;

    FLOAT_T * px = pSrcx;
    FLOAT_T * py = p_m_buff;
    FLOAT_T * pr = pSrcx;
    FLOAT_T * ph = p_n_buff;
    FLOAT_T * pg = p_m_buff;

    FLOAT_T yhdot;

    ext_vec_mat_mult(px, pSrcU, py); // y = x@U
    ext_mat_vec_mult_sub(pSrcU, py, pr); // r = x - U@y (r shares memory with x)
    ARM_MAT_VEC_MULT(pSrcPT, py, ph); // P.T = Q.T@y
    ARM_DOT_PROD(py, ph, m, &yhdot); // yhdot = y.T@h
    ARM_SCALE(ph, 1.0/(yhdot+BETA), pg, m); // g = h/(yhdot+BETA)

    ext_sub_outer_prod(pSrcPT, ph, pg, pSrcPT); // P.T = P.T - h@g.T

    ARM_MAT_SCALE(pSrcPT, 1.0/BETA, pSrcPT); // P.T = P.T/BETA

    ext_add_outer_prod(pSrcU, pr, pg, pSrcU); // U = U + r@g.T
}
