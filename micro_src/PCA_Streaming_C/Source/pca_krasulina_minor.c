#include "pca_streaming_c.h"

#define RATE(tau) 0.1

void pca_krasulina_minor(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                uint32_t tau)          /* step count */
{
    uint16_t n = pSrcU->numRows, m = pSrcU->numCols;
    FLOAT_T * px = pSrcx;
    FLOAT_T * py = p_mxm_buff;
    FLOAT_T * pp = p_nxm_buff;
    FLOAT_T * pr = pSrcx;
    FLOAT_T * pS = p_nxm_buff;
    FLOAT_T * pR = p_mxm_buff;

    ARM_MATRIX_INSTANCE S, R;
    ARM_MAT_INIT(&S, n, m, pS);
    ARM_MAT_INIT(&R, m, m, pR);

    ext_vec_mat_mult(px, pSrcU, py); // y = x@U
    ARM_MAT_VEC_MULT(pSrcU, py, pp); // p = U@y
    ARM_SUB(px, pp, pr, n); // r = x - p
    if(RATE(tau) != 1.0) ARM_SCALE(py, RATE(tau), py, m); // y = rate(tau)*y
    ext_sub_outer_prod(pSrcU, pr, py, &S); // S = U + r.T@y
    ext_mat_qr(&S, pSrcU, &R); // qr(S) -> U, R
}

void pca_krasulina_minor_lowmem(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                uint32_t tau)          /* step count */
{
    uint16_t n = pSrcU->numRows, m = pSrcU->numCols;
    FLOAT_T * pU = pSrcU->pData;
    FLOAT_T * px = pSrcx;
    FLOAT_T * py = p_mxm_buff;
    FLOAT_T * pr = pSrcx;
    FLOAT_T * pS = pU;
    FLOAT_T * pR = p_mxm_buff;

    ARM_MATRIX_INSTANCE S, R;
    ARM_MAT_INIT(&S, n, m, pS);
    ARM_MAT_INIT(&R, m, m, pR);

    ext_vec_mat_mult(px, pSrcU, py); // y = x@U
    ext_mat_vec_mult_sub(pSrcU, py, pr); // r = x-U@y (r and x -> same address)
    if(RATE(tau) != 1.0) ARM_SCALE(py, RATE(tau), py, m); // y = rate(tau)*y
    ext_sub_outer_prod(pSrcU, pr, py, &S); // S = U + r@y.T
    ext_mat_qr_overwrite(&S, &R, px); /* qr(S) -> U, R (S and U same address,
                                       * use x as buffer) */
}
