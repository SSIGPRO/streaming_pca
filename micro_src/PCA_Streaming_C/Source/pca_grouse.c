#include "pca_streaming_c.h"

#define ALPHA 0.2

void pca_grouse(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_n_buff,    /* n sized vector buffer */
                FLOAT_T * p_m_buff)    /* m sized vector buffer */
{
    uint16_t n = pSrcU->numRows, m = pSrcU->numCols;
    FLOAT_T * px = pSrcx;
    FLOAT_T * py = p_m_buff;
    FLOAT_T * pp = p_n_buff;
    FLOAT_T * pr = pSrcx;
    FLOAT_T * pz = pSrcx;
    FLOAT_T * pd = pSrcx;


    FLOAT_T ynorm, pnorm, rnorm, znorm, theta, costheta;


    ext_vec_mat_mult(px, pSrcU, py); // y = x@U
    ARM_MAT_VEC_MULT(pSrcU, py, pp); // p = U@y
    ARM_SUB(px, pp, pr, n); // r = x - p


    // Normalize p
    ARM_DOT_PROD(pp, pp, n, &pnorm);
    ARM_SQRT(pnorm, &pnorm);
    ARM_SCALE(pp, 1.0/pnorm, pp, n);


    // Normalize r
    ARM_DOT_PROD(pr, pr, n, &rnorm);
    ARM_SQRT(rnorm, &rnorm);
    ARM_SCALE(pr, 1.0/rnorm, pr, n);

    // Evaluate zeta
#ifndef __EXT_MATH_DOUBLE
    theta = atanf((1.0-ALPHA)*rnorm/pnorm);
    costheta = cosf(theta);
    ARM_SCALE(pp, costheta, pp, n);
    ARM_SCALE(pr, sinf(theta), pr, n);
#else
    theta = atan((1.0-ALPHA)*rnorm/pnorm);
    costheta = cos(theta);
    ARM_SCALE(pp, costheta, pp, n);
    ARM_SCALE(pr, sin(theta), pr, n);
#endif
    ARM_ADD(pp, pr, pz, n); // z = cos(theta)/norm(p)*p + sin(theta)/norm(r)*r
    ARM_SCALE(pp, 1.0/costheta, pp, n); // restore pp

    // Normalize z
    ARM_DOT_PROD(pz, pz, n, &znorm);
    ARM_SQRT(znorm, &znorm);
    ARM_SCALE(pz, 1.0/znorm, pz, n);


    // Normalize y
    ARM_DOT_PROD(py, py, m, &ynorm);
    ARM_SQRT(ynorm, &ynorm);
    ARM_SCALE(py, 1.0/ynorm, py, m);


    ARM_SUB(pz, pp, pd, n); // d = z/norm(z)-p/norm(p)
    ext_add_outer_prod(pSrcU, pd, py, pSrcU); // U = U + d@y/norm(y)
}
