#ifndef EXT_MATH_C
#define EXT_MATH_C



/* CMSIS DSP library is required **************************/

#include "arm_math.h"



/* Redefinition of CMSIS functions (32/64 bits selection) ******/

#ifndef __EXT_MATH_DOUBLE
    #define FLOAT_T float32_t
    #define ARM_MATRIX_INSTANCE arm_matrix_instance_f32
    #define ARM_MAT_INIT arm_mat_init_f32
    #define ARM_MAT_TRANS arm_mat_trans_f32
    #define ARM_MAT_MULT arm_mat_mult_f32
    #define ARM_MAT_VEC_MULT arm_mat_vec_mult_f32
    #define ARM_MAT_SCALE arm_mat_scale_f32
    #define ARM_MAT_ADD arm_mat_add_f32
    #define ARM_MAT_SUB arm_mat_sub_f32
    #define ARM_MAT_CHOLESKY arm_mat_cholesky_f32
    #define ARM_COPY arm_copy_f32
    #define ARM_ADD arm_add_f32
    #define ARM_SUB arm_sub_f32
    #define ARM_SCALE arm_scale_f32
    #define ARM_DOT_PROD arm_dot_prod_f32
    #define ARM_SQRT arm_sqrt_f32
    #define ARM_FILL arm_fill_f32
    #define ARM_MULT arm_mult_f32
#else
    #define FLOAT_T float64_t
    #define ARM_MATRIX_INSTANCE arm_matrix_instance_f64
    #define ARM_MAT_INIT arm_mat_init_f64
    #define ARM_MAT_TRANS arm_mat_trans_f64
    #define ARM_MAT_MULT arm_mat_mult_f64
    #define ARM_MAT_VEC_MULT arm_mat_vec_mult_f64
    #define ARM_MAT_SCALE arm_mat_scale_f64
    #define ARM_MAT_ADD arm_mat_add_f64
    #define ARM_MAT_SUB arm_mat_sub_f64
    #define ARM_MAT_CHOLESKY arm_mat_cholesky_f64
    #define ARM_COPY arm_copy_f64
    #define ARM_ADD arm_add_f64
    #define ARM_SUB arm_sub_f64
    #define ARM_SCALE arm_scale_f64
    #define ARM_DOT_PROD arm_dot_prod_f64
    #define ARM_SQRT arm_sqrt_f64
    #define ARM_FILL arm_fill_f64
    #define ARM_MULT arm_mult_f64
#endif



/* Macros *************************************************/

#define SIGN(x) ((x<0)?-1:1)
#define MIN(a, b) ((a<b)?a:b)



/* Matrix print and scan **********************************/

// Global pointers to print/scan functions
void (*_print_stream)(char*, unsigned int);
void (*_scan_stream)(char*, unsigned int);

// Function setters
void ext_set_print(void (*print_stream)(char*, unsigned int));
void ext_set_scan(void (*scan_stream)(char*, unsigned int));

// Actual print and scan
void ext_print(const FLOAT_T * pSrcv, uint32_t length);
void ext_scan(FLOAT_T * pSrcv, uint32_t length);
void ext_mat_print(const ARM_MATRIX_INSTANCE * pSrcA);
void ext_mat_scan(ARM_MATRIX_INSTANCE * pDstA);



/* Matrix fill functions **********************************/

// Fill with zeros
void ext_mat_zeros(ARM_MATRIX_INSTANCE * pDstA);

// Fill with zeros
void ext_mat_ones(ARM_MATRIX_INSTANCE * pDstA);

// Put ones on the diagonal
void ext_mat_eye(ARM_MATRIX_INSTANCE * pDstA);

// Fill lower triangular part with zeros
void ext_mat_triu(ARM_MATRIX_INSTANCE * pDstA);

// Fill upper triangular part with zeros
void ext_mat_tril(ARM_MATRIX_INSTANCE * pDstA);



/* Arithmetic operators ***********************************/

// Make X@Y (optimized with register blocking)
void ext_mat_mult(
    ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix */
    ARM_MATRIX_INSTANCE * pDstC); /* n x l output matrix */

// Make X@Y (optimized with register blocking, column first)
void ext_mat_mult2(
    ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix */
    ARM_MATRIX_INSTANCE * pDstC); /* n x l output matrix */

// Make X@Y (optimized with loop unroll)
void ext_mat_mult_simple(
    ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix */
    ARM_MATRIX_INSTANCE * pDstC); /* n x l output matrix */

// Make X.T@X product
void ext_mat_xtx_mult(
    const ARM_MATRIX_INSTANCE * pSrcX,
          ARM_MATRIX_INSTANCE * pDst);

// Make X@Y.T product, where Y is considered lower triangular
void ext_mat_xyt_mult_tril(
    const arm_matrix_instance_f32 * pSrcA,
    const arm_matrix_instance_f32 * pSrcB,
          arm_matrix_instance_f32 * pDst);

/* Make X@Y.T product, where Y is considered lower triangular; result overwrites
 * X (lower memory requirements) */
void ext_mat_xyt_mult_tril_overwriteA(
          ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix, result will be 
                                        * saved here (overwrite) */
    const ARM_MATRIX_INSTANCE * pSrcB, /* m x m square input matrix */
                      FLOAT_T * buff); /* m sized buffer */

/* Make X@Y product, where Y is square and the result overwrites X (lower memory
 * requirements) */
void ext_mat_mult_overwriteA(
    ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix, result will be stored
                                  * here */
    ARM_MATRIX_INSTANCE * pSrcB, /* m x m square input matrix */
                FLOAT_T * buff); /* m sized buffer */

void ext_mat_mult_overwriteB(
    ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix, result will be stored 
                                   * here (NOTE works iff n >= m !!!) */
                FLOAT_T * buff);  /* m sized buffer */

// Vector outer product
void ext_outer_prod(
                FLOAT_T * pSrcA,
                FLOAT_T * pSrcB,
    ARM_MATRIX_INSTANCE * pDstC);

// Matrix-vector multiplication
void ext_mat_vec_mult(
    ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix */
                FLOAT_T * pSrcB,  /* m x l input matrix */
                FLOAT_T * pDstC); /* n x l output matrix */

/* Matrix-vector multiplication and then subtraction from destionation vector
 * (r = r - Av) */
void ext_mat_vec_mult_sub(
    const arm_matrix_instance_f32 * pSrcMat,
                    const FLOAT_T * pVec,
                          FLOAT_T * pDst);

// Vector-matrix multiplication
void ext_vec_mat_mult(
                FLOAT_T * pSrcA,  /* m sized input vector */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix */
                FLOAT_T * pDstC); /* m sized output vector */

/* Vector-matrix multiplication and then subtraction from destination vector
 * (r = r - vA) */

// Vector-matrix multiplication
void ext_vec_mat_mult_sub(
                FLOAT_T * pSrcA,  /* m sized input vector */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix */
                FLOAT_T * pDstC); /* m sized output vector */

// Matrix + outer product addition
void ext_add_outer_prod(
    ARM_MATRIX_INSTANCE * pSrcA,  // matrix A
                FLOAT_T * pSrcv,  // vector v
                FLOAT_T * pSrca,  // vector a
    ARM_MATRIX_INSTANCE * pDstB); // result = A + outer(v, a)

// Matrix - outer product subtraction
void ext_sub_outer_prod(
    ARM_MATRIX_INSTANCE * pSrcA,  // matrix A
                FLOAT_T * pSrcv,  // vector v
                FLOAT_T * pSrca,  // vector a
    ARM_MATRIX_INSTANCE * pDstB); // result = A + outer(v, a)


/* Matrix inversion ***************************************/

// Evaluate the inverse of the upper triangular
// (clean matrix -> fill with zeros pDst before storing the inverse)
void ext_mat_inverse_triu(
    const ARM_MATRIX_INSTANCE * pSrcA,
          ARM_MATRIX_INSTANCE * pDstB);

// Evaluate the inverse of the upper triangular and store its transpose
void ext_mat_inverse_triu_t(
    const ARM_MATRIX_INSTANCE * pSrcA,
          ARM_MATRIX_INSTANCE * pDstB);

// Invert the values on the diagonal of the matrix
void ext_mat_inverse_diag(
    const ARM_MATRIX_INSTANCE * pSrcA,
          ARM_MATRIX_INSTANCE * pDstB);

/* Matrix transpositions **********************************/

void ext_mat_square_transpose_overwrite(
    ARM_MATRIX_INSTANCE * pSrcA); /* n x n square matrix, transpose will be
                                   * saved here */



/* Cholesky decomposition *********************************/

// Evaluate cholesky decomposition and output L.T
arm_status ext_mat_cholesky_t(
    const ARM_MATRIX_INSTANCE * pSrc,
          ARM_MATRIX_INSTANCE * pDst);



/* Choleksy QR decomposition ******************************/

// Evaluate QR decomposition using Cholesky
// (Fast, inaccurate for ill-conditioned input matrices)
void ext_mat_qr(
    const ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix A */
          ARM_MATRIX_INSTANCE * pDstQ,  /* n x m output Q matrix */
          ARM_MATRIX_INSTANCE * pDstR); /* m x m output R matrix */

// (Slightly slower, lower memory requirements)
void ext_mat_qr_overwrite(
    ARM_MATRIX_INSTANCE * pSrcQ, /* n x m input matrix A (will be overwritten by
                                  * Q) */
    ARM_MATRIX_INSTANCE * pDstR, /* m x m output matrix R */
                FLOAT_T * buff); /* m sized buffer */



/* SVD decomposition **************************************/

// Evaluate SVD decomposition using Golub-Reinsch method

// Keep original A, evaluate U, s and V
void ext_mat_svd(
    ARM_MATRIX_INSTANCE * pSrcA,     /* m x n input matrix */
    ARM_MATRIX_INSTANCE * pDstU,     /* m x n output matrix */
                FLOAT_T * pDsts,     /* m sized singular values vector */
    ARM_MATRIX_INSTANCE * pDstV,     /* m x m output matrix */
                FLOAT_T * p_m_buff); /* m sized vector buffer */

// Replace original matrix with U, evaluate U and s only
void ext_mat_svd_compact(
    ARM_MATRIX_INSTANCE * pSrcU,     /* m x n output matrix */
                FLOAT_T * pDsts,     /* m sized singular values vector */
                FLOAT_T * p_m_buff); /* m sized vector buffer */



#endif
