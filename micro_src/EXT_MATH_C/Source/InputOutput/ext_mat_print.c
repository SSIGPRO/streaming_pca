#include "ext_math_c.h"

void ext_mat_print(const ARM_MATRIX_INSTANCE * pSrcA)
{
    uint16_t rows = pSrcA->numRows, cols = pSrcA->numCols;
    uint32_t length = rows*cols;
    FLOAT_T *pA = pSrcA->pData;

    // send matrix content
    uint32_t total_print_size = sizeof(FLOAT_T)*length;
    uint32_t printed_size = 0;
    while(printed_size < total_print_size)
    {
        uint32_t print_size = total_print_size - printed_size;
        if(print_size > 0xFFFF)
        {
            print_size = 0xFFFF; //max size for uint16_t
        }
        _print_stream(((char*)pA)+printed_size, (uint16_t)print_size);
        printed_size += 0xFFFF;
    }
}
