#include "ext_math_c.h"

void ext_scan(FLOAT_T * pDstv, uint32_t length)
{
    FLOAT_T *pv = pDstv;

    // fill matrix content with input
    uint32_t total_scan_size = sizeof(FLOAT_T)*length;
    uint32_t scanned_size = 0U;
    while(scanned_size < total_scan_size)
    {
        uint32_t scan_size = total_scan_size - scanned_size;
        if(scan_size > 0xFFFF)
        {
            scan_size = 0xFFFF; //max size for uint16_t
        }
        _scan_stream(((char*)pv)+scanned_size, (uint16_t)scan_size);
        scanned_size += 0xFFFF;
    }
}
