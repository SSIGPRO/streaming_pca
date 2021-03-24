#include "ext_math_c.h"

void ext_set_scan(void (*scan_stream)(char*, unsigned int))
{
    _scan_stream = scan_stream;
}
