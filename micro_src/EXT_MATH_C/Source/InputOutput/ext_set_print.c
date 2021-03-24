#include "ext_math_c.h"

void ext_set_print(void (*print_stream)(char*, unsigned int))
{
    _print_stream = print_stream;
}
