#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        BOX_SPATIAL_DIMENSION=2\n"
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        GRAVITY_ANALYTIC\n"
"        ENERGY_ENTROPY_SWITCH_IS_ACTIVE\n"
"\n");
}
