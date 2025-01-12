#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        STARS\n"
"        WENDLAND_C2_KERNEL\n"
"        PASSIVE_SCALARS=1\n"
"        REGULARIZE_MESH_CM_DRIFT\n"
"        REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED\n"
"        REGULARIZE_MESH_FACE_ANGLE\n"
"        REFINEMENT_SPLIT_CELLS\n"
"        REFINEMENT_MERGE_CELLS\n"
"        NODEREFINE_BACKGROUND_GRID\n"
"        COOLING\n"
"        LOW_TEMP_COOLING\n"
"        USE_SFR\n"
"        SELFGRAVITY\n"
"        GRAVITY_NOT_PERIODIC\n"
"        MULTIPLE_NODE_SOFTENING\n"
"        ADAPTIVE_HYDRO_SOFTENING\n"
"        DOUBLEPRECISION=1\n"
"        OUTPUT_IN_DOUBLEPRECISION\n"
"        INPUT_IN_DOUBLEPRECISION\n"
"        HAVE_HDF5\n"
"        DEBUG\n"
"\n");
}
