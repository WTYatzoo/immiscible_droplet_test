#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} SWE_solver EXE=${exe} OUT_DIR=${out_dir} DT=${dt} DMETRIC=${dmetric} N_X=${n_x} N_Y=${n_y} FRAME=${frame} HEIGHT_MAX_ORI=${height_max_ori} HEIGHT_MIN_ORI=${height_min_ori} GRAVITY=${gravity} RIGID_BODY_DENSITY=${rigid_body_density} FLUID_DENSITY=${fluid_density}
