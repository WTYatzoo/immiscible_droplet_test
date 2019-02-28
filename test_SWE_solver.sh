#! /bin/bash

frame=200
dt=0.0001
dmetric=0.6
n_x=151
n_y=151
height_max_ori=15
height_min_ori=5
gravity=9.8
fluid_density=1000
rigid_body_density=300

project_dir=${HOME}/project/mytest/immiscible_droplet_test
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/immiscible_droplet
base_output_dir=${project_dir}/data
out_dir=${base_output_dir}/SWE_solver_result_13
mkdir -p ${out_dir}
source ${script_dir}/SWE_solver.sh
