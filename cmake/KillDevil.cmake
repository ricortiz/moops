
include(cmake/SubmitJob.cmake)

message(STATUS "Adding killdevil job submision target")
option(CLONE_PROJECT_KILLDEVIL "Clone entire project in killdevil" ON)
option(KILL_JOB_FIRST "Kill job in cluster" OFF)

set(SCRATCH_REMOTE_DIRECTORY /lustre/scr/r/i/ricortiz)
set(MOOPS_SOURCE_DIR /lustre/scr/r/i/ricortiz/moops)
set(MOOPS_BINARY_DIR ${MOOPS_SOURCE_DIR}/bin)
set(REMOTE_SERVER ricortiz@killdevil.unc.edu)
set(CMAKE_CACHE "-DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DUSE_CUDA=ON -DUSE_QT_GUI=OFF -DVTK_DIR=~/opt/lib/vtk-5.8 -DENABLE_EXAMPLES=ON")
set(KILLDEVIL_ENV "bash <<EOF \nmodule initclear")
set(KILLDEVIL_ENV "${KILLDEVIL_ENV} \nmodule initadd mvapich2_gcc/4.4.6 cuda cmake git")
set(QUEUE gpu)
set(NUM_PROCESSORS 1)
set(PROGRAM_PATH ${MOOPS_BINARY_DIR}/examples/valveless_heart)
set(PROGRAM_NAME valveless_heart)
set(JOB_NAME moops_${PROGRAM_NAME})


execute_process(COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} ${KILLDEVIL_ENV} RESULT_VARIABLE out)
#
if(CLONE_PROJECT_KILLDEVIL)
    clone_repository()
else()
    pull_repository()
endif()

set_build_commands(${CMAKE_CACHE})
add_data_path(${SCRATCH_REMOTE_DIRECTORY}/${PROGRAM_NAME}_data)


set_bsub_commands("${PROGRAM_PATH}/${PROGRAM_NAME} ${SCRATCH_REMOTE_DIRECTORY}/${PROGRAM_NAME}_data" ${JOB_NAME} ${QUEUE})
submit_job("SubmitValvelesHeartSolvers")
