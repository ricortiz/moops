
include(${PROJECT_SOURCE_DIR}/cmake/SubmitJob.cmake)
macro(SUBMIT_TO_KILLDEVIL PROGRAM_NAME)
    set(JOB_NAME ${PROGRAM_NAME}_job)
    set(SCRATCH_REMOTE_DIRECTORY /lustre/scr/r/i/ricortiz)
    set(REMOTE_SOURCE_DIR /lustre/scr/r/i/ricortiz/moops)
    set(REMOTE_BINARY_DIR ${REMOTE_SOURCE_DIR}/bin)
    set(REMOTE_SERVER ricortiz@killdevil.unc.edu)
    set(KILLDEVIL_ENV module initclear && module initadd mvapich2_gcc/4.4.6 cuda cmake git)
    set(CMAKE_CACHE "-DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DUSE_CUDA=ON -DUSE_QT_GUI=OFF -DVTK_DIR=~/opt/lib/vtk-5.8 -DENABLE_EXAMPLES=ON")

    add_custom_target(${PROGRAM_NAME}_killdevil_env COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} ${KILLDEVIL_ENV})

    copy_repository(${PROGRAM_NAME} ${CMAKE_SOURCE_DIR} ${REMOTE_SERVER}:${REMOTE_SOURCE_DIR})
    build_app(${PROGRAM_NAME} ${REMOTE_BINARY_DIR} ${CMAKE_CACHE} ${REMOTE_SOURCE_DIR})
    add_data_path(${PROGRAM_NAME} ${SCRATCH_REMOTE_DIRECTORY}/${PROGRAM_NAME}_data)
    set_bsub_commands("${PROGRAM_PATH}/${PROGRAM_NAME} ${SCRATCH_REMOTE_DIRECTORY}/${PROGRAM_NAME}_data" ${JOB_NAME} ${QUEUE} ${NUM_PROCESSORS} ${NUM_CORES_PER_PROCESSORS})
    submit_job(${PROGRAM_NAME}_to_killdevil killdevil_env)
    terminate_job(${PROGRAM_NAME}_at_killdevil ${JOB_NAME})
endmacro(SUBMIT_TO_KILLDEVIL)
#
