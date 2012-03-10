
include(cmake/SubmitJob.cmake)

message(STATUS "Adding killdevil job submision target")
option(CLONE_PROJECT_KILLDEVIL "Clone entire project in killdevil" ON)
option(PULL_PROJECT_KILLDEVIL "Pull master branch changes in killdevil" OFF)
option(KILL_JOB_FIRST "Kill job in cluster" OFF)

set(SCRATCH_REMOTE_DIRECTORY /lustre/scr/r/i/ricortiz)
set(MOOPS_REMOTE_DIRECTORY /lustre/scr/r/i/ricortiz/moops)
set(SERVER ricortiz@killdevil.unc.edu)
set(JOB_NAME moops)
set(PROGRAM_NAME TestFluidSolvers)

set(SET_ENV "module initclear")
set(SET_ENV "${SET_ENV} \nmodule initadd mvapich2_gcc/4.4.6 cuda cmake git")

set(GIT_COMMANDS "")
set_git_cmd(${SCRATCH_REMOTE_DIRECTORY})

set(BUILD_CMD "${BUILD_CMD} \ncd ${MOOPS_REMOTE_DIRECTORY}/bin")
set(BUILD_CMD "${BUILD_CMD} \ncmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DUSE_CUDA=ON -DUSE_QT_GUI=OFF -DVTK_DIR=~/opt/lib/vtk-5.8 ..")
set(BUILD_CMD "${BUILD_CMD} \nmake ${PROGRAM_NAME}\n")

set(BSUB_COMMANDS "")
set_bsub_cmd(testing/math/${PROGRAM_NAME} ${JOB_NAME} gpu)

file(WRITE ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh "#!/bin/bash\n\n")
file(APPEND ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh "${SET_ENV}\n")
file(APPEND ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh "${GIT_COMMANDS}\n")
file(APPEND ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh "${BUILD_CMD}\n")
file(APPEND ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh "${BSUB_COMMANDS}\n")

# CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/Submit.sh.in"
#                "${PROJECT_BINARY_DIR}/KillDevilSubmit.sh" @ONLY IMMEDIATE)

# if(KILL_JOB_FIRST)
#     execute_process(COMMAND ${SSH_COMMAND} ${SERVER} bkill `bjobs -u ricortiz | grep ${JOB_NAME} | cut -f1 -d\" \"`)
# endif()
add_custom_target(KillDevilSubmit ALL
        COMMAND ${SCP_COMMAND} ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh ${SERVER}:${SCRATCH_REMOTE_DIRECTORY}
        COMMAND ${SSH_COMMAND} ${SERVER} sh ${SCRATCH_REMOTE_DIRECTORY}/KillDevilSubmit.sh
        COMMAND ${SSH_COMMAND} ${SERVER} bjobs -J ${JOB_NAME}
        )

