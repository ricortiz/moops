
include(SubmitJob.cmake)
option(CLONE_PROJECT_KILLDEVIL "Clone entire project in killdevil" OFF)
option(PULL_PROJECT_KILLDEVIL "Pull master branch changes in killdevil" OFF)
set(SCRATCH_REMOTE_DIRECTORY /lustre/scr/r/i/ricortiz)
set(MOOPS_REMOTE_DIRECTORY /lustre/scr/r/i/ricortiz/moops)
set(SERVER ricortiz@killdevil.unc.edu)
message(STATUS "Adding killdevil submision")


set(SET_ENV "module initclear")
set(SET_ENV "${SET_ENV} \nmodule initadd mvapich2_gcc/4.4.6 cuda cmake git")
set(SET_ENV "${SET_ENV} \ncd ${SCRATCH_REMOTE_DIRECTORY}")
if(CLONE_PROJECT_KILLDEVIL)
    set(SET_ENV "${SET_ENV} \nif [ -x moops ]; then")
    set(SET_ENV "${SET_ENV} \nrm -Rf moops")
    set(SET_ENV "${SET_ENV} \nfi")
    set(GIT_CMD "git clone git@bitbucket.org:ricortiz/moops.git ${MOOPS_REMOTE_DIRECTORY}")
    set(EXEC "mkdir ${MOOPS_REMOTE_DIRECTORY}/bin")
elseif(PULL_PROJECT_KILLDEVIL)
    set(GIT_CMD "${GIT_CMD} \ncd ${MOOPS_REMOTE_DIRECTORY}")
    set(GIT_CMD "git pull")
else()
    set(GIT_CMD "")
endif()
set(EXEC "${EXEC} \ncd ${MOOPS_REMOTE_DIRECTORY}/bin")
set(EXEC "${EXEC} \ncmake -DUSE_MPI=ON -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DUSE_CUDA=ON -DUSE_QT_GUI=OFF -DVTK_DIR=~/opt/lib/vtk-5.8 ..")
set(EXEC "${EXEC} \nmake TestFluidSolvers")
set(BSUB "bsub <<EOF")
set(BSUB "${BSUB} \n####   Run program    ####")
set(BSUB "${BSUB} \n#BSUB -q gpu")
set(BSUB "${BSUB} \n#BSUB -J moops")
set(BSUB "${BSUB} \n#BSUB -n 2")
set(BSUB "${BSUB} \n#BSUB -o out.%J")
set(BSUB "${BSUB} \n#BSUB -e err.%J")
set(BSUB "${BSUB} \nmpirun testing/math/TestFluidSolvers stokesParallel")
set(BSUB "${BSUB} \n####   End run program    ####")
set(BSUB "${BSUB} \nEOF")

# 
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/Submit.sh.in"
               "${PROJECT_BINARY_DIR}/KillDevilSubmit.sh" @ONLY IMMEDIATE)

add_custom_target(KillDevilSubmit
                    COMMAND ${SCP_COMMAND} ${PROJECT_BINARY_DIR}/KillDevilSubmit.sh ${SERVER}:${SCRATCH_REMOTE_DIRECTORY}
                    COMMAND ${SSH_COMMAND} ${SERVER} sh ${SCRATCH_REMOTE_DIRECTORY}/KillDevilSubmit.sh
                    )