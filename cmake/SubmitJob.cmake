# find local commands
find_program(SSH_COMMAND NAMES ssh)
find_program(SSH_COMMAND NAMES scp)
find_program(RSYNC_COMMAND NAMES rsync)

# macro to clone repo into remote directory
macro(CLONE_REPOSITORY APP_NAME SOURCE_DIR)
    set(CMD "if [ -d ${SOURCE_DIR} ]\; then rm -Rf ${SOURCE_DIR} \; fi")
    set(CMD "${CMD} && git clone git@bitbucket.org:ricortiz/moops.git ${SOURCE_DIR} && mkdir ${SOURCE_DIR}/bin")
    add_custom_target(${APP_NAME}_clone_repository COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} ${CMD} VERBATIM)
endmacro(CLONE_REPOSITORY)

# macro to pull repo into remote directory
macro(PULL_REPOSITORY APP_NAME SOURCE_DIR)
    set(CMD "cd ${SOURCE_DIR} && git pull")
    add_custom_target(${APP_NAME}_pull_repository COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} ${CMD} VERBATIM)
endmacro(PULL_REPOSITORY)

# macro to rsync repo in to remote directory
macro(COPY_REPOSITORY APP_NAME FROM TO)
    set(CMD "${RSYNC_COMMAND} -av --rsh=ssh --exclude build* --exclude *.git* --exclude *docs* --exclude *~ ${FROM} ${TO}")
    add_custom_target(${APP_NAME}_copy_repository COMMAND ${CMD} VERBATIM)
endmacro(COPY_REPOSITORY)

# macro to build app in remote binary directory
macro(BUILD_APP APP_NAME BINARY_DIR CMAKE_CACHE SOURCE_DIR)
    set(CMD "if [ -f ${BINARY_DIR}/CMakeCache.txt ]\; then rm ${BINARY_DIR}/CMakeCache.txt \; fi")
    set(CMD "${CMD} && cd ${BINARY_DIR} && cmake ${CMAKE_CACHE} ${SOURCE_DIR} && make ${APP_NAME}")
    add_custom_target(${APP_NAME}_remote_build COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} ${CMD} VERBATIM)
endmacro(BUILD_APP)

# macro to set up BSUB_COMMANDS (bsub commands to be sent to the cluster)
macro(SET_BSUB_COMMANDS EXE JOB_NAME QUEUE NUM_MPI_PROCESSORS NUM_CORES_PROCESSORS)
    set(BSUB_COMMANDS "\nbsub <<EOF")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n####   Run program    ####")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -q ${QUEUE}")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -J ${JOB_NAME}")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -n ${NUM_MPI_PROCESSORS}")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -T ${NUM_CORES_PROCESSORS}")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -x")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -o out.%J")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n#BSUB -e err.%J")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \nmpirun ${EXE}")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \n####   End run program    ####")
    set(BSUB_COMMANDS "${BSUB_COMMANDS} \nEOF")
endmacro(SET_BSUB_COMMANDS)

# macro to create a remote data directory 
macro(ADD_DATA_PATH APP_NAME DATA_PATH)
    set(CMD "if [ ! -d ${DATA_PATH} ]\; then mkdir -p ${DATA_PATH} \; fi")
    add_custom_target(${APP_NAME}_add_data_path COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} ${CMD} VERBATIM)
endmacro()

# macro to create a target that submit the job to the cluster
macro(SUBMIT_JOB APP_NAME)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Submit${PROGRAM_NAME}Job.sh "#!/bin/bash \n\n${BSUB_COMMANDS}")
    execute_process(COMMAND ${SCP_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/Submit${PROGRAM_NAME}Job.sh ${REMOTE_SERVER}:${SCRATCH_REMOTE_DIRECTORY}/ RESULT_VARIABLE out OUTPUT_QUIET)
    add_custom_target(${APP_NAME} COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} sh ${SCRATCH_REMOTE_DIRECTORY}/Submit${PROGRAM_NAME}Job.sh VERBATIM)
    add_dependencies(${APP_NAME} ${APP_NAME}_clone_repository ${APP_NAME}_add_data_path ${APP_NAME}_remote_build)
endmacro()

# macro to create a target that terminates the job with specified job name
macro(TERMINATE_JOB APP_NAME)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Terminate${PROGRAM_NAME}Job.sh "#!/bin/bash \n\nbkill `bjobs -u ricortiz -J ${JOB_NAME} -l | sed -n '/Job </p' | sed -e 's/.*Job <//;s/>.*//'`")
    execute_process(COMMAND ${SCP_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/Terminate${PROGRAM_NAME}Job.sh ${REMOTE_SERVER}:${SCRATCH_REMOTE_DIRECTORY}/ RESULT_VARIABLE out OUTPUT_QUIET)
    add_custom_target(${APP_NAME} COMMAND ${SSH_COMMAND} ${REMOTE_SERVER} sh ${SCRATCH_REMOTE_DIRECTORY}/Terminate${PROGRAM_NAME}Job.sh VERBATIM)
endmacro()

