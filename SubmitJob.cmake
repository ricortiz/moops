

find_program(GIT_COMMAND NAMES git)
find_program(SSH_COMMAND NAMES ssh)
find_program(SCP_COMMAND NAMES scp)
set(MOOPS_CHECKOUT_COMMAND ${GIT_COMMAND} clone git@bitbucket.org:ricortiz/moops.git )
