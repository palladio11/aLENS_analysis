#!/usr/bin/bash

# Warn on undefined variable access
set -u

COPYTIME=20 # (minutes) time before job ends to start syncing data
GRACETIME=10 # (minutes) don't start more work if COPYTIME+GRACETIME < JOBTIMELEFT
SYNCTIME=60 # (minutes) how often to update trajectories from temp to global
SIMSOURCE="$(realpath $1)"
SCRIPTPATH="$(realpath $2)"

if [[ $# -ne 2 || ! -d "$1" || ! -f "$2" ]]; then
    echo "Error! Must provide argument to an existent path and a script to run"
    exit 1
fi

# Gets the time left on the job from the current time, converts it to minutes, subtracts
# COPYTIME, and returns if number is less than GRACETIME. Result is stored in DURATION
get_time_to_exit() {
    if [[ -z $SLURM_JOBID ]]; then return 1; fi
       
    TIME=$(squeue -j $SLURM_JOBID -h --format %L | tr '-' ':' | tr ':' '\n' | tac | tr '\n' ' ')    
    read -a TIMEARR <<< "$TIME"
    TIMETOEXIT="$((${TIMEARR[1]:-0} + ${TIMEARR[2]:-0} * 60 + ${TIMEARR[3]:-0} * 24 * 60 - $COPYTIME))"
    DURATION="${TIMETOEXIT}m"

    return $((${TIMETOEXIT:-0} < $GRACETIME))
}

# Copy SIMSOURCE to a temporary directory
copy_to_tmp() {
    tmpdir=$(mktemp -d -p /dev/shm)
    export tmpdir
    rm -f "$SIMSOURCE/copy.log"
    echo "${SIMSOURCE}: Syncing to $tmpdir" > $tmpdir/copy.log
    rsync -av  --exclude={'*.zip','*.log','*.out'} "$SIMSOURCE/." $tmpdir >> $tmpdir/copy.log
    echo "${SIMSOURCE}: Done syncing" >> $tmpdir/copy.log
}

# Sync temporary directory back to SIMSOURCE
update_zip() {
    echo "${SIMSOURCE}: Syncing back from $tmpdir"
    cd result
    zip -ur ${SIMSOURCE}/result.zip *
    cd $tmpdir
    # now that we've zipped things, we could clean up data from the tmpdir that's already in
    # the zip file to alleviate memory pressure, if we wanted
    return $?
}

# finalize. Copy all logs, etc if you want
finalize() {
    # stop our watch script
    kill -SIGTERM $WATCHPID
    rsync -rv --exclude=result $tmpdir/* $SIMSOURCE/ >> copy.log
    update_zip >> copy.log
    cp copy.log $SIMSOURCE
    cd $SIMSOURCE
    ls
    rm -rf $tmpdir
}

# If script receives a SIGUSR1 signal (kill -SIGUSR1), copy data back and quit
# not using yet because disbatch doesn't support forwarding signals from slurm.. yet
trap "finalize; exit" SIGUSR1

# export anything that update_zip needs
export -f update_zip
export SIMSOURCE

# Load time to exit and abort if not enough time
if ! get_time_to_exit; then
    echo "${SIMSOURCE}: Not enough time to copy sim! Aborting!"
    exit 1
fi

#######################################################################
#                   Actual script stuff starts here                   #
#######################################################################

copy_to_tmp
cd $tmpdir

#If there already exists a result.zip file in here, 
#  extract to tmpdir so that the run can continue.
#  Might want to be more selective about this in the future.
[[ -f ${SIMSOURCE}/result.zip ]] && { echo "Extracting result files";  
    unzip -u "${SIMSOURCE}/result.zip" -d result > /dev/null; 
    echo "Result files extracted."; }


# we just want to run the update_zip command every SYNCTIME minutes...
# the watch command annoyingly uses clearing control sequences, making logging crap complicated
# so we just run a loop in a subshell and kill that when we're ready
(while true; do
    sleep $((60 * $SYNCTIME)) || break
    update_zip >> copy.log || break
 done) &
WATCHPID=$!
disown # disown process so shell doesn't care when we kill this. otherwise ugly terminate text on finalize
#conda activate alens

# SIMULATION MAGIC TIME
timeout $DURATION $SCRIPTPATH

# wrap it up
finalize
