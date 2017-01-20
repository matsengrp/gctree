#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SCons.Action import ActionFactory
from SCons.Util import AddMethod
from SCons.Script import Environment

import SCons.Util
import os
import time
import subprocess


# Running commands on the cluster sometimes has the unfortunate side-effect of
# letting distributed filesystems get out of sync.  A file that is written on
# the cluster may not be visible on local machines for several seconds.  This
# doesn't happen all the time, and it can be difficult to demonstrate, but it
# often occurs when local and cluster commands manipulate the same file in
# quick succession.  The Wait() action is meant to wait for a file to become
# locally visible after running a command on the cluster (via srun or salloc).
#
# Wait() usage is typically,
# 	target='output.file'
# 	env.Command(target, 'source.file',
#   	        [ "srun some-command <${TARGET}",
#    			   Wait(target)
#				])
#
# This will cause the execution to pause after running 'some-command' until the target shows up on the local machine.
# The target will be polled on a 2-second interval, and the command will fail if the target does not show up within about 10 seconds.

def get_paths_str(dest):
    # If dest is a list, we need to manually call str() on each element
    if SCons.Util.is_List(dest):
        elem_strs = []
        for element in dest:
            elem_strs.append('"' + str(element) + '"')
        return '[' + ', '.join(elem_strs) + ']'
    else:
        return '"' + str(dest) + '"'

# https://github.com/azatoth/scons/blob/73f996e59902d03ec432cc662252aac5fb72f1f8/src/engine/SCons/Defaults.py 
def wait_func(dest):
    SCons.Node.FS.invalidate_node_memos(dest)
    if not SCons.Util.is_List(dest):
        dest = [dest]
    for entry in dest:
        count = 0
        limit = 30
        while not os.path.isfile(entry) or os.stat(entry).st_size == 0:
            print("waiting for {}...".format(entry))
            time.sleep(2)
            count = count + 1
            if count >limit:
                return 1
    return 0

Wait = ActionFactory(wait_func, lambda dir: 'Wait(%s)' % get_paths_str(dir))

            
# Define one of two versions of the SRun method,
# depending on whether the `srun` command is available or not.
exit_code = subprocess.call("type srun", shell=True, 
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
srun_exists = (exit_code == 0)
def SRun(env, target, source, action, **kwargs):
    if not hasattr(target, '__iter__'):
        target = [target]
    waitfor = target
    if 'chdir' in kwargs:
        waitfor = [os.path.basename(w) for w in waitfor]
    if srun_exists:
        action = ["srun sh -c ' "+action + " '", Wait(waitfor)]
    result = env.Command(target=target, source=source, action=action, **kwargs)
    return result

# Use the global AddMethod function to add a method to the global Environment class,
# so it will be used in all subsequently created environments.
AddMethod(Environment, SRun)



            
