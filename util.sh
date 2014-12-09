#!/bin/bash

# Define a function for printing error messages
# From: http://stackoverflow.com/a/2990533/1647819
to_stderr(){ cat <<< "$@" 1>&2; }
