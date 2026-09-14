#!/bin/bash -eu
# Linux only, does not work with BSD sed
fgrep -Rwl "$2" "$1" | xargs sed -i "s/\b$2\b/$3/g"
