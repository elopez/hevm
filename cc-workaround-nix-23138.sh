#!/bin/sh
if [ "$1" = "--print-file-name" ] && [ "$2" = "c++" ]; then
    echo c++
else
    exec cc "$@"
fi
