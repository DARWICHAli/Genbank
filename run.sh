#!/bin/bash

if [ $# = 0 ]; then
    cd project/src
    echo "compile"
    echo 'mv bin to bin'
    #mv X ../Bin
    exit 0

elif [[ $1 == "clean" ]]; then
    #rm -r bin
    exit 0
elif [[ $1 == "test" ]]; then
    echo "begin Testing"
    # lance sh qui lance des test unit && intégration
    sh test.sh
    echo "Done testing"
fi
exit 0
