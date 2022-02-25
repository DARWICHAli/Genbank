#!/bin/bash



if [ $# = 0 ]; then
    rm -r project/bin
    echo "compile"
    for file in project/src/*.java
    do
        javac $file
    done
    cd project
    mkdir bin
    pwd
    mv src/*.class bin/
    echo 'mv bin to bin'
    #mv X ../Bin
    exit 0

elif [[ $1 == "clean" ]]; then
    rm -r project/Bin
    exit 0
elif [[ $1 == "test" ]]; then
    echo "begin Testing"
    # lance sh qui lance des test unit && intégration
    sh test.sh
    echo "Done testing"
fi
exit 0
