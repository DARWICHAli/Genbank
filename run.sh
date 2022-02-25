#!/bin/bash



if [ $# = 0 ]; then
    rm -r project/bin
    echo "compile"
    for file in project/src/*.java
    do
        echo "compiling" $file
        javac $file
    done
    cd project
    mkdir bin
    mv src/*.class bin/
    exit 0

elif [[ $1 == "clean" ]]; then
    rm -r project/bin
    exit 0
elif [[ $1 == "test" ]]; then
    echo "begin Testing"
    # lance sh qui lance des test unit && intégration
    sh test.sh
    echo "Done testing"
fi
exit 0
