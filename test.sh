#!/usr/bin/env bash


echo "testing"

for file in project/*
do
    for res in project/test/test_u/*
    do
        echo $res $file

    done
done
