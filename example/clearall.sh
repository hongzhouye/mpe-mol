#!/bin/bash

find . -name "*.txt" | grep -v "xc_func" | xargs rm -f
find . -name "*pathtable*" | xargs rm -f
