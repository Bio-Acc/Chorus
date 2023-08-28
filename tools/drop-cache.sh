#!/bin/bash
param=${1:-"4"}
sync; echo ${param} > /proc/sys/vm/drop_caches 
sync; echo 4 > /proc/sys/vm/drop_caches
