#!/bin/sh -ex
echo "Wzy@12345" | sudo -S /usr/local/sbin/drop-cache.sh 3
echo "done"