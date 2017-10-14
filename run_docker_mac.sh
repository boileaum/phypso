#!/bin/bash

socat TCP-LISTEN:6001,reuseaddr,fork UNIX-CLIENT:\"$DISPLAY\" &
docker run -ti --rm  -e XAUTHORITY=/tmp/xauth \
       -v ~/.Xauthority:/tmp/xauth \
       -e DISPLAY=$(ifconfig en0 | grep 'inet '|awk '{print $2}'):1 \
       --net host \
       -v $(pwd):/home/euler/phypso boileaum/phypso-env

