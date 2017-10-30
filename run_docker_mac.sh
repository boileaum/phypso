#!/bin/bash

trap 'kill $BGPID; exit' SIGINT
socat TCP-LISTEN:6001,reuseaddr,fork UNIX-CLIENT:\"$DISPLAY\" &
BGPID=$!
active_interface_name=$(ifconfig | pcregrep -M -o '^[^\t:]+:([^\n]|\n\t)*status: active' | egrep -o -m 1 '^[^\t:]+')
docker run -ti --rm  -e XAUTHORITY=/tmp/xauth \
       -v ~/.Xauthority:/tmp/xauth \
       -e DISPLAY=$(ifconfig $active_interface_name | grep 'inet '|awk '{print $2}'):1 \
       --net host \
       -v $(pwd):/home/euler/phypso boileaum/phypso-env
kill $BGPID; exit
