#!/bin/sh
#
# This script can be used on Linux servers to detect httpd processes with tcp connections
# in CLOSE_WAIT state.   This can happen when a user closes their browser/tab or hits this stop
# button on the sample web interface.  This will not work with other non-Linux versions of netstat.
#

while true
do
  pid=`netstat -anp | grep CLOSE_WAIT | head -1 | awk '{print $7}' | cut -f 1 -d '/'`
  if [ -n "${pid}" ]
  then
    kill -9 ${pid}
  fi
  sleep 10
done
