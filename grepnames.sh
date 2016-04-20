#!/bin/bash
echo $1
echo $2
echo $3
path=$1
while read line
do
    var=$(LC_ALL=C fgrep -A3 "$line" "$2")
    echo "$var" >> $3
done < "$path"
