#!/bin/bash

awk '{print}' list.txt

awk '/henda.sh/ {print}' list.txt

awk '/1/ {print $2}' awkTest.txt

awk '$3 == "1" { print $0 }' points.txt

awk -F',' '$3 == "1" { print $0 }' points.txt


awk '$3==6 {print $0}' awkTest.txt

#original
awk '{$'"${1}"'=""}'"${1}"'' "$2".txt | awk '{$'"${1}"'=$'"${1}"'}'"${1}"'' > "$3".txt

#new
awk '{$'"${1}"'=""}'"${1}"'' points.txt | awk '{$'"${1}"'=$'"${1}"'}'"${1}"'' > ptKolbt.txt

awk 'FNR == 2 {print $1 $3 $2}' awkTest.txt

