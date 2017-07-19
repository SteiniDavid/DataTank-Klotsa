#!/bin/bash

#set the particle number
particleNum=31

mkdir particleFiles

for ((i=1; i<=$particleNum; i++)); do
	awk '$3=='"${i}"' {print $0}' points.txt > particleFiles/p$i.txt
	echo "Working"
done

touch particleFiles/positions.txt
echo "$particleNum" > particleFiles/positions.txt
echo "$particleNum"
echo " " >> particleFiles/positions.txt

numLines=$(wc -l < particleFiles/p1.txt)
echo "$numLines is number of lines"

for ((line=1; line<=$numLines; line++)); do
	for ((part=1; part<=$particleNum; part++)); do
		awk 'FNR == '"$line"' {print $3 " " $1 " " $2}' particleFiles/p$part.txt >> particleFiles/positions.txt
	done
		echo "$particleNum" >> particleFiles/positions.txt
		echo "  " >> particleFiles/positions.txt
done

return 0