#!/bin/bash

#  chmod +x henda.sh 
#  ./henda.sh 

# nice commands
# df - h (shows you the space on your computer)
# watch (runs the comand every 2 seconds)

echo "Hello World"

HELLO="Hello World"

echo $HELLO

######################
#If Statement
MYNUM=100

if [ $MYNUM -eq 200 ]
	then
		echo "Testing"
	else
		echo "100"
fi

####################
if [ -d ~/Klosta\ Group/Data\ Processing/Bash/MyFolder ]
	then
		echo "The Folder Exists."
	else
		echo "Nope"
fi



if [ -d ~/Klotsa\ Group/Data\ Processing/Bash/MyFolder ]
	then
		echo "The Folder Exists."
	else
		echo "Nope"
fi

#############

echo "Name"
read myName
echo "You entered: $myName"

######################

#piping
ls > list.txt
# or to append
ls >> list.txt

#####################
myVar=1

while [ $myVar -le 10 ]
do
	echo $myVar
	myVar=$(( $myVar + 1 ))
	sleep 0.5
done

#####################

echo "1,2,3"

read NUM;

case $NUM in
	1) echo "1 is okay";;
	2) echo "2 is good";;
    3) echo "3 is better";;
	*) echo "Wrong"
esac







