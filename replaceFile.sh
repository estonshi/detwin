# !/bin/bash

set -e
rootpath=`pwd`

# crystfel version
flag=0
while [ $flag -le 0 ]
do
	if [ ! -z $1 ] && [ $flag -eq 0 ]
	then
		crpath=$1
	else
		echo "CrystFEL package path (e.g '/home/myuser/Documents/crystfel-0.8.0'): "
		read crpath
	fi
	if [ -d "$crpath/libcrystfel" ]
	then
		flag=1
	else
		flag=-1
		echo "Invalid path !"
	fi
done

tmp=`find $crpath/relnotes-*`
tmp=(${tmp//relnotes-/ })
version=${tmp[1]}
if [ ! -d "${rootpath}/${version}" ]
then
	echo "CrystFEL-$version is not supported. Use version higher than 0.6.3."
	exit 1
fi

# get new file
detwinc=${rootpath}/${version}/detwin.c
detwinh=${rootpath}/${version}/detwin.h
if [ ${version}x = "0.8.0"x ] || [ ${version}x = "0.9.0"x ]
then
	makefile=${rootpath}/${version}/CMakeLists.txt
else
	makefile="${rootpath}/${version}/Makefile.in ${rootpath}/${version}/Makefile.am"
fi
manpage=${rootpath}/${version}/detwin.1

# replace
cp -f $detwinc $crpath/src
cp -f $detwinh $crpath/src
cp -f $makefile $crpath
cp -f $manpage $crpath/doc/man

echo "Done !"