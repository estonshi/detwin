# !/bin/bash

set -e
rootpath=`pwd`

# crystfel version
flag=0
while [ $flag -eq 0 ]
do
	echo "CrystFEL package path (e.g '/home/myuser/Documents/crystfel-0.8.0'): "
	read crpath
	if [ -d "$crpath/libcrystfel" ]
	then
		flag=1
	else
		echo "Invalid path !"
	fi
done

tmp=`find $crpath/relnotes-*`
tmp=(${tmp//relnotes-/ })
version=${tmp[1]}
if [ ! -d "${rootpath}/${version}" ]
then
	echo "CrystFEL version $version is not supported. Use (0.6.30 or (0.7.0) or (0.8.0) ."
	exit 1
fi

# get new file
detwinc=${rootpath}/${version}/detwin.c
detwinh=${rootpath}/${version}/detwin.h
cmakelists=${rootpath}/${version}/CMakeLists.txt
manpage=${rootpath}/${version}/detwin.1

# replace
cp -f $detwinc $crpath/src
cp -f $detwinh $crpath/src
cp -f $cmakelists $crpath
cp -f $manpage $crpath/doc/man

echo "Done !"