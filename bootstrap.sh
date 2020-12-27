#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mkdir -p $SCRIPTDIR/include
mkdir -p $SCRIPTDIR/lib

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # ...
	cpu_count=$(cat /proc/cpuinfo |grep processor| wc -l)
elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
	cpu_count=$(sysctl -n hw.ncpu)
fi

# Build & install sfml-widgets
cd $SCRIPTDIR/build
mkdir -p sfml-widgets/release && cd sfml-widgets/release
make -j$cpu_count
