#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

git submodule update --init --recursive

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
cd $SCRIPTDIR/third-parties/sfml-widgets/
make -j$cpu_count
cd src/Gui
rsync -a --include '*/' --include '*.hpp' --exclude '*' . $SCRIPTDIR/include
rsync -a --include '*/' --include '*.inl' --exclude '*' . $SCRIPTDIR/include
cd ../../lib
cp -a libsfml-widgets.a $SCRIPTDIR/lib/
echo 
echo
echo "Successfully built and installed sfml-widgets dependency!"

# Build conjugate-gradient
cd $SCRIPTDIR
mkdir -p build && cd build/
cmake ..
make -j$cpu_count
echo
echo
echo "Successfully built conjugate-gradient! Now run ./build/gradient to start the program!"
echo "Proudly made by Kiiiiiira. Enjoy!"
