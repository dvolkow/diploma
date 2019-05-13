#! /bin/bash

export CC=/usr/bin/clang-7
export ASAN_OPTIONS=fast_unwind_on_malloc=0
cmake .. -DCMAKE_C_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -g" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined" -DCMAKE_MODULE_LINKER_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined" && make

cp ../data/apogee-rc-DR14_hsoy_with_err.txt ./apogee-hsoy.txt
cp ../data/apogee-rc-DR14_ucac_with_err.txt ./apogee-ucac.txt
