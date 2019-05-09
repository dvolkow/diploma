#! /bin/bash

export CC=/usr/bin/clang-6.0
export ASAN_OPTIONS=fast_unwind_on_malloc=0
cmake .. -DCMAKE_C_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -g" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined" -DCMAKE_MODULE_LINKER_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined" && make
