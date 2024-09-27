build:
    mkdir -p build
    cd build && cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    cmake --build build
    build/ex_bhc3d

clean:
    rm -rf build
