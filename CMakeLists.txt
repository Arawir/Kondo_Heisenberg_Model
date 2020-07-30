project(Kondo_Heisenberg)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 17)

set(SOURCES
    main.cpp)
set(HEADERS
    kondo_heisenberg.h
    tdvp.h
    operators.h
    experiment1.h
)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC include)

target_link_libraries(${PROJECT_NAME} ${GTEST_LIBRARIES} -m64 -std=c++17 -fconcepts -fPIC -I'/home/maciej/itensor' -O2 -DNDEBUG -Wall -Wno-unknown-pragmas -Wno-unused-variable -L'/home/maciej/itensor/lib' -litensor -larmadillo -lpthread -L/usr/lib -lblas -llapack)
