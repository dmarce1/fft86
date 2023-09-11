cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)
project(fft86 CXX C ASM)
set(CMAKE_CXX_FLAGS "-march=native -ffast-math")

set(fft86_source_files
    ./src/apply_twiddles.S
    ./src/fft_1d_batch_dif.S
    ./src/fft_1d_batch_dit.S
 	./src/fft86.cpp
	./src/fftw.cpp
	./src/fourstep.cpp
	./src/scramble.cpp
	./src/scramble.S
	./src/mixed.cpp
	./src/real.cpp
	./src/transpose.cpp
	./src/transpose.S
	./src/transpose_prime.cpp
	./src/twiddles.cpp
	./src/pfa.cpp
)      

set(fft86_include_files
	./fft86/defs.hpp
	./fft86/fft.hpp
	./fft86/fftw.hpp
	./fft86/scramble.hpp
	./fft86/timer.hpp
	./fft86/transpose.hpp
	./fft86/twiddles.hpp
	./fft86/util.hpp
	./fft86/vec.hpp
)

add_executable(fft86 ${fft86_source_files} ${fft86_include_files})

target_link_directories(fft86 PUBLIC ${CMAKE_SOURCE_DIR}/lib)
target_include_directories(fft86 PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(fft86 PUBLIC fftw3 sfft)

