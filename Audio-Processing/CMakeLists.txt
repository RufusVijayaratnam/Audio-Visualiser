add_library(
    Audio STATIC
    src/Processing/SoundProcessing.hpp
    src/Processing/SoundProcessing.cpp
    src/Parse/AudioFile.h
)

target_include_directories(Audio PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")  