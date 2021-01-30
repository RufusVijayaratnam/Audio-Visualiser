add_executable(Visualise main.cpp)

find_package(GSL REQUIRED)
target_link_libraries(Visualise PRIVATE GSL::gsl GSL::gslcblas)
target_link_libraries(Visualise PRIVATE Audio)
