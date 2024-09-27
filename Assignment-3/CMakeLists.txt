find_package(Eigen3 CONFIG REQUIRED)

add_executable(main3 main.cpp)

target_link_libraries(main3 PUBLIC math PRIVATE Eigen3::Eigen)
