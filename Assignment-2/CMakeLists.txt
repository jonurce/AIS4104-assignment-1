find_package(Eigen3 CONFIG REQUIRED)

add_executable(main2 main.cpp)

target_link_libraries(main2 PUBLIC math PRIVATE Eigen3::Eigen)
