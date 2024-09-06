find_package(Eigen3 CONFIG REQUIRED)

add_executable(main1 main.cpp)

target_link_libraries(main1 PUBLIC math PRIVATE Eigen3::Eigen)
