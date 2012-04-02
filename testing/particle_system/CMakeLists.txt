
set(particle_system elastic_system.cpp spring_system.cpp spring.cpp)

make_tests("${particle_system}" "TestElasticSystem")
