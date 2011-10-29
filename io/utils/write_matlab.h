#ifndef WRITE_MATLAB
#define WRITE_MATLAB

#include <iostream>
#include <sstream>
#include <fstream>

template<typename system_type>
bool write_matlab_velocities(system_type &system, const std::string &filename, const std::string array_name = "a")
{
  typedef typename system_type::particle_iterator particle_iterator;
  
  std::stringstream file_stream;
  file_stream.precision(16);
  file_stream << array_name << " = [";
  particle_iterator p = system.particle_begin(), end = system.particle_end();
  for(;p != end; ++p)
  {
    file_stream << p->velocity(); 
    (p != end-1) ? file_stream << ";\n" : file_stream << "];\n";
  }
  std::ofstream file(filename.c_str());

  file << file_stream.str();
  file.flush();
  file.close();
  return 1;
}


#endif
