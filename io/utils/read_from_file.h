#ifndef READ_FROM_FILE_H
#define READ_FROM_FILE_H
//
// Author: Ricardo Ortiz <rortizro@caco>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
// $Id: read_from_file.h 31 2009-10-14 15:56:30Z rortiz $

#include <OpenTissue/core/containers/mesh/mesh.h>
namespace IO {

namespace polymesh = OpenTissue::polymesh;

template<typename surface_type>
bool read_from_file(surface_type &surface, const std::string &file_name)
{
  typedef typename surface_type::math_types           math_types;
  typedef typename math_types::vector3_type          vector3_type;
  typedef typename surface_type::particle_iterator	
particle_iterator;
  typedef typename surface_type::particle_type	    particle_type;
	
  // Read File "IBPoints.vertex"
  std::cout << "Read File " << file_name << std::endl;
  std::ifstream linkfile(file_name.c_str());
  if (!linkfile.is_open()){
    std::cout << "Error: can't access file " <<  file_name << std::endl;
  }

  int num_points;// number of IB points
  linkfile >> num_points;	
  double x, y, z;
  for (int i = 0; i < num_points; ++i)
    {
        linkfile >> x >> y >> z;
        vector3_type position(x, y, z);
		particle_iterator p = surface.create_particle(particle_type());
		p->position() = position;
//std::cout << p->position() << std::endl;

    }
  linkfile.close();

}
}
#endif
// kate: indent-mode cstyle; space-indent on; indent-width 4; 
