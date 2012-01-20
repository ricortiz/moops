#ifndef MESH_H
#define MESH_H
//
// OpenTissue Template Library
// - A generic toolbox for physics-based modeling and simulation.
// Copyright (C) 2008 Department of Computer Science, University of Copenhagen.
//
// OTTL is licensed under zlib: http://opensource.org/licenses/zlib-license.php
//
#include "mesh/polymesh/polymesh.h"
#include "mesh/polymesh/util/polymesh_util.h"

#include "mesh/trimesh/trimesh.h"
#include "mesh/trimesh/util/trimesh_util.h"

#include "mesh/common/util/mesh_compute_face_center.h"
#include "mesh/common/util/mesh_compute_mesh_center.h"
#include "mesh/common/util/mesh_compute_mesh_mean.h"
#include "mesh/common/util/mesh_compute_face_plane.h"
#include "mesh/common/util/mesh_compute_mesh_minimum_coord.h"
#include "mesh/common/util/mesh_compute_mesh_maximum_coord.h"
#include "mesh/common/util/mesh_compute_face_minimum_coord.h"
#include "mesh/common/util/mesh_compute_face_maximum_coord.h"

#include "mesh/common/util/mesh_profile_sweep.h"
#include "mesh/common/util/mesh_make_disk.h"
#include "mesh/common/util/mesh_make_cylinder.h"
#include "mesh/common/util/mesh_make_sphere.h"
#include "mesh/common/util/mesh_make_box.h"
#include "mesh/common/util/mesh_make_plane.h"

#include "mesh/common/util/mesh_deformation_modifiers.h"
#include "mesh/common/util/mesh_make_unit.h"
#include "mesh/common/util/mesh_compute_angle_weighted_vertex_normals.h"
#include "mesh/common/util/mesh_compute_mean_vertex_normals.h"

#include "mesh/common/util/mesh_plane_clipper.h"
#include "mesh/common/util/mesh_convex_hull.h"
#include "mesh/common/util/mesh_volume_integrator.h"
#include "mesh/common/util/mesh_compute_surface_covariance.h"

#include "mesh/common/util/mesh_convert.h"
#include "mesh/common/util/mesh_flip.h"
#include "mesh/common/util/mesh_isosurface.h"
#include "mesh/common/util/mesh_smooth_isosurface.h"

#include "mesh/common/io/mesh_default_read.h"
#include "mesh/common/io/mesh_default_write.h"
#include "mesh/common/io/mesh_tetgen_write.h"
#include "mesh/common/io/mesh_obj_read.h"
#include "mesh/common/io/mesh_obj_write.h"
#include "mesh/common/io/mesh_vrml_write.h"

#include "mesh/common/util/mesh_coordinate_iterator.h"

#include "mesh/common/util/mesh_clear_vertex_tags.h"
#include "mesh/common/util/mesh_clear_halfedge_tags.h"
#include "mesh/common/util/mesh_clear_face_tags.h"

#include "mesh/common/util/mesh_remove_redundant_vertices.h"
#include "mesh/common/util/mesh_compute_minmax_face_area.h"



// MESH_H
#endif
