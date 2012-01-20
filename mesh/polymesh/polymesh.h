#ifndef OPENTISSUE_CORE_CONTAINERS_MESH_POLYMESH_POLYMESH_H
#define OPENTISSUE_CORE_CONTAINERS_MESH_POLYMESH_POLYMESH_H
//
// OpenTissue Template Library
// - A generic toolbox for physics-based modeling and simulation.
// Copyright (C) 2008 Department of Computer Science, University of Copenhagen.
//
// OTTL is licensed under zlib: http://opensource.org/licenses/zlib-license.php
//

#include "mesh/polymesh/polymesh_mesh.h"
#include "mesh/mesh_default_traits.h"
#include "mesh/polymesh/kernels/polymesh_kernels.h"

namespace polymesh
{

    /**
    * This is basically just another name for polymesh::detail::mesh
    *
    * The first template argument is supposed to be a math types type binder.
    * OpenTissue provides a simple basic math type-binder in the math sub-library.
    * See OpenTissue::math::BasicMathTypes<real_type, size_type>
    * The next four template arguments are supposed to be vertex traits,
    * halfedge traits, edge traits, and face traits. The last template
    * argument is the polymesh kernel type that is supposed to be used.
    *
    */
    template <
    typename M
    , typename V = mesh::DefaultVertexTraits< M >
    , typename H = mesh::DefaultHalfEdgeTraits
    , typename E = mesh::DefaultEdgeTraits
    , typename F = mesh::DefaultFaceTraits
    , template <typename, typename, typename, typename> class K = PolyMeshListKernel
    >
    class PolyMesh
        : public detail::PMesh<M, V, H, E, F, K>
    {};

} // namespace polymesh

//OPENTISSUE_CORE_CONTAINERS_MESH_POLYMESH_POLYMESH_H
#endif
