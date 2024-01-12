//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AllSideSetsByNormalsGenerator.h"
#include "Parser.h"
#include "InputParameters.h"
#include "CastUniquePointer.h"

#include "libmesh/fe_base.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"

#include <typeinfo>
#include "MooseMeshUtils.h"
#include "MooseTypes.h"

registerMooseObject("MooseApp", AllSideSetsByNormalsGenerator);

InputParameters
AllSideSetsByNormalsGenerator::validParams()
{
  InputParameters params = SideSetsGeneratorBase::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addClassDescription("Adds sidesets to the entire mesh based on unique normals.");
  params.addParam<std::vector<std::string>>("boundary_name_prefix",
                               "If provided, prefix the built in boundary names with this string");
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The blocks around which to create sidesets");
  params.addParam<std::vector<std:: unsigned int>>("offset",
                                              "offset for sideset numbers");

  return params;
}



AllSideSetsByNormalsGenerator::AllSideSetsByNormalsGenerator(const InputParameters & parameters)
  : SideSetsGeneratorBase(parameters),
    _input(getMesh("input")),
    _boundary_name_prefix(getParam<std::vector<std::string>>("boundary_name_prefix")),
    _has_subdomain_ids(isParamValid("block")),
    _has_offset(isParamValid("offset")),
    _offset(getParam<std::vector<std:: unsigned int>>("offset")),
    _boundary_to_normal_map(
            declareMeshProperty<std::map<BoundaryID, RealVectorValue>>("boundary_normals"))
{

}

std::unique_ptr<MeshBase>
AllSideSetsByNormalsGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  std::vector<SubdomainName> block_names = getParam<std::vector<SubdomainName>>("block");
  // check that the blocks exist in the mesh
  for (const auto & name : block_names)
    if (!MooseMeshUtils::hasSubdomainName(*mesh, name))
      paramError("block", "The block '", name, "' was not found in the mesh");

  auto blocks = MooseMeshUtils::getSubdomainIDs(*mesh, block_names);
  std::set<subdomain_id_type> block_ids(blocks.begin(), blocks.end());

    // Create the boundary IDs from the list of names provided (the true flag creates ids from unknown
  // names)
  // std::vector<boundary_id_type> boundary_ids =
  //     MooseMeshUtils::getBoundaryIDs(*mesh, _boundary_names, true);


  if (!mesh->is_replicated())
    mooseError("AllSideSetsByNormalsGenerator is not implemented for distributed meshes");
  setup(*mesh);


  // Get the current list of boundaries so we can generate new ones that won't conflict
  _mesh_boundary_ids = mesh->get_boundary_info().get_boundary_ids();

  _visited.clear();



  // We'll need to loop over all of the elements to find ones that match this normal.
  // We can't rely on flood catching them all here...
  for (const auto & elem : mesh->element_ptr_range())
  {
    subdomain_id_type curr_subdomain = elem->subdomain_id();

    // // We only need to loop over elements in the source subdomain
    // if (_has_subdomain_ids && block_ids.count(curr_subdomain) == 0)
    //   continue;

    for (unsigned int side = 0; side < elem->n_sides(); ++side)
    {
      // auto neighbor = elem->neighbor_ptr(side);


      if (elem->neighbor_ptr(side) && !_has_subdomain_ids)
        continue;
      else if (_has_subdomain_ids && elem->neighbor_ptr(side) != NULL &&
               block_ids.count(elem->neighbor_ptr(side)->subdomain_id())== 0)
        continue;
      

      const std::vector<Point> & normals = _fe_face->get_normals();
      _fe_face->reinit(elem, side);

      {
        // See if we've seen this normal before (linear search)
        const std::map<BoundaryID, RealVectorValue>::value_type * item = nullptr;
        for (const auto & id_pair : _boundary_to_normal_map)
          if (std::abs(1.0 - id_pair.second * normals[0]) < 1e-5)
          {
            item = &id_pair;
            break;
          }

        if (item)
          flood(elem, normals[0], item->first, *mesh);
        else
        {
          boundary_id_type id = getNextBoundaryID();
          _boundary_to_normal_map[id] = normals[0];
          flood(elem, normals[0], id, *mesh);
        }
      }
    }
  }

  // Rename the sidesets
    BoundaryInfo & boundary_info = mesh->get_boundary_info();

    for (unsigned int i = 0; i < _block_names.size(); ++i)
    {
      // Copy, since we're modifying the container mid-iteration
      const auto mesh_boundary_ids = boundary_info.get_global_boundary_ids();
      for (auto rit = mesh_boundary_ids.rbegin(); rit != mesh_boundary_ids.rend(); ++rit)
      {
        const std::string old_sideset_name = boundary_info.sideset_name(*rit);
        const std::string old_nodeset_name = boundary_info.nodeset_name(*rit);

        boundary_info.sideset_name(*rit) = _boundary_name_prefix[i] + old_sideset_name;
        boundary_info.nodeset_name(*rit) = _boundary_name_prefix[i] + old_nodeset_name;
      }
    }


  finalize();

  mesh->set_isnt_prepared();
  return dynamic_pointer_cast<MeshBase>(mesh);
}

boundary_id_type
AllSideSetsByNormalsGenerator::getNextBoundaryID()
{
  std::set<boundary_id_type>::iterator it;
  boundary_id_type next_id = 1;

  while ((it = _mesh_boundary_ids.find(next_id)) != _mesh_boundary_ids.end())
    ++next_id;

  _mesh_boundary_ids.insert(next_id);

  return next_id;
}

// boundary_id_type
// AllSideSetsByNormalsGenerator::getOffsetBoundaryID()
// {
//   std::set<boundary_id_type>::iterator it;
//   boundary_id_type next_id = 1;

//   while ((it = _mesh_boundary_ids.find(next_id)) != _mesh_boundary_ids.end())
//     ++next_id;

//   _mesh_boundary_ids.insert(next_id+_offse);

//   return next_id;
// }
