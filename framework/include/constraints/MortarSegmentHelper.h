//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "libmesh/point.h"

#include <vector>

using libMesh::Point;
using libMesh::Real;

/**
 * This class supports defining mortar segment mesh elements in 3D by projecting secondary and
 * primary elements onto a linearized plane, computing the overlapping polygon formed by their
 * projections, and triangulating the resulting nodes.
 */

class MortarSegmentHelper
{
public:
  MortarSegmentHelper(const std::vector<Point> secondary_nodes,
                      const Point & center,
                      const Point & normal);

  /**
   * Computes the intersection between line segments defined by point pairs (p1,p2) and (q1,q2)
   * Also computes s, the ratio of distance between (p1,p2) that the intersection falls,
   * quantity s is useful in avoiding adding nearly degenerate nodes
   */
  Point getIntersection(
      const Point & p1, const Point & p2, const Point & q1, const Point & q2, Real & s) const;

  /**
   * Check that a point is inside the secondary polygon (for verification only)
   */
  bool isInsideSecondary(const Point & pt) const;

  /**
   * Checks whether polygons are disjoint for an easy out
   */
  bool isDisjoint(const std::vector<Point> & poly) const;

  /**
   * Clip secondary element (defined in instantiation) against given primary polygon
   * result is a set of 2D nodes defining clipped polygon
   */
  std::vector<Point> clipPoly(const std::vector<Point> & primary_nodes) const;

  /**
   * Triangulate a polygon (currently uses center of polygon to define triangulation)
   * @param poly_nodes List of 2D nodes defining polygon
   * @param offset Current size of 3D nodes array (not poly_nodes)
   * @return tri_map List of integer arrays defining which nodes belong to each triangle
   */
  void triangulatePoly(std::vector<Point> & poly_nodes,
                       const unsigned int offset,
                       std::vector<std::vector<unsigned int>> & tri_map) const;

  /**
   * Get mortar segments generated by a secondary and primary element pair
   * @param primary_nodes List of primary element 3D nodes
   * @return nodes List of 3D mortar segment nodes
   * @return tri_map List of integer arrays defining which nodes belong to each mortar segment
   */
  void getMortarSegments(const std::vector<Point> & primary_nodes,
                         std::vector<Point> & nodes,
                         std::vector<std::vector<unsigned int>> & elem_to_nodes);

  /**
   * Compute area of polygon
   */
  Real area(const std::vector<Point> & nodes) const;

  /**
   * Get center point of secondary element
   */
  const Point & center() const { return _center; }

  /**
   * Get area fraction remaining after clipping against primary elements
   */
  Real remainder() const { return _remaining_area_fraction; }

  /**
   * Get 3D position of node of linearized secondary element
   */
  Point point(unsigned int i) const
  {
    return (_secondary_poly[i](0) * _u) + (_secondary_poly[i](1) * _v) + _center;
  }

private:
  /**
   * Geometric center of secondary element
   */
  Point _center;

  /**
   * Normal at geometric center of secondary element
   */
  Point _normal;

  /**
   * Vectors orthogonal to normal that span the plane projection will be performed on.
   * These vectors are used to project the polygon clipping problem on a 2D plane,
   * they are defined so the nodes of the projected polygon are listed with positive orientation
   */
  Point _u, _v;

  /**
   * Area of projected secondary element
   */
  Real _secondary_area;

  /**
   * Fraction of area remaining after overlapping primary polygons clipped
   */
  Real _remaining_area_fraction;

  bool _debug;

  /**
   * Tolerance for intersection and clipping
   */
  Real _tolerance = 1e-8;

  /**
   * Tolerance times secondary area for dimensional consistency
   */
  Real _area_tol;

  /**
   * Tolerance times secondary area for dimensional consistency
   */
  Real _length_tol;

  /**
   * List of projected points on the linearized secondary element
   */
  std::vector<Point> _secondary_poly;
};
