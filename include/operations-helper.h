/***********************************************************************************
 * Copyright 2012, Oscar Efrain Ramos Ponce, LAAS-CNRS 
 *
 * This file is part of basic-collision
 * basic-collision is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * basic-collision is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Lesser Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License
 * along with basic-collision.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************************/

/*
 * Basic components and class to be used for the other classes
 *
 */


#ifndef __OPERATIONS_HELPER_H__
#define __OPERATIONS_HELPER_H__

#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>


/* Definition of a 'minimum value' to be used as a threshold for the distance
   Two points are considered to be the same if their distance is smaller than DIST_TOLERANCE*/
#define DIST_TOLERANCE 0.001

// Tolerance for the collinear distance (when the vectors are collinear)
#define tolerance_distance_collinear 0.00001

/* Tolerance distance to prune the points (if the distance is smaller than this, then, ignore one
   of the points */
#define tolerance_distance 0.001


class OperationsHelper
{

 public:

  void sort2(Eigen::Vector2d &vect, int &smallest);
  void pushVectorAsPoint(Eigen::Vector3d v, std::vector<Eigen::Vector3d>& points);
  void prunePoints(std::vector<Eigen::Vector3d>& points);
  int convexHull( std::vector<Eigen::Vector3d>& points );

 /* private: */
 /*  double tolerance_distance; */

};


#endif
