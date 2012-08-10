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
 * One Hexahedron/Many hexahedra intersection test:
 *   The collision between one hexahedron "B1" and many hexahedra "Bi" is tested. It
 *   decomposes the problem in many hexahedron/hexahedron problems. The hexahedra 
 *   Bi is supposed not to move after fixing its initial position, but B1 can move.
 *  
 */

#ifndef __COLLISION_ONE_MANY_BOXES_H__
#define __COLLISION_ONE_MANY_BOXES_H__

#include <math.h>
#include <vector>
#include <operations-helper.h>
#include <collision-two-boxes.h>


class CollisionOneManyBoxes:
public OperationsHelper
{
 public:
  /* -- Constructors -- */
  CollisionOneManyBoxes( void );

  /* -- Set values -- */
  void setTolerance(double toler);
  void setLengthB1(double lx1, double ly1, double lz1);
  void setTransformationB1(Eigen::Matrix3d R, Eigen::Vector3d T);
  int addBox(double lx, double ly, double lz, Eigen::Matrix3d R, Eigen::Vector3d T);

  /* -- Compute the collision detection -- */
  int computeIntersections( void );

  /* -- Get some values -- */
  double getTolerance( void );
  void getLengthB1(double & lx1, double & ly1, double & lz1);

  /* --- Print Different Information --- */
  void printLengths( void );
  void printBoxesVertices( void );
  void printBBsCollisionInformation( void );

  /* -- Collision Points -- */
  std::vector<Eigen::Vector3d> pointsBBs;
  std::vector<CollisionTwoBoxes> collisionBoxes;

 private:
  double Lx1, Ly1, Lz1;
  Eigen::Matrix3d R1; Eigen::Vector3d T1;
  double tolerance;    // Tolerance to average the close points in collision-two-boxes
  int id, collisionIndicator;

};



#endif
