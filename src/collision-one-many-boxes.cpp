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

#include <collision-one-many-boxes.h>
#include <iostream>

#include <Eigen/Dense>
using namespace Eigen;

/* -- Constructors -- */
CollisionOneManyBoxes::
CollisionOneManyBoxes( void )
{
  id = 0;
  // Initialize the values to an arbitrary default
  Lx1=0.0; Ly1=0.0; Lz1=0.0;
  R1.setIdentity();  T1.setZero();
  tolerance = 0.0005;         //TRI_EPSILON;
  collisionIndicator = 0;
}

/* -- Tolerance related -- */
void CollisionOneManyBoxes::
setTolerance(double toler)
{
  tolerance = toler;
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionBoxes[i].setTolerance(tolerance);
  }

}

double CollisionOneManyBoxes::
getTolerance( void )
{
  return tolerance;
}


/* -- Set lengths for the box 1 -- */
void CollisionOneManyBoxes::
setLengthB1(double lx1, double ly1, double lz1)
{
  Lx1=lx1; Ly1=ly1; Lz1=lz1;
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionBoxes[i].setLengthB1(Lx1, Ly1, Lz1);
  }

}

/* -- Get lengths for the box 1 -- */
void CollisionOneManyBoxes::
getLengthB1(double & lx1, double & ly1, double & lz1)
{
  lx1=Lx1; ly1=Ly1; lz1=Lz1;
}


/* -- Set transformation for box 1 -- */
void CollisionOneManyBoxes::
setTransformationB1(Matrix3d Rin, Vector3d Tin)
{
  R1 = Rin; T1= Tin;

  // Update the value of all the box-box pairs
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionBoxes[i].setTransformation1(R1, T1);
  }

}

// Once the object Bi is added, its transformation cannot be modified. It is
// thus assumed to be a fixed body forever.
int CollisionOneManyBoxes::
addBox(double lx, double ly, double lz, Matrix3d Rin, Vector3d Tin)
{
  CollisionTwoBoxes temporalPair(Lx1, Ly1, Lz1, lx, ly, lz);
  temporalPair.setTolerance(tolerance);
  temporalPair.setTransformation1(R1, T1);
  temporalPair.setTransformation2(Rin, Tin);
  collisionBoxes.push_back(temporalPair);
  
  return (id++);
}


/* -- Compute the intersections -- */
int CollisionOneManyBoxes::
computeIntersections( void )
{
  /* Clear the possible 'previous' contact points */
  pointsBBs.clear();

  // Update the tolerance of each pair of boxes
  // for (int i=0; i<collisionBoxes.size(); i++){
  //   collisionBoxes[i].setTolerance(tolerance);
  // }
  
  int collisionTemp; 
  collisionIndicator = 0;
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionTemp = collisionBoxes[i].computeBBintersections();
    collisionIndicator = collisionIndicator || collisionTemp;
  }

  //Point3d tempPoint;
  Vector3d tempPoint;
  if (collisionIndicator){
    for (int i=0; i<collisionBoxes.size(); i++){
      for (int k=0; k<collisionBoxes[i].pointsBB.size(); k++){
	tempPoint(0) = collisionBoxes[i].pointsBB[k](0);
	tempPoint(1) = collisionBoxes[i].pointsBB[k](1);
	tempPoint(2) = collisionBoxes[i].pointsBB[k](2);
	pointsBBs.push_back(tempPoint);
      }
    }
  }
  
  return collisionIndicator;

}


/* --- Print verbose information relative to the collision --- */
void CollisionOneManyBoxes::
printBBsCollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Boxes are intersecting" << std::endl;
    for (int i=0; i<pointsBBs.size(); i++){
      std::cout << "  (" << pointsBBs[i].transpose() << ")  ";
    }
    std::cout << std::endl;
  }
  else{
    std::cout << "Result: Boxes are NOT intersecting" << std::endl;
  }
  
}

void CollisionOneManyBoxes::
printLengths( void )
{
  double tempLx, tempLy, tempLz;
  for (int i=0; i<collisionBoxes.size(); i++){
    if (i==0) {
      collisionBoxes[i].getLengthB1(tempLx, tempLy, tempLz);
      std::cout << "Length box 1 = (" << Lx1 << ", " << Ly1 << ", " << Lz1 << " )" << std::endl;
    }
    collisionBoxes[i].getLengthB2(tempLx, tempLy, tempLz);
    std::cout << "Length box " << i+2 << " = (" << tempLx << ", " << tempLy
	      << ", " << tempLz << " )" << std::endl;
  }
}
 
void CollisionOneManyBoxes::
printBoxesVertices( void )
{
  Vector3d V1, V2, V3, V4, V5, V6, V7, V8;
  for (int i=0; i<collisionBoxes.size(); i++){
    if (i==0) {
      collisionBoxes[i].getVerticesB1(V1, V2, V3, V4, V5, V6, V7, V8);
      std::cout << "Vertices box 1: [" << V1.transpose() << "], [" << V2.transpose() << "], [" << V3.transpose() << "],[" << V4.transpose() << "],["
		<< V5.transpose() << "], [" << V6.transpose() << "], [" << V7.transpose() << "],[" << V8.transpose() << "]" << std::endl;
    }
    collisionBoxes[i].getVerticesB2(V1, V2, V3, V4, V5, V6, V7, V8);
    std::cout << "Vertices box "<< i+2 << ": [" << V1.transpose() << "], [" << V2.transpose() << "], [" << V3.transpose() << "],[" << V4.transpose() << "],["
	      << V5.transpose() << "], [" << V6.transpose() << "], [" << V7.transpose() << "],[" << V8.transpose() << "]" << std::endl;
  }

}

