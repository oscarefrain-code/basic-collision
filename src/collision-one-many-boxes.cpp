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


/* -- Constructors -- */
CollisionOneManyBoxes::
CollisionOneManyBoxes( void )
{
  id = 0;
  // Initialize the values to an arbitrary default
  Lx1=0.0; Ly1=0.0; Lz1=0.0;
  R1[0]=1.0; R1[1]=0.0; R1[2]=0.0;
  R1[3]=0.0; R1[4]=1.0; R1[5]=0.0;
  R1[6]=0.0; R1[7]=0.0; R1[8]=1.0;
  T1[0]=0.0; T1[1]=0.0; T1[2]=0.0;
  tolerance = TRI_EPSILON;
  collisionIndicator = 0;
}

/* -- Tolerance related -- */
void CollisionOneManyBoxes::
setTolerance(double toler)
{
  tolerance = toler;
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
setTransformationB1(double R[9], double T[3])
{
  R1[0]=R[0]; R1[1]=R[1]; R1[2]=R[2];
  R1[3]=R[3]; R1[4]=R[4]; R1[5]=R[5];
  R1[6]=R[6]; R1[7]=R[7]; R1[8]=R[8];
  T1[0]=T[0]; T1[1]=T[1]; T1[2]=T[2];

  // Update the value of all the box-box pairs
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionBoxes[i].setTransformation1(R1, T1);
  }

}

// Once the object Bi is added, its transformation cannot be modified. It is
// thus assumed to be a fixed body forever.
int CollisionOneManyBoxes::
addBox(double lx, double ly, double lz, double R[9], double T[3])
{
  CollisionTwoBoxes temporalPair(Lx1, Ly1, Lz1, lx, ly, lz);
  temporalPair.setTolerance(tolerance);
  temporalPair.setTransformation1(R1, T1);
  temporalPair.setTransformation2(R, T);
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
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionBoxes[i].setTolerance(tolerance);
  }
  
  int collisionTemp; 
  collisionIndicator = 0;
  for (int i=0; i<collisionBoxes.size(); i++){
    collisionTemp = collisionBoxes[i].computeBBintersections();
    collisionIndicator = collisionIndicator || collisionTemp;
  }

  Point3d tempPoint;
  if (collisionIndicator){
    for (int i=0; i<collisionBoxes.size(); i++){
      for (int k=0; k<collisionBoxes[i].pointsBB.size(); k++){
	tempPoint.x = collisionBoxes[i].pointsBB[k].x;
	tempPoint.y = collisionBoxes[i].pointsBB[k].y;
	tempPoint.z = collisionBoxes[i].pointsBB[k].z;
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
    std::cout << " - Intersections occur at: " << std::endl;
    for (int i=0; i<pointsBBs.size(); i++){
      std::cout << "     "; pointsBBs[i].print();
      std::cout << std::endl;
    }
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
      std::cout << "Lengths: Box 1 = (" << Lx1 << ", " << Ly1 << ", " << Lz1 << " )" << std::endl;
    }
    collisionBoxes[i].getLengthB2(tempLx, tempLy, tempLz);
    std::cout << "         Box " << i+2 << " = (" << tempLx << ", " << tempLy
	      << ", " << tempLz << " )" << std::endl;
  }
}
 
void CollisionOneManyBoxes::
printBoxesVertices( void )
{
  double V1[3], V2[3], V3[3], V4[3], V5[3], V6[3], V7[3], V8[3];
  std::cout << "Vertices: ";
  for (int i=0; i<collisionBoxes.size(); i++){
    if (i==0) {
      collisionBoxes[i].getVerticesB1(V1, V2, V3, V4, V5, V6, V7, V8);
      std::cout << "\n  - Box 1: "; printVector(V1); printVector(V2); printVector(V3); printVector(V4);
      std::cout << std::endl
		<< "           "; printVector(V5); printVector(V6); printVector(V7); printVector(V8);
    }
    collisionBoxes[i].getVerticesB2(V1, V2, V3, V4, V5, V6, V7, V8);
    std::cout << "\n  - Box " << i+2 << ": "; printVector(V1); printVector(V2); printVector(V3); printVector(V4);
    std::cout << std::endl
	      << "           "; printVector(V5); printVector(V6); printVector(V7); printVector(V8);
    std::cout << std::endl;
  }

}

