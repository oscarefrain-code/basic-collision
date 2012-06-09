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
 * Rectangle/rectangle intersection test:
 *   Each rectangle is decomposed into two triangles and the collision between 
 *   triangles is computed.   
 *
 *  
 */

#include <collision-two-rectangles.h>
#include <iostream>


/* --- Constructor that initializes all the vertices to zero --- */
CollisionTwoRectangles::
CollisionTwoRectangles( void )
{
  double zeros[3] = {0.0, 0.0, 0.0};
  set(V0, zeros); set(V1, zeros); set(V2, zeros); set(V3, zeros);
  set(U0, zeros); set(U1, zeros); set(U2, zeros); set(U3, zeros);
  collisionIndicator = 0;
}


/* --- Constructor that initializes the vertices of both rectangles --- */
CollisionTwoRectangles::
CollisionTwoRectangles(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3],
		       double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3])
{
  setVerticesAll(R1_V0, R1_V1, R1_V2, R1_V3, R2_V0, R2_V1, R2_V2, R2_V3);
  collisionIndicator = 0;
}


/* --- Set the vertices of both rectangles --- */
void CollisionTwoRectangles::
setVerticesAll(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3],
	       double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3])
{
  setVerticesR1(R1_V0, R1_V1, R1_V2, R1_V3);
  setVerticesR2(R2_V0, R2_V1, R2_V2, R2_V3);
}


/* --- Set the vertices of rectangle 1 --- */
void CollisionTwoRectangles::
setVerticesR1(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3])
{
  set(V0, R1_V0); set(V1, R1_V1); set(V2, R1_V2); set(V3, R1_V3);
}


/* --- Set the vertices of rectangle 2 --- */
void CollisionTwoRectangles::
setVerticesR2(double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3])
{
  set(U0, R2_V0); set(U1, R2_V1); set(U2, R2_V2); set(U3, R2_V3);
}


/* --- Get the vertices of rectangle 1 --- */
void CollisionTwoRectangles::
getVerticesR1(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3])
{
  set(R1_V0, V0); set(R1_V1, V1); set(R1_V2, V2); set(R1_V3, V3);
}


/* --- Get the vertices of rectangle 2 --- */
void CollisionTwoRectangles::
getVerticesR2(double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3])
{
  set(R2_V0, U0); set(R2_V1, U1); set(R2_V2, U2); set(R2_V3, U3);
}


/*================================================================================
                   COMPUTATION FUNCTIONS
 ================================================================================= */

int CollisionTwoRectangles::
computeRRintersections( void )
{
  int collisionTT[4];

  /* Clear the possible 'previous' contact points */
  pointsRR.clear();

  /* Divide the rectangle/rectangle problem in 4 triangle/triangle problem  */
  CollisionTwoTriangles scene1(V0, V1, V2, U0, U1, U2);
  CollisionTwoTriangles scene2(V0, V1, V2, U2, U3, U0);
  CollisionTwoTriangles scene3(V2, V3, V0, U0, U1, U2);
  CollisionTwoTriangles scene4(V2, V3, V0, U2, U3, U0);

  collisionTT[0] = scene1.computeTTintersections();
  collisionTT[1] = scene2.computeTTintersections();
  collisionTT[2] = scene3.computeTTintersections();
  collisionTT[3] = scene4.computeTTintersections();

  collisionIndicator = collisionTT[0] || collisionTT[1] ||
                       collisionTT[2] || collisionTT[3];

  std::cout << "Collisions : ";
  for (int i=0; i<4; i++)
    std::cout << collisionTT[i] << " ";
  std::cout << std::endl;

  for (int i=0; i<scene1.pointsTT.size(); i++)
    scene1.pointsTT[i].print();
  std::cout << std::endl;
  for (int i=0; i<scene2.pointsTT.size(); i++)
    scene2.pointsTT[i].print();
  std::cout << std::endl;
  for (int i=0; i<scene3.pointsTT.size(); i++)
    scene3.pointsTT[i].print();
  std::cout << std::endl;
  for (int i=0; i<scene4.pointsTT.size(); i++)
    scene4.pointsTT[i].print();
  std::cout << std::endl;


  /* Copy all the collision points in pointsRR  */
  if (collisionIndicator) {
    for (int i=0; i<scene1.pointsTT.size(); i++)
      pointsRR.push_back(scene1.pointsTT[i]);
    for (int i=0; i<scene2.pointsTT.size(); i++)
      pointsRR.push_back(scene2.pointsTT[i]);
    for (int i=0; i<scene3.pointsTT.size(); i++)
      pointsRR.push_back(scene3.pointsTT[i]);
    for (int i=0; i<scene4.pointsTT.size(); i++)
      pointsRR.push_back(scene4.pointsTT[i]);
  }

  prunePoints(pointsRR);
  convexHull(pointsRR);

  return collisionIndicator;
}


/* --- Print verbose information relative to the collision --- */
void CollisionTwoRectangles::
printRRcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Rectangles are intersecting" << std::endl;
    std::cout << " - Intersections occur at: " << std::endl;
    for (int i=0; i<pointsRR.size(); i++){
      std::cout << "     "; pointsRR[i].print();
      std::cout << std::endl;
    }
  }
  else{
    std::cout << "Result: Rectangles are NOT intersecting" << std::endl;
  }
}


/* --- Print vertices of the rectangles --- */
void CollisionTwoRectangles::
printRectanglesVertices( void )
{
  std::cout << "Vertices: ";
  std::cout << "\n  - Rectangle 1: "; printVector(V0); printVector(V1); printVector(V2); printVector(V2);
  std::cout << "\n  - Rectangle 2: "; printVector(U0); printVector(U1); printVector(U2); printVector(U3);
  std::cout << std::endl;

}

