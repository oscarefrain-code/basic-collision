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
  collisionIndicator = 0;
}


/* --- Constructor that initializes the vertices of both rectangles --- */
CollisionTwoRectangles::
CollisionTwoRectangles(Vector3d r1v1, Vector3d r1v2, Vector3d r1v3, Vector3d r1v4,
		       Vector3d r2v1, Vector3d r2v2, Vector3d r2v3, Vector3d r2v4)
{
  R1.setVertices(r1v1, r1v2, r1v3, r1v4);
  R2.setVertices(r2v1, r2v2, r2v3, r2v4);
  collisionIndicator = 0;
}


/* --- Constructor that initializes the vertices of both rectangles --- */
CollisionTwoRectangles::
CollisionTwoRectangles(Rectangle r1in, Rectangle r2in)
{
  R1=r1in; R2=r2in;
  collisionIndicator = 0;
}


/* --- Set the vertices of both rectangles --- */
void CollisionTwoRectangles::
setVerticesAll(Vector3d r1v1, Vector3d r1v2, Vector3d r1v3, Vector3d r1v4,
	       Vector3d r2v1, Vector3d r2v2, Vector3d r2v3, Vector3d r2v4)
{
  R1.setVertices(r1v1, r1v2, r1v3, r1v4);
  R2.setVertices(r2v1, r2v2, r2v3, r2v4);
}


/* --- Set the vertices of rectangle 1 --- */
void CollisionTwoRectangles::
setVerticesR1(Vector3d V1in, Vector3d V2in, Vector3d V3in, Vector3d V4in)
{
  R1.setVertices(V1in, V2in, V3in, V4in);
}


/* --- Set the vertices of rectangle 2 --- */
void CollisionTwoRectangles::
setVerticesR2(Vector3d V1in, Vector3d V2in, Vector3d V3in, Vector3d V4in)
{
  R2.setVertices(V1in, V2in, V3in, V4in);
}


/* --- Set both rectangles --- */
void CollisionTwoRectangles::
setRectangles(Rectangle r1in, Rectangle r2in)
{
  R1=r1in; R2=r2in;
}


/* --- Set rectangle 1 --- */
void CollisionTwoRectangles::
setRectangle1(Rectangle r1in)
{
  R1=r1in;
}


/* --- Set rectangle 2 --- */
void CollisionTwoRectangles::
setRectangle2(Rectangle r2in)
{
  R2=r2in;
}


/* --- Get the vertices of rectangle 1 --- */
void CollisionTwoRectangles::
getVerticesR1(Vector3d &V1out, Vector3d &V2out, Vector3d &V3out, Vector3d &V4out)
{
  R1.getVertices(V1out, V2out, V3out, V4out);
}


/* --- Get the vertices of rectangle 2 --- */
void CollisionTwoRectangles::
getVerticesR2(Vector3d &V1out, Vector3d &V2out, Vector3d &V3out, Vector3d &V4out)
{
  R2.getVertices(V1out, V2out, V3out, V4out);
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
  CollisionTwoTriangles scene1(R1.v1, R1.v2, R1.v3, R2.v1, R2.v2, R2.v3);
  CollisionTwoTriangles scene2(R1.v1, R1.v2, R1.v3, R2.v3, R2.v4, R2.v1);
  CollisionTwoTriangles scene3(R1.v3, R1.v4, R1.v1, R2.v1, R2.v2, R2.v3);
  CollisionTwoTriangles scene4(R1.v3, R1.v4, R1.v1, R2.v3, R2.v4, R2.v1);

  collisionTT[0] = scene1.computeTTintersections();
  collisionTT[1] = scene2.computeTTintersections();
  collisionTT[2] = scene3.computeTTintersections();
  collisionTT[3] = scene4.computeTTintersections();

  collisionIndicator = collisionTT[0] || collisionTT[1] ||
                       collisionTT[2] || collisionTT[3];

  // std::cout << "Collisions : ";
  // for (int i=0; i<4; i++)
  //   std::cout << collisionTT[i] << " ";
  // std::cout << std::endl;

  // for (int i=0; i<scene1.pointsTT.size(); i++)
  //   scene1.pointsTT[i].print();
  // std::cout << std::endl;
  // for (int i=0; i<scene2.pointsTT.size(); i++)
  //   scene2.pointsTT[i].print();
  // std::cout << std::endl;
  // for (int i=0; i<scene3.pointsTT.size(); i++)
  //   scene3.pointsTT[i].print();
  // std::cout << std::endl;
  // for (int i=0; i<scene4.pointsTT.size(); i++)
  //   scene4.pointsTT[i].print();
  // std::cout << std::endl;


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
    std::cout << "Result: Rectangles are intersecting at" << std::endl;
    for (int i=0; i<pointsRR.size(); i++){
      std::cout << "  (" << pointsRR[i].transpose() << ")  ";
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
  std::cout << "Rectangle 1: "; R1.print();
  std::cout << "Rectangle 2: "; R2.print();
}

