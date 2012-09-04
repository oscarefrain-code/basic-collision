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

#include <Eigen/Dense>
using namespace Eigen;


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


/* -- Set the coplanar tolerance (of the triangles) -- */
void CollisionTwoRectangles::
setCoplanarTolerance(const double tol)
{
  scene1.setCoplanarTolerance(tol);
  scene2.setCoplanarTolerance(tol);
  scene3.setCoplanarTolerance(tol);
  scene4.setCoplanarTolerance(tol);
}

double CollisionTwoRectangles::
getCoplanarTolerance( void )
{
  double t1, t2, t3, t4;
  t1 = scene1.getCoplanarTolerance();
  t2 = scene2.getCoplanarTolerance();
  t3 = scene3.getCoplanarTolerance();
  t4 = scene4.getCoplanarTolerance();
  if ((t1==t2) && (t2==t3) && (t3==t4))
    return (t1);
  else
    std::cerr << "Two Rectangles: Coplanar tolerance is not the same "
	      << "for all the triangles !!!" << std::endl;
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

  /* Initialization to zero of collisions indicators */
  for (int i=0; i<4; i++)
    collisionTT[i]=0;

  /* Get the vertices */
  Vector3d R1v1, R1v2, R1v3, R1v4, R2v1, R2v2, R2v3, R2v4;
  R1.getVertices(R1v1, R1v2, R1v3, R1v4);
  R2.getVertices(R2v1, R2v2, R2v3, R2v4);

  /* Clear the possible 'previous' contact points */
  pointsRR.clear();

  /* Divide the rectangle/rectangle problem in 4 triangle/triangle problem  */
  scene1.setVerticesAll(R1v1, R1v2, R1v3, R2v1, R2v2, R2v3);
  scene2.setVerticesAll(R1v1, R1v2, R1v3, R2v3, R2v4, R2v1);
  scene3.setVerticesAll(R1v3, R1v4, R1v1, R2v1, R2v2, R2v3);
  scene4.setVerticesAll(R1v3, R1v4, R1v1, R2v3, R2v4, R2v1);

  //scene1.printTrianglesVertices();
  collisionTT[0] = scene1.computeTTintersections();    // R1T1 - R2R1
  collisionTT[1] = scene2.computeTTintersections();    // R1T1 - R2T2
  collisionTT[2] = scene3.computeTTintersections();    // R1T2 - R2T1
  collisionTT[3] = scene4.computeTTintersections();    // R1T2 - R2T2

  collisionIndicator = collisionTT[0] || collisionTT[1] ||
                       collisionTT[2] || collisionTT[3];

  // std::cout << std::endl;
  // scene1.printTrianglesVertices();
  // scene1.printTTcollisionInformation();
  // std::cout << std::endl;
  // scene2.printTrianglesVertices();
  // scene2.printTTcollisionInformation();
  // std::cout << std::endl;
  // scene3.printTrianglesVertices();
  // scene3.printTTcollisionInformation();
  // std::cout << std::endl;
  // scene4.printTrianglesVertices();
  // scene4.printTTcollisionInformation();

  // printRectanglesVertices();

  // std::cout << "Collisions (triangle/triangle level):\n      ";
  // for (int i=0; i<4; i++)
  //   std::cout << collisionTT[i] << " ";
  // std::cout << std::endl;

  // for (int i=0; i<scene1.pointsTT.size(); i++)
  //   std::cout << "  (" << scene1.pointsTT[i].transpose() << ")  " << std::endl;
  // for (int i=0; i<scene2.pointsTT.size(); i++)
  //   std::cout << "  (" << scene2.pointsTT[i].transpose() << ")  " << std::endl;
  // for (int i=0; i<scene3.pointsTT.size(); i++)
  //   std::cout << "P" << i+1 << "= [" << scene3.pointsTT[i].transpose() << "]; " << std::endl;
  // for (int i=0; i<scene4.pointsTT.size(); i++)
  //   std::cout << "P" << i+1 << "= [" << scene4.pointsTT[i].transpose() << "]; " << std::endl;


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

  // printRRcollisionInformation();

  prunePoints(pointsRR);
  convexHull(pointsRR);

  // printRRcollisionInformation();

  return collisionIndicator;
}


/* --- Print verbose information relative to the collision --- */
void CollisionTwoRectangles::
printRRcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Rectangles are intersecting at" << std::endl;
    for (int i=0; i<pointsRR.size(); i++){
      //std::cout << "  [" << pointsRR[i].transpose() << "]  ";
      std::cout << "  P" << i+1 << "=[" << pointsRR[i](0) << ", " << pointsRR[i](1) << ", "
		<< pointsRR[i](2) << "];";
    }
    std::cout << std::endl;
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

