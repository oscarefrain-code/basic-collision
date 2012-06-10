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
 * Tests for Triangle/Triangle class
 *
 */

#include <collision-two-triangles.h>
#include <iostream>
#include <Eigen/Dense>

void testTriangle( void )
{
  Eigen::Vector3d V1(0.0, 0.0, 0.0);
  Eigen::Vector3d V2(1.0, 0.0, 0.0);
  Eigen::Vector3d V3(1.0, 1.0, 0.0);

  Triangle T1(V1, V2, V3);
  T1.print();

  Eigen::Vector3d u1, u2, u3;
  T1.getVertices(u1, u2, u3);

  std::cout << u2.transpose() << std::endl;

}

void evaluateTTCollision( CollisionTwoTriangles &scene1 )
{
  int collisionDetected;
  collisionDetected = scene1.computeTTintersections();
  scene1.printTrianglesVertices();
  scene1.printTTcollisionInformation();
  std::cout << std::endl;
}

void testTriangleCollision( void )
{
  Eigen::Vector3d V1, V2, V3;

  V1 << 0, 0, 0;  V2 << 1, 0, 0; V3 << 1, 1, 0;
  Triangle T1(V1, V2, V3);

  V1 << 2, 2, 0;  V2 << 3, 2, 1; V3 << 2, 3, 2;
  Triangle T2(V1, V2, V3);

  CollisionTwoTriangles scene;
  scene.setTriangles(T1,T2);
  evaluateTTCollision( scene );

  V1 << 0, 2, 0;  V2 << 1, 2, 0; V3 << 1, 3, 0;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  V1 << 0, 0.5, -0.5;  V2 << 1, 0.5, -0.5; V3 << 1, 0.5, 0.5;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

}


int main( void )
{

  std::cout << "Testing the class Triangle" << std::endl;
  testTriangle();

  std::cout << " --- Testing triangle collision --- " << std::endl << std::endl;
  testTriangleCollision();

  return 0;
}


