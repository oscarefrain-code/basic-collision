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


void evaluateTTCollision( CollisionTwoTriangles &scene1 );

int main( void )
{

  Eigen::Vector3d V1, V2, V3;

  V1 << 0, 0, 0;  V2 << 1, 0, 0; V3 << 1, 1, 0;
  Triangle T1(V1, V2, V3);

  V1 << 2, 2, 0;  V2 << 3, 2, 1; V3 << 2, 3, 2;
  Triangle T2(V1, V2, V3);

  CollisionTwoTriangles scene;
  scene.setTriangles(T1,T2);
  evaluateTTCollision( scene );

  /* No intersection (coplanar) */
  V1 << 0, 2, 0;  V2 << 1, 2, 0; V3 << 1, 3, 0;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  /* Intersection: triangle is perpendicular */
  /*     (0.5,0.5,0) and (1,0.5,0)           */ 
  V1 << 0, 0.5, -0.5;  V2 << 1, 0.5, -0.5; V3 << 1, 0.5, 0.5;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  /* Intersection: One point (perpendicular) */
  /*   (0.5, 0.25, 0.0)                      */
  V1 << 0.5, 0.25, 0;  V2 << 1, 0, 1; V3 << 0, 1, 1;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  /* Intersection: coplanar (one part intersecting) */
  /* ( 1, 0, 0), ( 1, 1, 0) and ( 0.5, 0.5, 0)      */
  V1 << 0, 1, 0;  V2 << 1, 0, 0; V3 << 1, 1, 0;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  /* Intersection: coplanar (one absolutely contained in the other) */
  /* ( 0, 0, 0), ( 1, 0, 0), ( 1, 1, 0)                             */
  V1 << -2, -1, 0;  V2 << 1, 2, 0; V3 << 3, -1, 0;
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  // /* Miscellaneous tests */
  // setVector(T1V1, 0.0, 3.0, 0.0);
  // setVector(T1V2, 3.0, 0.0, 0.0);
  // setVector(T1V3, 3.0, 3.0, 0.0);
  // setVector(T2V1, -1.0, 2.0, 0.0);
  // setVector(T2V2, 2.0, -1.0, 0.0);
  // setVector(T2V3, 2.0, 2.0, 0.0);
  // evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  std::cout << "Last Test --- " << std::endl;
  //V1 << 0, 0, 0;  V2 << 1, 0, 0; V3 << 0, 1, 0;
  V1 << -0.135861, -0.00388724, 0.08; V2 << 0.0941395, -0.00388724, 0.08; V3 << 0.0941395, -0.138887,0.08; 
  scene.setVerticesT1(V1, V2, V3);
  //V1 << 0.4, 0.4, 0;  V2 << 1, 1, 1; V3 << 1, 1, -1;
  //V1 << 1.1, 0, 0;  V2 << 0, 1.1, 0; V3 << 2, 2, 0.01;
  V1 << 0.125, -0.1875, 0.08; V2 << -0.065, -0.1875, 0.08; V3<< -0.065, -0.0725, 0.08001; 
  scene.setVerticesT2(V1, V2, V3);
  evaluateTTCollision( scene );

  return 0;
}


void evaluateTTCollision( CollisionTwoTriangles &scene1 )
{
  int collisionDetected;
  collisionDetected = scene1.computeTTintersections();
  scene1.printTrianglesVertices();
  scene1.printTTcollisionInformation();
  std::cout << std::endl;
}



