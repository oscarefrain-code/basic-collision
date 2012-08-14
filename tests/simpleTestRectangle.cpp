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
 * Tests for Rectangle/Rectangle class
 * 
 */


#include <collision-two-rectangles.h>
#include <iostream>


void evaluateRRCollision( CollisionTwoRectangles &scene1);

int main( void )
{
  //double R1V1[3], R1V2[3], R1V3[3], R1V4[3], R2V1[3], R2V2[3], R2V3[3], R2V4[3];
  Eigen::Vector3d V1, V2, V3, V4;

  V1 << 0, 0, 0;  V2 << 0, 1, 0; V3 << 2, 1, 0; V4 << 2, 0, 0;
  Rectangle R1(V1, V2, V3, V4);

  V1 << 0, 2, 0;  V2 << 0, 3, 0; V3 << 2, 3, 0; V4 << 2, 2, 0;
  Rectangle R2(V1, V2, V3, V4);

  CollisionTwoRectangles scene;
  scene.setRectangles(R1, R2);
  /* No intersection (coplanar) */
  evaluateRRCollision( scene );

  /* No intersection (not coplanar) */
  V1 << 0, 2, 1;  V2 << 0, 3, 1; V3 << 2, 3, -1; V4 << 2, 2, -1;
  scene.setVerticesR2(V1, V2, V3, V4);
  evaluateRRCollision( scene );

  /* Intersection (coplanar) */
  /*  ( 2, 0, 0), ( 2, 1, 0), ( 1, 1, 0), ( 1, 0, 0)  */
  V1 << 1, -1, 0;  V2 << 1, 2, 0; V3 << 3, 2, 0; V4 << 3, -1, 0;
  scene.setVerticesR2(V1, V2, V3, V4);
  evaluateRRCollision( scene );

  /* Intersection is a line with collinear points: it selects the extremes */
  V1 << 4, 0, 1;  V2 << 0, 0, 1; V3 << 0, 2, 1; V4 << 4, 2, 1;
  scene.setVerticesR1(V1, V2, V3, V4);
  V1 << 1, 1, 1;  V2 << 1, 1, 2; V3 << 1, 3, 2; V4 << 1, 3, 1;
  scene.setVerticesR2(V1, V2, V3, V4);
  evaluateRRCollision( scene );

  /* Intersection */
  V1 << 0, 0, 0;  V2 << 0, 5, 0; V3 << 4, 5, 0; V4 << 4, 0, 0;
  scene.setVerticesR1(V1, V2, V3, V4);
  V1 << -2, 1, 0;  V2 << 1, 4, 0; V3 << 3, 3, 0; V4 << 1, -3, 0;
  scene.setVerticesR2(V1, V2, V3, V4);
  evaluateRRCollision( scene );

  // setVector(R1V1, 0.0988299, -0.173027, 0.100004);  setVector(R1V2, 0.0988298, -0.173027, 0.0999938);
  // setVector(R1V3, 0.0987169, -0.0380271, 0.0999966);  setVector(R1V4, 0.0987169, -0.0380271, 0.0999866);
  // // setVector(R1V1, 0.0988299, -0.173027, 0.100004);  setVector(R1V2, 0.0988298, -0.173027, 0.0999938);
  // // setVector(R1V3, 0.0987169, -0.0380271, 0.0999966);  setVector(R1V4, 0.0987169, -0.0380271, 0.0999866);
  // setVector(R2V1, 0.125, -0.1875, 0.1); setVector(R2V2, -0.065, -0.1875, 0.1);
  // setVector(R2V3, -0.065, -0.0725, 0.1);  setVector(R2V4, 0.125, -0.0725, 0.1);
  // evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  std::cout << "\nTesting case: " << std::endl;

  V1 << 0.0941395, -0.00388724, 0.08; V2 << -0.135861, -0.00388724, 0.08; V3 << -0.135861, -0.138887, 0.08; V4 << 0.0941395, -0.138887, 0.08;
  scene.setVerticesR1(V1, V2, V3, V4);

  V1 << 0.125, -0.0725, 0.08; V2 << -0.065, -0.0725, 0.08; V3 << -0.065, -0.1875, 0.08; V4 << 0.125, -0.1875, 0.08;
  scene.setVerticesR2(V1, V2, V3, V4);

  evaluateRRCollision( scene );

  return 0;
}


void evaluateRRCollision( CollisionTwoRectangles &scene1)
{
  int collisionDetected;
  collisionDetected = scene1.computeRRintersections();
  scene1.printRectanglesVertices();
  scene1.printRRcollisionInformation();
  std::cout << std::endl;
}








  // if (collisionDetected){
  //   std::cout << "Rectangles are intersecting" << std::endl;
  //   std::cout << "Intersections occur at: " << std::endl;
  //   for (int i=0; i<scene.pointsRR.size(); i++){
  //     std::cout << "( " << scene.pointsRR[i].x << ", " << scene.pointsRR[i].y
  //     		<< ", " << scene.pointsRR[i].z << ")" << std::endl;
  //   }
  // }
  // else{
  //   std::cout << "Rectangles are not intersecting" << std::endl;
  // }
