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


void setVector(double v[3], double v0, double v1, double v2);
void evaluate(CollisionTwoTriangles & scene1,
	      double T1V1[3], double T1V2[3], double T1V3[3],
	      double T2V1[3], double T2V2[3], double T2V3[3]);


int main( void )
{
  double T1V1[3], T1V2[3], T1V3[3], T2V1[3], T2V2[3], T2V3[3];
  CollisionTwoTriangles scene;

  setVector(T1V1, 0.0, 0.0, 0.0);
  setVector(T1V2, 1.0, 0.0, 0.0);
  setVector(T1V3, 1.0, 1.0, 0.0);

  /* No intersection */
  setVector(T2V1, 2.0, 2.0, 0.0);
  setVector(T2V2, 3.0, 2.0, 1.0);
  setVector(T2V3, 2.0, 3.0, 2.0);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  /* No intersection (coplanar) */
  setVector(T2V1, 0.0, 2.0, 0.0);
  setVector(T2V2, 1.0, 2.0, 0.0);
  setVector(T2V3, 1.0, 3.0, 0.0);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  /* Intersection: triangle is perpendicular */
  /*     (0.5,0.5,0) and (1,0.5,0)           */ 
  setVector(T2V1, 0.0, 0.5, -0.5);
  setVector(T2V2, 1.0, 0.5, -0.5);
  setVector(T2V3, 1.0, 0.5, 0.5);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  /* Intersection: One point (perpendicular) */
  /*   (0.5, 0.25, 0.0)                      */
  setVector(T2V1, 0.5, 0.25, 0.0);
  setVector(T2V2, 1.0, 0.0, 1.0);
  setVector(T2V3, 0.0, 1.0, 1.0);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  /* Intersection: coplanar (one part intersecting) */
  /* ( 1, 0, 0), ( 1, 1, 0) and ( 0.5, 0.5, 0)      */
  setVector(T2V1, 0.0, 1.0, 0.0);
  setVector(T2V2, 1.0, 0.0, 0.0);
  setVector(T2V3, 1.0, 1.0, 0.0);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  /* Intersection: coplanar (one absolutely contained in the other) */
  /* ( 0, 0, 0), ( 1, 0, 0), ( 1, 1, 0)                             */
  setVector(T2V1, -2.0, -1.0, 0.0);
  setVector(T2V2, 1.0, 2.0, 0.0);
  setVector(T2V3, 3.0, -1.0, 0.0);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);


  /* Miscellaneous tests */
  setVector(T1V1, 0.0, 3.0, 0.0);
  setVector(T1V2, 3.0, 0.0, 0.0);
  setVector(T1V3, 3.0, 3.0, 0.0);
  setVector(T2V1, -1.0, 2.0, 0.0);
  setVector(T2V2, 2.0, -1.0, 0.0);
  setVector(T2V3, 2.0, 2.0, 0.0);
  evaluate(scene, T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);

  return 0;
}


void setVector(double v[3], double v0, double v1, double v2)
{
  v[0]=v0; v[1]=v1; v[2]=v2;
}


void evaluate(CollisionTwoTriangles & scene1,
	      double T1V1[3], double T1V2[3], double T1V3[3],
	      double T2V1[3], double T2V2[3], double T2V3[3])
{
  int collisionDetected;
  scene1.setVerticesAll(T1V1, T1V2, T1V3, T2V1, T2V2, T2V3);
  collisionDetected = scene1.computeTTintersections();
  std::cout << std::endl;
  scene1.printTrianglesVertices();
  scene1.printTTcollisionInformation();
}


// if (collisionDetected){
//    std::cout << "Triangles are intersecting" << std::endl;
//    std::cout << " - Intersections occur at: " << std::endl;
//    for (int i=0; i<scene.interPoints.size(); i++){
//      std::cout << "     "; scene.interPoints[i].print();
//      std::cout << std::endl;
//    }
//  }
//  else{
//    std::cout << "Triangles are NOT intersecting" << std::endl;
//  }

