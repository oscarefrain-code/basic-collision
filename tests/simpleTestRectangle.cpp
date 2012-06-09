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


void setVector(double v[3], double v0, double v1, double v2);
void evaluate(CollisionTwoRectangles & scene1,
	      double R1V1[3], double R1V2[3], double R1V3[3], double R1V4[3],
	      double R2V1[3], double R2V2[3], double R2V3[3], double R2V4[3]);


int main( void )
{
  double R1V1[3], R1V2[3], R1V3[3], R1V4[3], R2V1[3], R2V2[3], R2V3[3], R2V4[3];
  CollisionTwoRectangles scene;

  // setVector(R1V1, 0.0, 0.0, 0.0);  setVector(R1V2, 0.0, 1.0, 0.0);
  // setVector(R1V3, 2.0, 1.0, 0.0);  setVector(R1V4, 2.0, 0.0, 0.0);

  // /* No intersection (coplanar) */
  // setVector(R2V1, 0.0, 2.0, 0.0);  setVector(R2V2, 0.0, 3.0, 0.0);
  // setVector(R2V3, 2.0, 3.0, 0.0);  setVector(R2V4, 2.0, 2.0, 0.0);
  // evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  // /* No intersection (not coplanar) */
  // setVector(R2V1, 0.0, 2.0, 1.0);  setVector(R2V2, 0.0, 3.0, 1.0);
  // setVector(R2V3, 2.0, 3.0, -1.0); setVector(R2V4, 2.0, 2.0, -1.0);
  // evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  // /* Intersection (coplanar) */
  // /*  ( 2, 0, 0), ( 2, 1, 0), ( 1, 1, 0), ( 1, 0, 0)  */
  // setVector(R2V1, 1.0, -1.0, 0.0); setVector(R2V2, 1.0, 2.0, 0.0);
  // setVector(R2V3, 3.0, 2.0, 0.0);  setVector(R2V4, 3.0, -1.0, 0.0);
  // evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  // /* Intersection is a line with collinear points: it selects the extremes */
  // setVector(R1V1, 4.0, 0.0, 1.0);  setVector(R1V2, 0.0, 0.0, 1.0);
  // setVector(R1V3, 0.0, 2.0, 1.0);  setVector(R1V4, 4.0, 2.0, 1.0);
  // setVector(R2V1, 1.0, 1.0, 1.0);  setVector(R2V2, 1.0, 1.0, 2.0);
  // setVector(R2V3, 1.0, 3.0, 2.0);  setVector(R2V4, 1.0, 3.0, 1.0);
  // evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  // /* Intersection */
  // setVector(R1V1, 0.0, 0.0, 0.0);  setVector(R1V2, 0.0, 5.0, 0.0);
  // setVector(R1V3, 4.0, 5.0, 0.0);  setVector(R1V4, 4.0, 0.0, 0.0);
  // setVector(R2V1, -2.0, 1.0, 0.0); setVector(R2V2, 1.0, 4.0, 0.0);
  // setVector(R2V3, 3.0, 3.0, 0.0);  setVector(R2V4, 1.0, -3.0, 0.0);
  // evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  setVector(R1V1, 0.0988299, -0.173027, 0.100004);  setVector(R1V2, 0.0988298, -0.173027, 0.0999938);
  setVector(R1V3, 0.0987169, -0.0380271, 0.0999966);  setVector(R1V4, 0.0987169, -0.0380271, 0.0999866);
  // setVector(R1V1, 0.0988299, -0.173027, 0.100004);  setVector(R1V2, 0.0988298, -0.173027, 0.0999938);
  // setVector(R1V3, 0.0987169, -0.0380271, 0.0999966);  setVector(R1V4, 0.0987169, -0.0380271, 0.0999866);
  setVector(R2V1, 0.125, -0.1875, 0.1); setVector(R2V2, -0.065, -0.1875, 0.1);
  setVector(R2V3, -0.065, -0.0725, 0.1);  setVector(R2V4, 0.125, -0.0725, 0.1);
  evaluate(scene, R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);

  return 0;
}


void setVector(double v[3], double v0, double v1, double v2)
{
  v[0]=v0; v[1]=v1; v[2]=v2;
}


void evaluate(CollisionTwoRectangles & scene1,
	      double R1V1[3], double R1V2[3], double R1V3[3], double R1V4[3],
	      double R2V1[3], double R2V2[3], double R2V3[3], double R2V4[3])
{
  int collisionDetected;
  scene1.setVerticesAll(R1V1, R1V2, R1V3, R1V4, R2V1, R2V2, R2V3, R2V4);
  collisionDetected = scene1.computeRRintersections();
  std::cout << std::endl;
  scene1.printRectanglesVertices();
  scene1.printRRcollisionInformation();
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
