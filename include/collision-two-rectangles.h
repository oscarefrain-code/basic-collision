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
 */


#ifndef __COLLISION_TWO_RECTANGLES_H__
#define __COLLISION_TWO_RECTANGLES_H__

#include <math.h>
#include <vector>
#include <operations-helper.h>
#include <collision-two-triangles.h>


class CollisionTwoRectangles
: public OperationsHelper
{
 public:
   /* -- Constructors -- */
  CollisionTwoRectangles( void );
  CollisionTwoRectangles(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3],
			 double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3]);

  /* -- Set values for the vertices (in succesive order, clockwise or anticlockwise) -- */
  void setVerticesAll(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3],
		      double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3]);
  void setVerticesR1(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3]);
  void setVerticesR2(double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3]);

  /* -- Get the values of the vertices -- */
  void getVerticesR1(double R1_V0[3], double R1_V1[3], double R1_V2[3], double R1_V3[3]);
  void getVerticesR2(double R2_V0[3], double R2_V1[3], double R2_V2[3], double R2_V3[3]);


  /* -- Compute the collision detection -- */
  int computeRRintersections( void );

  /* -- Print additional information -- */
  void printRRcollisionInformation( void );
  void printRectanglesVertices( void );

  /* -- Collision Points -- */
  std::vector<Point3d> pointsRR;
 
 private:
  /* --- Variables --- */
  double V0[3], V1[3], V2[3], V3[3];   // Vertices of rectangle 1: V0, V1, V2, V3
  double U0[3], U1[3], U2[3], U3[3];   // Vertices of rectangle 2: U0, U1, U2, U3
  int collisionIndicator;

};


#endif


