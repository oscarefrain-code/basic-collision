/***********************************************************************************
 * Rectangle/rectangle intersection test:
 *   Each rectangle is decomposed into two triangles and the collision between 
 *   triangles is computed.   
 *
 *  
 *                                           Oscar E. Ramos Ponce, LAAS-CNRS, France
 ************************************************************************************/


#ifndef __COLLISION_TWO_RECTANGLES_H__
#define __COLLISION_TWO_RECTANGLES_H__

#include <math.h>
#include <vector>
#include <collision-helper.h>
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


