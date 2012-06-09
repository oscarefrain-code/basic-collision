/***********************************************************************************
 * Triangle/triangle intersection test
 * Based on the paper:
 *     Thomas Moller, "A Fast Triangle-Triangle Intersection Test",
 *     Journal of Graphics Tools, 2(2), 1997
 *
 *                                           Oscar Efrain Ramos Ponce, LAAS-CNRS
 ************************************************************************************/

#ifndef __COLLISION_TWO_TRIANGLES_H__
#define __COLLISION_TWO_TRIANGLES_H__

#include <math.h>
#include <vector>
#include <iostream>
#include <collision-helper.h>


class CollisionTwoTriangles
: public OperationsHelper
{
 public:

  /* -- Constructors -- */
  CollisionTwoTriangles();
  CollisionTwoTriangles(double T1_V0[3], double T1_V1[3], double T1_V2[3],
			double T2_V0[3], double T2_V1[3], double T2_V2[3]);

  /* -- Set values for the vertices -- */
  void setVerticesAll(double T1_V0[3], double T1_V1[3], double T1_V2[3],
		      double T2_V0[3], double T2_V1[3], double T2_V2[3]);
  void setVerticesT1(double T1_V0[3], double T1_V1[3], double T1_V2[3]);
  void setVerticesT2(double T2_V0[3], double T2_V1[3], double T2_V2[3]);

  /* -- Get the values of the vertices -- */
  void getVerticesT1(double T1_V0[3], double T1_V1[3], double T1_V2[3]);
  void getVerticesT2(double T2_V0[3], double T2_V1[3], double T2_V2[3]);

  /* -- Compute the collision detection algorithm -- */
  int computeTTintersectionsNoLine( void );
  int computeTTintersections( void );

  /* -- Print additional information -- */
  void printTTcollisionInformation( void );
  void printTrianglesVertices( void );

  /* -- Vector to store the collision points -- */
  std::vector<Point3d> pointsTT;

 private:

  /* --- Variables --- */
  double V0[3], V1[3], V2[3];   // Vertices of triangle 1: V0, V1, V2
  double U0[3], U1[3], U2[3];   // Vertices of triangle 2: U0, U1, U2
  double D[3];                  // Direction of the (intersecting) line
  int i0, i1, inot;             // Indexes
  int collisionIndicator;
  
  /* --- Internal Functions --- */
  void defaultInit( void );
  void computeIntervals(double VV[3], double DD[3], double D0D1, double D0D2, 
			double isect[2], bool &coplanar);
  void find_t(double VV0, double VV1, double VV2, double D0, double D1, double D2,
	      double isect[2]);
  bool coplanar_tri_tri(double N[3]);
  bool edge_against_tri_edges(double v0[3], double v1[3]);
  bool edge_edge_test(double v0[3], double u0[3], double u1[3], double Ax, double Ay);
  bool point_in_tri(double v0[3], double u0[3], double u1[3], double u2[3]);
  
  int compute_intervals_isectline(int id, double VV[3], double DD[3],
				  double D0D1,double D0D2,double *isect0, double *isect1,
				  double isectpoint0[3],double isectpoint1[3]);
  void isect2(double VTX0[3], double VTX1[3], double VTX2[3], double VV0, double VV1, double VV2,
	      double D0,double D1,double D2,double *isect0,double *isect1,double isectpoint0[3],double isectpoint1[3]);

  
};

#endif
