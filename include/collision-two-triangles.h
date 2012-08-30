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
 * Triangle/triangle intersection test
 * Based on the paper:
 *     Thomas Moller, "A Fast Triangle-Triangle Intersection Test",
 *     Journal of Graphics Tools, 2(2), 1997
 *
 */

#ifndef __COLLISION_TWO_TRIANGLES_H__
#define __COLLISION_TWO_TRIANGLES_H__

#include <math.h>
#include <vector>
#include <iostream>
#include <operations-helper.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>


/* Definition of a 'minimum value' to be used as a threshold for the triangle test*/
//#define TRI_EPSILON 0.000001

/* ------------------------------------------------------------------ */
class Triangle
{
 public:
  /* -- Vertices -- */
  Eigen::Vector3d v1;   
  Eigen::Vector3d v2;   
  Eigen::Vector3d v3;
  /* -- Constructors -- */
  Triangle()
    {
      v1 << 0, 0, 0;
      v2 << 0, 0, 0;
      v3 << 0, 0, 0;
    }
  Triangle(const Eigen::Vector3d & v1in, const Eigen::Vector3d & v2in, const Eigen::Vector3d & v3in)
    { 
      setVertices(v1in, v2in, v3in);
    }
  /* -- Functions -- */
  void setVertices(const Eigen::Vector3d & v1in, const Eigen::Vector3d & v2in, const Eigen::Vector3d & v3in)
  {  v1=v1in; v2=v2in; v3=v3in;
  }
  void getVertices(Eigen::Vector3d &v1out, Eigen::Vector3d &v2out, Eigen::Vector3d &v3out)
  {  v1out=v1; v2out=v2; v3out=v3;
  }
  void print()
  {
    /* Print as: V1=[x y z], V2=[x y z], V3=[x y z] */
    /* std::cout << "V1=["<< v1.transpose() << "]; V2=[" << v2.transpose() */
    /* 	      << "]; V3=[" << v3.transpose() << "]; " << std::endl; */
    std::cout << "V1=[" << v1(0) << ", " << v1(1) << ", " << v1(2) << "]; "
	      << "V2=[" << v2(0) << ", " << v2(1) << ", " << v2(2) << "]; "
	      << "V3=[" << v3(0) << ", " << v3(1) << ", " << v3(2) << "]; " << std::endl;
  }
};


/* ------------------------------------------------------------------ */
class CollisionTwoTriangles
: public OperationsHelper
{
 public:

  /* -- Constructors -- */
  CollisionTwoTriangles();
  CollisionTwoTriangles(Triangle t1, Triangle t2);
  CollisionTwoTriangles(Eigen::Vector3d t1v1, Eigen::Vector3d t1v2, Eigen::Vector3d t1v3,
  			Eigen::Vector3d t2v1, Eigen::Vector3d t2v2, Eigen::Vector3d t2v3);

  /* -- Set the triangles -- */
  void setTriangles(Triangle t1, Triangle t2);
  void setTriangle1(Triangle t1);
  void setTriangle2(Triangle t2);

  /* -- Set values for the vertices -- */
  void setVerticesAll(Eigen::Vector3d t1v1, Eigen::Vector3d t1v2, Eigen::Vector3d t1v3,
		      Eigen::Vector3d t2v1, Eigen::Vector3d t2v2, Eigen::Vector3d t2v3);
  void setVerticesT1(Eigen::Vector3d V1in, Eigen::Vector3d V2in, Eigen::Vector3d V3in);
  void setVerticesT2(Eigen::Vector3d V1in, Eigen::Vector3d V2in, Eigen::Vector3d V3in);

  /* -- Set and get the coplanar tolerance -- */
  void setCoplanarTolerance(double tol);
  double getCoplanarTolerance(void);

  /* -- Get the values of the vertices -- */
  void getVerticesT1(Eigen::Vector3d &V1out, Eigen::Vector3d &V2out, Eigen::Vector3d &V3out);
  void getVerticesT2(Eigen::Vector3d &V1out, Eigen::Vector3d &V2out, Eigen::Vector3d &V3out);

  /* -- Compute the collision detection algorithm -- */
  int computeTTintersections( void );

  /* -- Print additional information -- */
  void printTTcollisionInformation( void );

  void printTrianglesVertices( void );

  /* -- Vector to store the collision points -- */
  std::vector<Eigen::Vector3d> pointsTT;

 private:

  /* --- Variables --- */
  Triangle T1;
  Triangle T2;

  Eigen::Vector3d D;           // Direction of the intersecting line 
  int i0, i1, inot;             // Indexes
  int collisionIndicator;

  /* Value that indicates the admitted tolerance for coplanarity
     If the distance from a vertex of Ti to the plane of Tj is smaller than this value,
     it is assumed that the vertex lies in the plane  */
  double coplanar_tolerance;   

  /* --- Internal Functions --- */
  void defaultInit( void );
  bool coplanar_tri_tri(Eigen::Vector3d N);
  bool edge_against_tri_edges(Eigen::Vector3d v0, Eigen::Vector3d v1);
  bool edge_edge_test(Eigen::Vector3d v0, Eigen::Vector3d u0, Eigen::Vector3d u1, double Ax, double Ay);
  bool point_in_tri(Eigen::Vector3d v0, Eigen::Vector3d u0, Eigen::Vector3d u1, Eigen::Vector3d u2);
  int compute_intervals_isectline(int id, Eigen::Vector3d VV, Eigen::Vector3d DD,
				  double D0D1, double D0D2, Eigen::Vector2d &isect,
				  Eigen::Vector3d &isectpoint0, Eigen::Vector3d &isectpoint1);
  void isect2(Eigen::Vector3d VTX0, Eigen::Vector3d VTX1, Eigen::Vector3d VTX2, 
	      double VV1, double VV2, double VV3, double D0,double D1,double D2,
	      Eigen::Vector2d &isect, Eigen::Vector3d &isectpoint0, Eigen::Vector3d &isectpoint1); 

  double distPoint2Segment( const Eigen::Vector3d & P, const Eigen::Vector3d & P1, const Eigen::Vector3d & P2);
  double euclideanDistance(const Eigen::Vector3d & P1, const Eigen::Vector3d & P2);

  
};

#endif
