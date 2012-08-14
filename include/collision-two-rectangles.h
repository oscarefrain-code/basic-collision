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
#include <Eigen/Dense>

/* ------------------------------------------------------------------ */
class Rectangle
{
 public:
  /* -- Constructors -- */
  Rectangle()
    {
      v1 << 0, 0, 0; v2 << 0, 0, 0;
      v3 << 0, 0, 0; v4 << 0, 0, 0;
    }
  Rectangle(const Eigen::Vector3d v1in, const Eigen::Vector3d v2in, const Eigen::Vector3d v3in, const Eigen::Vector3d v4in)
    { 
      setVertices(v1in, v2in, v3in, v4in);
    }
  /* -- Functions -- */
  void setVertices(const Eigen::Vector3d v1in, const Eigen::Vector3d v2in, const Eigen::Vector3d v3in, const Eigen::Vector3d v4in)
  {  v1=v1in; v2=v2in; v3=v3in; v4=v4in;
  }
  void getVertices(Eigen::Vector3d &v1out, Eigen::Vector3d &v2out, Eigen::Vector3d &v3out, Eigen::Vector3d &v4out)
  {  v1out=v1; v2out=v2; v3out=v3; v4out=v4;
  }
  void print()
  {
    /* Print as: V1=[x y z], V2=[x y z], V3=[x y z], V4=[x y z] */
    std::cout << "V1=["<< v1.transpose() << "], V2=[" << v2.transpose() << "], V3=["
	      << v3.transpose() << "], V4=[" << v4.transpose() << "]" << std::endl;
  }
 private:
  /* -- Vertices -- */
  Eigen::Vector3d v1, v2, v3, v4;   
};


class CollisionTwoRectangles
: public OperationsHelper
{
 public:
   /* -- Constructors -- */
  CollisionTwoRectangles( void );
  CollisionTwoRectangles(Eigen::Vector3d R1_V0, Eigen::Vector3d R1_V1, Eigen::Vector3d R1_V2, Eigen::Vector3d R1_V3,
			 Eigen::Vector3d R2_V0, Eigen::Vector3d R2_V1, Eigen::Vector3d R2_V2, Eigen::Vector3d R2_V3);
  CollisionTwoRectangles(Rectangle r1in, Rectangle r2in);

  /* -- Set values for the vertices (in succesive order, clockwise or anticlockwise) -- */
  void setVerticesAll(Eigen::Vector3d r1v1, Eigen::Vector3d r1v2, Eigen::Vector3d r1v3, Eigen::Vector3d r1v4,
		      Eigen::Vector3d r2v1, Eigen::Vector3d r2v2, Eigen::Vector3d r2v3, Eigen::Vector3d r2v4);
  void setVerticesR1(Eigen::Vector3d V1in, Eigen::Vector3d V2in, Eigen::Vector3d V3in, Eigen::Vector3d V4in);
  void setVerticesR2(Eigen::Vector3d V1in, Eigen::Vector3d V2in, Eigen::Vector3d V3in, Eigen::Vector3d V4in);
  
  /* -- Set/get the coplanar tolerance (of the triangles) -- */
  void setCoplanarTolerance(const double tol);
  double getCoplanarTolerance( void );

  /* -- Set the rectangles -- */
  void setRectangles(Rectangle r1in, Rectangle r2in);
  void setRectangle1(Rectangle r1in);
  void setRectangle2(Rectangle r2in);

  /* -- Get the values of the vertices -- */
  void getVerticesR1(Eigen::Vector3d &V1out, Eigen::Vector3d &V2out, Eigen::Vector3d &V3out, Eigen::Vector3d &V4out);
  void getVerticesR2(Eigen::Vector3d &V1out, Eigen::Vector3d &V2out, Eigen::Vector3d &V3out, Eigen::Vector3d &V4out);

  /* -- Compute the collision detection -- */
  int computeRRintersections( void );

  /* -- Print additional information -- */
  void printRRcollisionInformation( void );
  void printRectanglesVertices( void );

  /* -- Collision Points -- */
  std::vector<Eigen::Vector3d> pointsRR;
  //std::vector<Point3d> pointsRR;
 
 private:
  /* --- Variables --- */
  Rectangle R1, R2;                        // Rectangles to be used
  int collisionIndicator;

  /* -- Objects that detect collision between triangles -- */
  CollisionTwoTriangles scene1;
  CollisionTwoTriangles scene2;
  CollisionTwoTriangles scene3;
  CollisionTwoTriangles scene4;


};


#endif


