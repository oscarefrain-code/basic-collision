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
 * Hexahedron/hexahedron intersection test:
 *   Each hexahedron is decomposed into a set of rectangles and the collision
 *   between the rectangles is tested.
 *  
 */

#ifndef __COLLISION_TWO_BOXES_H__
#define __COLLISION_TWO_BOXES_H__

#include <math.h>
#include <vector>
#include <operations-helper.h>
#include <collision-two-rectangles.h>
#include <iostream>


class Box
{
 public:
  /* -- Constructors -- */
  Box()
    {
      Lx=0.0; Ly=0.0; Lz=0.0;
      setLength(Lx, Ly, Lz);
      T.setZero(); R.setIdentity(); setTransformation(R,T);
    }

  Box(const double lx, const double ly, const double lz)
    {  
      Lx=lx, Ly=ly, Lz=lz;
      setLength(lx, ly, lz);
      T.setZero(); R.setIdentity(); setTransformation(R,T);
    }

  /* -- Functions to modify the properties of the box -- */
  void setLength(const double lx, const double ly, const double lz)
    {
      Lx=lx; Ly=ly; Lz=lz;
      vo1 <<  Lx/2,  Ly/2,  Lz/2;
      vo2 << -Lx/2,  Ly/2,  Lz/2;
      vo3 << -Lx/2,  Ly/2, -Lz/2;
      vo4 <<  Lx/2,  Ly/2, -Lz/2;
      vo5 <<  Lx/2, -Ly/2,  Lz/2;
      vo6 << -Lx/2, -Ly/2,  Lz/2;
      vo7 << -Lx/2, -Ly/2, -Lz/2;
      vo8 <<  Lx/2, -Ly/2, -Lz/2;
      setTransformation(R,T);
    }

  void setTransformation(const Eigen::Matrix3d & Rin, const Eigen::Vector3d & Tin)
  {
    R=Rin; T=Tin;
    v1=applyRT(vo1); v2=applyRT(vo2); v3=applyRT(vo3); v4=applyRT(vo4);
    v5=applyRT(vo5); v6=applyRT(vo6); v7=applyRT(vo7); v8=applyRT(vo8);
  }

  /* -- Getters -- */
  void getLength(double &lx, double &ly, double &lz)
  {
    lx=Lx; ly=Ly; lz=Lz; 
  }

  void getVertices(Eigen::Vector3d &v1out, Eigen::Vector3d &v2out, Eigen::Vector3d &v3out, Eigen::Vector3d &v4out,
		   Eigen::Vector3d &v5out, Eigen::Vector3d &v6out, Eigen::Vector3d &v7out, Eigen::Vector3d &v8out)
  { 
    v1out=v1; v2out=v2; v3out=v3; v4out=v4; v5out=v5; v6out=v6; v7out=v7; v8out=v8;
  }

  void getTransformation(Eigen::Matrix3d &Rout, Eigen::Vector3d &Tout)
  {
    Rout=R; Tout=T;
  }

  /* -- Print on screen -- */
  void printVertices()
  {
    /* std::cout << "V1=["<< v1.transpose() << "], V2=[" << v2.transpose() << "], V3=[" << v3.transpose() */
    /* 	      << "], V4=[" << v4.transpose() << "], V5=[" << v5.transpose() << "], V6=[" << v6.transpose() */
    /*   	      << "], V7=[" << v7.transpose() << "], V8=[" << v8.transpose() << "]" << std::endl; */
    std::cout << "V1="; printV(v1);
    std::cout << "; V2="; printV(v2);
    std::cout << "; V3="; printV(v3);
    std::cout << "; V4="; printV(v4);
    std::cout << "; V5="; printV(v5);
    std::cout << "; V6="; printV(v6);
    std::cout << "; V7="; printV(v7);
    std::cout << "; V8="; printV(v8);
    std::cout << std::endl;
  }

  void printLengths()
  {
    std::cout << "(Lx, Ly, Lz) = (" << Lx << ", " << Ly << ", " << Lz << ")" << std::endl;
  }

  void printTransformation()
  {
    std::cout << R.block<1,3>(0,0) << " " << T(0) << std::endl;
    std::cout << R.block<1,3>(1,0) << " " << T(1) << std::endl;
    std::cout << R.block<1,3>(2,0) << " " << T(2) << std::endl;
  }

 private:
  /* -- Vertices with respect to the world -- */
  Eigen::Vector3d v1, v2, v3, v4, v5, v6, v7, v8;

  /* -- Vertices with respect to the origin of the world (symmetric) -- */
  Eigen::Vector3d vo1, vo2, vo3, vo4, vo5, vo6, vo7, vo8;
  Eigen::Vector3d T; Eigen::Matrix3d R;
  double Lx, Ly, Lz;

  /* --- Transformation of a Point --- */
  Eigen::Vector3d applyRT(const Eigen::Vector3d & pin)
  {
    Eigen::Vector3d pout;
    pout(0) = R(0,0)*pin(0) + R(0,1)*pin(1) + R(0,2)*pin(2) + T(0);
    pout(1) = R(1,0)*pin(0) + R(1,1)*pin(1) + R(1,2)*pin(2) + T(1);
    pout(2) = R(2,0)*pin(0) + R(2,1)*pin(1) + R(2,2)*pin(2) + T(2);
    return pout;
  }

  void printV(const Eigen::Vector3d & v)
  {
    std::cout << "[ " << v(0) << ", " << v(1) << ", " << v(2) << " ]";
  }

};




class CollisionTwoBoxes
: public OperationsHelper
{
 public:
  /* -- Constructors -- */
  CollisionTwoBoxes( void );
  CollisionTwoBoxes(double lx1, double ly1, double lz1,
		    double lx2, double ly2, double lz2);

  /* -- Set values for the box -- */
  void setTolerance(double toler);
  void setLengthB1(double lx1, double ly1, double lz1);
  void setLengthB2(double lx2, double ly2, double lz2);
  void setTransformation1(const Eigen::Matrix3d & R1, const Eigen::Vector3d & T1);
  void setTransformation2(const Eigen::Matrix3d & R2, const Eigen::Vector3d & T2);

  /* -- Set/get the coplanar tolerance (of the low level triangles) -- 
        This tolerance is to decide when the triangles are coplanar    */
  void setCoplanarTolerance(const double coplanarToler);
  double getCoplanarTolerance( void );

  /* -- Compute the collision detection -- */
  int computeBBintersections( void );

  /* -- Get values from the box -- */
  double getTolerance( void );
  void getLengthB1(double & lx1, double & ly1, double & lz1);
  void getLengthB2(double & lx2, double & ly2, double & lz2);
  void getTransformation1(Eigen::Matrix3d &Rout, Eigen::Vector3d &Tout);
  void getTransformation2(Eigen::Matrix3d &Rout, Eigen::Vector3d &Tout);
  void getVerticesB1(Eigen::Vector3d &V1, Eigen::Vector3d &V2, Eigen::Vector3d &V3, Eigen::Vector3d &V4,
		     Eigen::Vector3d &V5, Eigen::Vector3d &V6, Eigen::Vector3d &V7, Eigen::Vector3d &V8);
  void getVerticesB2(Eigen::Vector3d &V1, Eigen::Vector3d &V2, Eigen::Vector3d &V3, Eigen::Vector3d &V4,
		     Eigen::Vector3d &V5, Eigen::Vector3d &V6, Eigen::Vector3d &V7, Eigen::Vector3d &V8);

  /* --- Print Different Information --- */
  void printLengths( void );
  void printBoxesVertices( void );
  void printBox1Vertices( void );
  void printBox2Vertices( void );
  void printBBcollisionInformation( void );
  void printTransformation( void );

  /* -- Collision Points -- */
  std::vector<Eigen::Vector3d> pointsBB;

 private:
  /* -- Variables -- */
  Box B1, B2;
  int collisionIndicator;
  double tolerance;   // Tolerance to average the close points

  /* -- Collision rect/rect objects -- */
  CollisionTwoRectangles To1_To2, To1_Bo2, To1_Le2, To1_Ri2, To1_Fr2, To1_Ba2,
    Bo1_To2, Bo1_Bo2, Bo1_Le2, Bo1_Ri2, Bo1_Fr2, Bo1_Ba2,
    Le1_To2, Le1_Bo2, Le1_Le2, Le1_Ri2, Le1_Fr2, Le1_Ba2,
    Ri1_To2, Ri1_Bo2, Ri1_Le2, Ri1_Ri2, Ri1_Fr2, Ri1_Ba2,
    Fr1_To2, Fr1_Bo2, Fr1_Le2, Fr1_Ri2, Fr1_Fr2, Fr1_Ba2,
    Ba1_To2, Ba1_Bo2, Ba1_Le2, Ba1_Ri2, Ba1_Fr2, Ba1_Ba2;

  /* -- Apply the rotation and translation -- */
  void averageSimilarValues( void );

  /* -- To initialize some values to zero -- */
  void defaultInit( void );

  /* -- Some functions used in debugging process -- */
  void debugPrintCollisionRectangles(int collisions[36]);
  
};

#endif
