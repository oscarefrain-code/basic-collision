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



class Box
{
 public:
  /* -- Vertices with respect to the world -- */
  Vector3d v1; Vector3d v2; Vector3d v3; Vector3d v4;
  Vector3d v5; Vector3d v6; Vector3d v7; Vector3d v8;

  /* -- Constructors -- */
  Box()
    {
      Lx=0.0; Ly=0.0; Lz=0.0;
      setLength(Lx, Ly, Lz);
      T.setZero(); R.setIdentity(); setTransformation(R,T);
    }

  Box(double lx, double ly, double lz)
    {  
      Lx=lx, Ly=ly, Lz=lz;
      setLength(lx, ly, lz);
      T.setZero(); R.setIdentity(); setTransformation(R,T);
    }

  /* -- Functions to modify the properties of the box -- */
  void setLength(double lx, double ly, double lz)
    {
      Lx=lx; Ly=ly; Lz=lz;
      vo1 <<  Lx/2, -Ly/2,  Lz/2;
      vo2 << -Lx/2, -Ly/2,  Lz/2;
      vo3 << -Lx/2, -Ly/2, -Lz/2;
      vo4 <<  Lx/2, -Ly/2, -Lz/2;
      vo5 <<  Lx/2,  Ly/2,  Lz/2;
      vo6 << -Lx/2,  Ly/2,  Lz/2;
      vo7 << -Lx/2,  Ly/2, -Lz/2;
      vo8 <<  Lx/2,  Ly/2, -Lz/2;
      setTransformation(R,T);
    }

  void setTransformation(Matrix3d Rin, Vector3d Tin)
  {
    R=Rin; T=Tin;
    applyRT(v1, vo1); applyRT(v2, vo2); applyRT(v3, vo3); applyRT(v4, vo4);
    applyRT(v5, vo5); applyRT(v6, vo6); applyRT(v7, vo7); applyRT(v8, vo8);
  }

  /* -- Getters -- */
  void getLength(double &lx, double &ly, double &lz)
  {
    lx=Lx; ly=Ly; lz=Lz; 
  }

  void getVertices(Vector3d &v1out, Vector3d &v2out, Vector3d &v3out, Vector3d &v4out,
		   Vector3d &v5out, Vector3d &v6out, Vector3d &v7out, Vector3d &v8out)
  { 
    v1out=v1; v2out=v2; v3out=v3; v4out=v4; v5out=v5; v6out=v6; v7out=v7; v8out=v8;
  }

  /* -- Print on screen -- */
  void printVertices()
  {
    /* Print as: V1=[x y z], V2=[x y z], V3=[x y z], V4=[x y z], ... */
    std::cout << "V1=["<< v1.transpose() << "], V2=[" << v2.transpose() << "], V3=[" << v3.transpose()
	      << "], V4=[" << v4.transpose() << "], V5=[" << v5.transpose() << "], V6=[" << v6.transpose()
      	      << "], V7=[" << v7.transpose() << "], V8=[" << v8.transpose() << "]" << std::endl;
  }
  void printLengths()
  {
    std::cout << "(Lx, Ly, Lz) = (" << Lx << ", " << Ly << ", " << Lz << ")" << std::endl;
  }

 private:
  /* -- Vertices with respect to the origin of the world (symmetric) -- */
  Vector3d vo1; Vector3d vo2; Vector3d vo3; Vector3d vo4;
  Vector3d vo5; Vector3d vo6; Vector3d vo7; Vector3d vo8;
  Vector3d T; Matrix3d R;
  double Lx, Ly, Lz;

  /* --- Transformation of a Point --- */
  void applyRT(Vector3d &pout, Vector3d pin)
  {
    pout(0) = R(0,0)*pin(0) + R(0,1)*pin(1) + R(0,2)*pin(2) + T(0);
    pout(1) = R(1,0)*pin(0) + R(1,1)*pin(1) + R(1,2)*pin(2) + T(1);
    pout(2) = R(2,0)*pin(0) + R(2,1)*pin(1) + R(2,2)*pin(2) + T(2);
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
  void setTransformation1(Matrix3d R1, Vector3d T1);
  void setTransformation2(Matrix3d R2, Vector3d T2);
  /* void setTransformation1(double R1[9], double T1[3]); */
  /* void setTransformation2(double R2[9], double T2[3]); */

  /* -- Compute the collision detection -- */
  int computeBBintersections( void );

  /* -- Get values from the box -- */
  double getTolerance( void );
  void getLengthB1(double & lx1, double & ly1, double & lz1);
  void getLengthB2(double & lx2, double & ly2, double & lz2);
  void getVerticesB1(double V1[3], double V2[3], double V3[3], double V4[3],
		     double V5[3], double V6[3], double V7[3], double V8[3]);
  void getVerticesB2(double V1[3], double V2[3], double V3[3], double V4[3],
		     double V5[3], double V6[3], double V7[3], double V8[3]);
  /* TODO */
  /* void getTransformation1(double R[9], double T[3]); */
  /* void getTransformation2(double R[9], double T[3]); */
  /* void printTransformation1( void ); */
  /* void printTransformation1( void ); */

  /* --- Print Different Information --- */
  void printLengths( void );
  void printBoxesVertices( void );
  void printBBcollisionInformation( void );

  /* -- Collision Points -- */
  std::vector<Vector3d> pointsBB;

 private:
  /* -- Variables -- */
  /* double Lx1, Ly1, Lz1, Lx2, Ly2, Lz2; */
  /* double Pc11[3], Pc12[3], Pc13[3], Pc14[3], Pc15[3], Pc16[3], Pc17[3], Pc18[3]; */
  /* double Pc21[3], Pc22[3], Pc23[3], Pc24[3], Pc25[3], Pc26[3], Pc27[3], Pc28[3]; */
  /* double Pw11[3], Pw12[3], Pw13[3], Pw14[3], Pw15[3], Pw16[3], Pw17[3], Pw18[3]; */
  /* double Pw21[3], Pw22[3], Pw23[3], Pw24[3], Pw25[3], Pw26[3], Pw27[3], Pw28[3]; */
  Box B1, B2;
  int collisionIndicator;
  double tolerance;

  /* -- Collision rect/rect objects -- */
  CollisionTwoRectangles To1_To2, To1_Bo2, To1_Le2, To1_Ri2, To1_Fr2, To1_Ba2,
    Bo1_To2, Bo1_Bo2, Bo1_Le2, Bo1_Ri2, Bo1_Fr2, Bo1_Ba2,
    Le1_To2, Le1_Bo2, Le1_Le2, Le1_Ri2, Le1_Fr2, Le1_Ba2,
    Ri1_To2, Ri1_Bo2, Ri1_Le2, Ri1_Ri2, Ri1_Fr2, Ri1_Ba2,
    Fr1_To2, Fr1_Bo2, Fr1_Le2, Fr1_Ri2, Fr1_Fr2, Fr1_Ba2,
    Ba1_To2, Ba1_Bo2, Ba1_Le2, Ba1_Ri2, Ba1_Fr2, Ba1_Ba2;

  /* -- Apply the rotation and translation -- */
  //void applyRT(double pout[3], double pin[3], double R[9], double T[3]);
  void averageSimilarValues( void );

  /* -- To initialize some values to zero -- */
  void defaultInit( void );

};

#endif
