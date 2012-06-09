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
  void setTransformation1(double R1[9], double T1[3]);
  void setTransformation2(double R2[9], double T2[3]);

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
  std::vector<Point3d> pointsBB;

 private:
  /* -- Variables -- */
  double Lx1, Ly1, Lz1, Lx2, Ly2, Lz2;
  double Pc11[3], Pc12[3], Pc13[3], Pc14[3], Pc15[3], Pc16[3], Pc17[3], Pc18[3];
  double Pc21[3], Pc22[3], Pc23[3], Pc24[3], Pc25[3], Pc26[3], Pc27[3], Pc28[3];
  double Pw11[3], Pw12[3], Pw13[3], Pw14[3], Pw15[3], Pw16[3], Pw17[3], Pw18[3];
  double Pw21[3], Pw22[3], Pw23[3], Pw24[3], Pw25[3], Pw26[3], Pw27[3], Pw28[3];
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
  void applyRT(double pout[3], double pin[3], double R[9], double T[3]);
  void averageSimilarValues( void );

  /* -- To initialize some values to zero -- */
  void defaultInit( void );

};

#endif
