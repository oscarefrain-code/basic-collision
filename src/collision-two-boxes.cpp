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


#include <collision-two-boxes.h>
#include <iostream>

#include <Eigen/Dense>
using namespace Eigen;

/* ------------------ Constructors ------------------------- */
CollisionTwoBoxes::
CollisionTwoBoxes( void )
{
  defaultInit();
}


CollisionTwoBoxes::
CollisionTwoBoxes(const double lx1, const double ly1, const double lz1,
		  const double lx2, const double ly2, const double lz2)
{
  B1.setLength(lx1, ly1, lz1);
  B2.setLength(lx2, ly2, lz2);
  defaultInit();
}


void CollisionTwoBoxes::
defaultInit( void )
{
  collisionIndicator = 0;
  tolerance = 0.0005; //TRI_EPSILON;
}


/* -- Set/get tolerances used to average close values (distance) -- */
void CollisionTwoBoxes::
setTolerance(const double toler)
{
  tolerance = toler;
}

double CollisionTwoBoxes::
getTolerance( void )
{
  return tolerance;
}


/* -- Set/get the coplanar tolerance (of the low level triangles) -- 
        This tolerance is to decide when the triangles are coplanar    */
void CollisionTwoBoxes::
setCoplanarTolerance(const double toler)
{
  To1_To2.setCoplanarTolerance(toler);  To1_Bo2.setCoplanarTolerance(toler);
  To1_Le2.setCoplanarTolerance(toler);  To1_Ri2.setCoplanarTolerance(toler);
  To1_Fr2.setCoplanarTolerance(toler);  To1_Ba2.setCoplanarTolerance(toler);
  Bo1_To2.setCoplanarTolerance(toler);  Bo1_Bo2.setCoplanarTolerance(toler);
  Bo1_Le2.setCoplanarTolerance(toler);  Bo1_Ri2.setCoplanarTolerance(toler);
  Bo1_Fr2.setCoplanarTolerance(toler);  Bo1_Ba2.setCoplanarTolerance(toler);
  Le1_To2.setCoplanarTolerance(toler);  Le1_Bo2.setCoplanarTolerance(toler);
  Le1_Le2.setCoplanarTolerance(toler);  Le1_Ri2.setCoplanarTolerance(toler);
  Le1_Fr2.setCoplanarTolerance(toler);  Le1_Ba2.setCoplanarTolerance(toler);
  Ri1_To2.setCoplanarTolerance(toler);  Ri1_Bo2.setCoplanarTolerance(toler);
  Ri1_Le2.setCoplanarTolerance(toler);  Ri1_Ri2.setCoplanarTolerance(toler);
  Ri1_Fr2.setCoplanarTolerance(toler);  Ri1_Ba2.setCoplanarTolerance(toler);
  Fr1_To2.setCoplanarTolerance(toler);  Fr1_Bo2.setCoplanarTolerance(toler);
  Fr1_Le2.setCoplanarTolerance(toler);  Fr1_Ri2.setCoplanarTolerance(toler);
  Fr1_Fr2.setCoplanarTolerance(toler);  Fr1_Ba2.setCoplanarTolerance(toler);
  Ba1_To2.setCoplanarTolerance(toler);  Ba1_Bo2.setCoplanarTolerance(toler);
  Ba1_Le2.setCoplanarTolerance(toler);  Ba1_Ri2.setCoplanarTolerance(toler);
  Ba1_Fr2.setCoplanarTolerance(toler);  Ba1_Ba2.setCoplanarTolerance(toler);
}

double CollisionTwoBoxes::
getCoplanarTolerance( void )
{
  double t[36];
  t[0] =To1_To2.getCoplanarTolerance();  t[1] =To1_Bo2.getCoplanarTolerance();
  t[2] =To1_Le2.getCoplanarTolerance();  t[3] =To1_Ri2.getCoplanarTolerance();
  t[4] =To1_Fr2.getCoplanarTolerance();  t[5] =To1_Ba2.getCoplanarTolerance();
  t[6] =Bo1_To2.getCoplanarTolerance();  t[7] =Bo1_Bo2.getCoplanarTolerance();
  t[8] =Bo1_Le2.getCoplanarTolerance();  t[9] =Bo1_Ri2.getCoplanarTolerance();
  t[10]=Bo1_Fr2.getCoplanarTolerance();  t[11]=Bo1_Ba2.getCoplanarTolerance();
  t[12]=Le1_To2.getCoplanarTolerance();  t[13]=Le1_Bo2.getCoplanarTolerance();
  t[14]=Le1_Le2.getCoplanarTolerance();  t[15]=Le1_Ri2.getCoplanarTolerance();
  t[16]=Le1_Fr2.getCoplanarTolerance();  t[17]=Le1_Ba2.getCoplanarTolerance();
  t[18]=Ri1_To2.getCoplanarTolerance();  t[19]=Ri1_Bo2.getCoplanarTolerance();
  t[20]=Ri1_Le2.getCoplanarTolerance();  t[21]=Ri1_Ri2.getCoplanarTolerance();
  t[22]=Ri1_Fr2.getCoplanarTolerance();  t[23]=Ri1_Ba2.getCoplanarTolerance();
  t[24]=Fr1_To2.getCoplanarTolerance();  t[25]=Fr1_Bo2.getCoplanarTolerance();
  t[26]=Fr1_Le2.getCoplanarTolerance();  t[27]=Fr1_Ri2.getCoplanarTolerance();
  t[28]=Fr1_Fr2.getCoplanarTolerance();  t[29]=Fr1_Ba2.getCoplanarTolerance();
  t[30]=Ba1_To2.getCoplanarTolerance();  t[31]=Ba1_Bo2.getCoplanarTolerance();
  t[32]=Ba1_Le2.getCoplanarTolerance();  t[33]=Ba1_Ri2.getCoplanarTolerance();
  t[34]=Ba1_Fr2.getCoplanarTolerance();  t[35]=Ba1_Ba2.getCoplanarTolerance();

  unsigned int sum=0;
  for(int i=0; i<35; i++){
    if (t[i]==t[i+1]){}
    else { sum++; }
  }

  if (sum==0)
    return (t[0]);
  else
    std::cerr << "TwoBoxes: Coplanar tolerance is not the same for all "
	      << "the triangles !!!" << std::endl;
}


/* ----------------- Set values for the box ---------------- */
void CollisionTwoBoxes::
setLengthB1(const double lx1, const double ly1, const double lz1)
{
  B1.setLength(lx1, ly1, lz1);
}


void CollisionTwoBoxes::
setLengthB2(const double lx2, const double ly2, const double lz2)
{
  B2.setLength(lx2, ly2, lz2);
}


/* ------------------ Get values from the box--------------- */

void CollisionTwoBoxes::
getLengthB1(double & lx1, double & ly1, double & lz1)
{
  B1.getLength(lx1, ly1, lz1);
}


void CollisionTwoBoxes::
getLengthB2(double & lx2, double & ly2, double & lz2)
{
  B2.getLength(lx2, ly2, lz2);
}

void CollisionTwoBoxes::
getTransformation1(Matrix3d &Rout, Vector3d &Tout)
{
  B1.getTransformation(Rout, Tout);
}
 

void CollisionTwoBoxes::
getTransformation2(Matrix3d &Rout, Vector3d &Tout)
{
  B2.getTransformation(Rout, Tout);
}


void CollisionTwoBoxes::
getVerticesB1(Vector3d &V1, Vector3d &V2, Vector3d &V3, Vector3d &V4,
	      Vector3d &V5, Vector3d &V6, Vector3d &V7, Vector3d &V8)
{
  B1.getVertices(V1, V2, V3, V4, V5, V6, V7, V8);
}


void CollisionTwoBoxes::
getVerticesB2(Vector3d &V1, Vector3d &V2, Vector3d &V3, Vector3d &V4,
	      Vector3d &V5, Vector3d &V6, Vector3d &V7, Vector3d &V8)
{
  B2.getVertices(V1, V2, V3, V4, V5, V6, V7, V8);
}


/* -------- Set the transformations for the boxes ---------- */

void CollisionTwoBoxes::
setTransformation1(const Matrix3d & R1, const Vector3d & T1)
{
  B1.setTransformation(R1, T1);

  /* Get the vertices */
  Vector3d B1v1, B1v2, B1v3, B1v4, B1v5, B1v6, B1v7, B1v8;
  B1.getVertices(B1v1, B1v2, B1v3, B1v4, B1v5, B1v6, B1v7, B1v8);

  // Create the Rectangles for the first box
  To1_To2.setVerticesR1(B1v1, B1v2, B1v6, B1v5);
  To1_Bo2.setVerticesR1(B1v1, B1v2, B1v6, B1v5);
  To1_Le2.setVerticesR1(B1v1, B1v2, B1v6, B1v5);
  To1_Ri2.setVerticesR1(B1v1, B1v2, B1v6, B1v5);
  To1_Fr2.setVerticesR1(B1v1, B1v2, B1v6, B1v5);
  To1_Ba2.setVerticesR1(B1v1, B1v2, B1v6, B1v5);

  Bo1_To2.setVerticesR1(B1v4, B1v3, B1v7, B1v8);
  Bo1_Bo2.setVerticesR1(B1v4, B1v3, B1v7, B1v8);
  Bo1_Le2.setVerticesR1(B1v4, B1v3, B1v7, B1v8);
  Bo1_Ri2.setVerticesR1(B1v4, B1v3, B1v7, B1v8);
  Bo1_Fr2.setVerticesR1(B1v4, B1v3, B1v7, B1v8);
  Bo1_Ba2.setVerticesR1(B1v4, B1v3, B1v7, B1v8);

  Le1_To2.setVerticesR1(B1v1, B1v2, B1v3, B1v4);
  Le1_Bo2.setVerticesR1(B1v1, B1v2, B1v3, B1v4);
  Le1_Le2.setVerticesR1(B1v1, B1v2, B1v3, B1v4);
  Le1_Ri2.setVerticesR1(B1v1, B1v2, B1v3, B1v4);
  Le1_Fr2.setVerticesR1(B1v1, B1v2, B1v3, B1v4);
  Le1_Ba2.setVerticesR1(B1v1, B1v2, B1v3, B1v4);

  Ri1_To2.setVerticesR1(B1v5, B1v6, B1v7, B1v8);
  Ri1_Bo2.setVerticesR1(B1v5, B1v6, B1v7, B1v8);
  Ri1_Le2.setVerticesR1(B1v5, B1v6, B1v7, B1v8);
  Ri1_Ri2.setVerticesR1(B1v5, B1v6, B1v7, B1v8);
  Ri1_Fr2.setVerticesR1(B1v5, B1v6, B1v7, B1v8);
  Ri1_Ba2.setVerticesR1(B1v5, B1v6, B1v7, B1v8);

  Fr1_To2.setVerticesR1(B1v1, B1v4, B1v8, B1v5);
  Fr1_Bo2.setVerticesR1(B1v1, B1v4, B1v8, B1v5);
  Fr1_Le2.setVerticesR1(B1v1, B1v4, B1v8, B1v5);
  Fr1_Ri2.setVerticesR1(B1v1, B1v4, B1v8, B1v5);
  Fr1_Fr2.setVerticesR1(B1v1, B1v4, B1v8, B1v5);
  Fr1_Ba2.setVerticesR1(B1v1, B1v4, B1v8, B1v5);

  Ba1_To2.setVerticesR1(B1v2, B1v3, B1v7, B1v6);
  Ba1_Bo2.setVerticesR1(B1v2, B1v3, B1v7, B1v6);
  Ba1_Le2.setVerticesR1(B1v2, B1v3, B1v7, B1v6);
  Ba1_Ri2.setVerticesR1(B1v2, B1v3, B1v7, B1v6);
  Ba1_Fr2.setVerticesR1(B1v2, B1v3, B1v7, B1v6);
  Ba1_Ba2.setVerticesR1(B1v2, B1v3, B1v7, B1v6);
}

void CollisionTwoBoxes::
setTransformation2( const Matrix3d & R2, const Vector3d & T2 )
{
  B2.setTransformation( R2, T2 );

  /* Get the vertices */
  Vector3d B2v1, B2v2, B2v3, B2v4, B2v5, B2v6, B2v7, B2v8;
  B2.getVertices(B2v1, B2v2, B2v3, B2v4, B2v5, B2v6, B2v7, B2v8);

  // Create the Rectangles for the second box
  To1_To2.setVerticesR2(B2v1, B2v2, B2v6, B2v5);
  To1_Bo2.setVerticesR2(B2v4, B2v3, B2v7, B2v8);
  To1_Le2.setVerticesR2(B2v1, B2v2, B2v3, B2v4);
  To1_Ri2.setVerticesR2(B2v5, B2v6, B2v7, B2v8);
  To1_Fr2.setVerticesR2(B2v1, B2v4, B2v8, B2v5);
  To1_Ba2.setVerticesR2(B2v2, B2v3, B2v7, B2v6);

  Bo1_To2.setVerticesR2(B2v1, B2v2, B2v6, B2v5);
  Bo1_Bo2.setVerticesR2(B2v4, B2v3, B2v7, B2v8);
  Bo1_Le2.setVerticesR2(B2v1, B2v2, B2v3, B2v4);
  Bo1_Ri2.setVerticesR2(B2v5, B2v6, B2v7, B2v8);
  Bo1_Fr2.setVerticesR2(B2v1, B2v4, B2v8, B2v5);
  Bo1_Ba2.setVerticesR2(B2v2, B2v3, B2v7, B2v6);

  Le1_To2.setVerticesR2(B2v1, B2v2, B2v6, B2v5);
  Le1_Bo2.setVerticesR2(B2v4, B2v3, B2v7, B2v8);
  Le1_Le2.setVerticesR2(B2v1, B2v2, B2v3, B2v4);
  Le1_Ri2.setVerticesR2(B2v5, B2v6, B2v7, B2v8);
  Le1_Fr2.setVerticesR2(B2v1, B2v4, B2v8, B2v5);
  Le1_Ba2.setVerticesR2(B2v2, B2v3, B2v7, B2v6);

  Ri1_To2.setVerticesR2(B2v1, B2v2, B2v6, B2v5);
  Ri1_Bo2.setVerticesR2(B2v4, B2v3, B2v7, B2v8);
  Ri1_Le2.setVerticesR2(B2v1, B2v2, B2v3, B2v4);
  Ri1_Ri2.setVerticesR2(B2v5, B2v6, B2v7, B2v8);
  Ri1_Fr2.setVerticesR2(B2v1, B2v4, B2v8, B2v5);
  Ri1_Ba2.setVerticesR2(B2v2, B2v3, B2v7, B2v6);

  Fr1_To2.setVerticesR2(B2v1, B2v2, B2v6, B2v5);
  Fr1_Bo2.setVerticesR2(B2v4, B2v3, B2v7, B2v8);
  Fr1_Le2.setVerticesR2(B2v1, B2v2, B2v3, B2v4);
  Fr1_Ri2.setVerticesR2(B2v5, B2v6, B2v7, B2v8);
  Fr1_Fr2.setVerticesR2(B2v1, B2v4, B2v8, B2v5);
  Fr1_Ba2.setVerticesR2(B2v2, B2v3, B2v7, B2v6);

  Ba1_To2.setVerticesR2(B2v1, B2v2, B2v6, B2v5);
  Ba1_Bo2.setVerticesR2(B2v4, B2v3, B2v7, B2v8);
  Ba1_Le2.setVerticesR2(B2v1, B2v2, B2v3, B2v4);
  Ba1_Ri2.setVerticesR2(B2v5, B2v6, B2v7, B2v8);
  Ba1_Fr2.setVerticesR2(B2v1, B2v4, B2v8, B2v5);
  Ba1_Ba2.setVerticesR2(B2v2, B2v3, B2v7, B2v6);

}


/* -------------- Calculate the collisions ----------------- */

int CollisionTwoBoxes::
computeBBintersections( void )
{
  int collisions[36];

  /* Initialization to zero of collisions indicators */
  for (int i=0; i<36; i++)
    collisions[i]=0;

  /* Clear the possible 'previous' contact points */
  pointsBB.clear();

  // collisions[0] = To1_To2.computeRRintersections(); // 
  // collisions[1] = To1_Bo2.computeRRintersections(); //
  // collisions[2] = To1_Le2.computeRRintersections(); // 
  // collisions[3] = To1_Ri2.computeRRintersections(); // 
  // collisions[4] = To1_Fr2.computeRRintersections(); // 
  // collisions[5] = To1_Ba2.computeRRintersections(); //

  collisions[6]  = Bo1_To2.computeRRintersections(); // y
  // collisions[7]  = Bo1_Bo2.computeRRintersections(); //
  // collisions[8]  = Bo1_Le2.computeRRintersections(); //
  // collisions[9]  = Bo1_Ri2.computeRRintersections(); // 
  // collisions[10] = Bo1_Fr2.computeRRintersections(); //
  // collisions[11] = Bo1_Ba2.computeRRintersections(); // 

  // collisions[12] = Le1_To2.computeRRintersections(); //
  // collisions[13] = Le1_Bo2.computeRRintersections(); //
  // collisions[14] = Le1_Le2.computeRRintersections(); //
  // collisions[15] = Le1_Ri2.computeRRintersections(); //
  // collisions[16] = Le1_Fr2.computeRRintersections(); //
  // collisions[17] = Le1_Ba2.computeRRintersections(); //

  // collisions[18] = Ri1_To2.computeRRintersections(); //
  // collisions[19] = Ri1_Bo2.computeRRintersections(); //
  // collisions[20] = Ri1_Le2.computeRRintersections(); //
  // collisions[21] = Ri1_Ri2.computeRRintersections(); //
  // collisions[22] = Ri1_Fr2.computeRRintersections(); //
  // collisions[23] = Ri1_Ba2.computeRRintersections(); //

  // collisions[24] = Fr1_To2.computeRRintersections(); // 
  // collisions[25] = Fr1_Bo2.computeRRintersections(); //
  // collisions[26] = Fr1_Le2.computeRRintersections(); //
  // collisions[27] = Fr1_Ri2.computeRRintersections(); //
  // collisions[28] = Fr1_Fr2.computeRRintersections(); //
  // collisions[29] = Fr1_Ba2.computeRRintersections(); //

  // collisions[30] = Ba1_To2.computeRRintersections(); //
  // collisions[31] = Ba1_Bo2.computeRRintersections(); //
  // collisions[32] = Ba1_Le2.computeRRintersections(); //
  // collisions[33] = Ba1_Ri2.computeRRintersections(); //
  // collisions[34] = Ba1_Fr2.computeRRintersections(); //
  // collisions[35] = Ba1_Ba2.computeRRintersections(); //

  // Collision detected? (yes:1, no:0)
  // collisionIndicator = collisions[0] || collisions[1] || collisions[2] ||
  //   collisions[3] || collisions[4] || collisions[5] || collisions[6];
  //collisionIndicator = collisions[6];

  collisionIndicator = collisions[0] || collisions[1] || collisions[2] ||
    collisions[3] || collisions[4] || collisions[5] || collisions[6] ||
    collisions[7] || collisions[8] || collisions[9] || collisions[10] ||
    collisions[11] || collisions[12] || collisions[13] || collisions[14] ||
    collisions[15] || collisions[16] || collisions[17] || collisions[18] ||
    collisions[19] || collisions[20] || collisions[21] || collisions[22] ||
    collisions[23] || collisions[24] || collisions[25] || collisions[26] ||
    collisions[27] || collisions[28] || collisions[29] || collisions[30] ||
    collisions[31] || collisions[32] || collisions[33] || collisions[34] ||
    collisions[35];

  std::cout << "Collisions : ";
  //for (int i=0; i<12; i++)
  for (int i=0; i<36; i++)
    {
      std::cout << collisions[i] << " ";
      if ((i+1)%6 == 0) std::cout << " - ";
    }
  std::cout << std::endl;

  // //printBoxesVertices();

  // for (int i=0; i<Bo1_To2.pointsRR.size(); i++)
  //   std::cout << " (" << Bo1_To2.pointsRR[i].transpose() << ") " << std::endl;

  // for (int i=0; i<To1_Ri2.pointsRR.size(); i++)
  //   std::cout << " (" << To1_Ri2.pointsRR[i].transpose() << ") " << std::endl;
      
  // Copy all the collision points in pointsRR
  if (collisionIndicator) {
    // for (int i=0; i<To1_To2.pointsRR.size(); i++)
    //   pointsBB.push_back(To1_To2.pointsRR[i]);
    // for (int i=0; i<To1_Bo2.pointsRR.size(); i++)
    //   pointsBB.push_back(To1_Bo2.pointsRR[i]);
    // for (int i=0; i<To1_Le2.pointsRR.size(); i++)
    //   pointsBB.push_back(To1_Le2.pointsRR[i]);
    // for (int i=0; i<To1_Ri2.pointsRR.size(); i++)
    //   pointsBB.push_back(To1_Ri2.pointsRR[i]);
    // for (int i=0; i<To1_Fr2.pointsRR.size(); i++)
    //   pointsBB.push_back(To1_Fr2.pointsRR[i]);
    // for (int i=0; i<To1_Ba2.pointsRR.size(); i++)
    //   pointsBB.push_back(To1_Ba2.pointsRR[i]);

    for (int i=0; i<Bo1_To2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_To2.pointsRR[i]);
    // for (int i=0; i<Bo1_Bo2.pointsRR.size(); i++)
    //   pointsBB.push_back(Bo1_Bo2.pointsRR[i]);
    // for (int i=0; i<Bo1_Le2.pointsRR.size(); i++)
    //   pointsBB.push_back(Bo1_Le2.pointsRR[i]);
    // for (int i=0; i<Bo1_Ri2.pointsRR.size(); i++)
    //   pointsBB.push_back(Bo1_Ri2.pointsRR[i]);
    // for (int i=0; i<Bo1_Fr2.pointsRR.size(); i++)
    //   pointsBB.push_back(Bo1_Fr2.pointsRR[i]);
    // for (int i=0; i<Bo1_Ba2.pointsRR.size(); i++)
    //   pointsBB.push_back(Bo1_Ba2.pointsRR[i]);

    // for (int i=0; i<Le1_To2.pointsRR.size(); i++)
    //   pointsBB.push_back(Le1_To2.pointsRR[i]);
    // for (int i=0; i<Le1_Bo2.pointsRR.size(); i++)
    //   pointsBB.push_back(Le1_Bo2.pointsRR[i]);
    // for (int i=0; i<Le1_Le2.pointsRR.size(); i++)
    //   pointsBB.push_back(Le1_Le2.pointsRR[i]);
    // for (int i=0; i<Le1_Ri2.pointsRR.size(); i++)
    //   pointsBB.push_back(Le1_Ri2.pointsRR[i]);
    // for (int i=0; i<Le1_Fr2.pointsRR.size(); i++)
    //   pointsBB.push_back(Le1_Fr2.pointsRR[i]);
    // for (int i=0; i<Le1_Ba2.pointsRR.size(); i++)
    //   pointsBB.push_back(Le1_Ba2.pointsRR[i]);

    // for (int i=0; i<Ri1_To2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ri1_To2.pointsRR[i]);
    // for (int i=0; i<Ri1_Bo2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ri1_Bo2.pointsRR[i]);
    // for (int i=0; i<Ri1_Le2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ri1_Le2.pointsRR[i]);
    // for (int i=0; i<Ri1_Ri2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ri1_Ri2.pointsRR[i]);
    // for (int i=0; i<Ri1_Fr2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ri1_Fr2.pointsRR[i]);
    // for (int i=0; i<Ri1_Ba2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ri1_Ba2.pointsRR[i]);

    // for (int i=0; i<Fr1_To2.pointsRR.size(); i++)
    //   pointsBB.push_back(Fr1_To2.pointsRR[i]);
    // for (int i=0; i<Fr1_Bo2.pointsRR.size(); i++)
    //   pointsBB.push_back(Fr1_Bo2.pointsRR[i]);
    // for (int i=0; i<Fr1_Le2.pointsRR.size(); i++)
    //   pointsBB.push_back(Fr1_Le2.pointsRR[i]);
    // for (int i=0; i<Fr1_Ri2.pointsRR.size(); i++)
    //   pointsBB.push_back(Fr1_Ri2.pointsRR[i]);
    // for (int i=0; i<Fr1_Fr2.pointsRR.size(); i++)
    //   pointsBB.push_back(Fr1_Fr2.pointsRR[i]);
    // for (int i=0; i<Fr1_Ba2.pointsRR.size(); i++)
    //   pointsBB.push_back(Fr1_Ba2.pointsRR[i]);

    // for (int i=0; i<Ba1_To2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ba1_To2.pointsRR[i]);
    // for (int i=0; i<Ba1_Bo2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ba1_Bo2.pointsRR[i]);
    // for (int i=0; i<Ba1_Le2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ba1_Le2.pointsRR[i]);
    // for (int i=0; i<Ba1_Ri2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ba1_Ri2.pointsRR[i]);
    // for (int i=0; i<Ba1_Fr2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ba1_Fr2.pointsRR[i]);
    // for (int i=0; i<Ba1_Ba2.pointsRR.size(); i++)
    //   pointsBB.push_back(Ba1_Ba2.pointsRR[i]);

  }

  prunePoints(pointsBB);
  averageSimilarValues();

  //printBBcollisionInformation();

  return collisionIndicator;
}


/* --- Get the average of very close points (to a certain tolerance) --- */
// void CollisionTwoBoxes::
// averageSimilarValues( void )
// {
//   double distance;
//   std::vector<int> indexToDiscard;
//   std::vector<Vector3d> result;
//   Vector3d sum; int Nsum;

//   for (int i=0; i<pointsBB.size(); i++){

//     /* Use i if it does not belong to the discarded indexes*/
//     if (find(indexToDiscard.begin(), indexToDiscard.end(), i) == indexToDiscard.end())
//       {
//   	/* Initialization of temporal variables*/
//   	sum = pointsBB[i]; Nsum = 1;

//   	for (int j=i+1; j<pointsBB.size(); j++){
//   	  /* Use j if it does not belong to the discarded indexes*/
//   	  if ( find(indexToDiscard.begin(), indexToDiscard.end(), j) == indexToDiscard.end())
//   	    {
//   	      distance = (pointsBB[i](0)-pointsBB[j](0))*(pointsBB[i](0)-pointsBB[j](0)) +
//   		(pointsBB[i](1)-pointsBB[j](1))*(pointsBB[i](1)-pointsBB[j](1)) +
//   		(pointsBB[i](2)-pointsBB[j](2))*(pointsBB[i](2)-pointsBB[j](2));
//   	      if ( distance < (tolerance*tolerance) )
//   		{
// 		  //std::cout << " ... distance smaller than tolerance -> averaging ..." << std::endl;
//   		  sum(0) += pointsBB[j](0); sum(1) += pointsBB[j](1); sum(2) += pointsBB[j](2);
//   		  Nsum++;
//   		  indexToDiscard.push_back(j);
//   		}
//   	    }
//   	}
// 	/* Get the average and push it back to 'result' */
//   	sum = (1.0/double(Nsum))*sum;
//  	result.push_back(sum);
// 	// std::cout << "sum: " << sum.transpose() << ", Nsum: " << Nsum << ", Div: " << average.transpose() << std::endl;
//       }
//   }

//   pointsBB.assign(result.begin(), result.end());

//   // std::cout << "Result: ";
//   // for (int i=0; i<result.size(); i++){ std::cout << "  (" << result[i].transpose() << ")  ";  }
//   // std::cout << std::endl;

// }



void CollisionTwoBoxes::
averageSimilarValues( void )
{
  double dist;
  for (int i=0; i<pointsBB.size(); i++){
    for (int j=i+1; j<pointsBB.size(); j++){
      dist = (pointsBB[i](0)-pointsBB[j](0))*(pointsBB[i](0)-pointsBB[j](0)) +
	     (pointsBB[i](1)-pointsBB[j](1))*(pointsBB[i](1)-pointsBB[j](1)) +
             (pointsBB[i](2)-pointsBB[j](2))*(pointsBB[i](2)-pointsBB[j](2));

      if ( dist < (tolerance*tolerance) ) {
	pointsBB[i](0) = (pointsBB[i](0)+pointsBB[j](0))/2;
	pointsBB[i](1) = (pointsBB[i](1)+pointsBB[j](1))/2;
	pointsBB[i](2) = (pointsBB[i](2)+pointsBB[j](2))/2;
	pointsBB.erase( pointsBB.begin()+j );
	j--;
      }
    }
  }
}


/* ----------------- PRINT INFORMATION --------------------- */

void CollisionTwoBoxes::
printLengths( void )
{
  std::cout << "Length box1: "; B1.printLengths();
  std::cout << "Length box2: "; B2.printLengths();
}


void CollisionTwoBoxes::
printBBcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Boxes are intersecting" << std::endl;
    for (int i=0; i<pointsBB.size(); i++){
      std::cout << "  (" << pointsBB[i].transpose() << ")  ";
    }
    std::cout << std::endl;
  }
  else{
    std::cout << "Result: Boxes are NOT intersecting" << std::endl;
  }
}


void CollisionTwoBoxes::
printBoxesVertices( void )
{
  std::cout << "Vertices Box 1: "; B1.printVertices();
  std::cout << "Vertices Box 2: "; B2.printVertices();

}


void CollisionTwoBoxes::
printTransformation( void )
{
  std::cout << "Transformation box 1: " << std::endl; B1.printTransformation();
  std::cout << "Transformation box 2: " << std::endl; B2.printTransformation();
}

