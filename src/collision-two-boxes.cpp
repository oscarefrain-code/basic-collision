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


/* --- Constructors --- */
CollisionTwoBoxes::
CollisionTwoBoxes( void )
{
  defaultInit();
}


CollisionTwoBoxes::
CollisionTwoBoxes(double lx1, double ly1, double lz1,
		  double lx2, double ly2, double lz2)
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

/* -- Tolerance related -- */
void CollisionTwoBoxes::
setTolerance(double toler)
{
  tolerance = toler;
}

double CollisionTwoBoxes::
getTolerance( void )
{
  return tolerance;
}


/* --- Set values for the box --- */
void CollisionTwoBoxes::
setLengthB1(double lx1, double ly1, double lz1)
{
  B1.setLength(lx1, ly1, lz1);
}


void CollisionTwoBoxes::
setLengthB2(double lx2, double ly2, double lz2)
{
  B2.setLength(lx2, ly2, lz2);
}

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


/* --- Set the transformations for the boxes --- */ 
void CollisionTwoBoxes::
setTransformation1(Matrix3d R1, Vector3d T1)
{
  B1.setTransformation(R1, T1);

  // Create the Rectangles for the first box
  To1_To2.setVerticesR1(B1.v1, B1.v2, B1.v6, B1.v5);
  To1_Bo2.setVerticesR1(B1.v1, B1.v2, B1.v6, B1.v5);
  To1_Le2.setVerticesR1(B1.v1, B1.v2, B1.v6, B1.v5);
  To1_Ri2.setVerticesR1(B1.v1, B1.v2, B1.v6, B1.v5);
  To1_Fr2.setVerticesR1(B1.v1, B1.v2, B1.v6, B1.v5);
  To1_Ba2.setVerticesR1(B1.v1, B1.v2, B1.v6, B1.v5);

  Bo1_To2.setVerticesR1(B1.v4, B1.v3, B1.v7, B1.v8);
  Bo1_Bo2.setVerticesR1(B1.v4, B1.v3, B1.v7, B1.v8);
  Bo1_Le2.setVerticesR1(B1.v4, B1.v3, B1.v7, B1.v8);
  Bo1_Ri2.setVerticesR1(B1.v4, B1.v3, B1.v7, B1.v8);
  Bo1_Fr2.setVerticesR1(B1.v4, B1.v3, B1.v7, B1.v8);
  Bo1_Ba2.setVerticesR1(B1.v4, B1.v3, B1.v7, B1.v8);

  Le1_To2.setVerticesR1(B1.v1, B1.v2, B1.v3, B1.v4);
  Le1_Bo2.setVerticesR1(B1.v1, B1.v2, B1.v3, B1.v4);
  Le1_Le2.setVerticesR1(B1.v1, B1.v2, B1.v3, B1.v4);
  Le1_Ri2.setVerticesR1(B1.v1, B1.v2, B1.v3, B1.v4);
  Le1_Fr2.setVerticesR1(B1.v1, B1.v2, B1.v3, B1.v4);
  Le1_Ba2.setVerticesR1(B1.v1, B1.v2, B1.v3, B1.v4);

  Ri1_To2.setVerticesR1(B1.v5, B1.v6, B1.v7, B1.v8);
  Ri1_Bo2.setVerticesR1(B1.v5, B1.v6, B1.v7, B1.v8);
  Ri1_Le2.setVerticesR1(B1.v5, B1.v6, B1.v7, B1.v8);
  Ri1_Ri2.setVerticesR1(B1.v5, B1.v6, B1.v7, B1.v8);
  Ri1_Fr2.setVerticesR1(B1.v5, B1.v6, B1.v7, B1.v8);
  Ri1_Ba2.setVerticesR1(B1.v5, B1.v6, B1.v7, B1.v8);

  Fr1_To2.setVerticesR1(B1.v1, B1.v4, B1.v8, B1.v5);
  Fr1_Bo2.setVerticesR1(B1.v1, B1.v4, B1.v8, B1.v5);
  Fr1_Le2.setVerticesR1(B1.v1, B1.v4, B1.v8, B1.v5);
  Fr1_Ri2.setVerticesR1(B1.v1, B1.v4, B1.v8, B1.v5);
  Fr1_Fr2.setVerticesR1(B1.v1, B1.v4, B1.v8, B1.v5);
  Fr1_Ba2.setVerticesR1(B1.v1, B1.v4, B1.v8, B1.v5);

  Ba1_To2.setVerticesR1(B1.v2, B1.v3, B1.v7, B1.v6);
  Ba1_Bo2.setVerticesR1(B1.v2, B1.v3, B1.v7, B1.v6);
  Ba1_Le2.setVerticesR1(B1.v2, B1.v3, B1.v7, B1.v6);
  Ba1_Ri2.setVerticesR1(B1.v2, B1.v3, B1.v7, B1.v6);
  Ba1_Fr2.setVerticesR1(B1.v2, B1.v3, B1.v7, B1.v6);
  Ba1_Ba2.setVerticesR1(B1.v2, B1.v3, B1.v7, B1.v6);
}

void CollisionTwoBoxes::
setTransformation2( Matrix3d R2, Vector3d T2 )
{
  B2.setTransformation( R2, T2 );

  // Create the Rectangles for the second box
  To1_To2.setVerticesR2(B2.v1, B2.v2, B2.v6, B2.v5);
  To1_Bo2.setVerticesR2(B2.v4, B2.v3, B2.v7, B2.v8);
  To1_Le2.setVerticesR2(B2.v1, B2.v2, B2.v3, B2.v4);
  To1_Ri2.setVerticesR2(B2.v5, B2.v6, B2.v7, B2.v8);
  To1_Fr2.setVerticesR2(B2.v1, B2.v4, B2.v8, B2.v5);
  To1_Ba2.setVerticesR2(B2.v2, B2.v3, B2.v7, B2.v6);

  Bo1_To2.setVerticesR2(B2.v1, B2.v2, B2.v6, B2.v5);
  Bo1_Bo2.setVerticesR2(B2.v4, B2.v3, B2.v7, B2.v8);
  Bo1_Le2.setVerticesR2(B2.v1, B2.v2, B2.v3, B2.v4);
  Bo1_Ri2.setVerticesR2(B2.v5, B2.v6, B2.v7, B2.v8);
  Bo1_Fr2.setVerticesR2(B2.v1, B2.v4, B2.v8, B2.v5);
  Bo1_Ba2.setVerticesR2(B2.v2, B2.v3, B2.v7, B2.v6);

  Le1_To2.setVerticesR2(B2.v1, B2.v2, B2.v6, B2.v5);
  Le1_Bo2.setVerticesR2(B2.v4, B2.v3, B2.v7, B2.v8);
  Le1_Le2.setVerticesR2(B2.v1, B2.v2, B2.v3, B2.v4);
  Le1_Ri2.setVerticesR2(B2.v5, B2.v6, B2.v7, B2.v8);
  Le1_Fr2.setVerticesR2(B2.v1, B2.v4, B2.v8, B2.v5);
  Le1_Ba2.setVerticesR2(B2.v2, B2.v3, B2.v7, B2.v6);

  Ri1_To2.setVerticesR2(B2.v1, B2.v2, B2.v6, B2.v5);
  Ri1_Bo2.setVerticesR2(B2.v4, B2.v3, B2.v7, B2.v8);
  Ri1_Le2.setVerticesR2(B2.v1, B2.v2, B2.v3, B2.v4);
  Ri1_Ri2.setVerticesR2(B2.v5, B2.v6, B2.v7, B2.v8);
  Ri1_Fr2.setVerticesR2(B2.v1, B2.v4, B2.v8, B2.v5);
  Ri1_Ba2.setVerticesR2(B2.v2, B2.v3, B2.v7, B2.v6);

  Fr1_To2.setVerticesR2(B2.v1, B2.v2, B2.v6, B2.v5);
  Fr1_Bo2.setVerticesR2(B2.v4, B2.v3, B2.v7, B2.v8);
  Fr1_Le2.setVerticesR2(B2.v1, B2.v2, B2.v3, B2.v4);
  Fr1_Ri2.setVerticesR2(B2.v5, B2.v6, B2.v7, B2.v8);
  Fr1_Fr2.setVerticesR2(B2.v1, B2.v4, B2.v8, B2.v5);
  Fr1_Ba2.setVerticesR2(B2.v2, B2.v3, B2.v7, B2.v6);

  Ba1_To2.setVerticesR2(B2.v1, B2.v2, B2.v6, B2.v5);
  Ba1_Bo2.setVerticesR2(B2.v4, B2.v3, B2.v7, B2.v8);
  Ba1_Le2.setVerticesR2(B2.v1, B2.v2, B2.v3, B2.v4);
  Ba1_Ri2.setVerticesR2(B2.v5, B2.v6, B2.v7, B2.v8);
  Ba1_Fr2.setVerticesR2(B2.v1, B2.v4, B2.v8, B2.v5);
  Ba1_Ba2.setVerticesR2(B2.v2, B2.v3, B2.v7, B2.v6);
}


/* -- Calculate the collisions --- */
int CollisionTwoBoxes::
computeBBintersections( void )
{
  int collisions[36];

  /* Clear the possible 'previous' contact points */
  pointsBB.clear();

  collisions[0] = To1_To2.computeRRintersections(); // Y
  collisions[1] = To1_Bo2.computeRRintersections(); //
  collisions[2] = To1_Le2.computeRRintersections(); // 
  collisions[3] = To1_Ri2.computeRRintersections(); // Y
  collisions[4] = To1_Fr2.computeRRintersections(); // 
  collisions[5] = To1_Ba2.computeRRintersections();

  collisions[6]  = Bo1_To2.computeRRintersections(); // Y
  collisions[7]  = Bo1_Bo2.computeRRintersections();
  collisions[8]  = Bo1_Le2.computeRRintersections();
  collisions[9]  = Bo1_Ri2.computeRRintersections(); // Y
  collisions[10] = Bo1_Fr2.computeRRintersections();
  collisions[11] = Bo1_Ba2.computeRRintersections();

  collisions[12] = Le1_To2.computeRRintersections();
  collisions[13] = Le1_Bo2.computeRRintersections();
  collisions[14] = Le1_Le2.computeRRintersections();
  collisions[15] = Le1_Ri2.computeRRintersections();
  collisions[16] = Le1_Fr2.computeRRintersections();
  collisions[17] = Le1_Ba2.computeRRintersections();

  collisions[18] = Ri1_To2.computeRRintersections();
  collisions[19] = Ri1_Bo2.computeRRintersections(); //
  collisions[20] = Ri1_Le2.computeRRintersections();
  collisions[21] = Ri1_Ri2.computeRRintersections();
  collisions[22] = Ri1_Fr2.computeRRintersections(); //
  collisions[23] = Ri1_Ba2.computeRRintersections();

  collisions[24] = Fr1_To2.computeRRintersections(); // Y
  collisions[25] = Fr1_Bo2.computeRRintersections();
  collisions[26] = Fr1_Le2.computeRRintersections();
  collisions[27] = Fr1_Ri2.computeRRintersections();
  collisions[28] = Fr1_Fr2.computeRRintersections();
  collisions[29] = Fr1_Ba2.computeRRintersections();

  collisions[30] = Ba1_To2.computeRRintersections();
  collisions[31] = Ba1_Bo2.computeRRintersections(); // 
  collisions[32] = Ba1_Le2.computeRRintersections(); // 
  collisions[33] = Ba1_Ri2.computeRRintersections();
  collisions[34] = Ba1_Fr2.computeRRintersections();
  collisions[35] = Ba1_Ba2.computeRRintersections();

  // Collision detected? (yes:1, no:0)
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

  // std::cout << "Collisions : ";
  // for (int i=0; i<36; i++)
  //   std::cout << collisions[i] << " ";
  // std::cout << std::endl;
  
  // for (int i=0; i<To1_To2.pointsRR.size(); i++)
  //   To1_To2.pointsRR[i].print();
  // std::cout << std::endl;

  // for (int i=0; i<To1_Ri2.pointsRR.size(); i++)
  //   To1_Ri2.pointsRR[i].print();
  // std::cout << std::endl;

  // for (int i=0; i<Bo1_To2.pointsRR.size(); i++)
  //   Bo1_To2.pointsRR[i].print();
  // std::cout << std::endl;

  // for (int i=0; i<Bo1_Ri2.pointsRR.size(); i++)
  //   Bo1_Ri2.pointsRR[i].print();
  // std::cout << std::endl;

  // for (int i=0; i<Fr1_To2.pointsRR.size(); i++)
  //   Fr1_To2.pointsRR[i].print();
  // std::cout << std::endl;

      
  // Copy all the collision points in pointsRR
  if (collisionIndicator) {
    for (int i=0; i<To1_To2.pointsRR.size(); i++)
      pointsBB.push_back(To1_To2.pointsRR[i]);
    for (int i=0; i<To1_Bo2.pointsRR.size(); i++)
      pointsBB.push_back(To1_Bo2.pointsRR[i]);
    for (int i=0; i<To1_Le2.pointsRR.size(); i++)
      pointsBB.push_back(To1_Le2.pointsRR[i]);
    for (int i=0; i<To1_Ri2.pointsRR.size(); i++)
      pointsBB.push_back(To1_Ri2.pointsRR[i]);
    for (int i=0; i<To1_Fr2.pointsRR.size(); i++)
      pointsBB.push_back(To1_Fr2.pointsRR[i]);
    for (int i=0; i<To1_Ba2.pointsRR.size(); i++)
      pointsBB.push_back(To1_Ba2.pointsRR[i]);

    for (int i=0; i<Bo1_To2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_To2.pointsRR[i]);
    for (int i=0; i<Bo1_Bo2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_Bo2.pointsRR[i]);
    for (int i=0; i<Bo1_Le2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_Le2.pointsRR[i]);
    for (int i=0; i<Bo1_Ri2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_Ri2.pointsRR[i]);
    for (int i=0; i<Bo1_Fr2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_Fr2.pointsRR[i]);
    for (int i=0; i<Bo1_Ba2.pointsRR.size(); i++)
      pointsBB.push_back(Bo1_Ba2.pointsRR[i]);

    for (int i=0; i<Le1_To2.pointsRR.size(); i++)
      pointsBB.push_back(Le1_To2.pointsRR[i]);
    for (int i=0; i<Le1_Bo2.pointsRR.size(); i++)
      pointsBB.push_back(Le1_Bo2.pointsRR[i]);
    for (int i=0; i<Le1_Le2.pointsRR.size(); i++)
      pointsBB.push_back(Le1_Le2.pointsRR[i]);
    for (int i=0; i<Le1_Ri2.pointsRR.size(); i++)
      pointsBB.push_back(Le1_Ri2.pointsRR[i]);
    for (int i=0; i<Le1_Fr2.pointsRR.size(); i++)
      pointsBB.push_back(Le1_Fr2.pointsRR[i]);
    for (int i=0; i<Le1_Ba2.pointsRR.size(); i++)
      pointsBB.push_back(Le1_Ba2.pointsRR[i]);

    for (int i=0; i<Ri1_To2.pointsRR.size(); i++)
      pointsBB.push_back(Ri1_To2.pointsRR[i]);
    for (int i=0; i<Ri1_Bo2.pointsRR.size(); i++)
      pointsBB.push_back(Ri1_Bo2.pointsRR[i]);
    for (int i=0; i<Ri1_Le2.pointsRR.size(); i++)
      pointsBB.push_back(Ri1_Le2.pointsRR[i]);
    for (int i=0; i<Ri1_Ri2.pointsRR.size(); i++)
      pointsBB.push_back(Ri1_Ri2.pointsRR[i]);
    for (int i=0; i<Ri1_Fr2.pointsRR.size(); i++)
      pointsBB.push_back(Ri1_Fr2.pointsRR[i]);
    for (int i=0; i<Ri1_Ba2.pointsRR.size(); i++)
      pointsBB.push_back(Ri1_Ba2.pointsRR[i]);

    for (int i=0; i<Fr1_To2.pointsRR.size(); i++)
      pointsBB.push_back(Fr1_To2.pointsRR[i]);
    for (int i=0; i<Fr1_Bo2.pointsRR.size(); i++)
      pointsBB.push_back(Fr1_Bo2.pointsRR[i]);
    for (int i=0; i<Fr1_Le2.pointsRR.size(); i++)
      pointsBB.push_back(Fr1_Le2.pointsRR[i]);
    for (int i=0; i<Fr1_Ri2.pointsRR.size(); i++)
      pointsBB.push_back(Fr1_Ri2.pointsRR[i]);
    for (int i=0; i<Fr1_Fr2.pointsRR.size(); i++)
      pointsBB.push_back(Fr1_Fr2.pointsRR[i]);
    for (int i=0; i<Fr1_Ba2.pointsRR.size(); i++)
      pointsBB.push_back(Fr1_Ba2.pointsRR[i]);

    for (int i=0; i<Ba1_To2.pointsRR.size(); i++)
      pointsBB.push_back(Ba1_To2.pointsRR[i]);
    for (int i=0; i<Ba1_Bo2.pointsRR.size(); i++)
      pointsBB.push_back(Ba1_Bo2.pointsRR[i]);
    for (int i=0; i<Ba1_Le2.pointsRR.size(); i++)
      pointsBB.push_back(Ba1_Le2.pointsRR[i]);
    for (int i=0; i<Ba1_Ri2.pointsRR.size(); i++)
      pointsBB.push_back(Ba1_Ri2.pointsRR[i]);
    for (int i=0; i<Ba1_Fr2.pointsRR.size(); i++)
      pointsBB.push_back(Ba1_Fr2.pointsRR[i]);
    for (int i=0; i<Ba1_Ba2.pointsRR.size(); i++)
      pointsBB.push_back(Ba1_Ba2.pointsRR[i]);

  }

  prunePoints(pointsBB);
  averageSimilarValues();

  return collisionIndicator;
}


/* --- Get the average of very close points (to a certain tolerance) --- */
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



/* --- Print Different Information --- */
void CollisionTwoBoxes::
printLengths( void )
{
  // std::cout << "Length: Box 1 = (" << B1.Lx << ", " << B1.Ly << ", " << B1.Lz << " ),   "
  // 	    << "Box 2 = (" << B2.Lx << ", " << B2.Ly << ", " << B2.Lz << " )" << std::endl;
  std::cout << "Length box1: "; B1.printLengths();
  std::cout << "Length box2: "; B2.printLengths();
}


/* --- Print verbose information relative to the collision --- */
void CollisionTwoBoxes::
printBBcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Boxes are intersecting" << std::endl;
    for (int i=0; i<pointsBB.size(); i++){
      std::cout << "  (" << pointsBB[i].transpose() << ")  ";
    }
  }
  else{
    std::cout << "Result: Boxes are NOT intersecting" << std::endl;
  }
}


/* --- Print vertices of the boxes  --- */
void CollisionTwoBoxes::
printBoxesVertices( void )
{
  std::cout << "Vertices Box 1: "; B1.printVertices();
  std::cout << "Vertices Box 2: "; B2.printVertices();

}


/* --- Get Vertices of the boxes --- */
// void CollisionTwoBoxes::
// getVerticesB1(double V1[3], double V2[3], double V3[3], double V4[3],
// 	      double V5[3], double V6[3], double V7[3], double V8[3])
// {
//   set(V1, Pw11); set(V2, Pw12); set(V3, Pw13); set(V4, Pw14);
//   set(V5, Pw15); set(V6, Pw16); set(V7, Pw17); set(V8, Pw18);
// }


// void CollisionTwoBoxes::
// getVerticesB2(double V1[3], double V2[3], double V3[3], double V4[3],
// 	      double V5[3], double V6[3], double V7[3], double V8[3])
// {
//   set(V1, Pw21); set(V2, Pw22); set(V3, Pw23); set(V4, Pw24);
//   set(V5, Pw25); set(V6, Pw26); set(V7, Pw27); set(V8, Pw28);
// }
