/***********************************************************************************
 * Hexahedron/hexahedron intersection test:
 *   Each hexahedron is decomposed into a set of rectangles and the collision
 *   between the rectangles is tested.
 *  
 *                                           Oscar E. Ramos Ponce, LAAS-CNRS, France
 ************************************************************************************/


#include <collision-two-boxes.h>
#include <iostream>


/* --- Constructors --- */
CollisionTwoBoxes::
CollisionTwoBoxes( void )
{
  setLengthB1(0.0, 0.0, 0.0);
  setLengthB2(0.0, 0.0, 0.0);
  defaultInit();
}


CollisionTwoBoxes::
CollisionTwoBoxes(double lx1, double ly1, double lz1,
		  double lx2, double ly2, double lz2)
{
  setLengthB1(lx1, ly1, lz1);
  setLengthB2(lx2, ly2, lz2);
  defaultInit();
}


void CollisionTwoBoxes::
defaultInit( void )
{
  double zeros[3] = {0.0, 0.0, 0.0};
  set(Pw11, zeros); set(Pw12, zeros); set(Pw13, zeros); set(Pw14, zeros);
  set(Pw15, zeros); set(Pw16, zeros); set(Pw17, zeros); set(Pw18, zeros);
  set(Pw21, zeros); set(Pw22, zeros); set(Pw23, zeros); set(Pw24, zeros);
  set(Pw25, zeros); set(Pw26, zeros); set(Pw27, zeros); set(Pw28, zeros);
  collisionIndicator = 0;
  tolerance = TRI_EPSILON;
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
  Lx1=lx1; Ly1=ly1; Lz1=lz1;
  // Vertices of the box with respect to its COM
  set(Pc11,  Lx1/2, -Ly1/2,  Lz1/2);
  set(Pc12, -Lx1/2, -Ly1/2,  Lz1/2);
  set(Pc13, -Lx1/2, -Ly1/2, -Lz1/2);
  set(Pc14,  Lx1/2, -Ly1/2, -Lz1/2);
  set(Pc15,  Lx1/2,  Ly1/2,  Lz1/2);
  set(Pc16, -Lx1/2,  Ly1/2,  Lz1/2);
  set(Pc17, -Lx1/2,  Ly1/2, -Lz1/2);
  set(Pc18,  Lx1/2,  Ly1/2, -Lz1/2);
}


void CollisionTwoBoxes::
setLengthB2(double lx2, double ly2, double lz2)
{
  Lx2=lx2; Ly2=ly2; Lz2=lz2;
  // Vertices of the box with respect to its CoM
  set(Pc21,  Lx2/2, -Ly2/2,  Lz2/2);
  set(Pc22, -Lx2/2, -Ly2/2,  Lz2/2);
  set(Pc23, -Lx2/2, -Ly2/2, -Lz2/2);
  set(Pc24,  Lx2/2, -Ly2/2, -Lz2/2);
  set(Pc25,  Lx2/2,  Ly2/2,  Lz2/2);
  set(Pc26, -Lx2/2,  Ly2/2,  Lz2/2);
  set(Pc27, -Lx2/2,  Ly2/2, -Lz2/2);
  set(Pc28,  Lx2/2,  Ly2/2, -Lz2/2);
}

void CollisionTwoBoxes::
getLengthB1(double & lx1, double & ly1, double & lz1)
{
  lx1=Lx1; ly1=Ly1; lz1=Lz1;
}

void CollisionTwoBoxes::
getLengthB2(double & lx2, double & ly2, double & lz2)
{
  lx2=Lx2; ly2=Ly2; lz2=Lz2;
}


/* --- Set the transformations for the boxes --- */ 
void CollisionTwoBoxes::
setTransformation1(double R1[9], double T1[3])
{
  // Vertices of the first box with respect to the world
  // Pw: wrt world, Pc: wrt CoM
  applyRT(Pw11, Pc11, R1, T1);
  applyRT(Pw12, Pc12, R1, T1);
  applyRT(Pw13, Pc13, R1, T1);
  applyRT(Pw14, Pc14, R1, T1);
  applyRT(Pw15, Pc15, R1, T1);
  applyRT(Pw16, Pc16, R1, T1);
  applyRT(Pw17, Pc17, R1, T1);
  applyRT(Pw18, Pc18, R1, T1);

  // Create the Rectangles for the first box
  To1_To2.setVerticesR1(Pw11, Pw12, Pw16, Pw15);
  To1_Bo2.setVerticesR1(Pw11, Pw12, Pw16, Pw15);
  To1_Le2.setVerticesR1(Pw11, Pw12, Pw16, Pw15);
  To1_Ri2.setVerticesR1(Pw11, Pw12, Pw16, Pw15);
  To1_Fr2.setVerticesR1(Pw11, Pw12, Pw16, Pw15);
  To1_Ba2.setVerticesR1(Pw11, Pw12, Pw16, Pw15);

  Bo1_To2.setVerticesR1(Pw14, Pw13, Pw17, Pw18);
  Bo1_Bo2.setVerticesR1(Pw14, Pw13, Pw17, Pw18);
  Bo1_Le2.setVerticesR1(Pw14, Pw13, Pw17, Pw18);
  Bo1_Ri2.setVerticesR1(Pw14, Pw13, Pw17, Pw18);
  Bo1_Fr2.setVerticesR1(Pw14, Pw13, Pw17, Pw18);
  Bo1_Ba2.setVerticesR1(Pw14, Pw13, Pw17, Pw18);

  Le1_To2.setVerticesR1(Pw11, Pw12, Pw13, Pw14);
  Le1_Bo2.setVerticesR1(Pw11, Pw12, Pw13, Pw14);
  Le1_Le2.setVerticesR1(Pw11, Pw12, Pw13, Pw14);
  Le1_Ri2.setVerticesR1(Pw11, Pw12, Pw13, Pw14);
  Le1_Fr2.setVerticesR1(Pw11, Pw12, Pw13, Pw14);
  Le1_Ba2.setVerticesR1(Pw11, Pw12, Pw13, Pw14);

  Ri1_To2.setVerticesR1(Pw15, Pw16, Pw17, Pw18);
  Ri1_Bo2.setVerticesR1(Pw15, Pw16, Pw17, Pw18);
  Ri1_Le2.setVerticesR1(Pw15, Pw16, Pw17, Pw18);
  Ri1_Ri2.setVerticesR1(Pw15, Pw16, Pw17, Pw18);
  Ri1_Fr2.setVerticesR1(Pw15, Pw16, Pw17, Pw18);
  Ri1_Ba2.setVerticesR1(Pw15, Pw16, Pw17, Pw18);

  Fr1_To2.setVerticesR1(Pw11, Pw14, Pw18, Pw15);
  Fr1_Bo2.setVerticesR1(Pw11, Pw14, Pw18, Pw15);
  Fr1_Le2.setVerticesR1(Pw11, Pw14, Pw18, Pw15);
  Fr1_Ri2.setVerticesR1(Pw11, Pw14, Pw18, Pw15);
  Fr1_Fr2.setVerticesR1(Pw11, Pw14, Pw18, Pw15);
  Fr1_Ba2.setVerticesR1(Pw11, Pw14, Pw18, Pw15);

  Ba1_To2.setVerticesR1(Pw12, Pw13, Pw17, Pw16);
  Ba1_Bo2.setVerticesR1(Pw12, Pw13, Pw17, Pw16);
  Ba1_Le2.setVerticesR1(Pw12, Pw13, Pw17, Pw16);
  Ba1_Ri2.setVerticesR1(Pw12, Pw13, Pw17, Pw16);
  Ba1_Fr2.setVerticesR1(Pw12, Pw13, Pw17, Pw16);
  Ba1_Ba2.setVerticesR1(Pw12, Pw13, Pw17, Pw16);
}

void CollisionTwoBoxes::
setTransformation2(double R2[9], double T2[3])
{
  // Vertices of the second box with respect to the world
  // Pw: wrt world, Pc: wrt CoM
  applyRT(Pw21, Pc21, R2, T2);
  applyRT(Pw22, Pc22, R2, T2);
  applyRT(Pw23, Pc23, R2, T2);
  applyRT(Pw24, Pc24, R2, T2);
  applyRT(Pw25, Pc25, R2, T2);
  applyRT(Pw26, Pc26, R2, T2);
  applyRT(Pw27, Pc27, R2, T2);
  applyRT(Pw28, Pc28, R2, T2);

  // std::cout<<std::endl;
  // printVector(Pw21);printVector(Pw22);printVector(Pw26);printVector(Pw25); std::cout<<std::endl;
  // printVector(Pw24);printVector(Pw23);printVector(Pw27);printVector(Pw28); std::cout<<std::endl;
  // printVector(Pw21);printVector(Pw22);printVector(Pw23);printVector(Pw24); std::cout<<std::endl;
  // printVector(Pw25);printVector(Pw26);printVector(Pw27);printVector(Pw28); std::cout<<std::endl;
  // printVector(Pw21);printVector(Pw24);printVector(Pw28);printVector(Pw25); std::cout<<std::endl;
  // printVector(Pw22);printVector(Pw23);printVector(Pw27);printVector(Pw26); std::cout<<std::endl;

  // Create the Rectangles for the second box
  To1_To2.setVerticesR2(Pw21, Pw22, Pw26, Pw25);
  To1_Bo2.setVerticesR2(Pw24, Pw23, Pw27, Pw28);
  To1_Le2.setVerticesR2(Pw21, Pw22, Pw23, Pw24);
  To1_Ri2.setVerticesR2(Pw25, Pw26, Pw27, Pw28);
  To1_Fr2.setVerticesR2(Pw21, Pw24, Pw28, Pw25);
  To1_Ba2.setVerticesR2(Pw22, Pw23, Pw27, Pw26);

  Bo1_To2.setVerticesR2(Pw21, Pw22, Pw26, Pw25);
  Bo1_Bo2.setVerticesR2(Pw24, Pw23, Pw27, Pw28);
  Bo1_Le2.setVerticesR2(Pw21, Pw22, Pw23, Pw24);
  Bo1_Ri2.setVerticesR2(Pw25, Pw26, Pw27, Pw28);
  Bo1_Fr2.setVerticesR2(Pw21, Pw24, Pw28, Pw25);
  Bo1_Ba2.setVerticesR2(Pw22, Pw23, Pw27, Pw26);

  Le1_To2.setVerticesR2(Pw21, Pw22, Pw26, Pw25);
  Le1_Bo2.setVerticesR2(Pw24, Pw23, Pw27, Pw28);
  Le1_Le2.setVerticesR2(Pw21, Pw22, Pw23, Pw24);
  Le1_Ri2.setVerticesR2(Pw25, Pw26, Pw27, Pw28);
  Le1_Fr2.setVerticesR2(Pw21, Pw24, Pw28, Pw25);
  Le1_Ba2.setVerticesR2(Pw22, Pw23, Pw27, Pw26);

  Ri1_To2.setVerticesR2(Pw21, Pw22, Pw26, Pw25);
  Ri1_Bo2.setVerticesR2(Pw24, Pw23, Pw27, Pw28);
  Ri1_Le2.setVerticesR2(Pw21, Pw22, Pw23, Pw24);
  Ri1_Ri2.setVerticesR2(Pw25, Pw26, Pw27, Pw28);
  Ri1_Fr2.setVerticesR2(Pw21, Pw24, Pw28, Pw25);
  Ri1_Ba2.setVerticesR2(Pw22, Pw23, Pw27, Pw26);

  Fr1_To2.setVerticesR2(Pw21, Pw22, Pw26, Pw25);
  Fr1_Bo2.setVerticesR2(Pw24, Pw23, Pw27, Pw28);
  Fr1_Le2.setVerticesR2(Pw21, Pw22, Pw23, Pw24);
  Fr1_Ri2.setVerticesR2(Pw25, Pw26, Pw27, Pw28);
  Fr1_Fr2.setVerticesR2(Pw21, Pw24, Pw28, Pw25);
  Fr1_Ba2.setVerticesR2(Pw22, Pw23, Pw27, Pw26);

  Ba1_To2.setVerticesR2(Pw21, Pw22, Pw26, Pw25);
  Ba1_Bo2.setVerticesR2(Pw24, Pw23, Pw27, Pw28);
  Ba1_Le2.setVerticesR2(Pw21, Pw22, Pw23, Pw24);
  Ba1_Ri2.setVerticesR2(Pw25, Pw26, Pw27, Pw28);
  Ba1_Fr2.setVerticesR2(Pw21, Pw24, Pw28, Pw25);
  Ba1_Ba2.setVerticesR2(Pw22, Pw23, Pw27, Pw26);
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

  std::cout << "Collisions : ";
  for (int i=0; i<36; i++)
    std::cout << collisions[i] << " ";
  std::cout << std::endl;
  
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

  for (int i=0; i<Fr1_To2.pointsRR.size(); i++)
    Fr1_To2.pointsRR[i].print();
  std::cout << std::endl;

      
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

/* --- Transformation of a Point --- */
void CollisionTwoBoxes::
applyRT(double pout[3], double pin[3], double R[9], double T[3])
{
  pout[0] = R[0]*pin[0] + R[1]*pin[1] + R[2]*pin[2] + T[0];
  pout[1] = R[3]*pin[0] + R[4]*pin[1] + R[5]*pin[2] + T[1];
  pout[2] = R[6]*pin[0] + R[7]*pin[1] + R[8]*pin[2] + T[2];
}


/* --- Get the average of very close points (to a certain tolerance) --- */
void CollisionTwoBoxes::
averageSimilarValues( void )
{
  double dist;
  for (int i=0; i<pointsBB.size(); i++){
    for (int j=i+1; j<pointsBB.size(); j++){
      dist = (pointsBB[i].x-pointsBB[j].x)*(pointsBB[i].x-pointsBB[j].x) +
	     (pointsBB[i].y-pointsBB[j].y)*(pointsBB[i].y-pointsBB[j].y) +
             (pointsBB[i].z-pointsBB[j].z)*(pointsBB[i].z-pointsBB[j].z);

      if ( dist < (tolerance*tolerance) ) {
	pointsBB[i].x = (pointsBB[i].x+pointsBB[j].x)/2;
	pointsBB[i].y = (pointsBB[i].y+pointsBB[j].y)/2;
	pointsBB[i].z = (pointsBB[i].z+pointsBB[j].z)/2;
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
  std::cout << "Lengths: Box 1 = (" << Lx1 << ", " << Ly1 << ", " << Lz1 << " ),   "
	    << "Box 2 = (" << Lx2 << ", " << Ly2 << ", " << Lz2 << " )" << std::endl;
}


/* --- Print verbose information relative to the collision --- */
void CollisionTwoBoxes::
printBBcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Boxes are intersecting" << std::endl;
    std::cout << " - Intersections occur at: " << std::endl;
    for (int i=0; i<pointsBB.size(); i++){
      std::cout << "     "; pointsBB[i].print();
      std::cout << std::endl;
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
  std::cout << "Vertices: ";
  std::cout << "\n  - Box 1: "; printVector(Pw11); printVector(Pw12); printVector(Pw13); printVector(Pw14);
  std::cout << std::endl
	    << "           "; printVector(Pw15); printVector(Pw16); printVector(Pw17); printVector(Pw18);
  std::cout << "\n  - Box 2: "; printVector(Pw21); printVector(Pw22); printVector(Pw23); printVector(Pw24);
  std::cout << std::endl
	    << "           "; printVector(Pw25); printVector(Pw26); printVector(Pw27); printVector(Pw28);
  std::cout << std::endl;

}


/* --- Get Vertices of the boxes --- */
void CollisionTwoBoxes::
getVerticesB1(double V1[3], double V2[3], double V3[3], double V4[3],
	      double V5[3], double V6[3], double V7[3], double V8[3])
{
  set(V1, Pw11); set(V2, Pw12); set(V3, Pw13); set(V4, Pw14);
  set(V5, Pw15); set(V6, Pw16); set(V7, Pw17); set(V8, Pw18);
}


void CollisionTwoBoxes::
getVerticesB2(double V1[3], double V2[3], double V3[3], double V4[3],
	      double V5[3], double V6[3], double V7[3], double V8[3])
{
  set(V1, Pw21); set(V2, Pw22); set(V3, Pw23); set(V4, Pw24);
  set(V5, Pw25); set(V6, Pw26); set(V7, Pw27); set(V8, Pw28);
}
