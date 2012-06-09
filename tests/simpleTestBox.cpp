/***********************************************************************************
 * Tests for Box/Box class
 *                                           Oscar E. Ramos Ponce, LAAS-CNRS, France
 ************************************************************************************/

#include <collision-two-boxes.h>
#include <iostream>
#include <math.h>

void rotMatrix(double R[9], double axis[3], double ang);


int main( void )
{
  int collisionDetected;
  double R1[9], R2[9], T1[3], T2[3], axis[3];

  //CollisionTwoBoxes scene(4.0, 2.0, 1.0, 4.0, 2.0, 1.0);
  //scene.setTolerance(0.001);

  // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // T1[0]=2.0; T1[1]=1.0; T1[2]=0.5;
  // //T1[0]=-1.0; T1[1]=2.0; T1[2]=1.5;
  // rotMatrix(R1, axis, 0.0);
  // scene.setTransformation1(R1, T1);

  // // Strange if T2[0]=-1 !!!! CHECK IT OUT !!!!x
  // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // //T2[0]=2.0; T2[1]=1.0; T2[2]=2.5;
  // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.5;
  // //T2[0]=-1.0; T2[1]=2.0; T2[2]=1.5;   
  // T2[0]=-0.5; T2[1]=2.0; T2[2]=1.5;   
  // //T2[0]=2.0; T2[1]=1.0; T2[2]=0.5;
  // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.4999;
  // rotMatrix(R2, axis, 0.0);
  // scene.setTransformation2(R2, T2);

  CollisionTwoBoxes scene;
  scene.setTolerance(0.01);

  R1[0]=0.9999860339278217 ; R1[1]=-0.00083655099522303; R1[2]=0.0052184415047020994;
  R1[3]=0.00083625863550401463; R1[4]=0.99999964864171165; R1[5]=5.8206082807325495e-05;
  R1[6]=-0.0052184883635159962; R1[7]=-5.3841303124812213e-05; R1[8]=0.99998638214743407;
  T1[0]=-0.016225030166884216; T1[1]=-0.10562326452894533; T1[2]=0.10059533757257406;
  scene.setLengthB1(0.23, 0.135, 0.00001);
  scene.setTransformation1(R1, T1);

  R2[0]=1.0; R2[1]=0.0; R2[2]=0.0;
  R2[3]=0.0; R2[4]=1.0; R2[5]=0.0;
  R2[6]=0.0; R2[7]=0.0; R2[8]=1.0;
  T2[0]=0.03;   T2[1]=-0.13;    T2[2]=0.05;
  scene.setLengthB2(0.19, 0.115, 0.10);
  scene.setTransformation2(R2, T2);

  collisionDetected = scene.computeBBintersections();
  scene.printBoxesVertices();
  scene.printBBcollisionInformation();

  return 0;
}




void rotMatrix(double R[9], double axis[3], double angle)
{
  double ca, sa, va, norm_axis;

  ca=cos(angle); sa=sin(angle); va=1-cos(angle);
  norm_axis = sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
  for (int i=0; i<3; i++)
    axis[i] = axis[i]/norm_axis;

  R[0] = 1 - va*(axis[1]*axis[1]+axis[2]*axis[2]);
  R[1] = axis[0]*axis[1]*va - axis[2]*sa;
  R[2] = axis[0]*axis[2]*va + axis[1]*sa;
  R[3] = axis[0]*axis[1]*va + axis[2]*sa;
  R[4] = 1 - va*(axis[0]*axis[0]+axis[2]*axis[2]);
  R[5] = axis[1]*axis[2]*va - axis[0]*sa;
  R[6] = axis[0]*axis[2]*va - axis[1]*sa;
  R[7] = axis[1]*axis[2]*va + axis[0]*sa;
  R[8] = 1 - va*(axis[0]*axis[0]+axis[1]*axis[1]);

  if (false) {
    std::cout << R[0] << " " << R[1] << " " << R[2] << std::endl
	      << R[3] << " " << R[4] << " " << R[5] << std::endl
	      << R[6] << " " << R[7] << " " << R[8] << std::endl;
  }
}

