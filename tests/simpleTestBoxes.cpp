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

#include <collision-one-many-boxes.h>
#include <iostream>
#include <math.h>

#include <Eigen/Dense>
using namespace Eigen;

void rotMatrix(Matrix3d &R, Vector3d axis, double angle);

int main( void )
{
  int collisionDetected;
  //double R1[9], R2[9], T1[3], T2[3], axis[3];
  Matrix3d R1, R2;
  Vector3d T1, T2, axis;

  CollisionOneManyBoxes scene;
  //scene.setTolerance(0.001);

  // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // T1[0]=2.0; T1[1]=1.0; T1[2]=0.5;
  // rotMatrix(R1, axis, 0.0);
  // scene.setLengthB1(4.0, 2.0, 1.0);
  // scene.setTransformationB1(R1, T1);

  // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // //T2[0]=2.0; T2[1]=1.0; T2[2]=2.5;
  // T2[0]=4.0; T2[1]=2.0; T2[2]=1.5;
  // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.4999;
  // rotMatrix(R2, axis, 0.0);
  // scene.addBox(4.0, 2.0, 1.0, R2, T2);

  // // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // // T2[0]=5.0; T2[1]=2.0; T2[2]=1.5;
  // // rotMatrix(R2, axis, 0.0);
  // // scene.addBox(4.0, 2.0, 1.0, R2, T2);

  // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // //T2[0]=-1; T2[1]=2.0; T2[2]=1.5;
  // T2[0]=-1.5; T2[1]=1.5; T2[2]=1.5;
  // rotMatrix(R2, axis, 0.0);
  // scene.addBox(4.0, 2.0, 1.0, R2, T2);


  // R1[0]=0.99986589763917377;    R1[1]=-0.00067914783993300103; R1[2]=0.016362319408345875;
  // R1[3]=0.00075979084842649745; R1[4]=0.99998759396560366;     R1[5]=-0.004922868345;
  // R1[6]=-0.016358773061445196;  R1[7]=0.0049346401180313396;   R1[8]=0.99985400927876944;
  // T1[0]=0.020394159554286626;   T1[1]=-0.10441503659450251;    T1[2]=0.087178068672960629;
  // scene.setLengthB1(0.23, 0.135, 0.001);
  // scene.setTransformationB1(R1, T1);
  // R2[0]=1.0; R2[1]=0.0; R2[2]=0.0;
  // R2[3]=0.0; R2[4]=1.0; R2[5]=0.0;
  // R2[6]=0.0; R2[7]=0.0; R2[8]=1.0;
  // T2[0]=0.03;   T2[1]=-0.13;    T2[2]=0.05;
  // scene.addBox(0.19, 0.115, 0.10, R2, T2);

  // R1[0]=; R1[1]=; R1[2]=;
  // R1[3]=; R1[4]=; R1[5]=;
  // R1[6]=; R1[7]=; R1[8]=;
  // T1[0]=; T1[1]=; T1[2]=;
  // scene.setLengthB1(0.23, 0.135, 0.001);
  // scene.setTransformationB1(R1, T1);
  // R2[0]=; R2[1]=; R2[2]=;
  // R2[3]=; R2[4]=; R2[5]=;
  // R2[6]=; R2[7]=; R2[8]=;
  // T2[0]=; T2[1]=; T2[2]=;
  // scene.addBox(0.19, 0.115, 0.10, R2, T2);


  // R2[0]=1.0; R2[1]=0.0; R2[2]=0.0;
  // R2[3]=0.0; R2[4]=1.0; R2[5]=0.0;
  // R2[6]=0.0; R2[7]=0.0; R2[8]=1.0;
  // T2[0]=0.03;   T2[1]=-0.13;    T2[2]=0.05;
  // scene.addBox(0.19, 0.115, 0.10, R2, T2);

  // R1[0]=0.93040071; R1[1]= 0.09043101; R1[2]=0.35521366;
  // R1[3]=-0.00390309; R1[4]=0.971478; R1[5]=-0.23709758;
  // R1[6]=-0.36652323; R1[7]=0.21920933; R1[8]=0.90421678;
  // T1[0]=-0.02154648; T1[1]=-0.0820861; T1[2]=0.14945323;
  // scene.setLengthB1(0.23, 0.135, 0.001);
  // scene.setTransformationB1(R1, T1);


  // scene.setTolerance(0.01);

  // R2[0]=1.0; R2[1]=0.0; R2[2]=0.0;
  // R2[3]=0.0; R2[4]=1.0; R2[5]=0.0;
  // R2[6]=0.0; R2[7]=0.0; R2[8]=1.0;
  // T2[0]=0.03;   T2[1]=-0.13;    T2[2]=0.05;
  // scene.addBox(0.19, 0.115, 0.10, R2, T2);

  // R1[0]=0.9999860339278217 ; R1[1]=-0.00083655099522303; R1[2]=0.0052184415047020994;
  // R1[3]=0.00083625863550401463; R1[4]=0.99999964864171165; R1[5]=5.8206082807325495e-05;
  // R1[6]=-0.0052184883635159962; R1[7]=-5.3841303124812213e-05; R1[8]=0.99998638214743407;
  // T1[0]=-0.016225030166884216; T1[1]=-0.10562326452894533; T1[2]=0.10059533757257406;
  // // R1[0]=1.0; R1[1]=0.0; R1[2]=0.0;
  // // R1[3]=0.0; R1[4]=1.0; R1[5]=0.0;
  // // R1[6]=0.0; R1[7]=0.0; R1[8]=1.0;
  // // T1[0]=-0.016; T1[1]=-0.105; T1[2]=0.095; //0.1004;
  // scene.setLengthB1(0.23, 0.135, 0.00001);
  // //scene.setLengthB1(0.23, 0.135, 0.02);
  // scene.setTransformationB1(R1, T1);


  scene.setTolerance(0.005);

  R2 << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
  T2 << 0.03, -0.13, 0.03;
  scene.addBox(0.19, 0.115, 0.10, R2, T2);


  // R1 << 0.99999999999999978, 3.738479408160609e-09, 2.7267014987878495e-08,
  //   -3.7384784942671144e-09, 0.99999999999999944, -3.3435310033084201e-08,
  //   -2.7267015140377106e-08, 3.3435309937393195e-08, 0.99999999999999922;
  // T1 << -0.020860537574661286, -0.071387244448416889, 0.09200000535012888;

  // R1 << 0.99999971973071911, 0.0001195694301036515, 0.00073908161579726815,
  //   -0.00011892518681012387, 0.99999961302983009, -0.00087166334682232146,
  //   -0.00073918555408437696, 0.00087157520710265065, 0.99999934698047421;
  // T1 << -0.021355335987184074, -0.071746781540336341, 0.092457378181053249;

  R1 << 0.99999970911564173, 0.00013513698057147076, 0.00075067078594650367,
    -0.00013447496703785513, 0.99999960210513528, -0.00088187643953444104,
    -0.00075078966137770561, 0.00088177523658116157, 0.99999932939343339;
  T1 << -0.021355971343321688, -0.07174552364551072, 0.0924592662811861;

  scene.setLengthB1(0.23, 0.135, 0.024);
  scene.setTransformationB1(R1, T1);


  collisionDetected = scene.computeIntersections();
  scene.printLengths();
  scene.printBoxesVertices();
  scene.printBBsCollisionInformation();

  return 0;
}




void rotMatrix(Matrix3d &R, Vector3d axis, double angle)
{
  double ca, sa, va, norm_axis;

  ca=cos(angle); sa=sin(angle); va=1-cos(angle);
  norm_axis = sqrt(axis(0)*axis(0)+axis(1)*axis(1)+axis(2)*axis(2));
  for (int i=0; i<3; i++)
    axis(i) = axis(i)/norm_axis;

  R(0,0) = 1 - va*(axis(1)*axis(1)+axis(2)*axis(2));
  R(0,1) = axis(0)*axis(1)*va - axis(2)*sa;
  R(0,2) = axis(0)*axis(2)*va + axis(1)*sa;
  R(1,0) = axis(0)*axis(1)*va + axis(2)*sa;
  R(1,1) = 1 - va*(axis(0)*axis(0)+axis(2)*axis(2));
  R(1,2) = axis(1)*axis(2)*va - axis(0)*sa;
  R(2,0) = axis(0)*axis(2)*va - axis(1)*sa;
  R(2,1) = axis(1)*axis(2)*va + axis(0)*sa;
  R(2,2) = 1 - va*(axis(0)*axis(0)+axis(1)*axis(1));

  if (false) {
    std::cout << R << std::endl;
  }
}


  // if (collisionDetected){
  //   std::cout << "Boxes are intersecting" << std::endl;
  //   std::cout << "Intersections occur at: " << std::endl;

  //   for (int i=0; i<scene.pointsBBs.size(); i++){
  //     std::cout << "( " << scene.pointsBBs[i].x << ", " << scene.pointsBBs[i].y
  //     		<< ", " << scene.pointsBBs[i].z << ")" << std::endl;
  //   }
  // }
  // else{
  //   std::cout << "Boxes are not intersecting" << std::endl;
  // }
