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
 * Tests for Box/Box class
 *
 */

#include <collision-two-boxes.h>
#include <iostream>
#include <math.h>

#include <Eigen/Dense>
using namespace Eigen;

Matrix3d rotMatrix(Vector3d axis, double angle);


int main( void )
{
  int collisionDetected;
  Matrix3d R1, R2;
  Vector3d T1, T2, axis;

  // CollisionTwoBoxes scene(4.0, 2.0, 1.0, 4.0, 2.0, 1.0);
  // scene.setTolerance(0.001);

  // axis << 1.0, 0.0, 0.0;
  // T1 << 2.0, 1.0, 0.5;
  // //T1[0]=-1.0; T1[1]=2.0; T1[2]=1.5;
  // rotMatrix(R1, axis, 0.0);
  // scene.setTransformation1(R1, T1);

  // // Strange if T2[0]=-1 !!!! CHECK IT OUT !!!!
  // axis << 1.0, 0.0, 0.0;
  
  // //T2[0]=2.0; T2[1]=1.0; T2[2]=2.5;
  // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.5;
  // //T2[0]=-1.0; T2[1]=2.0; T2[2]=1.5;   
  // T2 << -0.5, 2.0, 1.5;   
  // //T2[0]=2.0; T2[1]=1.0; T2[2]=0.5;
  // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.4999;
  // rotMatrix(R2, axis, 0.0);
  // scene.setTransformation2(R2, T2);


  // //CollisionTwoBoxes scene(4.0, 2.0, 1.0, 4.0, 2.0, 1.0);
  // //scene.setTolerance(0.001);

  // // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // // T1[0]=2.0; T1[1]=1.0; T1[2]=0.5;
  // // //T1[0]=-1.0; T1[1]=2.0; T1[2]=1.5;
  // // rotMatrix(R1, axis, 0.0);
  // // scene.setTransformation1(R1, T1);

  // // // Strange if T2[0]=-1 !!!! CHECK IT OUT !!!!x
  // // axis[0]=1.0; axis[1]=0.0; axis[2]=0.0;
  // // //T2[0]=2.0; T2[1]=1.0; T2[2]=2.5;
  // // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.5;
  // // //T2[0]=-1.0; T2[1]=2.0; T2[2]=1.5;   
  // // T2[0]=-0.5; T2[1]=2.0; T2[2]=1.5;   
  // // //T2[0]=2.0; T2[1]=1.0; T2[2]=0.5;
  // // //T2[0]=4.0; T2[1]=2.0; T2[2]=1.4999;
  // // rotMatrix(R2, axis, 0.0);
  // // scene.setTransformation2(R2, T2);


  // // CollisionTwoBoxes scene;
  // // scene.setTolerance(0.01);

  // // R1[0]=0.9999860339278217 ; R1[1]=-0.00083655099522303; R1[2]=0.0052184415047020994;
  // // R1[3]=0.00083625863550401463; R1[4]=0.99999964864171165; R1[5]=5.8206082807325495e-05;
  // // R1[6]=-0.0052184883635159962; R1[7]=-5.3841303124812213e-05; R1[8]=0.99998638214743407;
  // // T1[0]=-0.016225030166884216; T1[1]=-0.10562326452894533; T1[2]=0.10059533757257406;
  // // scene.setLengthB1(0.23, 0.135, 0.00001);
  // // scene.setTransformation1(R1, T1);

  // // R2[0]=1.0; R2[1]=0.0; R2[2]=0.0;
  // // R2[3]=0.0; R2[4]=1.0; R2[5]=0.0;
  // // R2[6]=0.0; R2[7]=0.0; R2[8]=1.0;
  // // T2[0]=0.03;   T2[1]=-0.13;    T2[2]=0.05;
  // // scene.setLengthB2(0.19, 0.115, 0.10);
  // // scene.setTransformation2(R2, T2);

  CollisionTwoBoxes scene(0.23, 0.135, 0.024, 0.19, 0.115, 0.10);
  scene.setTolerance(0.005);

  R2 << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
  T2 << 0.03, -0.13, 0.03;

  // R1 << 0.99999999999999978, 3.738479408160609e-09, 2.7267014987878495e-08,
  //   -3.7384784942671144e-09, 0.99999999999999944, -3.3435310033084201e-08,
  //   -2.7267015140377106e-08, 3.3435309937393195e-08, 0.99999999999999922;
  // T1 << -0.020860537574661286, -0.071387244448416889, 0.09200000535012888;

  // R1 << 0.99999971973071911, 0.0001195694301036515, 0.00073908161579726815,
  //   -0.00011892518681012387, 0.99999961302983009, -0.00087166334682232146,
  //   -0.00073918555408437696, 0.00087157520710265065, 0.99999934698047421;
  // T1 << -0.021355335987184074, -0.071746781540336341, 0.092457378181053249;

  R1 << 0.99999971973071911, 0.0001195694301036515, 0.00073908161579726815,
    -0.00011892518681012387, 0.99999961302983009, -0.00087166334682232146,
    -0.00073918555408437696, 0.00087157520710265065, 0.99999934698047421;
  T1 << -0.021355335987184074, -0.071746781540336341, 0.092457378181053249;

  scene.setTransformation1(R1, T1);
  scene.setTransformation2(R2, T2);

  collisionDetected = scene.computeBBintersections();
  scene.printBoxesVertices();
  // scene.printLengths();
  // scene.printTransformation();

  scene.printBBcollisionInformation();
  // std::cout << std::endl;

  return 0;
}


Matrix3d rotMatrix(Vector3d axis, double angle)
{
  double ca, sa, va, norm_axis;
  Matrix3d R;

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

  return R;
}


