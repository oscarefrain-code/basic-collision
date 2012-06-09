/***********************************************************************************
 *
 * Basic components and class to be used for the other classes
 *                                   Oscar E. Ramos Ponce
 *
 ********************************************************************************** */


#ifndef __COLLISION_HELPER_H__
#define __COLLISION_HELPER_H__

#include <math.h>
#include <vector>
#include <iostream>


/* Definition of a 'minimum value' to be used as a threshold */
#define TRI_EPSILON 0.000001
//#define TRI_EPSILON 0.00000001


struct Point3d
{
  double x;
  double y;
  double z;

  double operator[](int i){
    if (i==0) return x;
    if (i==1) return y;
    if (i==2) return z;
  }
  void print(){
    std::cout << " ( " << x << ", " << y << ", " << z << ") ";
  }
};


class OperationsHelper
{
 
 public:
  /* --- Basic Vector and Scalar operations --- */

  void add(double result[3], double v1[3], double v2[3]);
  void sub(double result[3], double v1[3], double v2[3]);
  void mult(double result[3], double v1[3], double factor);

  double dot(double v1[3], double v2[3]);
  void cross(double result[3], double v1[3], double v2[3]);

  void sort(double &a, double &b);
  void sort2(double &a, double &b, int &smallest);

  void printVector(double v[3]);
  void printVector(double v1, double v2, double v3);

  void set(double result[3], double source[3]);
  void set(double v[3], double v0, double v1, double v2);
 
  void pushVectorAsPoint(double v[3], std::vector<Point3d>& points);
  void prunePoints(std::vector<Point3d>& points);
  int convexHull( std::vector<Point3d>& points );


};


#endif
