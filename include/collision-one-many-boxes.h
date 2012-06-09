/***********************************************************************************
 * One Hexahedron/Many hexahedra intersection test:
 *   The collision between one hexahedron "B1" and many hexahedra "Bi" is tested. It
 *   decomposes the problem in many hexahedron/hexahedron problems. The hexahedra 
 *   Bi is supposed not to move after fixing its initial position, but B1 can move.
 *  
 *                                           Oscar E. Ramos Ponce, LAAS-CNRS, France
 ************************************************************************************/

#ifndef __COLLISION_ONE_MANY_BOXES_H__
#define __COLLISION_ONE_MANY_BOXES_H__

#include <math.h>
#include <vector>
#include <collision-helper.h>
#include <collision-two-boxes.h>


class CollisionOneManyBoxes:
public OperationsHelper
{
 public:
  /* -- Constructors -- */
  CollisionOneManyBoxes( void );

  /* -- Set values -- */
  void setTolerance(double toler);
  void setLengthB1(double lx1, double ly1, double lz1);
  void setTransformationB1(double R[9], double T[3]);
  int addBox(double lx, double ly, double lz, double R[9], double T[3]);

  /* -- Compute the collision detection -- */
  int computeIntersections( void );

  /* -- Get some values -- */
  double getTolerance( void );
  void getLengthB1(double & lx1, double & ly1, double & lz1);

  /* --- Print Different Information --- */
  void printBBsCollisionInformation( void );
  void printLengths( void );
  void printBoxesVertices( void );

  /* -- Collision Points -- */
  std::vector<Point3d> pointsBBs;
  std::vector<CollisionTwoBoxes> collisionBoxes;

 private:
  double Lx1, Ly1, Lz1;
  double R1[9], T1[3];
  double tolerance;
  int id, collisionIndicator;

};



#endif
