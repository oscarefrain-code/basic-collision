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
 * Basic components and class to be used for the other classes
 *
 */


#include <operations-helper.h>


/* ---  Sort a,b so that a < = b  --- */
/* smallest: initial smallest value  */
void OperationsHelper::
sort2(Vector2d &vect, int &smallest)
{
  if(vect(0)>vect(1)) {
    float c;
    c=vect(0); vect(0)=vect(1); vect(1)=c;
    smallest=1;      
  }
  else smallest=0;
}


/* --- Push (add) the value of the array (double[3]) as element of the vector of points --- */

void OperationsHelper::
pushVectorAsPoint(Vector3d v, std::vector<Vector3d>& points)
{
  Vector3d tempoPoint;
  tempoPoint = v;
  // tempoPoint(0) = fabs(v(0)) > TRI_EPSILON ? v(0) : 0.0;
  // tempoPoint(1) = fabs(v(1)) > TRI_EPSILON ? v(1) : 0.0;
  // tempoPoint(2) = fabs(v(2)) > TRI_EPSILON ? v(2) : 0.0;
  points.push_back(tempoPoint);
}


/* --- Eliminate those points that are very similar to each other
       keeping just one of them (prune the points)                --- */
void OperationsHelper::
prunePoints( std::vector<Vector3d>& points )
{
  for (int i=0; i<points.size(); i++){
    for (int j=i+1; j<points.size(); j++)
      {
	double norm2 = ( (points[i](0)-points[j](0)) * (points[i](0)-points[j](0)) ) +
	  ( (points[i](1)-points[j](1)) * (points[i](1)-points[j](1)) ) + 
	  ( (points[i](2)-points[j](2)) * (points[i](2)-points[j](2)) );
	if (norm2 < (tolerance_distance*tolerance_distance) )
	  {
	    points.erase( points.begin()+j );
	    j--;
	  }
      }
  }
}
// void OperationsHelper::
// prunePoints( std::vector<Point3d>& points )
// {
//   for (int i=0; i<points.size(); i++){
//     for (int j=i+1; j<points.size(); j++){
//       if ( fabs(points[i].x-points[j].x) < DIST_TOLERANCE ){
// 	if ( fabs(points[i].y-points[j].y) < DIST_TOLERANCE ){
// 	  if ( fabs(points[i].z-points[j].z) < DIST_TOLERANCE ){
// 	    points.erase( points.begin()+j );
// 	    j--;
// 	  }
// 	}
//       }
//     }
//   }
// }



/* --- Find the convex hull (on the plane that maximizes the area of the points) of a set of points --- */
int OperationsHelper::
convexHull( std::vector<Vector3d>& points )
{
  // Discard the verification if there are only 2 or 1 contact points
  if (points.size() < 3)
    { 
      //std::cout << "Less than 3 points" << std::endl;
      return 0;
    }

  // Find the normal to the collision plane
  /* Let the points be p1, p2, p3, p4, ...
     - Use (p1, p2) to get the first vector 'v1'. 
     - Use (p1, p3) to get the second vector 'v2'
     - Find the 'normal' as the cross product of v1 x v2. If is different to zero (v1,v2 are not collinear),
       then use the value as the normal. Otherwise, compute v2 using (p1,p4), then (p1,p5), ... until the
       'normal' is not zero (if possible). If the normal is always zero, then, all the points lie on the same
       line and the extremes of the line have to be found
  */
  // double p1[3] = {points[0].x, points[0].y, points[0].z};
  // double p2[3] = {points[1].x, points[1].y, points[1].z};
  // Vector3d p1 = points[0];
  // Vector3d p2 = points[1];

  //double p3[3], v1[3], v2[3], normal[3];
  Vector3d p3, v1, v2, normal;
  bool normal_not_zero;
  //v1 = p1-p2;
  v1 = points[0]-points[1];

  for (int i=2; i<points.size(); i++) {
    // p3[0] = points[i].x; p3[1] = points[i].y; p3[2] = points[i].z;
    // sub(v2, p1, p3);
    // cross(normal, v1, v2);
    //p3 = points[i];
    //v2 = p1 - p3;
    v2 = points[0] - points[i];
    normal = v1.cross(v2);
    normal_not_zero = fabs(normal(0))>DIST_TOLERANCE || fabs(normal(1))>DIST_TOLERANCE || fabs(normal(2))>DIST_TOLERANCE;
    if (normal_not_zero)
      break;
  }
  
  //std::cout << "Normal not zero (if 0: >= 3 colinear points): " << normal_not_zero << std::endl;

  /* All the points are collinear (the intersection is a line), then, find the extremes */

  /* The equation of the line passing through the points is P = p1 + t*v1, then, the extremes
     wil be those points for which 't' is minimum and maximum. Find t as: t=(P-p1)/v1, where
     only the first component (or second [or third], if the first [or second] is equal to zero) 
     of the vectors is considered. We can just compare one component, since we know that the points
     are collinear, and thus, t must be the same for any component different to zero      
  */
  if (!normal_not_zero)
    {
      // std::cout << "Some points are collinear (convex hull function)" << std::endl;
      int ind=0;                // Index for v1 (and the other vectors)
      int imax=0, imin=0;       // Indexes corresponding to the max and min of 't'
      //double tmax=0, tmin=0;    // Maximum and minimum values of 't' (initially 0 since P=p1)
      double tmax=-1e10, tmin=1e10;    // Maximum and minimum values of 't' (initially 0 since P=p1)
      double ttemp;             // Temporal value to store 't'

      // std::cout << v1.transpose() << std::endl;
      if (fabs(v1(ind))<tolerance_distance_collinear){
	ind=1;
	if (fabs(v1(ind))<tolerance_distance_collinear){
	  ind=2;
	}
      }
 
      // Find the maximum and minimum t
      for (int j=0; j<points.size(); j++) {
	ttemp = ( points[j](ind) - points[1](ind) )/v1(ind);
	// std::cout << "* (" <<  points[j](ind) << " - " << points[0](ind) << ") / " << v1(ind)
	// 	  << " = " << ttemp << std::endl;
	if (ttemp>tmax){
	  tmax=ttemp; imax=j;
	}
	else if(ttemp<tmin){
	  tmin=ttemp; imin=j;
	}
      }
      
      // std::cout << "- Min: " << tmin << " (j=" << imin << "), Max: " << tmax
      // 		<< " (j=" << imax << ")" << std::endl;

      // Create a copy with the same elements as points (collision points)
      std::vector<Vector3d> pointsTemp;
      pointsTemp = points;
      points.clear();
      // Keep only the points corresponding to the extremes of the line
      points.push_back(pointsTemp[imin]);
      points.push_back(pointsTemp[imax]);

      return 0;
    }
  

  /* Project onto an axis-aligned plane that maximizs the area of the polygon: compute indexes */
  // Discard the largest normal component
  //double absN[3];
  Vector3d absN;
  int x, y, inot;
  absN(0) = fabs(normal(0)); absN(1)=fabs(normal(1)); absN(2)=fabs(normal(2));
  if( absN(0)>absN(1) ) {
    if( absN(0)>absN(2) ) {
      x=1; y=2;  inot=0;     /* absN[0] is greatest */
    }
    else {
      x=0; y=1;  inot=2;     /* absN[2] is greatest */
    }
  }
  else {                       /* absN[0]<=absN[1] */
    if( absN(2)>absN(1) ) {
      x=0; y=1;  inot=2;     /* absN[2] is greatest */
    }
    else {
      x=0; y=2;  inot=1;     /* absN[1] is greatest */
    }
  }               
  
  /*  Find the minimum point in y direction */
  double ymin=points[0](y), imin=0;
  
  for (int i=1; i<points.size(); i++) {
    if (points[i](y) < ymin ){
      ymin = points[i](y);
      imin = i;
    }
  }

  //std::cout << "y: " << y << ", imin: " << imin << std::endl;

  /* Compute the angle as an estimator of the desired edge */
  /*    Starting with the minimum point in the 'x' direction, and assuming that
	we want to go through the convex hull in counterclockwise sense, then, the 
	initial 'vector' is horizontal (1,0). Then, find the angles, the minimum 
	is the next point in the convex hull. 
  */

  int i = imin, k;  // k: index of point with smallest theta
  double e1[2]={1.0, 0.0};    // First edge
  double e2[2];          // Second edge
  double p0[2] = {points[i](x), points[i](y)};   // First Point (lowest)
  double pi[2];          // Second Point
  double mod_e1=1, mod_e2, th, thmin;   // |e1|, |e2|
  std::vector<Vector3d> interPointsRectCH;
  
  //std::cout << "p0 = " << p0[0] << ", " << p0[1] <<std::endl;
  int counter=0;
  do{
    //std::cout << "e1 = (" << e1[0] << ", " << e1[1] << ")\n";
    thmin = 7.0;
    for (int j=0; j<points.size(); j++){
      if (j==i) continue;    // for each j!=i
      pi[0]=points[j](x); pi[1]=points[j](y);    // Point to test
      //std::cout << "pi = (" << pi[0] << ", " << pi[1] << ")  ";

      e2[0]=pi[0]-p0[0]; e2[1]=pi[1]-p0[1];                        // edge from p0 to pi
      //std::cout << "e2 = (" << e2[0] << ", " << e2[1] << ") ";

      mod_e2 = sqrt(e2[0]*e2[0]+e2[1]*e2[1]);                      // |e2|
      th = acos( (e2[0]*e1[0]+e2[1]*e1[1])/(mod_e2*mod_e1) );      // theta
      //std::cout << "  th = " << th << "  num: " << e2[0]*e1[0]+e2[1]*e1[1] << "  den: " << mod_e2*mod_e1 <<std::endl;

      // If the angle of the previous edge with the next edge in the convex hull is zero, the current point lies
      // in a line (is part of the segment), then, eliminate the last point added to the convex hull since it lied
      // in the line. 
      if (counter!=0 && th==0)
	  interPointsRectCH.pop_back();

      if (th<thmin){
	thmin = th;
	k = j;
      }
    }
    interPointsRectCH.push_back(points[k]);
    i=k;
    e1[0]=points[k](x)-p0[0]; e1[1]=points[k](y)-p0[1];
    mod_e1 = sqrt(e1[0]*e1[0]+e1[1]*e1[1]);                      // |e2|
    p0[0]=points[i](x); p0[1]=points[i](y);
    
    counter++;
    } while(i!=imin);

  // for (int ii=0; ii<interPointsRectCH.size(); ii++){
  //   std::cout << "( " << interPointsRectCH[ii].x << ", " << interPointsRectCH[ii].y
  // 	      << ", " << interPointsRectCH[ii].z << ")" << std::endl;
  // }

  points = interPointsRectCH;
  return 0;

}



/* --- Find the convex hull (on the plane that maximizes the area of the points) of a set of points --- */
// int OperationsHelper::
// convexHull( std::vector<Point3d>& points )
// {
//   // Discard the verification if there are only 2 or 1 contact points
//   if (points.size() < 3)
//     return 0;

//   // Find the normal to the collision plane
//   /* Let the points be p1, p2, p3, p4, ...
//      - Use (p1, p2) to get the first vector 'v1'. 
//      - Use (p1, p3) to get the second vector 'v2'
//      - Find the 'normal' as the cross product of v1 x v2. If is different to zero (v1,v2 are not collinear),
//        then use the value as the normal. Otherwise, compute v2 using (p1,p4), then (p1,p5), ... until the
//        'normal' is not zero (if possible). If the normal is always zero, then, all the points lie on the same
//        line and the extremes of the line have to be found
//   */
//   double p1[3] = {points[0].x, points[0].y, points[0].z};
//   double p2[3] = {points[1].x, points[1].y, points[1].z};
//   double p3[3], v1[3], v2[3], normal[3];
//   bool normal_not_zero;
//   sub(v1, p1, p2); 

//   for (int i=2; i<points.size(); i++) {
//     p3[0] = points[i].x; p3[1] = points[i].y; p3[2] = points[i].z;
//     sub(v2, p1, p3);
//     cross(normal, v1, v2);
//     normal_not_zero = fabs(normal[0])>DIST_TOLERANCE || fabs(normal[1])>DIST_TOLERANCE || fabs(normal[2])>DIST_TOLERANCE;
//     if (normal_not_zero)
//       break;
//   }
  
//   /* All the points are collinear (the intersection is a line), then, find the extremes */

//   /* The equation of the line passing through the points is P = p1 + t*v1, then, the extremes
//      wil be those points for which 't' is minimum and maximum. Find t as: t=(P-p1)/v1, where
//      only the first component (or second [or third], if the first [or second] is equal to zero) 
//      of the vectors is considered. We can just compare one component, since we know that the points
//      are collinear, and thus, t must be the same for any component different to zero      
//   */
//   if (!normal_not_zero){
//     // std::cout << "Some points are collinear (convex hull function)" << std::endl;
//     int ind=0;             // Index for v1 (and the other vectors)
//     int imax=0, imin=0;    // Indexes corresponding to the max and min of 't'
//     double tmax=0, tmin=0;    // Maximum and minimum values of 't' (initially 0 since P=p1)
//     double ttemp;             // Temporal value to store 't'

//     if (fabs(v1[ind])<DIST_TOLERANCE){
//       ind=1;
//       if (fabs(v1[ind])<DIST_TOLERANCE){
// 	ind=2;
//       }
//     }
 
//     // Find the maximum and minimum t
//     for (int j=1; j<points.size(); j++) {
//       ttemp = ( points[j][ind] - points[0][ind])/v1[ind];
//       // std::cout << "tnum: " <<  points[j][ind] - points[0][ind] << "  v1: " << v1[ind]
//       // 		<< "  t= " << ttemp << std::endl;
//       if (ttemp>tmax){
// 	tmax=ttemp; imax=j;
//       }
//       else if(ttemp<tmin){
// 	tmin=ttemp; imin=j;
//       }
//     }
    
//     // std::cout << "- Extreme Values of t: " << tmin << "  " << tmax << std::endl;
//     // std::cout << "- Extreme Values of the index:  " << imin << "  " << imax << std::endl;

//     // Eliminate the points that are not the extremes
//     for (int j=0; j<points.size(); j++) {
//       if (j==imax || j==imin)
// 	continue;
//       else {
// 	// std::cout << "index erased: " << j << std::endl;
// 	points.erase( points.begin()+j );
//       }
//     }
    
//     return 0;
//    }
    

//   /* Project onto an axis-aligned plane that maximizs the area of the polygon: compute indexes */
//   // Discard the largest normal component
//   double absN[3];
//   int x, y, inot;
//   absN[0] = fabs(normal[0]); absN[1]=fabs(normal[1]); absN[2]=fabs(normal[2]);
//   if( absN[0]>absN[1] ) {
//     if( absN[0]>absN[2] ) {
//       x=1; y=2;  inot=0;     /* absN[0] is greatest */
//     }
//     else {
//       x=0; y=1;  inot=2;     /* absN[2] is greatest */
//     }
//   }
//   else {                       /* absN[0]<=absN[1] */
//     if( absN[2]>absN[1] ) {
//       x=0; y=1;  inot=2;     /* absN[2] is greatest */
//     }
//     else {
//       x=0; y=2;  inot=1;     /* absN[1] is greatest */
//     }
//   }               
  
//   /*  Find the minimum point in y direction */
//   double ymin=points[0][y], imin=0;
  
//   for (int i=1; i<points.size(); i++) {
//     if (points[i][y] < ymin ){
//       ymin = points[i][y];
//       imin = i;
//     }
//   }

//   //std::cout << "y: " << y << ", imin: " << imin << std::endl;

//   /* Compute the angle as an estimator of the desired edge */
//   /*    Starting with the minimum point in the 'x' direction, and assuming that
// 	we want to go through the convex hull in counterclockwise sense, then, the 
// 	initial 'vector' is horizontal (1,0). Then, find the angles, the minimum 
// 	is the next point in the convex hull. 
//   */

//   int i = imin, k;  // k: index of point with smallest theta
//   double e1[2]={1.0, 0.0};    // First edge
//   double e2[2];          // Second edge
//   double p0[2] = {points[i][x], points[i][y]};   // First Point (lowest)
//   double pi[2];          // Second Point
//   double mod_e1=1, mod_e2, th, thmin;   // |e1|, |e2|
//   std::vector<Point3d> interPointsRectCH;
  
//   //std::cout << "p0 = " << p0[0] << ", " << p0[1] <<std::endl;
//   int counter=0;
//   do{
//     //std::cout << "e1 = (" << e1[0] << ", " << e1[1] << ")\n";
//     thmin = 7.0;
//     for (int j=0; j<points.size(); j++){
//       if (j==i) continue;    // for each j!=i
//       pi[0]=points[j][x]; pi[1]=points[j][y];    // Point to test
//       //std::cout << "pi = (" << pi[0] << ", " << pi[1] << ")  ";

//       e2[0]=pi[0]-p0[0]; e2[1]=pi[1]-p0[1];                        // edge from p0 to pi
//       //std::cout << "e2 = (" << e2[0] << ", " << e2[1] << ") ";

//       mod_e2 = sqrt(e2[0]*e2[0]+e2[1]*e2[1]);                      // |e2|
//       th = acos( (e2[0]*e1[0]+e2[1]*e1[1])/(mod_e2*mod_e1) );      // theta
//       //std::cout << "  th = " << th << "  num: " << e2[0]*e1[0]+e2[1]*e1[1] << "  den: " << mod_e2*mod_e1 <<std::endl;

//       // If the angle of the previous edge with the next edge in the convex hull is zero, the current point lies
//       // in a line (is part of the segment), then, eliminate the last point added to the convex hull since it lied
//       // in the line. 
//       if (counter!=0 && th==0)
// 	  interPointsRectCH.pop_back();

//       if (th<thmin){
// 	thmin = th;
// 	k = j;
//       }
//     }
//     interPointsRectCH.push_back(points[k]);
//     i=k;
//     e1[0]=points[k][x]-p0[0]; e1[1]=points[k][y]-p0[1];
//     mod_e1 = sqrt(e1[0]*e1[0]+e1[1]*e1[1]);                      // |e2|
//     p0[0]=points[i][x]; p0[1]=points[i][y];
    
//     counter++;
//     } while(i!=imin);

//   // for (int ii=0; ii<interPointsRectCH.size(); ii++){
//   //   std::cout << "( " << interPointsRectCH[ii].x << ", " << interPointsRectCH[ii].y
//   // 	      << ", " << interPointsRectCH[ii].z << ")" << std::endl;
//   // }

//   points = interPointsRectCH;
//   return 0;

// }

