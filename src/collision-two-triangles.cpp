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
 * Triangle/triangle intersection test
 * Based on the paper:
 *     Thomas Moller, "A Fast Triangle-Triangle Intersection Test",
 *     Journal of Graphics Tools, 2(2), 1997
 *
 */

#include <collision-two-triangles.h>
#include <iostream>



/*==============================================================================
          INITIALIZATION, SETTERS AND GETTERS   
/*==============================================================================


/* --- Default Constructor that initializes all the vertices to zero --- */
CollisionTwoTriangles::
CollisionTwoTriangles( void )
{
  defaultInit();
}


/* --- Constructor that initializes the vertices of both triangles --- */
CollisionTwoTriangles::
CollisionTwoTriangles(Vector3d t1v1, Vector3d t1v2, Vector3d t1v3,
		      Vector3d t2v1, Vector3d t2v2, Vector3d t2v3)
{
  T1.setVertices(t1v1, t1v2, t1v3);
  T2.setVertices(t2v1, t2v2, t2v3);
  defaultInit();
}


/* --- Constructor that initializes both triangles --- */
CollisionTwoTriangles::
CollisionTwoTriangles(Triangle t1in, Triangle t2in)
{
  T1=t1in; T2 = t2in;
  defaultInit();
}


/* --- Initialize some values by default --- */
void CollisionTwoTriangles::
defaultInit( void )
{
  D << 0, 0, 0;
  i0=0, i1=0, inot=0; collisionIndicator=0;
  coplanar_tolerance = 0.0005;   // 0.5mm
}


/* --- Set the vertices of both triangles --- */
void CollisionTwoTriangles::
setVerticesAll(Vector3d t1v1, Vector3d t1v2, Vector3d t1v3,
	       Vector3d t2v1, Vector3d t2v2, Vector3d t2v3)
{
  T1.setVertices(t1v1, t1v2, t1v3);
  T2.setVertices(t2v1, t2v2, t2v3);
}


/* --- Set the vertices of triangle 1 --- */
void CollisionTwoTriangles::
setVerticesT1(Vector3d V1in, Vector3d V2in, Vector3d V3in)
{
  T1.setVertices(V1in, V2in, V3in);
  //V1=V1in; V2=V2in; V3=V3in;
}


/* --- Set the vertices of triangle 2 --- */
void CollisionTwoTriangles::
setVerticesT2(Vector3d V1in, Vector3d V2in, Vector3d V3in)
{
  T2.setVertices(V1in, V2in, V3in);
}


/* --- Set both triangles --- */
void CollisionTwoTriangles::
setTriangles(Triangle t1in, Triangle t2in)
{
  T1=t1in; T2=t2in;
}


/* --- Set triangle 1 --- */
void CollisionTwoTriangles::
setTriangle1(Triangle t1in)
{
  T1=t1in;
}


/* --- Set triangle 2 --- */
void CollisionTwoTriangles::
setTriangle2(Triangle t2in)
{
  T2=t2in;
}


/* --- Get the vertices of triangle 1 --- */
void CollisionTwoTriangles::
getVerticesT1(Vector3d &V1out, Vector3d &V2out, Vector3d &V3out)
{
  T1.getVertices(V1out, V2out, V3out);
}


/* --- Get the vertices of triangle 2 --- */
void CollisionTwoTriangles::
getVerticesT2(Vector3d &V1out, Vector3d &V2out, Vector3d &V3out)
{
  T2.getVertices(V1out, V2out, V3out);
}



/*==============================================================================
         --------- MAIN COMPUTATION OF THE INTERSECTION -------------------
  ============================================================================== */

/* --- Detect if both triangles are colliding --- */
/*  It returns 1 if collision, 0 if no collision  */
/*  The intersection points are also computed     */

int CollisionTwoTriangles::
computeTTintersections( void )
{
  Vector2d isect1, isect2;
  double d1, d2;
  Vector3d distT2, distT1;    // distances from the vertices to the plane of the other triangle
  Vector3d N1, N2;            // Normals
  
  double du0du1, du0du2, dv0dv1, dv0dv2;
  Vector3d vp, up;
  int index;

  /* Clear the possible 'previous' contact points */
  pointsTT.clear();

  /* Compute the plane equation of triangle T1 (V1,V2,V3) */
  /*   plane equation 1: N1.X+d1=0                        */
  N1 = (T1.v2 - T1.v1).cross(T1.v3 - T1.v1);
  d1 = -N1.dot(T1.v1);     //d1 = -dot(N1,V1);   

  /* Put the vertices of T2 into the plane of T1 to compute the signed distances to the plane of T1
     distT2: distance from the vertices of T2 to the plane of T1 */
  distT2(0) = N1.dot(T2.v1) + d1;
  distT2(1) = N1.dot(T2.v2) + d1;
  distT2(2) = N1.dot(T2.v3) + d1;

  /* Coplanarity robustness check
     If the distance is smaller than coplanar_tolerance, the vertex is assumed to be in the plane */
  if( fabs(distT2(0) ) < coplanar_tolerance ) distT2(0) = 0.0;
  if( fabs(distT2(1) ) < coplanar_tolerance ) distT2(1) = 0.0;
  if( fabs(distT2(2) ) < coplanar_tolerance ) distT2(2) = 0.0;
  du0du1 = distT2(0)*distT2(1); du0du2 = distT2(0)*distT2(2);          // To get the sign

  /* If the distances of all the vertices of T2 to the plane of T1 have the same sign,
     they are on one side of T1, and theere is no intersection  */
  if(du0du1>0.0f && du0du2>0.0f) 
    {
      collisionIndicator = 0;
      return 0;                    /* no intersection occurs */
    }

  /* Compute the plane equation of triangle T2 (U1,U2,U3) */
  /*   plane equation 2: N2.X+d2=0                        */
  N2 = (T2.v2-T2.v1).cross(T2.v3-T2.v1);
  d2 = -N2.dot(T2.v1);

  /* put V1,V2,V3 into plane equation 2 */
  /* Put the vertices of T1 into the plane of T2 to compute the signed distances to the plane of T2
     distT1: signed distance from the vertices of T1 to the plane of T2 */
  distT1(0) = N2.dot(T1.v1) + d2;
  distT1(1) = N2.dot(T1.v2) + d2;
  distT1(2) = N2.dot(T1.v3) + d2;

  /* Coplanarity robustness check */
  if( fabs(distT1(0) ) < coplanar_tolerance) distT1(0) = 0.0;
  if( fabs(distT1(1) ) < coplanar_tolerance) distT1(1) = 0.0;
  if( fabs(distT1(2) ) < coplanar_tolerance) distT1(2) = 0.0;
  dv0dv1 = distT1(0)*distT1(1); dv0dv2 = distT1(0)*distT1(2);

  /* If the distances of all the vertices of T1 to the plane of T2 have the same sign,
     they are on one side of T2, and theere is no intersection  */
  if(dv0dv1>0.0f && dv0dv2>0.0f) {
    collisionIndicator = 0;
    return 0;                    /* no intersection occurs */
  }


  /* Compute the direction of the intersection line */
  D = N1.cross(N2);

  /* Compute and index to the largest component of D */
  // TODO ------- CREATE A FUNCTION TO FIND THE MAX AND ITS INDEX
  double b, c, max;
  max = fabs( D(0) );
  index = 0;
  b = fabs( D(1) );
  c = fabs( D(2) );
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* Simplified projection onto L*/
  vp(0)=T1.v1(index); vp(1)=T1.v2(index); vp(2)=T1.v3(index);
  up(0)=T2.v1(index); up(1)=T2.v2(index); up(2)=T2.v3(index);
  
  Vector3d isectpointA1, isectpointA2, isectpointB1, isectpointB2;
  Vector3d isectpt1, isectpt2;
  int smallest1,smallest2, coplanar;
  bool coplanarCollision;
  
  /* compute interval for triangle 1 */
  coplanar = compute_intervals_isectline(1, vp, distT1, dv0dv1, dv0dv2,
  					 isect1, isectpointA1, isectpointA2);

  /* Coplanar Planes? */
  if(coplanar){
    //std::cout << "\nTriangles are coplanar" << std::endl;
    coplanarCollision = coplanar_tri_tri(N1);
    //std::cout << "Coplanar Collision: " << coplanarCollision << std::endl;
    if (coplanarCollision)
      {
  	//printTTcollisionInformation( );
  	prunePoints(pointsTT);             	// To remove repeated elements
  	collisionIndicator = 1;
  	return 1;
      }
    collisionIndicator = 0;
    return 0;
  }

  /* compute interval for triangle 2 */
  compute_intervals_isectline(2, up, distT2, du0du1, du0du2,
  			      isect2, isectpointB1, isectpointB2);

  sort2(isect1, smallest1);
  sort2(isect2, smallest2);
  
  if(isect1(1)<isect2(0) || isect2(1)<isect1(0)) {
    collisionIndicator = 0;
    return 0;
  }

  /* at this point, we know that the triangles intersect */

  if(isect2(0)<isect1(0))
  {
    if(smallest1==0) {
      isectpt1 = isectpointA1; 
    }
    else {
      isectpt1 = isectpointA2;
    }
    if(isect2(1)<isect1(1)) {
      if(smallest2==0) {
  	isectpt2 = isectpointB2;
      }
      else {
  	isectpt2 = isectpointB1;
      }
    }
    else {
      if(smallest1==0) {
  	isectpt2 = isectpointA2;
      }
      else {
  	isectpt2 = isectpointA1; 
      }
    }
  }
  else  {
    if(smallest2==0) {
      isectpt1 = isectpointB1;
    }
    else {
      isectpt1 = isectpointB2; 
    }
    if(isect2(1)>isect1(1))
    {
      if(smallest1==0) {
  	isectpt2 = isectpointA2; 
      }
      else {
  	isectpt2 = isectpointA1; 
      }      
    }
    else {
      if(smallest2==0) {
  	isectpt2 = isectpointB2;
      }
      else {
  	isectpt2 = isectpointB1;
      } 
    }
  }
    
  // std::cout << "( " << isectpt1[0] << ", " << isectpt1[1] << ", " << isectpt1[2] << ")\n"
  // 	    << "( " << isectpt2[0] << ", " << isectpt2[1] << ", " << isectpt2[2] << ")\n";
  pushVectorAsPoint(isectpt1, pointsTT);
  pushVectorAsPoint(isectpt2, pointsTT);

  // To remove repeated elements
  prunePoints(pointsTT);
  
  collisionIndicator = 1;
  return 1;
 }





/*==============================================================================
  ---------------- INTERNAL FUNCTIONS -------------------------------------
  ============================================================================== */

// Returns 1 if the triangles are coplanar, 0 if they are not coplanar
int CollisionTwoTriangles::
compute_intervals_isectline(int id, Vector3d VV, Vector3d DD,
			    double D0D1, double D0D2, Vector2d &isect,
			    Vector3d &isectpoint0, Vector3d &isectpoint1)
{
  Vector3d Vert1, Vert2, Vert3;
  // For triangle 1
  if (id==1) { Vert1=T1.v1; Vert2=T1.v2; Vert3=T1.v3; }
  // For triangle 2
  else if (id==2){ Vert1=T2.v1; Vert2=T2.v2; Vert3=T2.v3; }

  if(D0D1>0.0f)                                        
  {                                                    
    /* here we know that D0D2<=0.0 */                  
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(Vert3, Vert1, Vert2, VV(2), VV(0), VV(1), DD(2), DD(0), DD(1),
	   isect, isectpoint0, isectpoint1);
  } 
  else if(D0D2>0.0f)                                   
    {                                                   
    /* here we know that d0d1<=0.0 */             
    isect2(Vert2, Vert1, Vert3, VV(1), VV(0), VV(2), DD(1), DD(0), DD(2),
	   isect, isectpoint0, isectpoint1);
  }                                                  
  else if(DD(1)*DD(2)>0.0f || DD(0)!=0.0f)   
  {                                   
    /* here we know that d0d1<=0.0 or that D0!=0.0 */
    isect2(Vert1, Vert2, Vert3, VV(0), VV(1), VV(2), DD(0), DD(1), DD(2),
	   isect, isectpoint0, isectpoint1);   
  }                                                  
  else if(DD(1)!=0.0f)                                  
  {                                               
    isect2(Vert2, Vert1, Vert3, VV(1), VV(0), VV(2), DD(1), DD(0), DD(2), 
	   isect, isectpoint0, isectpoint1); 
  }                                         
  else if(DD(2)!=0.0f)                                  
  {                                                   
    isect2(Vert3, Vert1, Vert2, VV(2), VV(0), VV(1), DD(2), DD(0), DD(1),
	   isect, isectpoint0, isectpoint1);     
  }                                                 
  else                                               
  {                                                   
    /* triangles are coplanar */    
    return 1;
  }
  return 0;
}


void CollisionTwoTriangles::
isect2(Vector3d VTX0, Vector3d VTX1, Vector3d VTX2, 
       double VV1, double VV2, double VV3, double D0,double D1,double D2,
       Vector2d &isect, Vector3d &isectpoint0, Vector3d &isectpoint1) 
{
  double tmp;
  Vector3d diff;

  tmp = D0/(D0-D1);
  isect(0) = VV1 + (VV2-VV1)*tmp;
  isectpoint0 =  tmp*(VTX1-VTX0) + VTX0;

  tmp = D0/(D0-D2);           
  isect(1) = VV1 + (VV3-VV1)*tmp;
  isectpoint1 = tmp*(VTX2-VTX0) + VTX0;
}


bool CollisionTwoTriangles::
coplanar_tri_tri(Vector3d N)
{
  Vector3d absN;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
  //absN[0]=fabs(N[0]);  absN[1]=fabs(N[1]);  absN[2]=fabs(N[2]);
  absN = N.cwiseAbs();
  if( absN(0)>absN(1) ) {
    if( absN(0)>absN(2) ) {
       i0=1; i1=2;  inot=0;     /* absN[0] is greatest */
     }
     else {
       i0=0; i1=1;  inot=2;     /* absN[2] is greatest */
     }
   }
   else {                       /* absN[0]<=absN[1] */
     if( absN(2)>absN(1) ) {
       i0=0; i1=1;  inot=2;     /* absN[2] is greatest */
     }
     else {
       i0=0; i1=2;  inot=1;     /* absN[1] is greatest */
     }
   }               

   bool collision[9];
   /* Test all edges of triangle 1 against the edges of triangle 2 */
   collision[0] = edge_against_tri_edges(T1.v1, T1.v2);
   collision[1] = edge_against_tri_edges(T1.v2, T1.v3);
   collision[2] = edge_against_tri_edges(T1.v3, T1.v1);
   
   /* Finally, test if tri1 is totally contained in tri2 or vice versa */
   /*  check if one point is contained and if so, assume all of it is contained */
   // collision[3] = point_in_tri(V1, U1, U2, U3);
   // collision[4] = point_in_tri(U1, V1, V2, V3);
   collision[3] = point_in_tri(T1.v1, T2.v1, T2.v2, T2.v3);
   collision[4] = point_in_tri(T1.v2, T2.v1, T2.v2, T2.v3);
   collision[5] = point_in_tri(T1.v3, T2.v1, T2.v2, T2.v3);
   collision[6] = point_in_tri(T2.v1, T1.v1, T1.v2, T1.v3);
   collision[7] = point_in_tri(T2.v2, T1.v1, T1.v2, T1.v3);
   collision[8] = point_in_tri(T2.v3, T1.v1, T1.v2, T1.v3);
   
   return (collision[0] || collision[1] || collision[2] || collision[3] || collision[6]);
}


/* --- Test if an edge of T1 intersects the edges of T2 -- */
bool CollisionTwoTriangles::
edge_against_tri_edges(Vector3d v0, Vector3d v1)
{
  double Ax, Ay;
  bool collision[3];       // True if there is collision
  // Components of the edge in T1
  Ax = v1(i0) - v0(i0);    
  Ay = v1(i1) - v0(i1);    
  /* test edge U1,U2 against v0,v1 */
  //std::cout << "\nEdge: "; printVector(v0); std::cout << " , "; printVector(v1); std::cout << "with:";

  //std::cout << "\n - edge "; printVector(U1); std::cout << " , "; printVector(U2);
  collision[0] = edge_edge_test(v0, T2.v1, T2.v2, Ax, Ay);
  //std::cout << ": result is " << collision[0] << std::endl;

  /* test edge U2,U3 against v0,v1 */
  //std::cout << "\n - edge "; printVector(U2); std::cout  << " , "; printVector(U3);
  collision[1] = edge_edge_test(v0, T2.v2, T2.v3, Ax, Ay);
  //std::cout << ": result is " << collision[1] << std::endl;

  /* test edge U3,U2 against v0,v1 */
  //std::cout << "\n - edge "; printVector(U3); std::cout  << " , "; printVector(U1);
  collision[2] = edge_edge_test(v0, T2.v3, T2.v1, Ax, Ay);
  //std::cout << ": result is " << collision[2] << std::endl;

  return (collision[0] || collision[1] || collision[2]);
}


/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
bool CollisionTwoTriangles::
edge_edge_test(Vector3d v0, Vector3d u0, Vector3d u1, double Ax, double Ay)
{
  double Bx, By, Cx, Cy, numu, numv, den;
  double alpha, beta;
  Vector3d Pintu, Pintv;
  Bx = u0(i0) - u1(i0);
  By = u0(i1) - u1(i1);
  Cx = v0(i0) - u0(i0);
  Cy = v0(i1) - u0(i1);
  den  = Ay*Bx - Ax*By;   // denominator
  numv  = By*Cx - Bx*Cy;   // numerator 1
  if((den>0 && numv>=0 && numv<=den) || (den<0 && numv<=0 && numv>=den)) {
    numu = Ax*Cy-Ay*Cx;      // numerator 2
    if(den>0) {
      if(numu>=0 && numu<=den) {
	beta = numu/den;
	Pintu(i0) = u0(i0) + beta*(u1(i0)-u0(i0));
	Pintu(i1) = u0(i1) + beta*(u1(i1)-u0(i1));
	Pintu(inot) = u0(inot) + beta*(u1(inot)-u0(inot));
	//std::cout << " Intersection at: ( " << Pintu[0) << ", " << Pintu[1) << ", " << Pintu[2) << ")";
	pushVectorAsPoint(Pintu, pointsTT);
	return (true);
      }
    }
    else {
      if(numu<=0 && numu>=den){
	beta = numu/den;
	Pintu(i0) = u0(i0) + beta*(u1(i0)-u0(i0));
	Pintu(i1) = u0(i1) + beta*(u1(i1)-u0(i1));
	Pintu(inot) = u0(inot) + beta*(u1(inot)-u0(inot));
	//std::cout << " Intersection at: ( " << Pintu[0] << ", " << Pintu[1] << ", " << Pintu[2] << ")";
	pushVectorAsPoint(Pintu, pointsTT);
	return (true);
      }
    }
  }
  return (false);
}


/* --- Check if v0 is inside tri(u0,u1,u2) --- */
bool CollisionTwoTriangles::
point_in_tri(Vector3d v0, Vector3d u0, Vector3d u1, Vector3d u2)
{
  double a,b,c,d0,d1,d2;
  /* is T1 completly inside T2? */
  /* check if v0 is inside tri(u0,u1,u2) */
  a  =  u1(i1) - u0(i1);
  b  = -(u1(i0)-u0(i0));
  c  = -a*u0(i0) - b*u0(i1);
  d0 = a*v0(i0) + b*v0(i1) + c;

  a  = u2(i1) - u1(i1);
  b  = -(u2(i0)-u1(i0));
  c  = -a*u1(i0) - b*u1(i1);
  d1 = a*v0(i0) + b*v0(i1) + c;

  a  = u0(i1) - u2(i1);
  b  = -(u0(i0)-u2(i0));
  c  = -a*u2(i0) - b*u2(i1);
  d2 = a*v0(i0) + b*v0(i1) + c;

  if(d0*d1>0.0) {
    if(d0*d2>0.0)
      {
	//std::cout << "\n Vertex "; printVector(v0); std::cout<< " inside.";
	// Point v0 inside the other triangle (extreme)
	pushVectorAsPoint(v0, pointsTT);
	return (true);
      }
  }
  return (false);
}


/* --- Print vertices of the triangles --- */
void CollisionTwoTriangles::
printTrianglesVertices( void )
{
  std::cout << "Triangle 1: "; T1.print();
  std::cout << "Triangle 2: "; T2.print();
}


/* --- Print verbose information relative to the collision --- */
void CollisionTwoTriangles::
printTTcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Triangles are intersecting at" << std::endl;
    for (int i=0; i<pointsTT.size(); i++){
      std::cout << "  (" << pointsTT[i].transpose() << ")  ";
    }
  }
  else{
    std::cout << "Result: Triangles are NOT intersecting" << std::endl;
  }
}
