/***********************************************************************************
 * Triangle/triangle intersection test
 * Based on the paper:
 *     Thomas Moller, "A Fast Triangle-Triangle Intersection Test",
 *     Journal of Graphics Tools, 2(2), 1997
 *
 *                                           Oscar Efrain Ramos Ponce, LAAS-CNRS
 ************************************************************************************/

#include <collision-two-triangles.h>
#include <iostream>


/* --- Constructor that initializes all the vertices to zero --- */
CollisionTwoTriangles::
CollisionTwoTriangles( void )
{
  double zeros[3] = {0.0, 0.0, 0.0};
  set(V0, zeros); set(V1, zeros); set(V2, zeros);
  set(U0, zeros); set(U1, zeros); set(U2, zeros);
  defaultInit();
}


/* --- Constructor that initializes the vertices of both triangles --- */
CollisionTwoTriangles::
CollisionTwoTriangles(double T1_V0[3], double T1_V1[3], double T1_V2[3],
		      double T2_V0[3], double T2_V1[3], double T2_V2[3])
{
  setVerticesAll(T1_V0, T1_V1, T1_V2, T2_V0, T2_V1, T2_V2);
  defaultInit();
}


/* --- Initialize value by default --- */
void CollisionTwoTriangles::
defaultInit( void )
{
  double zeros[3] = {0.0, 0.0, 0.0};
  set(D, zeros);
  i0=0, i1=0, inot=0; collisionIndicator=0;
}


/* --- Set the vertices of both triangles --- */
void CollisionTwoTriangles::
setVerticesAll(double T1_V0[3], double T1_V1[3], double T1_V2[3],
	       double T2_V0[3], double T2_V1[3], double T2_V2[3])
{
  setVerticesT1(T1_V0, T1_V1, T1_V2);
  setVerticesT2(T2_V0, T2_V1, T2_V2);
}


/* --- Set the vertices of triangle 1 --- */
void CollisionTwoTriangles::
setVerticesT1(double T1_V0[3], double T1_V1[3], double T1_V2[3])
{
  set(V0, T1_V0); set(V1, T1_V1); set(V2, T1_V2);
}


/* --- Set the vertices of triangle 2 --- */
void CollisionTwoTriangles::
setVerticesT2(double T2_V0[3], double T2_V1[3], double T2_V2[3])
{
  set(U0, T2_V0); set(U1, T2_V1); set(U2, T2_V2);
}


/* --- Get the vertices of triangle 1 --- */
void CollisionTwoTriangles::
getVerticesT1(double T1_V0[3], double T1_V1[3], double T1_V2[3])
{
  set(T1_V0, V0); set(T1_V1, V1); set(T1_V2, V2);
}


/* --- Get the vertices of triangle 2 --- */
void CollisionTwoTriangles::
getVerticesT2(double T2_V0[3], double T2_V1[3], double T2_V2[3])
{
  set(T2_V0, U0); set(T2_V1, U1); set(T2_V2, U2);
}



/*==============================================================================

         --------- MAIN COMPUTATION OF THE INTERSECTION -------------------
 
  ============================================================================== */

/* --- Detect if both triangles are colliding --- */
/*  It returns 1 if collision, 0 if no collision  */
/*  The intersection points are NOT computed      */

int CollisionTwoTriangles::
computeTTintersectionsNoLine( void )
{
  double E1[3], E2[3], N1[3],N2[3], isect1[2], isect2[2];
  double d1, d2, du[3], dv[3];
  double du0du1, du0du2, dv0dv1, dv0dv2;
  double vp[3], up[3], b, c, max;
  int index;

  /* Clear the possible 'previous' contact points */
  pointsTT.clear();

  /* Compute the plane equation of triangle T1 (V0,V1,V2) */
  /*   plane equation 1: N1.X+d1=0                        */
  sub(E1,V1,V0); sub(E2,V2,V0);
  cross(N1,E1,E2);
  d1 = -dot(N1,V0);   

  /* Put U0,U1,U2 into plane equation 1 to compute the signed distances to the plane*/
  du[0] = dot(N1,U0) + d1;
  du[1] = dot(N1,U1) + d1;
  du[2] = dot(N1,U2) + d1;

  /* Coplanarity robustness check */
  if( fabs(du[0]) < TRI_EPSILON ) du[0] = 0.0;
  if( fabs(du[1]) < TRI_EPSILON ) du[1] = 0.0;
  if( fabs(du[2]) < TRI_EPSILON ) du[2] = 0.0;
  du0du1 = du[0]*du[1];  du0du2 = du[0]*du[2];

  if(du0du1>0.0f && du0du2>0.0f) { /* same sign on all of them + not equal 0 ? */
    collisionIndicator = 0;
    return 0;                      /* no intersection occurs */
  }

  /* Compute the plane equation of triangle T2 (U0,U1,U2) */
  /*   plane equation 2: N2.X+d2=0                        */
  sub(E1,U1,U0); sub(E2,U2,U0);
  cross(N2,E1,E2);
  d2 = -dot(N2,U0);

  /* put V0,V1,V2 into plane equation 2 */
  dv[0] = dot(N2,V0) + d2;
  dv[1] = dot(N2,V1) + d2;
  dv[2] = dot(N2,V2) + d2;

  /* Coplanarity robustness check */
  if( fabs(dv[0]) < TRI_EPSILON) dv[0] = 0.0;
  if( fabs(dv[1]) < TRI_EPSILON) dv[1] = 0.0;
  if( fabs(dv[2]) < TRI_EPSILON) dv[2] = 0.0;
  dv0dv1 = dv[0]*dv[1];  dv0dv2 = dv[0]*dv[2];
        
  if(dv0dv1>0.0f && dv0dv2>0.0f) { /* same sign on all of them + not equal 0 ? */
    collisionIndicator = 0;
    return 0;                      /* no intersection occurs */
  }

  /* Compute the direction of the intersection line */
  cross(D,N1,N2);

  /* Compute and index to the largest component of D */
  max = fabs(D[0]);
  index = 0;
  b = fabs(D[1]);
  c = fabs(D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* Simplified projection onto L*/
  vp[0]=V0[index]; vp[1]=V1[index]; vp[2]=V2[index];
  up[0]=U0[index]; up[1]=U1[index]; up[2]=U2[index];

  
  /* Compute interval for the triangles */
  bool coplanar, coplanarCollision;
  computeIntervals(vp, dv, dv0dv1, dv0dv2, isect1, coplanar); // For T1
  computeIntervals(up, du, du0du1, du0du2, isect2, coplanar); // For T2
  
  if (coplanar) {
    //std::cout << "Triangles are coplanar" << std::endl;
    coplanarCollision = coplanar_tri_tri(N1);
    if (coplanarCollision) {
      collisionIndicator = 1;
      return 1;
    }
  }
  else {
    //std::cout << "Triangles are not coplanar" << std::endl;
    sort(isect1[0],isect1[1]);
    sort(isect2[0],isect2[1]);
    if(isect1[1]<isect2[0] || isect2[1]<isect1[0]){
      collisionIndicator = 0;
      return 0;
    }
    collisionIndicator = 1;
    return 1;
  }
}


/* --- Detect if both triangles are colliding --- */
/*  It returns 1 if collision, 0 if no collision  */
/*  The intersection points are also computed     */

int CollisionTwoTriangles::
computeTTintersections( void )
{
  double E1[3], E2[3], N1[3],N2[3], isect1[2], isect2[2];
  double d1, d2, du[3], dv[3];
  double du0du1, du0du2, dv0dv1, dv0dv2;
  double vp[3], up[3], b, c, max;
  int index;

  /* Clear the possible 'previous' contact points */
  pointsTT.clear();

  /* Compute the plane equation of triangle T1 (V0,V1,V2) */
  /*   plane equation 1: N1.X+d1=0                        */
  sub(E1,V1,V0); sub(E2,V2,V0);
  cross(N1,E1,E2);
  d1 = -dot(N1,V0);   

  /* Put U0,U1,U2 into plane equation 1 to compute the signed distances to the plane*/
  du[0] = dot(N1,U0) + d1;
  du[1] = dot(N1,U1) + d1;
  du[2] = dot(N1,U2) + d1;

  /* Coplanarity robustness check */
  if( fabs(du[0]) < TRI_EPSILON ) du[0] = 0.0;
  if( fabs(du[1]) < TRI_EPSILON ) du[1] = 0.0;
  if( fabs(du[2]) < TRI_EPSILON ) du[2] = 0.0;
  du0du1 = du[0]*du[1]; du0du2 = du[0]*du[2];

  if(du0du1>0.0f && du0du2>0.0f) { /* same sign on all of them + not equal 0 ? */
    collisionIndicator = 0;
    return 0;                    /* no intersection occurs */
  }

  /* Compute the plane equation of triangle T2 (U0,U1,U2) */
  /*   plane equation 2: N2.X+d2=0                        */
  sub(E1,U1,U0); sub(E2,U2,U0);
  cross(N2,E1,E2);
  d2 = -dot(N2,U0);

  /* put V0,V1,V2 into plane equation 2 */
  dv[0] = dot(N2,V0) + d2;
  dv[1] = dot(N2,V1) + d2;
  dv[2] = dot(N2,V2) + d2;

  /* Coplanarity robustness check */
  if( fabs(dv[0]) < TRI_EPSILON) dv[0] = 0.0;
  if( fabs(dv[1]) < TRI_EPSILON) dv[1] = 0.0;
  if( fabs(dv[2]) < TRI_EPSILON) dv[2] = 0.0;
  dv0dv1 = dv[0]*dv[1]; dv0dv2 = dv[0]*dv[2];
        
  if(dv0dv1>0.0f && dv0dv2>0.0f) { /* same sign on all of them + not equal 0 ? */
    collisionIndicator = 0;
    return 0;                    /* no intersection occurs */
  }

  /* Compute the direction of the intersection line */
  cross(D,N1,N2);

  /* Compute and index to the largest component of D */
  max = fabs(D[0]);
  index = 0;
  b = fabs(D[1]);
  c = fabs(D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* Simplified projection onto L*/
  vp[0]=V0[index]; vp[1]=V1[index]; vp[2]=V2[index];
  up[0]=U0[index]; up[1]=U1[index]; up[2]=U2[index];
  
  double isectpointA1[3],isectpointA2[3], isectpointB1[3],isectpointB2[3];
  double isectpt1[3], isectpt2[3];
  int smallest1,smallest2, coplanar;
  bool coplanarCollision;
  
  /* compute interval for triangle 1 */
  coplanar = compute_intervals_isectline(1,vp,dv,dv0dv1,dv0dv2,
					 &isect1[0],&isect1[1],isectpointA1,isectpointA2);

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
  compute_intervals_isectline(2,up,du,du0du1,du0du2,
  			      &isect2[0],&isect2[1],isectpointB1,isectpointB2);

  sort2(isect1[0], isect1[1], smallest1);
  sort2(isect2[0], isect2[1], smallest2);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) {
    collisionIndicator = 0;
    return 0;
  }

  /* at this point, we know that the triangles intersect */

  if(isect2[0]<isect1[0])
  {
    if(smallest1==0) {
      set(isectpt1,isectpointA1); 
    }
    else {
      set(isectpt1,isectpointA2); 
    }
    if(isect2[1]<isect1[1]) {
      if(smallest2==0) {
  	set(isectpt2,isectpointB2);
      }
      else {
  	set(isectpt2,isectpointB1);
      }
    }
    else {
      if(smallest1==0) {
  	set(isectpt2,isectpointA2);
      }
      else {
  	set(isectpt2,isectpointA1); 
      }
    }
  }
  else  {
    if(smallest2==0) {
      set(isectpt1,isectpointB1);
    }
    else {
      set(isectpt1,isectpointB2); 
    }
    if(isect2[1]>isect1[1])
    {
      if(smallest1==0) {
  	set(isectpt2,isectpointA2); 
      }
      else {
  	set(isectpt2,isectpointA1); 
      }      
    }
    else {
      if(smallest2==0) {
  	set(isectpt2,isectpointB2);
      }
      else {
  	set(isectpt2,isectpointB1);
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


int CollisionTwoTriangles::
compute_intervals_isectline(int id, double VV[3], double DD[3],
			    double D0D1,double D0D2,double *isect0, double *isect1,
			    double isectpoint0[3],double isectpoint1[3])
{
  double VERT0[3], VERT1[3], VERT2[3];
  if (id==1){
    set(VERT0,V0); set(VERT1,V1); set(VERT2,V2);
  }
  else if (id==2){
    set(VERT0,U0); set(VERT1,U1); set(VERT2,U2);
  }

  if(D0D1>0.0f)                                        
  {                                                    
    /* here we know that D0D2<=0.0 */                  
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(VERT2, VERT0, VERT1, VV[2], VV[0], VV[1], DD[2], DD[0], DD[1],
	   isect0, isect1, isectpoint0, isectpoint1);
  } 
  else if(D0D2>0.0f)                                   
    {                                                   
    /* here we know that d0d1<=0.0 */             
    isect2(VERT1,VERT0,VERT2,VV[1],VV[0],VV[2],DD[1],DD[0],DD[2],isect0,isect1,isectpoint0,isectpoint1);
  }                                                  
  else if(DD[1]*DD[2]>0.0f || DD[0]!=0.0f)   
  {                                   
    /* here we know that d0d1<=0.0 or that D0!=0.0 */
    isect2(VERT0,VERT1,VERT2,VV[0],VV[1],VV[2],DD[0],DD[1],DD[2],isect0,isect1,isectpoint0,isectpoint1);   
  }                                                  
  else if(DD[1]!=0.0f)                                  
  {                                               
    isect2(VERT1,VERT0,VERT2,VV[1],VV[0],VV[2],DD[1],DD[0],DD[2],isect0,isect1,isectpoint0,isectpoint1); 
  }                                         
  else if(DD[2]!=0.0f)                                  
  {                                                   
    isect2(VERT2,VERT0,VERT1,VV[2],VV[0],VV[1],DD[2],DD[0],DD[1],isect0,isect1,isectpoint0,isectpoint1);     
  }                                                 
  else                                               
  {                                                   
    /* triangles are coplanar */    
    return 1;
  }
  return 0;
}

void CollisionTwoTriangles::
isect2(double VTX0[3], double VTX1[3], double VTX2[3], double VV0, double VV1, double VV2,
       double D0,double D1,double D2,double *isect0,double *isect1,double isectpoint0[3],double isectpoint1[3]) 
{
  double tmp;
  double diff[3];

  tmp = D0/(D0-D1);
  *isect0 = VV0 + (VV1-VV0)*tmp;
  sub(diff,VTX1,VTX0);
  mult(diff,diff,tmp);
  add(isectpoint0,diff,VTX0);

  tmp = D0/(D0-D2);           
  *isect1 = VV0 + (VV2-VV0)*tmp;
  sub(diff,VTX2,VTX0);          
  mult(diff,diff,tmp);
  add(isectpoint1,VTX0,diff);
}


void CollisionTwoTriangles::
computeIntervals(double VV[3], double DD[3], double D0D1, double D0D2,
		 double t[2], bool &coplanar)
{
  coplanar = false;
  if(D0D1 > 0.0f) {
    /* Here we know that D0D2<=0.0 
       that is D0, D1 are on the same side, D2 on the other or on the plane */
    find_t(VV[2], VV[0], VV[1], DD[2], DD[0], DD[1], t);
  }
  else if(D0D2 > 0.0f) {
    /* Here we know that d0d1<=0.0: D0,D2 on the same side */
    find_t(VV[1], VV[0], VV[2], DD[1], DD[0], DD[2], t);
  }
  else if( DD[1]*DD[2]>0.0f || DD[0]!=0.0f ) {
    /* Here we know that d0d1<=0.0 or that D0!=0.0 */
    find_t(VV[0], VV[1], VV[2], DD[0], DD[1], DD[2], t);
  }
  else if(D[1] != 0.0f) {
    find_t(VV[1], VV[0], VV[2], DD[1], DD[0], DD[2], t);
  }
  else if(D[2] != 0.0f) {
    find_t(VV[2], VV[0], VV[1], D[2], D[0], D[1], t);
  }
  else {
    /* Triangles are coplanar */
    coplanar = true;
  }
}

/* --- Find the parameters 't' in the line that represent the intersection
       extremes of L with the triangle: L = t*D                            --- */
void CollisionTwoTriangles::
find_t(double VV0, double VV1, double VV2, double D0, double D1, double D2, double t[2])
{    
  t[0] = VV0 + (VV1-VV0)*D0/(D0-D1);
  t[1] = VV0 + (VV2-VV0)*D0/(D0-D2);
}


bool CollisionTwoTriangles::
coplanar_tri_tri(double N[3])
{
   double absN[3];
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
   absN[0]=fabs(N[0]);  absN[1]=fabs(N[1]);  absN[2]=fabs(N[2]);
   if( absN[0]>absN[1] ) {
     if( absN[0]>absN[2] ) {
       i0=1; i1=2;  inot=0;     /* absN[0] is greatest */
     }
     else {
       i0=0; i1=1;  inot=2;     /* absN[2] is greatest */
     }
   }
   else {                       /* absN[0]<=absN[1] */
     if( absN[2]>absN[1] ) {
       i0=0; i1=1;  inot=2;     /* absN[2] is greatest */
     }
     else {
       i0=0; i1=2;  inot=1;     /* absN[1] is greatest */
     }
   }               

   bool collision[9];
   /* Test all edges of triangle 1 against the edges of triangle 2 */
   collision[0] = edge_against_tri_edges(V0, V1);
   collision[1] = edge_against_tri_edges(V1, V2);
   collision[2] = edge_against_tri_edges(V2, V0);
   
   /* Finally, test if tri1 is totally contained in tri2 or vice versa */
   /*  check if one point is contained and if so, assume all of it is contained */
   // collision[3] = point_in_tri(V0, U0, U1, U2);
   // collision[4] = point_in_tri(U0, V0, V1, V2);
   collision[3] = point_in_tri(V0, U0, U1, U2);
   collision[4] = point_in_tri(V1, U0, U1, U2);
   collision[5] = point_in_tri(V2, U0, U1, U2);
   collision[6] = point_in_tri(U0, V0, V1, V2);
   collision[7] = point_in_tri(U1, V0, V1, V2);
   collision[8] = point_in_tri(U2, V0, V1, V2);
   
   return (collision[0] || collision[1] || collision[2] || collision[3] || collision[6]);
}


/* --- Test if an edge of T1 intersects the edges of T2 -- */
bool CollisionTwoTriangles::
edge_against_tri_edges(double v0[3], double v1[3])
{
  double Ax, Ay;
  bool collision[3];       // True if there is collision
  // Components of the edge in T1
  Ax = v1[i0] - v0[i0];    
  Ay = v1[i1] - v0[i1];    
  /* test edge U0,U1 against v0,v1 */
  //std::cout << "\nEdge: "; printVector(v0); std::cout << " , "; printVector(v1); std::cout << "with:";

  //std::cout << "\n - edge "; printVector(U0); std::cout << " , "; printVector(U1);
  collision[0] = edge_edge_test(v0, U0, U1, Ax, Ay);
  //std::cout << ": result is " << collision[0] << std::endl;

  /* test edge U1,U2 against v0,v1 */
  //std::cout << "\n - edge "; printVector(U1); std::cout  << " , "; printVector(U2);
  collision[1] = edge_edge_test(v0, U1, U2, Ax, Ay);
  //std::cout << ": result is " << collision[1] << std::endl;

  /* test edge U2,U1 against v0,v1 */
  //std::cout << "\n - edge "; printVector(U2); std::cout  << " , "; printVector(U0);
  collision[2] = edge_edge_test(v0, U2, U0, Ax, Ay);
  //std::cout << ": result is " << collision[2] << std::endl;

  return (collision[0] || collision[1] || collision[2]);
}

/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
bool CollisionTwoTriangles::
edge_edge_test(double v0[3], double u0[3], double u1[3], double Ax, double Ay)
{
  double Bx, By, Cx, Cy, numu, numv, den;
  double alpha, beta, Pintu[3], Pintv[3];
  Bx = u0[i0] - u1[i0];
  By = u0[i1] - u1[i1];
  Cx = v0[i0] - u0[i0];
  Cy = v0[i1] - u0[i1];
  den  = Ay*Bx - Ax*By;   // denominator
  numv  = By*Cx - Bx*Cy;   // numerator 1
  if((den>0 && numv>=0 && numv<=den) || (den<0 && numv<=0 && numv>=den)) {
    numu = Ax*Cy-Ay*Cx;      // numerator 2
    if(den>0) {
      if(numu>=0 && numu<=den) {
	beta = numu/den;
	Pintu[i0] = u0[i0] + beta*(u1[i0]-u0[i0]);
	Pintu[i1] = u0[i1] + beta*(u1[i1]-u0[i1]);
	Pintu[inot] = u0[inot] + beta*(u1[inot]-u0[inot]);
	//std::cout << " Intersection at: ( " << Pintu[0] << ", " << Pintu[1] << ", " << Pintu[2] << ")";
	pushVectorAsPoint(Pintu, pointsTT);
	return (true);
      }
    }
    else {
      if(numu<=0 && numu>=den){
	beta = numu/den;
	Pintu[i0] = u0[i0] + beta*(u1[i0]-u0[i0]);
	Pintu[i1] = u0[i1] + beta*(u1[i1]-u0[i1]);
	Pintu[inot] = u0[inot] + beta*(u1[inot]-u0[inot]);
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
point_in_tri(double v0[3], double u0[3], double u1[3], double u2[3])
{
  double a,b,c,d0,d1,d2;
  /* is T1 completly inside T2? */
  /* check if v0 is inside tri(u0,u1,u2) */
  a  = u1[i1] - u0[i1];
  b  = -(u1[i0]-u0[i0]);
  c  = -a*u0[i0] - b*u0[i1];
  d0 = a*v0[i0] + b*v0[i1] + c;

  a  = u2[i1] - u1[i1];
  b  = -(u2[i0]-u1[i0]);
  c  = -a*u1[i0] - b*u1[i1];
  d1 = a*v0[i0] + b*v0[i1] + c;

  a  = u0[i1] - u2[i1];
  b  = -(u0[i0]-u2[i0]);
  c  = -a*u2[i0] - b*u2[i1];
  d2 = a*v0[i0] + b*v0[i1] + c;

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
  std::cout << "Vertices: ";
  std::cout << "\n  - Triangle 1: "; printVector(V0); printVector(V1); printVector(V2);
  std::cout << "\n  - Triangle 2: "; printVector(U0); printVector(U1); printVector(U2);
  std::cout << std::endl;
}

/* --- Print verbose information relative to the collision --- */
void CollisionTwoTriangles::
printTTcollisionInformation( void )
{
  if (collisionIndicator) {
    std::cout << "Result: Triangles are intersecting" << std::endl;
    std::cout << " - Intersections occur at: " << std::endl;
    for (int i=0; i<pointsTT.size(); i++){
      std::cout << "     "; pointsTT[i].print();
      std::cout << std::endl;
    }
  }
  else{
    std::cout << "Result: Triangles are NOT intersecting" << std::endl;
  }
}
