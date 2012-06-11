/*
 * Copyright 2012, Oscar E. Ramos Ponce, LAAS-CNRS
 *
 * This file is part of sot-dyninv.
 * sot-dyninv is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * sot-dyninv is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.  You should
 * have received a copy of the GNU Lesser General Public License along
 * with sot-dyninv.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <sot-dyninv/box-box-collision-detector.h>
#include <sot-dyninv/commands-helper.h>
#include <sot/core/debug.hh>
#include <dynamic-graph/factory.h>

#include <vector>

namespace dynamicgraph
{
  namespace sot
  {
    namespace dyninv
    {

      namespace dg = ::dynamicgraph;
      using namespace dg;
      using dg::SignalBase;

      /* --- DG FACTORY ------------------------------------------------------- */
      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN(BoxBoxCollisionDetector,"BoxBoxCollisionDetector");

      /* ---------------------------------------------------------------------- */
      /* --- CONSTRUCTION ----------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      BoxBoxCollisionDetector::
      BoxBoxCollisionDetector( const std::string & name )
	: Entity(name)
	  ,scene()
	,CONSTRUCT_SIGNAL_OUT(tolerance, double, sotNOSIGNAL)
	,CONSTRUCT_SIGNAL_OUT(inContact, int, sotNOSIGNAL)
	,CONSTRUCT_SIGNAL_OUT(contactPoints, ml::Matrix,
			      inContactSOUT)
	,CONSTRUCT_SIGNAL_OUT(transformationOut1, ml::Matrix, sotNOSIGNAL)
	,CONSTRUCT_SIGNAL_OUT(transformationOut2, ml::Matrix, sotNOSIGNAL)	  
	,CONSTRUCT_SIGNAL_OUT(verticesB1, ml::Matrix, sotNOSIGNAL)
	,CONSTRUCT_SIGNAL_OUT(verticesB2, ml::Matrix, sotNOSIGNAL)

      {
	signalRegistration(toleranceSOUT << inContactSOUT << contactPointsSOUT
			   << transformationOut1SOUT << transformationOut2SOUT
			   << verticesB1SOUT << verticesB2SOUT ); 
	toleranceSOUT.setNeedUpdateFromAllChildren( true );
	inContactSOUT.setNeedUpdateFromAllChildren( true );
	transformationOut1SOUT.setNeedUpdateFromAllChildren( true );
	transformationOut2SOUT.setNeedUpdateFromAllChildren( true );
	verticesB1SOUT.setNeedUpdateFromAllChildren( true );
	verticesB2SOUT.setNeedUpdateFromAllChildren( true );
	// contactPointsSOUT.setNeedUpdateFromAllChildren( true );
	initCommands();
      }

      void BoxBoxCollisionDetector::
      initCommands( void )
      {
	using namespace dynamicgraph::command;

	addCommand("setTolerance",
		   makeCommandVoid1(*this, &BoxBoxCollisionDetector::setTolerance,
				    docCommandVoid1("Set the tolerance used to mix closed points.",
						    "tolerance (double)")));

	addCommand("addBox1",
		   makeCommandVoid2(*this, &BoxBoxCollisionDetector::addBox1,
				    docCommandVoid2("Add a fixed box that cannot move (many per scene).",
						    "Lenghts in x, y, z (vector)",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("addBox2",
		   makeCommandVoid2(*this, &BoxBoxCollisionDetector::addBox2,
				    docCommandVoid2("Add a box that can move (note: only ONE per scene!).",
						    "Lenghts in x, y, z (vector)",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("setTransformationB1",
		   makeCommandVoid1(*this, &BoxBoxCollisionDetector::setTransformationB1,
				    docCommandVoid1("Set the position and orientation of the moving box.",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("setTransformationB2",
		   makeCommandVoid1(*this, &BoxBoxCollisionDetector::setTransformationB2,
				    docCommandVoid1("Set the position and orientation of the moving box.",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("printInformation",
		   makeCommandVoid0(*this, &BoxBoxCollisionDetector::cmd_printInformation,
				    docCommandVoid0("Print information about the contacts.")));

      }
      
      /* ---------------------------------------------------------------------- */
      /* --- COMMANDS --------------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      void BoxBoxCollisionDetector::
      setTolerance( const double & toleranceValue )
      {
	scene.setTolerance(toleranceValue);
      }

      void BoxBoxCollisionDetector::
      addBox1( const ml::Vector & len,
	       const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;
	setRTfromM(R, T, M);
	scene.setLengthB1(len(0), len(1), len(2));
	scene.setTransformation1(R, T);

	if(false){
	  printRotTrans(R, T);
	  printLengths(len);
	}
      }

      void BoxBoxCollisionDetector::
      addBox2( const ml::Vector & len,
	       const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;
	setRTfromM(R, T, M);
	scene.setLengthB2(len(0), len(1), len(2));
	scene.setTransformation2(R, T);

	if(false){
	  printRotTrans(R, T);
	  printLengths(len);
	}
      }

      void BoxBoxCollisionDetector::
      setTransformationB1( const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;
	setRTfromM(R, T, M);
	scene.setTransformation1(R, T);

	if(false)
	  printRotTrans(R, T);
      }

      void BoxBoxCollisionDetector::
      setTransformationB2( const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;
	setRTfromM(R, T, M);
	scene.setTransformation2(R, T);

	if(false)
	  printRotTrans(R, T);
      }

      void BoxBoxCollisionDetector::
      cmd_printInformation( void )
      {
	scene.printLengths();
	scene.printBoxesVertices();
      }


      /* ---------------------------------------------------------------------- */
      /* --- SIGNALS ---------------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      double& BoxBoxCollisionDetector::
      toleranceSOUT_function( double& tol, int t)
      {
	tol = scene.getTolerance();
	return tol;
      }


      int& BoxBoxCollisionDetector::
      inContactSOUT_function( int & collisionDetected, int t )
      {
	collisionDetected = scene.computeBBintersections();
      	return collisionDetected;
      }


      ml::Matrix& BoxBoxCollisionDetector::
      contactPointsSOUT_function( ml::Matrix& mlpoints, int t )
      {
	const int & collision = inContactSOUT(t);
	int numContactPoints = scene.pointsBB.size();
	mlpoints.resize(numContactPoints, 3);

	if (collision) {
	  for (int i=0; i<numContactPoints; i++) {
	    mlpoints(i,0) = scene.pointsBB[i](0);
	    mlpoints(i,1) = scene.pointsBB[i](1);
	    mlpoints(i,2) = scene.pointsBB[i](2);
	  }
	}
	return mlpoints;
      }


      ml::Matrix& BoxBoxCollisionDetector::
      transformationOut1SOUT_function( ml::Matrix& mlM1out, int t )
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;

	mlM1out.resize(4,4);
	scene.getTransformation1(R, T);
	setMfromRT(mlM1out, R, T);
	return mlM1out;
      }


      ml::Matrix& BoxBoxCollisionDetector::
      transformationOut2SOUT_function( ml::Matrix& mlM2out, int t )
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;
	scene.getTransformation2(R, T);

	mlM2out.resize(4,4);
	setMfromRT(mlM2out, R, T);
	return mlM2out;
      }


      ml::Matrix& BoxBoxCollisionDetector::
      verticesB1SOUT_function( ml::Matrix& mlVert, int t )
      {
	Eigen::Vector3d v1, v2, v3, v4, v5, v6, v7, v8;
	scene.getVerticesB1(v1, v2, v3, v4, v5, v6, v7, v8);

	mlVert.resize(8,3);
	for (int i=0; i<3; i++)
	  {
	    mlVert(0,i) = v1(i);
	    mlVert(1,i) = v2(i);
	    mlVert(2,i) = v3(i);
	    mlVert(3,i) = v4(i);
	    mlVert(4,i) = v5(i);
	    mlVert(5,i) = v6(i);
	    mlVert(6,i) = v7(i);
	    mlVert(7,i) = v8(i);
	  }
	return mlVert;
      }


      ml::Matrix& BoxBoxCollisionDetector::
      verticesB2SOUT_function( ml::Matrix& mlVert, int t )
      {
	Eigen::Vector3d v1, v2, v3, v4, v5, v6, v7, v8;
	scene.getVerticesB2(v1, v2, v3, v4, v5, v6, v7, v8);

	mlVert.resize(8,3);
	for (int i=0; i<3; i++)
	  {
	    mlVert(0,i) = v1(i);
	    mlVert(1,i) = v2(i);
	    mlVert(2,i) = v3(i);
	    mlVert(3,i) = v4(i);
	    mlVert(4,i) = v5(i);
	    mlVert(5,i) = v6(i);
	    mlVert(6,i) = v7(i);
	    mlVert(7,i) = v8(i);
	  }
	return mlVert;
      }


      /* ---------------------------------------------------------------------- */
      /* --- ENTITY ----------------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      void BoxBoxCollisionDetector::
      display( std::ostream& os ) const
      {
	os << "BoxBoxCollisionDetector "<<getName() << "." << std::endl;
      }

      void BoxBoxCollisionDetector::
      commandLine( const std::string& cmdLine,
		   std::istringstream& cmdArgs,
		   std::ostream& os )
      {
	if( cmdLine == "help" )
	  {
	    os << "BoxBoxCollisionDetector:\n"
	       << "\t- ." << std::endl;
	    Entity::commandLine( cmdLine,cmdArgs,os );
	  }
	else
	  {
	    Entity::commandLine( cmdLine,cmdArgs,os );
	  }
      }

      /* ---------------------------------------------------------------------- */
      /* --- HELPER FUNCTIONS ------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      void BoxBoxCollisionDetector::
      setRTfromM(Eigen::Matrix3d &R, Eigen::Vector3d &T, const ml::Matrix & M)
      {
	R(0,0) = M(0,0); R(0,1) = M(0,1); R(0,2) = M(0,2); 
	R(1,0) = M(1,0); R(1,1) = M(1,1); R(1,2) = M(1,2);
	R(2,0) = M(2,0); R(2,1) = M(2,1); R(2,2) = M(2,2);
	T(0) = M(0,3); T(1) = M(1,3); T(2) = M(2,3);
      }

      void BoxBoxCollisionDetector::
      setMfromRT(ml::Matrix & mlM, Eigen::Matrix3d R, Eigen::Vector3d T)
      {
	mlM(0,0) = R(0,0); mlM(0,1) = R(0,1); mlM(0,2) = R(0,2); mlM(0,3) = T(0);
	mlM(1,0) = R(1,0); mlM(1,1) = R(1,1); mlM(1,2) = R(1,2); mlM(1,3) = T(1);
	mlM(2,0) = R(2,0); mlM(2,1) = R(2,1); mlM(2,2) = R(2,2); mlM(2,3) = T(2);
	mlM(3,0) =    0.0; mlM(3,1) =    0.0; mlM(3,2) =    0.0; mlM(3,3) =  1.0;
      }

      void BoxBoxCollisionDetector::
      printRotTrans(Eigen::Matrix3d R, Eigen::Vector3d T )
      {
	std::cout << "Rotation: " << std::endl;
	std::cout << R << std::endl;
	std::cout << "Translation: " << std::endl;
	std::cout << T.transpose() << std::endl;
      }
     
      void BoxBoxCollisionDetector::
      printLengths(const ml::Vector & length)
      {
	std::cout << "Lengths: "
	   	  << "Lx="   << length(0) << ", Ly=" << length(1)
		  << ", Lz=" << length(2) << std::endl;
      }

    } // namespace dyninv
  } // namespace sot
} // namespace dynamicgraph

