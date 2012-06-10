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
      {
	signalRegistration(toleranceSOUT << inContactSOUT << contactPointsSOUT ); 
	toleranceSOUT.setNeedUpdateFromAllChildren( true );
	inContactSOUT.setNeedUpdateFromAllChildren( true );
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
	//double R[9], T[3];   // Rotation and translation matrices
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
	//double R[9], T[3];   // Rotation and translation matrices
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
	//double R[9], T[3];   // Rotation and translation matrices
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
	//double R[9], T[3];   // Rotation and translation matrices
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

