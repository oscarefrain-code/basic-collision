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


#include <sot-dyninv/boxes-collision-detector.h>
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
      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN(BoxesCollisionDetector,"BoxesCollisionDetector");

      /* ---------------------------------------------------------------------- */
      /* --- CONSTRUCTION ----------------------------------------------------- */
      /* ---------------------------------------------------------------------- */
      BoxesCollisionDetector::
      BoxesCollisionDetector( const std::string & name )
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
	initCommands();
      }

      void BoxesCollisionDetector::
      initCommands( void )
      {
	using namespace dynamicgraph::command;

	addCommand("setTolerance",
		   makeCommandVoid1(*this, &BoxesCollisionDetector::setTolerance,
				    docCommandVoid1("Set the tolerance used to mix closed points.",
						    "tolerance (double)")));

	addCommand("addFixedBox",
		   makeCommandVoid2(*this, &BoxesCollisionDetector::addFixedBox,
				    docCommandVoid2("Add a fixed box that cannot move (many per scene).",
						    "Lenghts in x, y, z (vector)",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("addMovingBox",
		   makeCommandVoid2(*this, &BoxesCollisionDetector::addMovingBox,
				    docCommandVoid2("Add a box that can move (note: only ONE per scene!).",
						    "Lenghts in x, y, z (vector)",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("poseMovingBox",
		   makeCommandVoid1(*this, &BoxesCollisionDetector::poseMovingBox,
				    docCommandVoid1("Set the position and orientation of the moving box.",
						    "Homogeneous transformation matrix (4x4)")));

	addCommand("printInformation",
		   makeCommandVoid0(*this, &BoxesCollisionDetector::cmd_printInformation,
				    docCommandVoid0("Print information about the contacts.")));

      }
      
      /* ---------------------------------------------------------------------- */
      /* --- COMMANDS --------------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      void BoxesCollisionDetector::
      setTolerance( const double & toleranceValue )
      {
	scene.setTolerance(toleranceValue);
      }

      void BoxesCollisionDetector::
      addFixedBox( const ml::Vector & len,
		   const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;   // Rotation and translation matrices
	setRTfromM(R, T, M);
	scene.addBox(len(0), len(1), len(2), R, T);

	if(false){
	  printRotTrans(R, T);
	  printLengths(len);
	}
      }

      void BoxesCollisionDetector::
      addMovingBox( const ml::Vector & len,
		    const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;   // Rotation and translation matrices
	setRTfromM(R, T, M);
	scene.setLengthB1(len(0), len(1), len(2));
	scene.setTransformationB1(R, T);

	if(false){
	  printRotTrans(R, T);
	  printLengths(len);
	}
      }

      void BoxesCollisionDetector::
      poseMovingBox( const ml::Matrix & M)
      {
	Eigen::Matrix3d R;
	Eigen::Vector3d T;   // Rotation and translation matrices
	setRTfromM(R, T, M);
	scene.setTransformationB1(R, T);

	if(false)
	  printRotTrans(R, T);
      }

      void BoxesCollisionDetector::
      cmd_printInformation( void )
      {
	scene.printLengths();
	scene.printBoxesVertices();
      }


      /* ---------------------------------------------------------------------- */
      /* --- SIGNALS ---------------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      double& BoxesCollisionDetector::
      toleranceSOUT_function( double& tol, int t)
      {
	tol = scene.getTolerance();
	return tol;
      }


      int& BoxesCollisionDetector::
      inContactSOUT_function( int & collisionDetected, int t )
      {
	collisionDetected = scene.computeIntersections();
      	return collisionDetected;
      }


      ml::Matrix& BoxesCollisionDetector::
      contactPointsSOUT_function( ml::Matrix& mlpoints, int t )
      {
	const int & collision = inContactSOUT(t);
	int numContactPoints = scene.pointsBBs.size();
	mlpoints.resize(numContactPoints, 3);

	if (collision) {
	  for (int i=0; i<numContactPoints; i++) {
	    mlpoints(i,0) = scene.pointsBBs[i](0);
	    mlpoints(i,1) = scene.pointsBBs[i](1);
	    mlpoints(i,2) = scene.pointsBBs[i](2);
	  }
	}
      	return mlpoints;
      }


      /* ---------------------------------------------------------------------- */
      /* --- ENTITY ----------------------------------------------------------- */
      /* ---------------------------------------------------------------------- */

      void BoxesCollisionDetector::
      display( std::ostream& os ) const
      {
	os << "BoxesCollisionDetector "<<getName() << "." << std::endl;
      }

      void BoxesCollisionDetector::
      commandLine( const std::string& cmdLine,
		   std::istringstream& cmdArgs,
		   std::ostream& os )
      {
	if( cmdLine == "help" )
	  {
	    os << "BoxesCollisionDetector:\n"
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

      void BoxesCollisionDetector::
      setRTfromM(Eigen::Matrix3d &R, Eigen::Vector3d &T, const ml::Matrix & M)
      {
	R(0,0) = M(0,0); R(0,1) = M(0,1); R(0,2) = M(0,2); 
	R(1,0) = M(1,0); R(1,1) = M(1,1); R(1,2) = M(1,2);
	R(2,0) = M(2,0); R(2,1) = M(2,1); R(2,2) = M(2,2);
	T(0) = M(0,3); T(1) = M(1,3); T(2) = M(2,3);
      }

      void BoxesCollisionDetector::
      printRotTrans(Eigen::Matrix3d R, Eigen::Vector3d T )
      {
	std::cout << "Rotation: " << std::endl;
	std::cout << R << std::endl;
	std::cout << "Translation: " << std::endl;
	std::cout << T.transpose() << std::endl;
      }
     
      void BoxesCollisionDetector::
      printLengths(const ml::Vector & length)
      {
	std::cout << "Lengths: "
	   	  << "Lx="   << length(0) << ", Ly=" << length(1)
		  << ", Lz=" << length(2) << std::endl;
      }

    } // namespace dyninv
  } // namespace sot
} // namespace dynamicgraph

