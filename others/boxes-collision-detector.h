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

#ifndef __sot_dyninv_BoxesCollisionDetector_H__
#define __sot_dyninv_BoxesCollisionDetector_H__
/* --------------------------------------------------------------------- */
/* --- API ------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

#if defined (WIN32)
#  if defined (dynamic_interpretor_EXPORTS)
#    define SOTBOXESCOLLISIONDETECTOR_EXPORT __declspec(dllexport)
#  else
#    define SOTBOXESCOLLISIONDETECTOR_EXPORT __declspec(dllimport)
#  endif
#else
#  define SOTBOXESCOLLISIONDETECTOR_EXPORT
#endif

/* ------------------------------------------------------------------
--- INCLUDE ---------------------------------------------------------
--------------------------------------------------------------------- */

// SOT
#include <sot-dyninv/signal-helper.h>
#include <sot-dyninv/entity-helper.h>
#include <collision-one-many-boxes.h>
#include <Eigen/Dense>


namespace dynamicgraph {
  namespace sot {
    namespace dyninv {

      /* ---------------------------------------------------------------------
         --- CLASS -----------------------------------------------------------
         --------------------------------------------------------------------- */

      class SOTBOXESCOLLISIONDETECTOR_EXPORT BoxesCollisionDetector
	:public ::dynamicgraph::Entity
	,public ::dynamicgraph::EntityHelper<BoxesCollisionDetector>
	{
	  
	public:
	  /* --- CONSTRUCTOR ---- */
	  BoxesCollisionDetector( const std::string & name );

	  /* --- ENTITY INHERITANCE --- */
	  
	  static const std::string CLASS_NAME;
	  virtual void display( std::ostream& os ) const;
	  virtual const std::string& getClassName( void ) const {return CLASS_NAME;}
	  virtual void commandLine( const std::string& cmdLine,
				    std::istringstream& cmdArgs,
				    std::ostream& os );
	  void initCommands( void );

	  /* --- SIGNALS --- */
	  DECLARE_SIGNAL_OUT(tolerance, double);
	  DECLARE_SIGNAL_OUT(inContact, int);
	  DECLARE_SIGNAL_OUT(contactPoints, ml::Matrix);

	  /* --- COMMANDS --- */
	  void setTolerance( const double & toleranceValue);
	  void addFixedBox( const ml::Vector & lengths,
			    const ml::Matrix & M);
	  void addMovingBox( const ml::Vector & lengths,
			     const ml::Matrix & M);
	  void poseMovingBox( const ml::Matrix & M);
	  void cmd_printInformation( void );

	private:
	  /* --- Variables --- */
	  CollisionOneManyBoxes scene;     // Scene that checks collision from one box to many boxes */
	  
	  /* --- Helper Printing functions --- */
	  void printRotTrans(Eigen::Matrix3d R, Eigen::Vector3d T);
	  void printLengths(const ml::Vector & length);
	  void setRTfromM(Eigen::Matrix3d &R, Eigen::Vector3d &T, const ml::Matrix & M);

	}; // class BoxesCollisionDetector

    } // namespace dyninv
  } // namespace sot
} // namespace dynamicgraph



#endif // #ifndef __sot_dyninv_BoxesCollisionDetector_H__
