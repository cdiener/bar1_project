/*
 * instationary.cc
 * This file is part of ExCeCo
 *
 * Copyright (C) 2011 - Christian Diener, Wolfgang Giese
 *
 * ExCeCo is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * ExCeCo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ExCeCo. If not, see <http://www.gnu.org/licenses/>.
 */
 


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>

#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>
#include<dune/istl/superlu.hh>
#include<dune/istl/supermatrix.hh>


#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/gridoperatorspace/instationarygridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
//#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#define VTK true
#define USE_SUPER_LU true
#define ORDER  1  //in case =0 we have a finite volume method (special case of DG-Method)
#define DG_METHOD  1  //choice of the method, e.g. 0 for OBB, 1 for NIPG or 2 for SIPG

const int BLOCKSIZE2DP1 = 3;
const int BLOCKSIZE2DP2 = 6;
const int BLOCKSIZE2DP3 = 10;
const int BLOCKSIZE3DP1 = 4;
const int BLOCKSIZE3DP2 = 10;
const int BLOCKSIZE3DP3 = 20;

//***************************************
// Standard units used:
// 		Concentrations: nanomolar [nm]
//		Time:		seconds [s]
//		Space:		micrometer [um]
//**************************************/

/** \var  	N_ALPHA    
 *  \brief 	number of MATalpha cells
 */
/**  \var  	N_CHEAT    	
 *  \brief 	number of cheaters
 */ 
/**  \var  	N_MATA     	
 *  \brief 	number of MATa cells
 */

static unsigned int N_ALPHA = 1;
static unsigned int N_MATA = 2;
static unsigned int N_CHEAT = 0;

/** \var D0 
 * \brief diffusion rate for alpha-factor
 *
 * \var D1 
 * \brief diffusion rate for Bar1 
 *
 * \var KBAR0 
 * \brief basic secretion of Bar1; independent of excitment by alpha-factor
 *
 * \var KBAR1
 * \brief maximal secretion of Bar1 dependent of excitment by alpha-factor
 *
 * \var BAR_DIR 
 * \brief Bar1 level on the outer boundary, that is defined by a Dirichlet boundary condition
 *
 * \var ALPHA_SEC
 * \brief secretion rate of alpha-factor by MATalpha cells in nmol/(l*s)
 *
 * \var K_M 
 * \brief Michelis-Menten constant that regulates the secretion of Bar1 dependet on alpha-factor stimulation
 *
 * \var N_H
 * \brief The Hill coefficient used for the Bar1 induction
 *
 * \var ALPHA_SOURCE 
 * \brief basic generation of alpha-factor in the computational domain 
 *
 * \var absorbtionAlpha 
 * \brief for each time step the total integrated alpha-factor-absorption of each MATa call is stored in this vector
 */

std::vector< std::vector<double> > coordsOfCells;

/*static double center1_x = 0.0;
static double center1_y = 0.0;
static double radius1 = 4.0;

static double center2_x = 15.0;
static double center2_y = 0.0;
static double radius2 = 3.5;

static double center3_x = 0.0;
static double center3_y = 25.0;
static double radius3 = 4.5;*/

static double max_radius = 5.0;
static double dish_radius = 170.0;

static double elementSize = 1.0;
static int elements2d = 100;

static int Integration_steps = 10;
static double Integration_height = 20.0;
static int QORDER = 10;

static double D0 = 324.03; 
static double D1 = 93.38;
static double ALPHA_SEC = 100.0;
static double KBAR0 = 5.6;
static double KBARI = 1.4;
static double K_M = 30.41;
static double N_H = 1.62;
static double ALPHA_SOURCE = 0.0;
static double BAR1_DIR = 0.28;
static double KAPPA = 1.0;
static std::vector<std::vector<double>> absorptionAlpha;
std::ofstream datei;



/** \var T_DELAY
 * \brief time delay for Bar1-Response
 *
 */

static double T_DELAY = 300.0; 



#include"alphaAbsorption.hh"
#include"bar1Boundary.hh"
#include"bc/alpha_bc.hh"
#include"bc/bar1p_bc.hh"
#include"bc/diffusion_tensor_alpha.hh"
#include"bc/diffusion_tensor_bar1p.hh"
#include"diffusiondg.hh"
#include"toperator.hh"
#include"driver_dg.hh"
#include"concentrationChange.hh"

//===============================================================
// Main program with grid setup
//===============================================================

 /**
 * \brief 	This function creates the mesh and starts the driver for the dg approximation
 *
 * \fn driver
 * 
 * \param	meshfile	name of the grid file
 * \param	level		initial refinement of the grid
 * \param	order		order of polynomials for approximation
 * \param	dg_method	choice of the method, e.g. 0 for OBB, 1 for NIPG or 2 for SIPG
 * \param	dt		lenght of time steps
 * \param	tend		end time
 * 
 */

template<typename GridType>
void driver(std::string const& meshfile, std::string const& output, int level, int order, int dg_method, double dt, double tend) {
   
   //UG Grid
   int dim = GridType::dimension;
   GridType grid(400);

   std::vector<int> boundaryIndexToPhysicalEntity;
   std::vector<int> elementIndexToPhysicalEntity;

   //read a gmsh file                                                                                                                                                                                                                              
   Dune::GmshReader<GridType> gmshreader;
   gmshreader.read(grid, meshfile, boundaryIndexToPhysicalEntity, elementIndexToPhysicalEntity, true, false);

   if(level>0)
     grid.globalRefine(level);
   
   if(!Dune::MPIHelper::isFake)
     grid.loadBalance();
    
 
   // get view
   typedef typename GridType::LeafGridView GV;
   const GV& gv = grid.leafView();

   // call driver
   switch(order){
     case 1:
       if(dim==2) driver_dg<GV,1,BLOCKSIZE2DP1>(gv, output, boundaryIndexToPhysicalEntity, dg_method, dt, tend, VTK);
       if(dim==3) driver_dg<GV,1,BLOCKSIZE3DP1>(gv, output, boundaryIndexToPhysicalEntity, dg_method, dt, tend, VTK);
       break;
     case 2:
       if(dim==2) driver_dg<GV,2,BLOCKSIZE2DP2>(gv, output, boundaryIndexToPhysicalEntity, dg_method, dt, tend, VTK);
       if(dim==3) driver_dg<GV,2,BLOCKSIZE3DP2>(gv, output, boundaryIndexToPhysicalEntity, dg_method, dt, tend, VTK);
       break;
     case 3:
       if(dim==2) driver_dg<GV,3,BLOCKSIZE2DP3>(gv, output, boundaryIndexToPhysicalEntity, dg_method, dt, tend, VTK);
       if(dim==3) driver_dg<GV,3,BLOCKSIZE3DP3>(gv, output, boundaryIndexToPhysicalEntity, dg_method, dt, tend, VTK);
       break;
   default:
     std::cout << "\nunsupported polynomial degree\n";
   }
 

}

 /**
 * \brief 	main function
 *
 * \fn main
 * 
 * \param argc number of arguments that is given by the user when the programm is started
 * \param argv arguments that are given by the user when the programm is started
 * 
 */


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      {
        if(helper.rank()==0)
          std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
      }

    if (argc!=5)
      {
        if(helper.rank()==0)
          std::cout << "usage: ./instationary <mesh file> <config file> <time step> <time max>" << std::endl;
        return 1;
      }
    
    int level = 0;
    int refs;
    int dim = 2;
    double dt, tend;
    std::string meshfile2d(argv[1]);
    std::string meshfile3d(argv[2]);
    /*unsigned int c_read = read_config(argv[2]);
	if(c_read>0) std::cout<<"Read config file ("<<c_read<<"MATa cells)."<<std::endl;
	else return 1;*/
    
    sscanf(argv[3],"%lg",&dt);
    sscanf(argv[4],"%lg",&tend);
    
    //datei.open("../output/bar1Boundary.txt"); 
    //datei << "Bdry-No. \t x \t y \t bar1-level \n";
	
    typedef Dune::UGGrid<3> GridType3D;
    //typedef Dune::UGGrid<3> GridType3;
    std::string output3d("3dtest");

    driver<GridType3D>(meshfile3d,output3d,level,ORDER,DG_METHOD,dt,tend);
    
    //typedef Dune::UGGrid<2> GridType2D;
//    std::string output2d("2dtest");
    
  //  driver<GridType2D>(meshfile2d,output2d,level,ORDER,DG_METHOD,dt,tend);
    
	/*
	switch(dim){
    case 2:
      driver<GridType2>(meshfile,level,ORDER,DG_METHOD,dt,tend);
      break;
    case 3:
      driver<GridType3>(meshfile,level,ORDER,DG_METHOD,dt,tend);
      break;
    }
	*/
    
    datei.close();
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}

