/*
 * driver_dg.hh
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



/**
 * \brief Reads the parameters for the system from the config file. Just that...
 * Will be called in the main file.
 *
 * @param[in] config_name Path to the config file.
 */
unsigned int read_config(const char* config_name, int cheat = 0)
{
	std::ifstream cf(config_name);
	if( !cf.is_open() ) return 0;
	
	cf >> N_ALPHA;
	cf >> N_MATA;
	if(cheat) cf >> N_CHEAT;	
	cf >> ALPHA_SEC;
	cf >> KBAR0;
	cf >> KBARI;
	cf >> BAR1_DIR;
	cf >> T_DELAY;
	
	cf.close();
	
	return N_MATA;
}

/**
 * \brief Heart and soul of the program. Basically does everything, sets upt the spaces, intitializes the solvers
 * and cycles through solving the first and second equation until the solution converges 
 * (time-shifted Picard's method).
 *
 * @param gv 		the grid view allows reading acces to the grid
 * @param pg 		a vector that assigns a boundary element the corresponding cell identifier
 * @param dt		lenght of time steps
 * @param tend		end time
 * @param dg_method	choice of the method, e.g. 0 for OBB, 1 for NIPG or 2 for SIPG
 * @param vtk 		Should VTK output be generated?
 */


template<class GV, unsigned int DEG, unsigned int BLOCKSIZE, typename PGMap>
void driver_dg (const GV& gv, std::string const& output, const PGMap& pg, int dg_method, double dt, double tend, bool vtk)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  Real time = 0.0;                                              // make a time variable

  typedef typename Dune::PDELab::MonomLocalFiniteElementMap<Coord,Real,dim,DEG> FEM0;
  FEM0 fem0(Dune::GeometryType::simplex);
  typedef typename Dune::PDELab::MonomLocalFiniteElementMap<Coord,Real,dim,DEG> FEM1;
  FEM1 fem1(Dune::GeometryType::simplex);
 
  // <<<2a>>> Make grid function space
  typedef Dune::PDELab::ISTLVectorBackend<BLOCKSIZE> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,Dune::PDELab::NoConstraints,VBE> GFS0;
  GFS0 gfs0(gv,fem0);
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM1,Dune::PDELab::NoConstraints,VBE> GFS1;
  GFS1 gfs1(gv,fem1);
  typedef BCType<GV,PGMap> B;
  B b(gv,pg);
  b.setTime(time);                                              
  typedef BCType1<GV,PGMap> B1;
  B1 b1(gv,pg);
  b1.setTime(time);                                              

  
  // <<<3>>> Make FE function with initial value / Dirichlet b.c.
  typedef typename GFS0::template VectorContainer<Real>::Type U0;
  U0 u0old(gfs0,0.0), nullvector(gfs0,0.0), u0new(gfs0,0.0);                                           
  typedef BCExtension<GV,Real,PGMap> G;                               // defines boundary condition,
  G g(gv,pg);                                                    
  g.setTime(time);                                              
  Dune::PDELab::interpolate(g,gfs0,u0old);
  typedef Initial<GV,Real,PGMap> Init;                                                                                                                         
  Init initial(gv,pg);                                                                                                                                                    
                                                                                                                                             
  typedef Neumann<GV,Real,PGMap> N;
  N neumann(gv,pg);
  neumann.setTime(time); 

  typedef Source<GV,Real> F;
  F f(gv);
  f.setTime(time);



  typedef K<GV,Real> K;
  K k(gv);
  k.setTime(time);


  typedef typename GFS1::template VectorContainer<Real>::Type U1;
  U1 u1old(gfs1,0.0), u1new(gfs1,0.0);                                                
  typedef BCExtension1<GV,Real,PGMap> G1;                               // defines boundary condition,
  G1 g1(gv,pg);                                                      
  g1.setTime(time);                                             
  typedef Initial1<GV,Real,PGMap> Init1;                                                                                                                         
  Init1 initial1(gv,pg);                                                                                                                                                    
  Dune::PDELab::interpolate(initial1,gfs1,u1old);

  typedef Neumann1<GV,Real,PGMap> N1;
  N1 neumann1(gv,pg);
  neumann1.setTime(time); 

  typedef Source1<GV,Real> F1;
  F1 f1(gv);
  f1.setTime(time);

  typedef Degradation1<GV,Real> Degradation1;
  Degradation1 a1(gv);
  a1.setTime(time);

  typedef K1<GV,Real> K1;
  K1 k1(gv);
  k1.setTime(time);


  // <<<4>>> Make instationary grid operator space
  typedef Dune::PDELab::DiscreteGridFunction<GFS0,U0> DGF;
  DGF nullfunction(gfs0,nullvector), u0dgf(gfs0,u0old);
  typedef Dune::PDELab::DiscreteGridFunction<GFS1,U1> DGF1;
  DGF1 u1dgf(gfs1,u1old);
  typedef Dune::PDELab::DiffusionDG<K,F,B,G,N,DGF1> LOP;
  LOP lop(k,f,b,g,neumann,u1dgf,dg_method,4);
  typedef HeatTimeLocalOperator TLOP; 
  TLOP tlop(4);                                                 
  typedef Dune::PDELab::DiffusionDG<K1,F1,B1,G1,N1,DGF> LOP1;
  LOP1 lop1(k1,f1,b1,g1,neumann1,nullfunction,dg_method,4);
  typedef HeatTimeLocalOperator TLOP1; 
  TLOP1 tlop1(4);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<BLOCKSIZE,BLOCKSIZE> MBE;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,U0,GFS0,GFS0,LOP,TLOP,
    Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,
    MBE > IGOS;
  IGOS igos(gfs0,gfs0,lop,tlop);

  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,U1,GFS1,GFS1,LOP1,TLOP1,
    Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,
    MBE > IGOS1;
  IGOS1 igos1(gfs1,gfs1,lop1,tlop1);

  // <<<5>>> Select a linear solver backend
  #ifdef USE_SUPER_LU // use lu decomposition as solver
  #if HAVE_SUPERLU
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(0);
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS1;
    LS1 ls1(0);
  #else error No superLU support, please install and configure it.
  #endif
  #else // Use iterative solver
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    LS ls(5000,false);
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS1;
    LS1 ls1(5000,false);
  #endif
  

  // <<<6>>> Solver for linear problem per stage
  
  typedef Dune::PDELab::Newton<IGOS,LS,U0> PDESOLVER;
  PDESOLVER pdesolver(igos,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setReduction(1e-10);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);
  
  typedef Dune::PDELab::Newton<IGOS1,LS1,U1> PDESOLVER1;
  PDESOLVER1 pdesolver1(igos1,ls1);
  pdesolver1.setReassembleThreshold(0.0);
  pdesolver1.setVerbosityLevel(2);
  pdesolver1.setReduction(1e-10);
  pdesolver1.setMinLinearReduction(1e-4);
  pdesolver1.setMaxIterations(25);
  pdesolver1.setLineSearchMaxIterations(10);


  // <<<7>>> time-stepper
  Dune::PDELab::Alexander3Parameter<Real> method;               // time stepping scheme
  Dune::PDELab::Alexander3Parameter<Real> method1;               // time stepping scheme
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,U0,U0> osm(method,igos,pdesolver);
  osm.setVerbosityLevel(1);                                     
  Dune::PDELab::OneStepMethod<Real,IGOS1,PDESOLVER1,U1,U1> osm1(method1,igos1,pdesolver1);
  osm1.setVerbosityLevel(1);                                   

  
  // <<<8>>> graphics for initial guess
  
  /*if (dim == 2) {
    Dune::PDELab::FilenameHelper fn("../output/2dtest");              // append number to file name
  }
  else {
    Dune::PDELab::FilenameHelper fn("../output/3dtest"); 
  }*/
  Dune::PDELab::FilenameHelper fn(output); 
  
  
  if(vtk)
  {
    //use the subsampling-VTK-writer for high quality output
    //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3); 
    Dune::VTKWriter<GV> vtkwriter(gv);
    
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(u0dgf,"alpha"));

    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(u1dgf,"bar1"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::ascii);
    fn.increment();                                             // increase file number
    
  }

  // <<<9>>> time loop
  while (time<tend-1e-8) {
    
      // do time step
      b.setTime(time+dt);
      f.setTime(time+dt);
      neumann.setTime(time+dt);
      g.setTime(time+dt);
      
      lop.setTime(time+dt);
      lop.set_dt(dt);
      //lop.update(u1new);
      
      osm.apply(time,dt,u0old,g,u0new);                           // do one time step

      // Calculation of absorption
      //alphaAbsorption(gv,gfs0,pg,absorptionAlpha,u0new,dt);
    
      // do time step
      b1.setTime(time+dt);
      f1.setTime(time+dt);
      neumann1.setTime(time+dt);
      neumann1.setDt(dt);
      g1.setTime(time+dt);
      g1.setDt(dt);
      
      lop1.setTime(time+dt);
      lop1.set_dt(dt);
      //lop1.update(nullvector);
      
      osm1.apply(time,dt,u1old,g1,u1new);                           // do one time step

      if(vtk) {
        // graphics

	//use the subsampling-VTK-writer for high quality output
        //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3); 
	Dune::VTKWriter<GV> vtkwriter(gv);
        
	vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(u0dgf,"alpha"));
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(u1dgf,"bar1"));
	vtkwriter.write(fn.getName(),Dune::VTKOptions::ascii);
		
        fn.increment();
      }

      u0old = u0new;  
      u1old = u1new;
                                          
      time += dt;		// advance time step
    }
  
  //bar1Boundary(gv,gfs1,pg,u1new);

    
}


