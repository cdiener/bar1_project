/*
 * yeast_dg_fix.hh
 * This file is part of ExCeCo
 *
 * Copyright (C) 2011 - Christian Diener, Wolfgang Chris
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
	
	
	cf.close();
	
	return N_MATA;
}

/**
 * \brief Heart and soul of the program. Basically does everything, sets upt the spaces, intitializes the solvers
 * and cycles through solving the first and second equation until the solution converges 
 *
 * @param gv 		the grid view allows reading acces to the grid
 * @param pg 		a vector that assigns a boundary element the corresponding cell identifier
 * @param dg_method	choice of the method, e.g. 0 for OBB, 1 for NIPG or 2 for SIPG
 * @param vtk 		Should VTK output be generated?
 */

template<class GV, unsigned int DEG, unsigned int BLOCKSIZE, typename PGMap>
void yeast_dg (const GV& gv, std::string const& output, const PGMap& pg, int dg_method, bool vtk)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  typedef Dune::PDELab::MonomLocalFiniteElementMap<Coord,Real,dim,DEG> FEM0;
  FEM0 fem0(Dune::GeometryType::simplex);
  typedef Dune::PDELab::MonomLocalFiniteElementMap<Coord,Real,dim,DEG> FEM1;
  FEM1 fem1(Dune::GeometryType::simplex);
 
  
  // <<<2a>>> Make grid function space
  typedef Dune::PDELab::ISTLVectorBackend<BLOCKSIZE> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,Dune::PDELab::NoConstraints,VBE> GFS0;
  GFS0 gfs0(gv,fem0);
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM1,Dune::PDELab::NoConstraints,VBE> GFS1;
  GFS1 gfs1(gv,fem1);
  typedef BCType<GV,PGMap> B;
  B b(gv,pg);
                                       
  typedef BCType1<GV,PGMap> B1;
  B1 b1(gv,pg);
                                             
  typedef typename GFS0::template ConstraintsContainer<Real>::Type CC0;
  CC0 cc0;
  Dune::PDELab::constraints(b,gfs0,cc0);
  typedef typename GFS1::template ConstraintsContainer<Real>::Type CC1;
  CC1 cc1;
  Dune::PDELab::constraints(b1,gfs1,cc1);

  
  // <<<3>>> Make FE function for Dirichlet/ Neumann b.c.
  typedef typename GFS0::template VectorContainer<Real>::Type U0;
  U0 u0old(gfs0,0.0), nullvector(gfs0,0.0); 
  U0 u0(gfs0,0.0),x0(gfs0,0.0);                                          
  typedef BCExtension<GV,Real,PGMap> G;                           
  G g(gv,pg);                                                 
                                   
  Dune::PDELab::interpolate(g,gfs0,u0);
                                                                                                                                               
                                                                                                                                   
  typedef Neumann<GV,Real,PGMap> N;
  N neumann(gv,pg);

  typedef Source<GV,Real> F;
  F f(gv);

  typedef K<GV,Real> K;
  K k(gv);

  // set absorptionAlpha to some value, so that the dirichlet can be initialized
  std::vector<double> temp(N_MATA);
  for(int i=0; i < temp.size(); i++)
    temp[i] = 0.0;
    
  absorptionAlpha.push_back(temp);
  
  //alphaAbsorption(gv,gfs0,pg,absorptionAlpha,u0,1.0);

  typedef typename GFS1::template VectorContainer<Real>::Type U1;
  U1 u1old(gfs1,0.0);   
  //U1 u1(gfs1,BAR1_DIR);  
  U1 u1(gfs1,0.0);  
  typedef BCExtension1<GV,Real,PGMap> G1;                           
  G1 g1(gv,pg);                                                   
 
  Dune::PDELab::interpolate(g1,gfs1,u1);                                          
                                                                                                                                                                                                                                                                   

  typedef Neumann1<GV,Real,PGMap> N1;
  N1 neumann1(gv,pg);

  typedef Source1<GV,Real> F1;
  F1 f1(gv);

  typedef K1<GV,Real> K1;
  K1 k1(gv);
  

  // <<<4>>> Make instationary grid operator space
  typedef Dune::PDELab::DiscreteGridFunction<GFS0,U0> DGF;
  DGF nullfunction(gfs0,nullvector), u0dgf(gfs0,u0);
  
  //typedef Dune::PDELab::DiffusionDG<K1,F1,B1,G1,N1,U0> LOP1;
  typedef Dune::PDELab::DiffusionDG<K1,F1,B1,G1,N1,DGF> LOP1;
  //LOP1 lop1(k1,f1,b1,g1,neumann1,nullvector,dg_method,4);
  LOP1 lop1(k1,f1,b1,g1,neumann1,nullfunction,dg_method,4);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<BLOCKSIZE,BLOCKSIZE> MBE;

  typedef Dune::PDELab::GridOperatorSpace<GFS1,GFS1,LOP1,Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,MBE> GOS1;
  GOS1 gos1(gfs1,gfs1,lop1);

  // <<<5>>> Select a linear solver backend
  #ifdef USE_SUPER_LU // use lu decomposition as solver
  #if HAVE_SUPERLU
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS1;
    LS1 ls1(1);
  #else error No superLU support, please install and configure it.
  #endif
  #else
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS1;
    LS1 ls1(5000,false);
  #endif
  
  /*typedef Dune::PDELab::Newton<GOS1,LS1,U1> PDESOLVER1;
  PDESOLVER1 pdesolver1(gos1,ls1);
  pdesolver1.setReassembleThreshold(0.0);
  pdesolver1.setVerbosityLevel(2);
  pdesolver1.setReduction(1e-10);
  pdesolver1.setMinLinearReduction(1e-4);
  pdesolver1.setMaxIterations(25);
  pdesolver1.setLineSearchMaxIterations(10);*/

  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS1,LS1,U1> PDESOLVER1;
  PDESOLVER1 pdesolver1(gos1,u1,ls1,1e-6);
  
  
  typedef Dune::PDELab::DiscreteGridFunction<GFS1,U1> DGF1;
  DGF1 u1dgf(gfs1,u1);
  //typedef Dune::PDELab::DiffusionDG<K,F,B,G,N,U1> LOP;
  typedef Dune::PDELab::DiffusionDG<K,F,B,G,N,DGF1> LOP;
  //LOP lop(k,f,b,g,neumann,u1,dg_method,4);
  LOP lop(k,f,b,g,neumann,u1dgf,dg_method,4);
  
  typedef Dune::PDELab::GridOperatorSpace<GFS0,GFS0,LOP,Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,MBE> GOS;
  GOS gos(gfs0,gfs0,lop);
  
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  //LS ls(5000,false);
  
  //typedef Dune::PDELab::ISTLBackend_SEQ_ILU0<Dune::LoopSolver> LS;
  //LS ls(5000,false);
  #ifdef USE_SUPER_LU // use lu decomposition as solver
  #if HAVE_SUPERLU
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(1);
  #else error No superLU support, please install and configure it.
  #endif
  #else // Use iterative solver*/
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    LS ls(5000,false);
  #endif
    
  /*typedef Dune::PDELab::Newton<GOS,LS,U0> PDESOLVER;
  PDESOLVER pdesolver(gos,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setReduction(1e-10);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);*/
  
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U0> PDESOLVER;
  PDESOLVER pdesolver(gos,u0,ls,1e-6);
  
  // for better initial guess
  /*std::vector<double> temp;
  
  for (int i = 0; i < N_MATA; i++) temp.push_back(1.0);
  
  absorptionAlpha.push_back(temp);
  
  pdesolver1.apply();*/
    
  double norm;
  
  
  
  // <<<6>>> iterate between the to equation until the solution converges
  
  for(unsigned int i=1; i<=100; i++)
  {  std::cout << "-----------------------" << std::endl ;
     std::cout << "Fixed point step: " << i << ".\n" ;
     std::cout << "-----------------------" << std::endl ;
    
     // Solve for bar1
     pdesolver1.apply();
     
     
   
     // Solve for alpha
     pdesolver.apply();
     alphaAbsorption(gv,gfs0,pg,absorptionAlpha,u0);
     
     typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS0,U0> DGFG;
     DGF u0dgfg(gfs0,u0);
     max_gradient(gv,gfs0,u0dgfg);
     
     //DGF udgf(gfs0,u0);
     //DGF1 udgf1(gfs1,u1);
     
     std::string s;
     std::stringstream out;
     out << i;
     s = out.str();
     
     if(FIXEDPOINT_VTK){
	//Dune::VTKWriter<GV> vtkwriter(gv);
	Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
	vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(u0dgf,"alpha"));
	vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(u1dgf,"bar1"));
	vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<G>(g,"boundary-function"));
	
	//vtkwriter.write("Pheromone-FixedpointStep-" + s,Dune::VTKOptions::binaryappended);
	vtkwriter.write(output + "-" + s,Dune::VTKOptions::ascii);
     
	//Dune::VTKWriter<GV> vtkwriter1(gv);
	//Dune::SubsamplingVTKWriter<GV> vtkwriter1(gv,3);
	//vtkwriter1.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(u1dgf,"bar1"));
	//vtkwriter1.write("bar1-FixedpointStep-"+ s ,Dune::VTKOptions::binaryappended);
	//vtkwriter1.write("../output/bar1-FP-"+ s ,Dune::VTKOptions::ascii);
     }   

     
     //lop.update(u1);
     //lop.update(udgf1);
     //lop.update(u1dgf);
     //lop1.update(nullvector);
      
     u0old -= u0;
     norm = u0old.one_norm()/u0.one_norm();
     std::cout << "Fixed point error: eps = " << norm << std::endl;
     if(norm < 1.0e-6) break;
     u0old = u0;
  }

  // <<<7>>> graphics for initial guess
  if(vtk)
  {
    //typedef Dune::PDELab::DiscreteGridFunction<GFS0,U0> DGF;
    //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    Dune::VTKWriter<GV> vtkwriter(gv);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(u0dgf,"alpha"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(u1dgf,"bar1"));
    //vtkwriter.write("PheromoneDistribution",Dune::VTKOptions::binaryappended);
    vtkwriter.write(output,Dune::VTKOptions::ascii);
    
    std::cout << "Final output!" << std::endl;
    
    // Output of the concentration values on the surface of the aCells
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    int physgroup_index =0; 
	
    for(int i= (N_ALPHA +1); i < (N_ALPHA + N_MATA + 1 ); i ++ ) {
    
      std::cout << absorptionAlpha.back()[i - N_ALPHA -1]/SIGMA << std::endl; //correct for adsorption parameters
    
    }

    alphaValues(gv,gfs0,pg,u0);
    
    typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS0,U0> DGFG;
    DGF u0dgfg(gfs0,u0);
    max_gradient(gv,gfs0,u0dgfg);
    
    //typedef Dune::PDELab::DiscreteGridFunction<GFS1,U1> DGF1;
    //Dune::SubsamplingVTKWriter<GV> vtkwriter1(gv,3);
    //Dune::VTKWriter<GV> vtkwriter1(gv);
    //vtkwriter1.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(u1dgf,"bar1"));
    //vtkwriter1.write("bar1",Dune::VTKOptions::binaryappended);
    //vtkwriter1.write("bar1",Dune::VTKOptions::ascii);

  }


    
}

