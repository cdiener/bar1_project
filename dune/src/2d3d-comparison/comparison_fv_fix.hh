/*
 * yeast_fv_fix.hh
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
 * \brief Heart and soul of the program. Basically does everything, sets upt the spaces, intitializes the solvers
 * and cycles through solving the first and second equation until the solution converges 
 *
 * @param gv 		the grid view allows reading acces to the grid
 * @param pg 		a vector that assigns a boundary element the corresponding cell identifier
 * @param vtk 		Should VTK output be generated?
 */

template<class GV2D, class GV3D, typename PGMap2D, typename PGMap3D>
void comparison_fv (const GV2D& gv2d, const GV3D& gv3d, std::string const& output2d, std::string const& output3d, const PGMap2D& pg2d, const PGMap3D& pg3d, bool vtk)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV3D::Grid::ctype Coord;
  typedef double Real;
  const int dim = 3;

  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM0;
  FEM0 fem0(Dune::GeometryType::simplex);
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM1;
  FEM1 fem1(Dune::GeometryType::simplex);
 
  // <<<2a>>> Make grid function space
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV3D,FEM0,Dune::PDELab::NoConstraints,VBE> GFS0;
  GFS0 gfs0(gv3d,fem0);
  typedef Dune::PDELab::GridFunctionSpace<GV3D,FEM1,Dune::PDELab::NoConstraints,VBE> GFS1;
  GFS1 gfs1(gv3d,fem1);
  typedef BCType<GV3D,PGMap3D> B;
  B b(gv3d,pg3d);
                                       
  typedef BCType1<GV3D,PGMap3D> B1;
  B1 b1(gv3d,pg3d);
                                             
  typedef typename GFS0::template ConstraintsContainer<Real>::Type CC0;
  CC0 cc0;
  Dune::PDELab::constraints(b,gfs0,cc0);
  typedef typename GFS1::template ConstraintsContainer<Real>::Type CC1;
  CC1 cc1;
  Dune::PDELab::constraints(b1,gfs1,cc1);

  
  std::vector<double> temp(N_MATA);
  for(int i=0; i < temp.size(); i++)
    temp[i] = 0.0;
    
  absorptionAlpha.push_back(temp);
  
  // <<<3>>> Make FE function for Dirichlet/ Neumann b.c.
  typedef typename GFS0::template VectorContainer<Real>::Type U0;
  U0 u0old(gfs0,0.0), nullvector(gfs0,0.0); 
  U0 u0(gfs0,0.0),x0(gfs0,0.0);                                          
  typedef BCExtension<GV3D,Real,PGMap3D> G;                           
  G g(gv3d,pg3d);                                                 
                                   
  Dune::PDELab::interpolate(g,gfs0,u0);
                                                                                                                                               
                                                                                                                                   
  typedef Neumann<GV3D,Real,PGMap3D> N;
  N neumann(gv3d,pg3d);

  typedef Source<GV3D,Real> F;
  F f(gv3d);

  typedef K<GV3D,Real> K;
  K k(gv3d);


  typedef typename GFS1::template VectorContainer<Real>::Type U1;
  U1 u1old(gfs1,0.0);   
  U1 u1(gfs1,BAR1_DIR);                                         
  typedef BCExtension1<GV3D,Real,PGMap3D> G1;                           
  G1 g1(gv3d,pg3d);                                                   
 
  Dune::PDELab::interpolate(g1,gfs0,u1);                                          
                                                                                                                                                                                                                                                                   

  typedef Neumann1<GV3D,Real,PGMap3D> N1;
  N1 neumann1(gv3d,pg3d);

  typedef Source1<GV3D,Real> F1;
  F1 f1(gv3d);

  typedef K1<GV3D,Real> K1;
  K1 k1(gv3d);


  // <<<4>>> Make instationary grid operator space
  typedef Dune::PDELab::DiffusionCCFV<K1,U0,F1,B1,N1,G1> LOP1;
  LOP1 lop1(k1,nullvector,f1,b1,neumann1,g1,D1);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

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
  #else // Use iterative solver*/
    // make ISTL solver
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS1;
    LS1 ls1(5000,false);
  #endif
  

  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS1,LS1,U1> PDESOLVER1;
  PDESOLVER1 pdesolver1(gos1,u1,ls1,1e-12);
  
  
  
  typedef Dune::PDELab::DiffusionCCFV<K,U1,F,B,N,G> LOP;
  LOP lop(k,u1,f,b,neumann,g,D0);
  
  typedef Dune::PDELab::GridOperatorSpace<GFS0,GFS0,LOP,Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,MBE> GOS;
  GOS gos(gfs0,gfs0,lop);
  
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
  
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U0> PDESOLVER;
  PDESOLVER pdesolver(gos,u0,ls,1e-12);
    
  double norm;
  
  // <<<6>>> iterate between the to equation until the solution converges
  
  std::vector<double> c(gv2d.size(0));
  std::fill(c.begin(),c.end(),0.0); 
  
  typedef typename GFS0::LocalFunctionSpace LFS0;
  typedef typename LFS0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  RangeType max_value = 0.0;
  
  for(unsigned int i=0; i<100; i++)
  {
     // Solve for alpha
     pdesolver.apply();
     
     
     //alphaAbsorption(gv3d,gfs0,pg3d,absorptionAlpha,u0);
     
     typedef Dune::PDELab::DiscreteGridFunction<GFS0,U0> DGF;
     DGF udgf(gfs0,u0);
     
     integrateOverHeight(gv2d,gv3d,udgf,c,max_value);
     alphaAbsorptionOverHeight(gv2d,absorptionAlpha,c);
  
     // Solve for bar1
     pdesolver1.apply();
     lop.update(u1);
      
      
     u0old -= u0;
     norm = u0old.two_norm()/u0.one_norm();
     std::cout << "eps = " << norm << std::endl;
     if(norm < 1.0e-4) break;
     u0old = u0;
  }

  // <<<7>>> graphics for initial guess
  if(vtk)
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS0,U0> DGF;
    DGF udgf(gfs0,u0);
    typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS0,U0> DGFG;
    DGF u0dgfg(gfs0,u0);
    max_gradient(gv3d,gfs0,u0dgfg);
    
    typedef Dune::PDELab::DiscreteGridFunction<GFS1,U1> DGF1;
    DGF1 udgf1(gfs1,u1);
    //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    Dune::VTKWriter<GV3D> vtkwriter(gv3d);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"alpha"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(udgf1,"bar1"));
    vtkwriter.write(output3d,Dune::VTKOptions::binaryappended);
    
    // Output of the concentration values on the surface of the aCells
    typedef typename GV3D::template Codim<0>::Iterator ElementIterator;
    typedef typename GV3D::IntersectionIterator IntersectionIterator;
    int physgroup_index =0; 
	
    for(int i= (N_ALPHA +1); i < (N_ALPHA + N_MATA + 1 ); i ++ ) {
    
      std::cout << absorptionAlpha.back()[i - N_ALPHA -1]/SIGMA << std::endl; //correct for adsorption parameters
    
    }
    
    double value = 0.0;
    Dune::FieldVector<Coord,3> CoordVec;
    CoordVec[0] = 0.0;
    CoordVec[1] = 0.0;
    CoordVec[2] = 6.0;
    
    
    //value = evaluateGlobal(gv3d,gfs0,udgf,CoordVec);
    
    //std::cout << "Value in (" <<  CoordVec[0] << "," << CoordVec[1] << "," << CoordVec[2] << ") is " << value << std::endl;
    
    typedef typename GV2D::template Codim<0>::Iterator ElementIterator2D;
    
    typedef typename GV3D::template Codim<0>::Iterator ElementIterator3D;
    
    integrateOverHeight(gv2d,gv3d,udgf,c,max_value);
    
    RangeType integrated_max_value = 0.0;
    
    for ( int i = 0; i < c.size(); i++)
    {
      if (integrated_max_value < c[i]) integrated_max_value = c[i];
    }
    std::cout << "\nMaximal pheromone concentration 3D: " << max_value << std::endl;
    std::cout << "\nMaximal integrated pheromone 3D: " << integrated_max_value << std::endl;
    
    //std::cout << "\nMaximal pheromone concentration(Integrated over h): " << max_value*Integration_height << std::endl;
    
     
    /*for (ElementIterator2D eit = gv2d.template begin<0>();
	   eit!= gv2d.template end<0>(); ++eit)
    {
      int index = gv2d.indexSet().index(*eit);
      Dune::FieldVector<Coord,2> globalcenter, distance;
      globalcenter = eit->geometry().center();
      
      CoordVec[0] = globalcenter[0] - 170.0;
      CoordVec[1] = globalcenter[1] - 170.0;
      CoordVec[2] = 0.1;
      
      double delta_h = Integration_height/((double)Integration_steps);
      
      for (int i = 0; i < Integration_steps + 1; i++) {
	CoordVec[2] = (double)i*delta_h;
	c[index] += delta_h*evaluateGlobal(gv3d,gfs0,udgf,CoordVec);
	//std::cout << "Integration up to height " << (double)i*delta_h << "yields: " << c[index] << std::endl;
      }
      
      //std::cout << "Pheromonkonzentration in : " << c[index] << std::endl;
      
      distance[0] = globalcenter[0] - 170.0;
      distance[1] = globalcenter[1] - 170.0;
      double temp = distance[0]*distance[0] + distance[1]*distance[1];
           
      if (temp > 170.0*170.0) c[index] = 0.0;
      
      distance[0] = globalcenter[0] - 170.0 - center1_x;
      distance[1] = globalcenter[1] - 170.0 - center1_y;
      
      temp = distance[0]*distance[0] + distance[1]*distance[1];
      
      if (temp < radius1*radius1) c[index] = 0.0;
      
      distance[0] = globalcenter[0] - 170.0 - center2_x;
      distance[1] = globalcenter[1] - 170.0 - center2_y;
      
      temp = distance[0]*distance[0] + distance[1]*distance[1];
      
      if (temp < radius2*radius2) c[index] = 0.0;
      
      distance[0] = globalcenter[0] - 170.0 - center3_x;
      distance[1] = globalcenter[1] - 170.0 - center3_y;
      
      temp = distance[0]*distance[0] + distance[1]*distance[1];
      
      if (temp < radius3*radius3) c[index] = 0.0;
      
    }*/
    
    //Dune::VTKWriter<GV2D> vtkwriter2d(gv2d);
    
    alphaAbsorptionOverHeight(gv2d,absorptionAlpha,c);
    
    Dune::VTKWriter<GV2D> vtkwriter2d(gv2d);
    //Dune::SubsamplingVTKWriter<GV2D> vtkwriter2d(gv2d,0);
    vtkwriter2d.addCellData (c , "integrated alpha" , 1);
    vtkwriter2d.write( output2d, Dune::VTKOptions::binaryappended  );

    
  }


    
}

