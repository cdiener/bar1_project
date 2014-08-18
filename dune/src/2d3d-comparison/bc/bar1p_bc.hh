/*
 * bar1_bc.hh
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

// Definition of boundary conditions, source term, and initial concentration for the Bar1 distribution                                                           

/** \brief boundary grid function selecting boundary conditions;                                                                                     
 * 0 means Neumann, 1 means Dirichlet 
 * 
 * 
 * \class BCType1 
 * 
 * \param time this function is dependend on the time
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */

template<typename GV,typename PGMap>
class BCType1 : public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,BCType1<GV,PGMap> >
{
  double time;
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view                                                                                                                      
  BCType1 (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection                                                                                                       
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
    typedef typename GV::Grid::ctype ct;
    Dune::FieldVector<ct,GV::dimension> x = i.geometry().global(xlocal), distance;

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];

//    y = 0; //Neumann everywhere

    //Optional Dirichlet
    y = 0;
    
    if( physgroup_index > N_ALPHA + N_CHEAT ) y = 1;
    
    /*distance[0] = center2_x;
    distance[1] = center2_y;
    distance[2] = 0.0;
    distance -= x;
    
    if ( distance.two_norm() < radius2 + 3.0*elementSize)
    {
      //std::cout << "Element-Coord:" << x[0] << "\t" << x[1] << "\t" << x[2] << std::endl;
      //std::cout << "Distance:" << distance[0] << "\t" << distance[1] << "\t" << distance[2] << std::endl;
      if (x[2] > 0.0)
	y = 1;
    }
    
    distance[0] = center3_x;
    distance[1] = center3_y;
    distance[2] = 0.0;
    distance -= x;
    
    if ( distance.two_norm() < radius3 + 3.0*elementSize)
    {
      //std::cout << "Element-Coord:" << x[0] << "\t" << x[1] << "\t" << x[2] << std::endl;
      //std::cout << "Distance:" << distance[0] << "\t" << distance[1] << "\t" << distance[2] << std::endl;
      if (x[2] > 0.0)
	y = 1;
    }*/
    
    
    //if( physgroup_index > 0 ) y = 0;
    
    /*y = 0;
    
    if ( x[0] > 0.0)
      y = 1;*/
    
                                                                                                                                       
    return;
  }

  //! get a reference to the grid view                                                                                                              
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation                                                                                                            
  void setTime (double t) {time = t;}
  private:

    const GV&  gv;
    const PGMap& pg;

};

/** \brief boundary grid function for the Neumann boundary conditions                                                                                    
 *  \class Neumann1
 * 
 * \param time this function is dependend on the time
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */                                                                                                                                                

template<typename GV,typename RF,typename PGMap>
class Neumann1 : public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,Neumann1<GV,RF,PGMap> >
{
  RF time;
  RF dt;
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                                                                                                   

  Neumann1 (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection                                                                                                                                                                                    

  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
    typedef typename GV::Grid::ctype ct;
    Dune::FieldVector<ct,GV::dimension> x = i.geometry().global(xlocal);
    typename Traits::RangeType norm = x.two_norm2();

                                                                                                                                                    
    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
      
    y = 0.0;
      
/*    if( physgroup_index > N_ALPHA + N_CHEAT ) {
       //RF integral = integrateAbsorption(dt,time - T_RECOV- T_DELAY, time - T_DELAY, physgroup_index - N_ALPHA -1);
       RF integral = integrateAlpha(dt, time - T_DELAY, physgroup_index - N_ALPHA -1);
       y = - KBAR0 - (double)KBARI*pow(integral,N_H)/(pow(integral,N_H)+pow(K_M,N_H));
    }*/

    return;
                                                                                                                                                                                 
  }

  //! get a reference to the grid view                                                                                                                                                                                           
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation                                                                                                                                                                                         
  void setTime (double t) {time = t;}

  void setDt (double Dt) {dt = Dt;}

private:

  const GV&  gv;
  const PGMap& pg;

};

/** \brief grid function for the Dirichlet boundary conditions                                                                                    
 *  \class BCExtension1
 * 
 * \param time this function is dependend on the time
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */      

template<typename GV, typename RF,typename PGMap>
class BCExtension1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension1<GV,RF,PGMap> >
{
  RF time;
  RF dt;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                                                                  
  BCExtension1 (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! evaluate extended function on element                                                                                                                                                     
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal), distance;
    
    
    
    y = BAR1_DIR;
    
    for (int i = N_ALPHA + N_CHEAT + 1; i < N_MATA + N_ALPHA + N_CHEAT + 1; i++)
    {
      distance[0] = coordsOfCells[i][0];
      distance[1] = coordsOfCells[i][1];
      distance[2] = 0.0;
      distance -= x;
      
      if ( distance.two_norm() < max_radius + 0.3*elementSize)
      {
	y = KBAR0; 
      }
    }

    typedef typename GV::IntersectionIterator IntersectionIterator;
    
    for (IntersectionIterator is = gv.ibegin(e); is!= gv.iend(e); ++is)                                                                                                                        
      {
        if(is->boundary()){
          int physgroup_index = pg[is->boundarySegmentIndex()];
     	  if( physgroup_index > N_ALPHA + N_CHEAT) {
	    RF integral = absorptionAlpha.back()[physgroup_index - N_ALPHA - N_CHEAT-1];
	    y = KBAR0 + (double)KBARI*pow(integral,N_H)/(pow(integral,N_H)+pow(K_M,N_H));
	  }
          //return;
        }
      }
    
    return;
  }

  //! get a reference to the grid view                                                                                                                                                          
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation                                                                                                                                                       
  void setTime (double t) {time = t;}
  void setDt (double Dt) {dt = Dt;}

private:

  const GV&  gv;
  const PGMap& pg;

};

/** \brief source function on the right hand side in the equation, for possibly produced or injected BAR1 in the intercellular domain                                                                                     
 *  \class Source1
 * 
 * \param time this function is dependend on the time
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */


template<typename GV, typename RF>
class Source1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Source1<GV,RF> >
{
  const GV& gv;
  RF time;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                            
  Source1 (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element                                                                                                               
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);                                                                                                                                              
    typename Traits::RangeType norm = x.two_norm2();

    y = 0.0;                                                                                                                                             
    return;
  }

  //! get a reference to the grid view                                                                                                                    
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation                                                                                                                  
  void setTime (double t) {time = t;}
};

/** \brief degradation of Bar1  ( actually set to zero in this model )                                                                                
 *  \class Degradation1
 * 
 * \param time this function is dependend on the time
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */

template<typename GV, typename RF>
class Degradation1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Degradation1<GV,RF> >
{
  const GV& gv;
  RF time;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                                                                  
  Degradation1 (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element                                                                                                                                                     
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);                                                                                                                              

    typename Traits::RangeType norm = x.two_norm2();

    y = 0.0;
    return;
  }

  //! get a reference to the grid view                                                                                                                                                          
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation                                                                                                                                                        
  void setTime (double t) {time = t;}
};

/** \brief initial concentration of Bar1 at time = 0                                                                                 
 *  \class Initial1
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */

template<typename GV, typename RF, typename PGMap>
class Initial1  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Initial1<GV,RF,PGMap> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                                                                 
                                                                                                                                                        
                                                                                                                                                                                               
Initial1 (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! evaluate extended function on element                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                               
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal), distance;

    //typename Traits::RangeType norm = x.two_norm2();

    typedef typename GV::IntersectionIterator IntersectionIterator;

    y = BAR1_DIR;
    
    for (int i = N_ALPHA + N_CHEAT + 1; i < N_MATA + N_ALPHA + N_CHEAT + 1; i++)
    {
      distance[0] = coordsOfCells[i][0];
      distance[1] = coordsOfCells[i][1];
      distance[2] = 0.0;
      distance -= x;
      
      if ( distance.two_norm() < max_radius + 0.3*elementSize)
      {
	y = KBAR0; 
      }
    }

    typedef typename GV::IntersectionIterator IntersectionIterator;
    
    for (IntersectionIterator is = gv.ibegin(e); is!= gv.iend(e); ++is)                                                                                                                        
      {
        if(is->boundary()){
          int physgroup_index = pg[is->boundarySegmentIndex()];
     	  if( physgroup_index > N_ALPHA + N_CHEAT) {
	    RF integral = absorptionAlpha.back()[physgroup_index - N_ALPHA - N_CHEAT-1];
	    double r = coordsOfCells[physgroup_index][2];
	    y = KBAR0 + (double)KBARI*pow(integral,N_H)/(pow(integral,N_H)+pow(r*K_M,N_H));
	  }
          //return;
        }
      }
    
    



    return;
  }

  //! get a reference to the grid view                                                                                                                                                         \ 
                                                                                                                                                                                               

  inline const GV& getGridView () {return gv;}

  private:
  const PGMap& pg;
};
