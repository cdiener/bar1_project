/*
 * bar1p_bc.hh
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



//  Definition of boundary conditions, source term, and initial concentration for the Bar1 distribution                                                                

/** \brief boundary grid function selecting boundary conditions;                                                                                     
 * 0 means Neumann, 1 means Dirichlet 
 * 
 * 
 * \class BCType1 
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */

template<typename GV,typename PGMap>
class BCType1 : public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,BCType1<GV,PGMap> >
{
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view                                                                                                                      
  BCType1 (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection                                                                                                       
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
    typedef typename GV::Grid::ctype ct;
    Dune::FieldVector<ct,GV::dimension> x = i.geometry().global(xlocal);

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];

    y = 0;
    //if( (physgroup_index < N_ALPHA) && (physgroup_index > 0) ) y = 0;
    if(physgroup_index > N_ALPHA + N_CHEAT) y = 1;
                                                                                                                                       
    return;
  }

  //! get a reference to the grid view                                                                                                              
  inline const GV& getGridView () {return gv;}

  private:

    const GV&  gv;
    const PGMap& pg;

};

/** \brief boundary grid function for the Neumann boundary conditions                                                                                    
 *  \class Neumann1
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */   
                                                                                                                                             

template<typename GV,typename RF,typename PGMap>
class Neumann1 : public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,Neumann1<GV,RF,PGMap> >
{
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
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */     

template<typename GV, typename RF,typename PGMap>
class BCExtension1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension1<GV,RF,PGMap> >
{

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
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
                                                                                                                                                                  
    typename Traits::RangeType norm = x.two_norm2();
                                                                                                                                                                                                                                                                                                    
    y = 0.0;
 
    typedef typename GV::IntersectionIterator IntersectionIterator;
    
    for (IntersectionIterator is = gv.ibegin(e); is!= gv.iend(e); ++is)                                                                                                                        
      {
        if(is->boundary()){
          int physgroup_index = pg[is->boundarySegmentIndex()];
     	  if( physgroup_index > N_ALPHA + N_CHEAT) {
	      RF integral = absorptionAlpha.back()[physgroup_index - N_ALPHA - N_CHEAT-1];
	      y = KBAR0 + (double)KBARI*pow(integral,N_H)/(pow(integral,N_H)+pow(K_M,N_H));
	  }
        }
      }
      
      
    return;
  }

  //! get a reference to the grid view                                                                                                                                                          
  inline const GV& getGridView () {return gv;}

private:

  const GV&  gv;
  const PGMap& pg;

};

/** \brief source function on the right hand side in the equation, for possibly produced or injected BAR1 in the intercellular domain                                                                                     
 *  \class Source1
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */


template<typename GV, typename RF>
class Source1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Source1<GV,RF> >
{
  const GV& gv;
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

};



