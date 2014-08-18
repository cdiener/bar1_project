/*
 * alpha_bc.hh
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



//	Definition of boundary conditions, source term, and initial concentration for the alpha-factor distribution                                                              

/** \brief boundary grid function selecting boundary conditions; 0 means Neumann, 1 means Dirichlet 
 * 
 * \class BCType 
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */

template<typename GV,typename PGMap>
class BCType : public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,BCType<GV,PGMap> >
{

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view                                                                                                                      
  BCType (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection                                                                                                       
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
    typedef typename GV::Grid::ctype ct;
    Dune::FieldVector<ct,GV::dimension> x = i.geometry().global(xlocal);

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];


    //Neumann default
    y = 0;
    
    //y = 1;
    //if( (physgroup_index < N_ALPHA) && (physgroup_index > 0) ) y = 0;
    //if( physgroup_index > 0 ) y = 0;
	//if(physgroup_index == 0) y = 1;
                                                                                                                                       
    return;
  }

  //! get a reference to the grid view                                                                                                              
  inline const GV& getGridView () {return gv;}

  private:

    const GV&  gv;
    const PGMap& pg;

};

/** \brief grid function for the Dirichlet boundary conditions                                                                                    
 *  \class BCExtension
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */                                                                             


template<typename GV, typename RF,typename PGMap>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension<GV,RF,PGMap> >
{

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                                                                  
  BCExtension (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

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


private:

  const GV&  gv;
  const PGMap& pg;

};

/** \brief boundary grid function for the Neumann boundary conditions                                                                                    
 *  \class Neumann
 * 
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */   
  

template<typename GV,typename RF,typename PGMap>
class Neumann : public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,Neumann<GV,RF,PGMap> >
{

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                                                                                                   

  Neumann (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection                                                                                                                                                                                    

  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {
    typedef typename GV::Grid::ctype ct;
    Dune::FieldVector<ct,GV::dimension> x = i.geometry().global(xlocal);
    typename Traits::RangeType norm = x.two_norm2();

    
    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    
    y = 0.0;
      
    if( (physgroup_index <= N_ALPHA) && (physgroup_index > 0) ) y = -ALPHA_SEC ; //Inflow
    
    /*else if(physgroup_index > 0 && !absorptionAlpha.empty())
    {   
    	RF integral = absorptionAlpha.back()[physgroup_index - N_ALPHA -1];
    	y = SIGMA*integral;
    }*/

    return;
                                                                                                                                                                                 
  }

  //! get a reference to the grid view                                                                                                                                                                                           
  inline const GV& getGridView () {return gv;}

private:

  const GV&  gv;
  const PGMap& pg;

};

/** \brief source function on the right hand side in the equation, for possibly produced or injected alpha-factor in the intercellular domain                                                                                     
 *  \class Source
 * 
 * \param time this function is dependend on the time
 * \param pg a vector that assigns a boundary element the corresponding cell identifier
 * \param gv grid view allows read only access to all grid elements
 */                                                                    

template<typename GV, typename RF>
class Source
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
					  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Source<GV,RF> >
{
  const GV& gv;
  RF time;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view                                                                                                                            
  Source (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element                                                                                                               
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);                                                                                                                                              
    typename Traits::RangeType norm = x.two_norm2();

    y = ALPHA_SOURCE;                                                                                                                                             
    return;
  }

  //! get a reference to the grid view                                                                                                                    
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation                                                                                                                  
  void setTime (double t) {time = t;}
};



