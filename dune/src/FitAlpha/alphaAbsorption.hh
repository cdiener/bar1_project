/*
 * alphaAbsorption.hh
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
 * \brief this function calculates the integral of the alpha-factor absorption of a given MATa-Cell from 0 to each time step
 *
 * \fn alphaAbsorption
 * 
 * \param[in]  	dt    		length of time step
 * \param[in]  	t_start    	start point for integration
 * \param[in]  	t_end    	end point for integration
 * \param[in]	cellNumber	identifier for the corresponding MATa-cell
 * \param[out]	absorbtion  	vector that stores the value of the integral for each time step
 * 
 * \return	value of the integral
 */


template<class GV,class GFS, class T, class PGMap, class X>
void alphaAbsorption(const GV& gv, const GFS& gfs,const PGMap& pg, T &absorption, const X& x)
{


  typedef typename GFS::LocalFunctionSpace::Traits::SizeType size_type;
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);

  //typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS,X> DGFG;
  //DGFG dgfg(gfs,x);

  typedef Dune::PDELab::DiscreteGridFunction<GFS,X> DGF;
  DGF dgf(gfs,x);

  //
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::Grid::ctype Coord;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
  typedef double Real;
  const int dim = GV::dimension;
  int physgroup_index =0; 
  std::vector<double> temp, circumference;

  for (int i = 0; i < N_MATA; i++) {
    temp.push_back(0.0);
    circumference.push_back(0.0);
  }  

  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    for (IntersectionIterator is = gv.ibegin(*eit); is!= gv.iend(*eit); ++is)
    {
      if(is->boundary())
      {
	physgroup_index = pg[is->boundarySegmentIndex()];
        if(physgroup_index > N_ALPHA) {
   	  Dune::FieldVector<DF, dim> globalcenter, localcenter;
	  globalcenter = is->geometry().center();
	  localcenter = eit->geometry().local(globalcenter);
          double length = is->geometry().volume();
          circumference[physgroup_index - N_ALPHA - 1] += length;
          
          RangeType alpha;
          dgf.evaluate(*eit,localcenter,alpha);
          temp[physgroup_index - N_ALPHA -1] += (double)alpha*length; //*SIGMA
        }
      }
    }
  }

  for (int i = 0; i < N_MATA; i++) {
    	temp[i] = temp[i]/circumference[i];	
  }  

  
  absorption.push_back(temp);
}

