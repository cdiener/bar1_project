/*
 * bar1Boundary.hh
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
 * \brief this function calculates boundary values on the "world" boundary for Bar1
 *
 * \fn bar1Boundary
 */


template<class GV,class GFS, class PGMap, class X>
void bar1Boundary(const GV& gv, const GFS& gfs,const PGMap& pg, const X& x)
{


  typedef typename GFS::LocalFunctionSpace::Traits::SizeType size_type;
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);

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

  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    for (IntersectionIterator is = gv.ibegin(*eit); is!= gv.iend(*eit); ++is)
    {
      if(is->boundary())
      {
	physgroup_index = pg[is->boundarySegmentIndex()];

	if(physgroup_index == 0) {
	
	  Dune::FieldVector<DF, dim> globalcenter, localcenter;
	  globalcenter = is->geometry().center();
	  localcenter = eit->geometry().local(globalcenter);
          double length = is->geometry().volume();
          RangeType bar1;
          dgf.evaluate(*eit,localcenter,bar1);
	  
          datei << physgroup_index << "\t" 
				<< globalcenter[0] << "\t"
				<< globalcenter[1] << "\t"
				<< bar1 << std::endl;
         
	}
      }
    }
  }
}



