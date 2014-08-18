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
        if(physgroup_index > N_ALPHA + N_CHEAT) {
   	  Dune::FieldVector<DF, dim> globalcenter, localcenter;
	  globalcenter = is->geometry().center();
	  localcenter = eit->geometry().local(globalcenter);
          double length = is->geometry().volume();
          circumference[physgroup_index - N_ALPHA - N_CHEAT- 1] += length;
          
          RangeType alpha;
          dgf.evaluate(*eit,localcenter,alpha);
          temp[physgroup_index - N_ALPHA - N_CHEAT - 1] += alpha*length; //*SIGMA
        }
      }
    }
  }

  for (int i = 0; i < N_MATA; i++) {
    temp[i] = temp[i]/circumference[i];	
    //test << i << "\t" << circumference[i] << std::endl;
    std::cout << "Cell Nr.:" << i << "\tSurface area:" << circumference[i] << std::endl;
    if (dim == 2) {
      std::cout << "\t\t2D-Absorption:" << temp[i] << std::endl;
      std::cout << "\t\t2D-Absorption/h:" << temp[i]/Integration_height << std::endl;
      //std::cout << "\t\t2D-Absorption/R_i:" << temp[i]/coordsOfCells[i + N_ALPHA + N_CHEAT][2] << std::endl;      
    }
    else {
      std::cout << "\t\t3D-Absorption:" << temp[i] << std::endl;
      std::cout << "\t\t3D-Absorption * h:" << temp[i]*Integration_height << std::endl;
      std::cout << "\t\t3D-Absorption * R_i:" << temp[i]*coordsOfCells[i + N_ALPHA + N_CHEAT][2] << std::endl;
    }
    
    
  }  
  
  
  
  absorption.push_back(temp);
}

template<class GV,class GFS, class PGMap, class X>
void alphaValues(const GV& gv, const GFS& gfs,const PGMap& pg, const X& x)
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

  datei << "CellId\tx\ty\talpha" << std::endl;
  
  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    for (IntersectionIterator is = gv.ibegin(*eit); is!= gv.iend(*eit); ++is)
    {
      if(is->boundary())
      {
	physgroup_index = pg[is->boundarySegmentIndex()];
        if(physgroup_index > N_ALPHA + N_CHEAT) {
   	  Dune::FieldVector<DF, dim> globalcenter, localcenter;
	  globalcenter = is->geometry().center();
	  localcenter = eit->geometry().local(globalcenter);
          
          RangeType alpha;
	  
          dgf.evaluate(*eit,localcenter,alpha);
	  
	  
	  datei << physgroup_index << "\t" << globalcenter[0] << "\t" << globalcenter[1] << "\t" << alpha << std::endl;
	  

        }
      }
    }
  }

}

template<class GV, class GFS, class DGFG>
void max_gradient(const GV& gv, const GFS& gfs, const DGFG& dgfg)
{


  typedef typename GFS::LocalFunctionSpace::Traits::SizeType size_type;
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);

  //
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::Grid::ctype Coord;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RT;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
  typedef double Real;
  const int dim = GV::dimension;
 
  
  
  for (int j = 1; j < 11; j++)
  {
    double Jz_max = 0.0;
    double Jxy_max = 0.0;
    double flux_max = 0.0;
    double slice = (double)j;
    
    for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
    {
      Dune::FieldVector<DF, dim> globalcenter, localcenter;
      globalcenter = eit->geometry().center();
      localcenter = eit->geometry().local(globalcenter);
      //Dune::FieldVector<RT, dim> gradient;
      RangeType gradient;
      dgfg.evaluate(*eit,localcenter,gradient);
      
      //std::cout << "Flux in z = " << globalcenter[2] << ":\t"<< gradient[0] << ":\t"<< gradient[1] << ":\t"<< gradient[2] << std::endl;
      
      double Jx = gradient[0];
      double Jy = gradient[1];
      double Jz = gradient[2];
      
      double flux = gradient[0];
      
      if( Jx < 0.0 ) Jx = - Jx;
      if( Jy < 0.0 ) Jy = - Jy;
      if( Jz < 0.0 ) Jz = - Jz;
      
      if( flux < 0.0 ) flux = - flux;
      
      flux = D0*flux/ALPHA_SEC;
      
      //std::cout << "...and absolute values " << Jx << ":\t"<< Jy << ":\t"<< Jz << std::endl;
      
      
      if ( (17.0*(slice-1.0) < globalcenter[2] ) && (17.0*slice > globalcenter[2] ) ){
	//if ( Jz > Jz_max  )   Jz_max = Jz;
	//if ( Jx + Jy > Jxy_max  )   Jxy_max = Jx + Jy;
	if ( flux > flux_max  )   flux_max = flux;
	
      }
    }
    
    //std::cout << "Maximal z-Flux in layer " << j << ":\t"<< Jz_max << std::endl;
    //std::cout << "Maximal xy-Flux in layer " << j << ":\t"<< Jxy_max << std::endl;
    std::cout << "Maximal flux in layer " << j << " at " << 8.5 + 17.0*(slice-1.0) << ":\t"<< flux_max << std::endl;
  }
  
}  

  
template<class GV, class GFS, class DGF, typename CoordVec>
double evaluateGlobal(const GV& gv, const GFS& gfs, const DGF& dgf, const CoordVec& coord)
{


  typedef typename GFS::LocalFunctionSpace::Traits::SizeType size_type;
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);
  //
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::Grid::ctype Coord;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RT;
  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
  typedef double Real;
  const int dim = GV::dimension;
 
  double value = 0.0;
  
  // look for nearest mesh element
  double distance = 170.0;
  
  int index = 0;
  
  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    int temp_index = gv.indexSet().index(*eit);
    //Dune::FieldVector<DF, dim> globalcenter, localcenter;
    //CoordVec globalcenter;
    Dune::FieldVector<DF, dim> globalcenter;
    globalcenter = eit->geometry().center();
    //globalcenter -= coord;
    globalcenter[0] -= coord[0];
    globalcenter[1] -= coord[1];
    globalcenter[2] -= coord[2];
    
    if ( distance > globalcenter.two_norm() ){
      distance = globalcenter.two_norm();
      index = temp_index;    
      //std::cout << "Index: " << index << " und die Distanz: " << distance << std::endl;
    }
    
  }
  
  //std::cout << "Index: " << index << " und die Distanz: " << distance << std::endl;
  
  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    int temp_index = gv.indexSet().index(*eit);
    //Dune::FieldVector<DF, dim> globalcenter, localcenter;
    if ( temp_index == index ) {
      CoordVec globalcenter, localcenter;
      globalcenter = eit->geometry().center();
      localcenter = eit->geometry().local(globalcenter);
      //Dune::FieldVector<RT, dim> gradient;
    
      RangeType temp_value = 0.0;
      dgf.evaluate(*eit,localcenter,temp_value);
      value = temp_value;
      break;
    }
    
  }
  
  //std::cout << "Index: " << index << " und die Distanz: " << distance << " und der Wert: " << value << std::endl;
  
  return value;
  
  
}

template<class GV, typename PGMap>
void determine_centers(const GV& gv, const PGMap& pg)
{
  typedef typename GV::Grid::ctype Coord;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  std::vector< int > anz_coords(N_ALPHA + N_MATA + N_CHEAT + 1,0);
  
  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    typedef typename GV::IntersectionIterator IntersectionIterator;
	
    for (IntersectionIterator is = gv.ibegin(*eit); is!= gv.iend(*eit); ++is)
    {
      if(is->boundary()) 
      {
	int physgroup_index = pg[is->boundarySegmentIndex()];
	//if ( physgroup_index > 0) 
	//{
	  anz_coords[physgroup_index] = anz_coords[physgroup_index] + 1;
	  Dune::FieldVector<Coord,3> global_face;
	  global_face = is->geometry().center();
	  coordsOfCells[physgroup_index][0] += global_face[0];
	  coordsOfCells[physgroup_index][1] += global_face[1];
	//}
      }
    }
  }
  
  
  for (int i = 0; i < N_ALPHA + N_MATA + N_CHEAT + 1; i++ )
  {
    coordsOfCells[i][0] = coordsOfCells[i][0]/((double)anz_coords[i]);
    coordsOfCells[i][1] = coordsOfCells[i][1]/((double)anz_coords[i]);
      
    std::cout << "Center Coords: (" << coordsOfCells[i][0] << "," << coordsOfCells[i][1] << ") of cell number " << i << std::endl;
  }
  
  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    typedef typename GV::IntersectionIterator IntersectionIterator;
	
    for (IntersectionIterator is = gv.ibegin(*eit); is!= gv.iend(*eit); ++is)
    {
      if(is->boundary()) 
      {
	int physgroup_index = pg[is->boundarySegmentIndex()];
	if ( physgroup_index > 0) 
	{
	  anz_coords[physgroup_index] = anz_coords[physgroup_index] + 1;
	  Dune::FieldVector<Coord,3> global_face, distance;
	  distance[0] = coordsOfCells[physgroup_index][0];
	  distance[1] = coordsOfCells[physgroup_index][1];
	  distance[2] = 0.0;
	  global_face = is->geometry().center();
	  distance -= global_face;
	  coordsOfCells[physgroup_index][2] = distance.two_norm();

	}
      }
    }
  }
  
  for (int i = 0; i < N_ALPHA + N_MATA + N_CHEAT + 1; i++ )
  {
      
    std::cout << "Radius: " << coordsOfCells[i][2] << " of cell number " << i << std::endl;
  }
}

template<class GV2D, class GV3D, class UDGF, typename Sol2D, typename RangeType>
void integrateOverHeight(const GV2D& gv2d, const GV3D& gv3d, const UDGF& udgf, Sol2D& c, RangeType& max_value)
{
    int qorder = QORDER;
  
    typedef typename GV3D::Grid::ctype Coord;
  
    typedef typename GV2D::template Codim<0>::Iterator ElementIterator2D;
    
    typedef typename GV3D::template Codim<0>::Iterator ElementIterator3D;
  
    std::vector< std::vector<int> > coordsToIndex(elements2d, std::vector<int>(elements2d,0));  
    std::vector< std::vector<double> > IndexToHeightAndValue(gv2d.size(0), std::vector<double>(Integration_steps,0.0));  
    double delta_x = 2.0*170.0/((double)elements2d), delta_y = 2.0*170.0/((double)elements2d), delta_z = Integration_height/((double)Integration_steps);
    
    std::cout << "\ndelta_x:" << delta_x;
    std::cout << "\ndelta_y:" << delta_y;
    std::cout << "\ndelta_z:" << delta_z << std::endl;
    
    for (ElementIterator2D eit = gv2d.template begin<0>();
	   eit!= gv2d.template end<0>(); ++eit)
    {
      int index = gv2d.indexSet().index(*eit);
      Dune::FieldVector<Coord,2> globalcenter;
      globalcenter = eit->geometry().center();

      int n_x = (int)(globalcenter[0]/delta_x);
      int n_y = (int)(globalcenter[1]/delta_y);
      coordsToIndex[n_x][n_y] = index;
    }
    
    for (ElementIterator3D eit = gv3d.template begin<0>();
	   eit!= gv3d.template end<0>(); ++eit)
    {
        Dune::GeometryType gt = eit->geometry().type();
        const Dune::QuadratureRule<Coord,3>& rule = Dune::QuadratureRules<Coord,3>::rule(gt,qorder);
	
	for (typename Dune::QuadratureRule<Coord,3>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
	  Dune::FieldVector<Coord,3> global,local;
	  local = it->position();
	  global = eit->geometry().global(local);
	  
	  if ( Integration_height < global[2] ) continue;
	
	  int n_x = (int)(global[0]/delta_x);
	  int n_y = (int)(global[1]/delta_y);
	  int n_z = (int)(global[2]/delta_z);
	
	  if ( (n_x >= elements2d) || (n_y >= elements2d) || (n_z >= Integration_steps) ) continue;
	  if ( (n_x < 0) || (n_y < 0) || (n_z < 0) ) continue;
	
	  RangeType value = 0.0;
	  udgf.evaluate(*eit,local,value);
	
	  if (max_value < value) max_value = value;
	
	  IndexToHeightAndValue[coordsToIndex[n_x][n_y]][n_z] = value; 
	}
	
    }
   
    /*for (ElementIterator3D eit = gv3d.template begin<0>();
	   eit!= gv3d.template end<0>(); ++eit)
    {
	Dune::FieldVector<Coord,3> global, local;
	global = eit->geometry().center();
	local = eit->geometry().local(global);
	
	if ( Integration_height < global[2] ) continue;
	
	//global[0] = global[0] + 170.0;
        //global[1] = global[1] + 170.0;
	int n_x = (int)(global[0]/delta_x);
        int n_y = (int)(global[1]/delta_y);
	int n_z = (int)(global[2]/delta_z);
	
	if ( (n_x > elements2d) || (n_y > elements2d) || (n_z > Integration_steps) ) continue;
	if ( (n_x < 0) || (n_y < 0) || (n_z < 0) ) continue;
	
	RangeType value = 0.0;
        udgf.evaluate(*eit,local,value);
	
	if (max_value < value) max_value = value;
	
	IndexToHeightAndValue[coordsToIndex[n_x][n_y]][n_z] = value; 
	
	for (unsigned int i = 0; i < 4; i++)
	{
	  global = eit->geometry().corner(i);
	  local = eit->geometry().local(global);
	
	  //global[0] = global[0] + 170.0;
	  //global[1] = global[1] + 170.0;
	  n_x = (int)(global[0]/delta_x);
	  n_y = (int)(global[1]/delta_y);
	  n_z = (int)(global[2]/delta_z);
	
	  //RangeType value = 0.0;
	  udgf.evaluate(*eit,local,value);
	  
	  if ( (n_x >= elements2d) || (n_y >= elements2d) || (n_z >= Integration_steps) ) continue;
	  if ( (n_x < 0) || (n_y < 0) || (n_z < 0) ) continue;
	
	  IndexToHeightAndValue[coordsToIndex[n_x][n_y]][n_z] = value; 
	  
	}
	
    }*/

    
    
    int not_count = 0;
    int elements_inside = 0;
    
    for (ElementIterator2D eit = gv2d.template begin<0>();
	   eit!= gv2d.template end<0>(); ++eit)
    { 
        int index = gv2d.indexSet().index(*eit);
	Dune::FieldVector<Coord,2> globalcenter, distance;
        globalcenter = eit->geometry().center();
	bool outside = false;
	
	//distance[0] = globalcenter[0] - 170.0;
	//distance[1] = globalcenter[1] - 170.0;
	
	for (int i = 1; i < N_ALPHA + N_MATA + N_CHEAT + 1; i++ )
	{
	  distance[0] = globalcenter[0] - coordsOfCells[i][0];
	  distance[1] = globalcenter[1] - coordsOfCells[i][1];
	  
	  if ( distance.two_norm() < coordsOfCells[i][2]) outside = true;
	  
	}
	
	distance[0] = globalcenter[0] - 170.0;
	distance[1] = globalcenter[1] - 170.0;
	
	if ( distance.two_norm() > 170.0 ) outside = true;
	
	c[index] = 0.0;
	
	if (outside)
	{  
	  c[index] = -0.000001;
	}
	else {
	  for ( int j = 0; j < Integration_steps; j++)
	  {
	    c[index] += delta_z*IndexToHeightAndValue[index][j];
	    elements_inside++;
	    if ( IndexToHeightAndValue[index][j] == 0.0 ) not_count++; 
	  }
	}
    }  
    
    std::cout << "\nThere are " << not_count << " coord without assigned value."
	      << "\nIf this value is high in comparions to "  <<  elements_inside//elements2d*elements2d*Integration_steps
	      << "\nthen improvement of algorithm/mesh is required!" 
	      << "\nThis equates to: " << (double)not_count*100.0/((double)elements_inside) << " % "<<std::endl;
}

template<class GV, class T,class C>
void alphaAbsorptionOverHeight(const GV& gv, T &absorption, const C& c)
{


  //
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::Grid::ctype Coord;

  typedef double Real;
  const int dim = GV::dimension;
  int physgroup_index =0; 
  std::vector<double> temp, circumference;
  std::vector<int> counter(N_MATA,0);
  std::vector<double> max_value(N_MATA,0.0), min_value(N_MATA,1.0);
  
  for (int i = 0; i < N_MATA; i++) {
    temp.push_back(0.0);
    double radius = coordsOfCells[i + N_ALPHA + N_CHEAT + 1][2];
    circumference.push_back(2.0*3.14*radius);
  }  

  for (ElementIterator eit = gv.template begin<0>();
	   eit!= gv.template end<0>(); ++eit)
  {
    int index = gv.indexSet().index(*eit);
    Dune::FieldVector<Coord,2> globalcenter, distance;
    globalcenter = eit->geometry().center();
    
    for (int i = 0; i < N_MATA; i++ )
    {
      distance[0] = globalcenter[0] - coordsOfCells[i + N_ALPHA + N_CHEAT + 1][0];
      distance[1] = globalcenter[1] - coordsOfCells[i + N_ALPHA + N_CHEAT + 1][1];
	  
      //bool outside = false;
      bool inside = false;
      //bool boundary = false;
      
      if ( distance.two_norm() < coordsOfCells[i + N_ALPHA + N_CHEAT + 1][2] + 0.3) inside = true;
      //if ( distance.two_norm() > coordsOfCells[i + N_ALPHA + N_CHEAT + 1][2]) outside = true;
      //if ( outside && inside ) boundary = true;
      
      if ( inside )
      {	
	if (c[index] > 0.0) {
	  counter[i]++;
	  temp[i] += c[index];
	  if (max_value[i] < c[index] ) max_value[i] = c[index];
	  if (min_value[i] > c[index] ) min_value[i] = c[index];
	}
      }
      
    }
	
  }
  
  for (int i = 0; i < N_MATA; i++) {
    temp[i] = temp[i]/((double)counter[i]);
  }  

  for (int i = 0; i < N_MATA; i++) {
    //temp[i] = temp[i]/circumference[i];	
    //test << i << "\t" << circumference[i] << std::endl;
    std::cout << "Cell Nr.:" << i << "\t with surface area:" << circumference[i] << std::endl;
    std::cout << "\t\tAbsorbtion(Integrated over h):" << temp[i] << std::endl;
    std::cout << "\t\tMax value:" << max_value[i] << std::endl;
    std::cout << "\t\tMin value:" << min_value[i] << std::endl;
    //std::cout << "\t\tAverage value:" << (max_value[i] + min_value[i])/2.0 << std::endl;
  } 
   
  absorption.push_back(temp);
}
