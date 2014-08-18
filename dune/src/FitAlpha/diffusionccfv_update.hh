/*
 * diffusionccfv_update.hh
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



// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_DIFFUSIONCCFV_HH
#define DUNE_PDELAB_DIFFUSIONCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/genericreferenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>

/*#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"pattern.hh"
#include"flags.hh"*/
#include"dune/pdelab/common/geometrywrapper.hh"
#include"dune/pdelab/gridoperatorspace/gridoperatorspace.hh"
#include"dune/pdelab/localoperator/pattern.hh"
#include"dune/pdelab/localoperator/flags.hh"
#include<dune/pdelab/localoperator/idefault.hh>
//#include "diffusionparam.hh"


namespace Dune {
  namespace PDELab {
    
 // \f{align*} {
 //            &- div (k(x) grad u) + a0*u = f \text{ in }  \Omega, \\
 //            &                u = g \text{ on }  \partial \Omega_D \\
 //           & - k(x) grad u * \nu = j \text{ on } \partial \Omega_N 
 // \f}


/** \brief the diffusion operator for calculating the stiffness matrix with cell centered finite volumes
 * 
 * 
 * \class DiffusionCCFV
 * 
 * \tparam K diffusion tensor
 * \tparam F function for the source term
 * \tparam B function for selecting the boundary condition, e.g. Dirichlet or Neumann 
 * \tparam G function the Dirichlet boundary condition 
 * \tparam J function the Neumann
 * \tparam A0 function defining the degradation 
 * 
 */

    template<typename K, typename A0, typename F, typename B, typename J, typename G>
	class DiffusionCCFV : public NumericalJacobianApplySkeleton<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianApplyBoundary<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianApplyVolume<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianSkeleton<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianBoundary<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianVolume<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public FullSkeletonPattern, 
                          public FullVolumePattern,
                          public LocalOperatorDefaultFlags

	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

	  // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = false };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = false };

      DiffusionCCFV (const K& k_, A0& a0_, const F& f_, const B& b_, const J& j_, const G& g_, const double Diffusion_)
        : k(k_), a0(a0_), f(f_), b(b_), j(j_), g(g_), Diffusion(Diffusion_)
      {}


	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate Helmholtz term
        //typename A0::Traits::RangeType a0value;
        //a0.evaluate(eg.entity(),inside_local,a0value);

        // evaluate u
        std::vector<RangeType> degradation_local(lfsu.size());
        lfsu.vread(a0, degradation_local);
        RF y=0.0;
        //for (size_t i=0; i<lfsu.size(); i++)
           y = degradation_local[0];

        r[0] += y*x[0]*eg.geometry().volume();
	//std::cout << "Degradation: " << degradation_local[0] << "Wert y: "<< y << std::endl;
	
        // evaluate source term
        typename F::Traits::RangeType fvalue;
        f.evaluate(eg.entity(),inside_local,fvalue);

        r[0] -= fvalue*eg.geometry().volume();
      }

	  // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n) const
	  {
 		// domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
 		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
 		  Traits::LocalBasisType::Traits::RangeFieldType RF;

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();

        // cell centers in references elements
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // evaluate diffusion coefficient
        //typename K::Traits::RangeType k_inside, k_outside;
        //k.evaluate(*(ig.inside()),inside_local,k_inside);
        //k.evaluate(*(ig.outside()),outside_local,k_outside);
        //typename K::Traits::RangeType k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.outside()->geometry().global(outside_local);

        // distance between the two cell centers
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
        
        // contribution to residual on inside element, other residual is computed by symmetric call
        //r_s[0] += k_avg*(x_s[0]-x_n[0])*face_volume/distance;
        //r_n[0] -= k_avg*(x_s[0]-x_n[0])*face_volume/distance;
        r_s[0] += Diffusion*(x_s[0]-x_n[0])*face_volume/distance;
        r_n[0] -= Diffusion*(x_s[0]-x_n[0])*face_volume/distance;
	  }

	  // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save som e geometry evaluations
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
	  {
 		// domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
 		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
 		  Traits::LocalBasisType::Traits::RangeFieldType RF;

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();
        
        // cell center in reference element
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);

        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,face_local,bctype);

        // distance between cell center and face center
        Dune::FieldVector<DF,IG::dimension> 
              inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension> 
              outside_global = ig.geometry().global(face_local);
              inside_global -= outside_global;
        RF distance = inside_global.two_norm();
        const Dune::FieldVector<DF,IG::dimension-1>& face_center =
                  Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        if ( bctype == 1)
          {
            // Dirichlet boundary

            
            // evaluate diffusion coefficient
            //typename K::Traits::RangeType k_inside;
            //k.evaluate(*(ig.inside()),inside_local,k_inside);
            
            // evaluate boundary condition function
            typename G::Traits::DomainType x = ig.geometryInInside().global(face_local);
            typename G::Traits::RangeType y;
            g.evaluate(*(ig.inside()),x,y);

    
            // contribution to residual on inside element
            //r_s[0] += k_inside*(x_s[0]-y[0])*face_volume/distance;
            r_s[0] += Diffusion*(x_s[0]-y[0])*face_volume/distance;
          }
        else // if (DiffusionBoundaryCondition::isNeumann(bctype))
          {
            // Neumann boundary
            // evaluate flux boundary condition
            typename J::Traits::RangeType jvalue;
            //j.evaluate(*(ig.inside()),inside_local,jvalue);
            j.evaluate(ig,face_center,jvalue);
	    // int physgroup_index = pg[ig.intersection().boundarySegmentIndex()];
            // contribution to residual on inside element
            r_s[0] += jvalue*face_volume;
          }
	  }
	  
    void update (A0& up) {a0 = up;}

    private:
      const K& k;
      A0& a0;
      const F& f;
      const B& b;
      const J& j;
      const G& g;
      const double Diffusion;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
