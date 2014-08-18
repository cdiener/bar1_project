#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

/**
 * \brief function that calculates the absolute value of the concentration change integrated over the whole computaional domain 
 * from the recent time step to the preceding time step
 * 
 * @param[in] unew 	solution from the recent time step
 * @param[in] uold 	solution from the preceding time step
 * @param[in] GFS  	Grid Function Space
 * @param[in] qorder  	order of the quadrature formular used for integration
 * 
 * \return absolute value of the concentration change integrated over the whole computaional domain 
 * 
 */


template<class U, class GFS> 
double concentrationChange (const U& unew, const GFS& gfs, U& uold, int qorder) {

  // constants and types
  typedef typename GFS::Traits::GridViewType GV;
  const int dim = GV::Grid::dimension;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GFS::Traits::LocalFiniteElementType::
	Traits::LocalBasisType::Traits FTraits;
  typedef typename FTraits::DomainFieldType D;
  typedef typename FTraits::RangeFieldType R;
  typedef typename FTraits::RangeType RangeType;
  
  // Make Discrete Grid Function from coefficent vector 
  typedef typename GFS::template VectorContainer<R>::Type V;
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS,V> DGFD;
  
  DGF dgfold(gfs,unew);
  DGF dgfnew(gfs,uold);

  // make local function space
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);                            /*@\label{l2int:lfs}@*/
  std::vector<R> xlnew(lfs.maxSize()),xlold(lfs.maxSize());        // local coefficients
  std::vector<RangeType> b(lfs.maxSize()); // shape function values
  
  // loop over grid view
  double sumAbs = 0.0;

  //typename DU::Traits::RangeType sumGrad2;
  std::vector<double> sumGrad2(dim);

  for (ElementIterator eit = gfs.gridview().template begin<0>();
	   eit!=gfs.gridview().template end<0>(); ++eit)
  {
    lfs.bind(*eit);      // bind local function space to element /*@\label{l2int:bind}@*/
    lfs.vread(unew,xlnew);     // get local coefficients
    lfs.vread(uold,xlold); 
    Dune::GeometryType gt = eit->geometry().type();
    const Dune::QuadratureRule<D,dim>&  /*@\label{l2int:quad}@*/
    rule = Dune::QuadratureRules<D,dim>::rule(gt,qorder);

    for (typename Dune::QuadratureRule<D,dim>::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
    {

      // evaluate the given grid function at integration point
      RangeType u_new,u_old;
      dgfnew.evaluate(*eit,qit->position(),u_new);
      dgfold.evaluate(*eit,qit->position(),u_old);

      //accumulate error
      u_new -= u_old;
     
      //Du_approx -= Du_given;
      sumAbs+= u_new.one_norm()*qit->weight()*eit->geometry().integrationElement(qit->position());
       
     }
  }
  
  return sumAbs;
}
