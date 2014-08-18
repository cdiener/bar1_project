
/** \brief function for defining the diffusion tensor for Bar1
 *  \class K1
 * 
 *  \tparam GV grid view allows read only access to all grid elements
 *  \tparam RF range field type
 */ 

template<typename GV, typename RF>
class K1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K1<GV,RF> >
{
  const GV& gv;
  RF time;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K1<GV,RF> > BaseT;

  K1 (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
          y[i][i] = D1;
        else
          y[i][j] = 0.0;
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }

  //! set time for subsequent evaluation
  void setTime (double t) {time = t;}
};

