namespace contrib{
  namespace lawlerem {
    template<class Type>
    struct sqrtm_typedefs{
      typedef Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > SAES;
    };

    TMB_ATOMIC_VECTOR_FUNCTION(
      // ATOMIC_NAME
      sqrtm
      ,
      // OUTPUT_DIM
      tx.size()
      ,
      // ATOMIC_DOUBLE
      typedef sqrtm_typedefs<double>::SAES SAES_t;
      int n=sqrt((double)tx.size());
      matrix<double> X=atomic::vec2mat(tx, n, n);
      SAES_t saes(X);
      matrix<double> sqrtX = saes.operatorSqrt();
      for(int i = 0; i < n * n; i++) ty[i]=sqrtX(i);
      ,
      // ATOMIC_REVERSE ( vec(f'(X)) =  (f(X)^T %kronecker sum% f(X))^-1 * vec(X') )
      int n = sqrt((double)ty.size());
      matrix<Type> Y = atomic::vec2mat(ty, n, n); // f(X)
      matrix<Type> Yt = Y.transpose(); // f(x)^T
      matrix<Type> I(Y.rows(), Y.cols());
      I.setIdentity();
      matrix<Type> kronSum = kronecker(Yt, I) + kronecker(I, Y);
      matrix<Type> kronSumInv = atomic::matinv(kronSum);
      px = kronSumInv * vector<Type>(py);
    )

    template<class Type>
    matrix<Type> sqrtm(matrix<Type> x){
      int n=x.rows();
      return atomic::vec2mat(sqrtm(atomic::mat2vec(x)), n, n);
    }
  }
}