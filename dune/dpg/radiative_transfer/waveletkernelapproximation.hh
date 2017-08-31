// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH

#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

namespace Dune {

namespace ScatteringKernelApproximation {

  namespace AlpertWavelet {

    void add(Eigen::VectorXd& target,
             const Eigen::VectorXd& data,
             const int pos0)
    {
      for(int i=0; i<data.size(); i++){
        target(pos0+i)=data(i);
      }
    }

    Eigen::VectorXd extract(const Eigen::VectorXd& data,
                            size_t pos0,
                            size_t pos1)
    {
      Eigen::VectorXd result(pos1-pos0);
      for (size_t i = pos0; i < pos1; i++)
      {
          result(i-pos0) = data(i);
      }
      return result;
    }

    double ip(const Eigen::VectorXd& f,
              const Eigen::VectorXd& g,
              const Eigen::VectorXd& quadWeight,
              double xmin, double xmax)
    {
      return 0.5*(xmax-xmin)*(quadWeight.array()*f.array()*g.array()).sum();
    }

    double l2norm(const Eigen::VectorXd& f,
                  const Eigen::VectorXd& quadWeight,
                  double xmin, double xmax)
    {
      double s = 0.5*(xmax-xmin)*(quadWeight.array()*pow(f.array(),2)).sum();
      return std::sqrt(s);
    }

    // Get Gauss-Legendre quadrature rule of order quadOrder in [xmin,xmax]
    void getLegendreQuad(Eigen::VectorXd& quadPos,
                         Eigen::VectorXd& quadWeight,
                         const size_t quadOrder,
                         double xmin,
                         double xmax)
    {
      // Get Gauss-Legendre quadrature rule of order quadOrder in [0,1]
      const Dune::QuadratureRule<double, 1>& quad =
            Dune::QuadratureRules<double, 1>::rule(Dune::GeometryType::simplex,
                                                     quadOrder,
                                                     QuadratureType::GaussLegendre);
      quadPos.resize(quad.size());
      quadWeight.resize(quad.size());

      for ( size_t iQuad=0; iQuad < quad.size(); iQuad++ ) {
        quadPos(iQuad)=quad[iQuad].position()*(xmax-xmin)+xmin;
        quadWeight(iQuad) = quad[iQuad].weight()*(xmax-xmin);
      }
      return;
    }

    Eigen::VectorXd getLagrangePoly(const Eigen::VectorXd& quadPos,
                                    const Eigen::VectorXd& xInterp,
                                    int index)
    {
      Eigen::VectorXd P=Eigen::VectorXd::Ones(quadPos.size());
      for(int i=0; i<xInterp.size(); i++){
        if(i==index) continue;
        // else P = P*(quadPos-xInterp(i))/(xInterp(index)-xInterp(i));
        else P = P.array()*(quadPos.array()-xInterp(i))/(xInterp(index)-xInterp(i));
      }
      return P;
    }

    // Legendre Poly orthonormalized in L2([-1,1])
    Eigen::VectorXd getLegendrePoly(const Eigen::VectorXd& quadPos,
                                    const size_t degree)
    {
      Eigen::VectorXd P=Eigen::VectorXd::Zero(quadPos.size());

      if(degree==0) P=Eigen::VectorXd::Ones(quadPos.size());
      else{
        if(degree==1) P=quadPos;
        else{
          Eigen::VectorXd P_pp=Eigen::VectorXd::Ones(quadPos.size());
          Eigen::VectorXd P_p=quadPos;

          for(size_t deg=2; deg<=degree; deg++){
            P=((2*deg-1)*quadPos.array()*P_p.array()-(deg-1)*P_pp.array())/deg;
            P_pp= P_p;
            P_p= P;
          }
        }
      }
      double l2NormP=std::sqrt(2./(2.*degree+1.));
      return P/l2NormP;
    }

    bool sort_criterion(const std::pair<Eigen::VectorXd,double>& left,
                   const std::pair<Eigen::VectorXd,double>& right)
    {
      return left.second > right.second;
    }

    void alpert_sort(std::vector<Eigen::VectorXd>& F,
                    const Eigen::VectorXd& P,
                    const Eigen::VectorXd& quadWeight,
                    double xmin, double xmax)
    {
      std::vector<std::pair<Eigen::VectorXd,double>> Fpair;
      for(size_t l=0; l<F.size(); l++){
        double cF = abs(ip(F[l],P,quadWeight,xmin,xmax));
        Fpair.push_back(std::pair<Eigen::VectorXd,double>(F[l],cF));
      }
      std::sort(Fpair.begin(),Fpair.end(),sort_criterion);
      for(size_t l=0; l<F.size(); l++){
        F[l]=Fpair[l].first;
      }
    }

    std::vector<Eigen::VectorXd> gram_schmidt(
      const std::vector<Eigen::VectorXd>& F,
      const Eigen::VectorXd& quadWeight,
      double xmin, double xmax)
    {
      std::vector<Eigen::VectorXd> F_ortho(F.size());
      F_ortho[0]=F[0]/l2norm(F[0],quadWeight,xmin,xmax);
      if(F.size()>1){
        for(size_t i=1; i<F.size(); i++){
          F_ortho[i]=F[i];
          for(size_t j=0; j<i; j++){
            F_ortho[i]-=ip(F_ortho[i],F_ortho[j],quadWeight,xmin,xmax)*F_ortho[j];
          }
          F_ortho[i]=F_ortho[i]/l2norm(F_ortho[i],quadWeight,xmin,xmax);
        }
      }
      return F_ortho;
    }

    std::vector<Eigen::VectorXd> orthogonalize_wrt_high_order_monomials(
      std::vector<Eigen::VectorXd>& F,
      std::vector<Eigen::VectorXd>& P,
      const Eigen::VectorXd& quadWeight,
      double xmin, double xmax)
    {
      if(F.size()==1) return F;
      else{
        alpert_sort(F,P[0],quadWeight,xmin,xmax);

        std::vector<Eigen::VectorXd> Fsubset(F.size()-1);
        std::vector<Eigen::VectorXd> Psubset(F.size()-1);
        for(size_t l=1;l<F.size();l++){
          double c=ip(F[l],P[0],quadWeight,xmin,xmax)/ip(F[0],P[0],quadWeight,xmin,xmax);
          Fsubset[l-1]=F[l]-c*F[0];
          Psubset[l-1]=P[l];
        }
        Fsubset=orthogonalize_wrt_high_order_monomials(Fsubset,Psubset,quadWeight,xmin,xmax);
        for(size_t l=1;l<F.size();l++){
          F[l]=Fsubset[l-1];
        }
        return F;
      }
    }

    Eigen::VectorXd orthogonalize_wrt_space_poly(
      const Eigen::VectorXd& f,
      size_t L,
      const Eigen::VectorXd& quadPos,
      const Eigen::VectorXd& quadWeight,
      double xmin, double xmax)
    {
      Eigen::VectorXd f_ortho=f;
      for(size_t l=0; l<L; l++){
        Eigen::VectorXd legendrePolyNormalized
          = getLegendrePoly(quadPos,l)
          /l2norm(getLegendrePoly(quadPos,l),quadWeight,xmin,xmax);
        f_ortho-=ip(f,legendrePolyNormalized,quadWeight,xmin,xmax)*legendrePolyNormalized;
      }
      return f_ortho/l2norm(f_ortho,quadWeight,xmin,xmax);
    }

    // Computes Alpert wlt with L vanishing moments in the interval [xmin,xmax]
    // The wlt functions are evaluated at Gauss-Legendre quad point in [xmin,xmax]
    // quadPos, quadWeight is a quad rule in [-1,1].
    std::vector<Eigen::VectorXd> getAlpertWlt(
      size_t L,double xmin,double xmax,
      const Eigen::VectorXd& quadPos,
      const Eigen::VectorXd& quadWeight)
    {
      // Translate x values from the interval [-1, 1] to [xmin, xmax]
      Eigen::VectorXd x=0.5*(xmax-xmin)*(quadPos.array()+1)+xmin;

      // Step 0: Create Functions F=[f_1^1,...,f_k^1] and P=[1.,x,...x^{K-1}] in [xmin,xmax]
      std::vector<Eigen::VectorXd> F(L);
      std::vector<Eigen::VectorXd> P_higher(L);

      for(size_t l=0; l<L; l++){
        P_higher[l]=pow(x.array(),l+L);
        F[l]=pow(x.array(),l);
        for(Eigen::Index i=0;i<x.size();i++){
          if(x(i)<= (xmin+xmax)/2.) F[l](i)=-pow(x(i),l);
        }
      }

      // Step 1: Orthogonalize [f_1^1,...,f_k^1] wrt P=[1,x,..,x^{k-1}].
      // This yields [f_1^2,...,f_k^2] (store in F1)
      std::vector<Eigen::VectorXd> F1(L);
      for(size_t l=0; l<L; l++){
        F1[l]= orthogonalize_wrt_space_poly(F[l],L,quadPos,quadWeight,xmin,xmax);
      }

      // Step 2: Transform [f_1^2,...,f_k^2] --> [f_1^2,...,f_k^{k+1}]
      // such that < f_j^{j+1},x^i >=0 for i <= j+k-2
      // Remark: To have wlt with a certain number of vanishing moments, this step is in theory not required. Alpert added it to make his construction unique. It is in theory not required. I did not observe any significant difference in the MRA of the kernel if it is not done.
      std::vector<Eigen::VectorXd> F2=orthogonalize_wrt_high_order_monomials(F1,P_higher,quadWeight,xmin,xmax);

      // Step 3: Gram-Schmidt orthonormalization of F2.
      std::reverse(F2.begin(),F2.end());

      return gram_schmidt(F2,quadWeight,xmin,xmax);
    }

    std::vector<Eigen::MatrixXd> get_alpert_transform_matrices(
        const size_t L,
        const size_t quadOrder)
    {
      Eigen::MatrixXd F1 = Eigen::MatrixXd::Zero(L,L);
      Eigen::MatrixXd F2 = Eigen::MatrixXd::Zero(L,L);
      Eigen::MatrixXd G1 = Eigen::MatrixXd::Zero(L,L);
      Eigen::MatrixXd G2 = Eigen::MatrixXd::Zero(L,L);

      // Get Gauss-Legendre quadrature rule of order quadOrder in [-1,0]
      Eigen::VectorXd quadPos0,quadWeight0;
      getLegendreQuad(quadPos0,quadWeight0,quadOrder,-1.,0.);
      // Get Gauss-Legendre quadrature rule of order quadOrder in [0,1]
      Eigen::VectorXd quadPos1,quadWeight1;
      getLegendreQuad(quadPos1,quadWeight1,quadOrder,0.,1.);
      // Get Gauss-Legendre quadrature rule of order quadOrder in [-1,1]
      Eigen::VectorXd quadPos2,quadWeight2;
      getLegendreQuad(quadPos2,quadWeight2,quadOrder,-1.,1.);
      // Concatenate quadratures 0 and 1 to get quadrature in [-1,1]
      Eigen::VectorXd quadPos(quadPos0.size()+quadPos1.size()),
                      quadWeight(quadPos0.size()+quadPos1.size());
      quadPos << quadPos0, quadPos1;
      quadWeight << quadWeight0, quadWeight1;
      // Get Alpert wavelets
      std::vector<Eigen::VectorXd> alpertWlt =
        getAlpertWlt(L,-1.,1.,quadPos,quadWeight);

      for(size_t li=0; li<L; li++){
        Eigen::VectorXd alpert_shift0=alpertWlt[li].head(quadPos0.size());
        Eigen::VectorXd alpert_shift1=alpertWlt[li].tail(quadPos1.size());

        Eigen::VectorXd legendre_shift0=getLegendrePoly(quadPos0,li);
        Eigen::VectorXd legendre_shift1=getLegendrePoly(quadPos1,li);

        for(size_t lj=0; lj<L; lj++){
          Eigen::VectorXd legendre_t=getLegendrePoly(quadPos2,lj);
          F1(li,lj)=ip(legendre_shift0,legendre_t,quadWeight2,-1.,1.);
          F2(li,lj)=ip(legendre_shift1,legendre_t,quadWeight2,-1.,1.);
          G1(li,lj)=ip(alpert_shift0,legendre_t,quadWeight2,-1.,1.);
          G2(li,lj)=ip(alpert_shift1,legendre_t,quadWeight2,-1.,1.);
        }
      }
      return {F1/std::sqrt(2.),F2/std::sqrt(2.),G1/std::sqrt(2.),G2/std::sqrt(2.)};
    }

    Eigen::MatrixXd
    get_lagrange_to_legendre_matrix(const Eigen::VectorXd& xInterp,
        double xmin, double xmax, double quadOrder)
    {

      int L=xInterp.size();
      Eigen::MatrixXd A(L,L);

      // Get Gauss-Legendre quadrature rule of order quadOrder in [-1,1]
      Eigen::VectorXd quadPos,quadWeight;
      getLegendreQuad(quadPos,quadWeight,quadOrder,-1.,1.);

      // Translate quadPos values from the interval [-1, 1] to [xmin, xmax]
      Eigen::VectorXd quadPosShift=0.5*(xmax-xmin)*(quadPos.array()+1)+xmin;

      for(int lx=0; lx<L; lx++){
        Eigen::VectorXd legendrePoly=getLegendrePoly(quadPos,lx);
        legendrePoly*=std::sqrt(2./(xmax-xmin)); // normalize in [xmin,xmax]
        for(int ly=0; ly<L; ly++){
          Eigen::VectorXd lagrangePoly=getLagrangePoly(quadPosShift,xInterp,ly);
          A(lx,ly)=ip(legendrePoly,lagrangePoly,quadWeight,xmin,xmax);
        }
      }
      return A;
    }

    template<class lambdaExpr> Eigen::VectorXd
    ProjectOntoVJ(const lambdaExpr& f,
                  double r, int J, int L, int quadOrder)
    {
      Eigen::VectorXd data=Eigen::VectorXd::Zero(static_cast<int>(L*pow(2,J)));
      double xmin;
      double xmax;
      for(int k=0; k<static_cast<int>(pow(2,J)); k++){
        xmin=-r+r*k*pow(2.,-J+1);
        xmax=-r+r*(k+1)*pow(2.,-J+1);
        Eigen::VectorXd x=Eigen::VectorXd::LinSpaced(L,xmin,xmax);
        Eigen::VectorXd fx(L);
        for(int l=0; l<L; l++) fx(l)=f(x(l),xmin,xmax);
        Eigen::MatrixXd A
          = get_lagrange_to_legendre_matrix(x,xmin,xmax,quadOrder);
        Eigen::VectorXd datak=A*fx;
        data.segment(L*k,L)=datak;
      }
      return data;
    }

    Eigen::VectorXd
    lagrange_to_legendre_VJ(Eigen::VectorXd& f, int maxLevel, int waveletOrder)
    {
      // Get Gauss-Legendre quadrature in [0,1]
      int quadOrder=2*waveletOrder+1;
      GeometryType type;
      type.makeLine();
      const size_t dim = 1;
      const Dune::QuadratureRule<double, dim>& quad =
        Dune::QuadratureRules<double, dim>::rule(type,
          quadOrder, QuadratureType::GaussLegendre);

      // Auxiliary variables
      double xmin;
      double xmax;
      Eigen::VectorXd angle(waveletOrder+1);
      Eigen::VectorXd data((waveletOrder+1) << maxLevel);
      using namespace boost::math::constants;
      double r = pi<double>();
      for(int k = 0; k < (1<<maxLevel); ++k)
      {
        for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
          // Angle in [-pi,pi]
          angle(pt) = 2*r*(k+quad[pt].position())/(1<<maxLevel) - r;
        }
        xmin = r *   k   * std::exp2(-(double)maxLevel+1) - r;
        xmax = r *  (k+1)* std::exp2(-(double)maxLevel+1) - r;
        Eigen::MatrixXd A
          = get_lagrange_to_legendre_matrix(angle,xmin,xmax,quadOrder);
        Eigen::VectorXd datak=
          A*(f.segment((waveletOrder+1)*k,(waveletOrder+1)));
        data.segment((waveletOrder+1)*k,waveletOrder+1)=datak;
      }
      return data;
    }

    Eigen::VectorXd
    legendre_to_lagrange_VJ(Eigen::VectorXd& f, int maxLevel, int waveletOrder)
    {
      // Get Gauss-Legendre quadrature in [0,1]
      int quadOrder=2*waveletOrder+1;
      GeometryType type;
      type.makeLine();
      const size_t dim = 1;
      const Dune::QuadratureRule<double, dim>& quad =
        Dune::QuadratureRules<double, dim>::rule(type,
          quadOrder, QuadratureType::GaussLegendre);

      // Auxiliary variables
      double xmin;
      double xmax;
      Eigen::VectorXd angle(waveletOrder+1);
      Eigen::VectorXd data((waveletOrder+1) << maxLevel);
      using namespace boost::math::constants;
      double r = pi<double>();
      for(int k = 0; k < (1<<maxLevel); ++k)
      {
        for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
          // Angle in [-pi,pi]
          angle(pt) = 2*r*(k+quad[pt].position())/(1<<maxLevel) - r;
        }
        xmin = r *   k   * std::exp2(-(double)maxLevel+1) - r;
        xmax = r *  (k+1)* std::exp2(-(double)maxLevel+1) - r;
        Eigen::MatrixXd A
          = get_lagrange_to_legendre_matrix(angle,xmin,xmax,quadOrder);
        Eigen::VectorXd datak =
          A.fullPivLu().solve(f.segment((waveletOrder+1)*k,(waveletOrder+1)));
        data.segment((waveletOrder+1)*k,waveletOrder+1)=datak;
      }
      return data;
    }

    // Computes direct DWT of function vJ with coefs of scaling functions
    // expressed in an ONB basis of V_J. For us, the basis is built with
    // Legendre polynomials normalized in the intervals I_{J,k}
    // The function returns a pair for which:
    // pair.first     --> components in V_0
    // pair.second[j] --> Components in W_j (0 <= j <= J)
    std::pair<Eigen::VectorXd,std::vector<Eigen::VectorXd>>
    DWT(const Eigen::VectorXd& data,
        size_t L,
        size_t J,
        size_t quadOrder)
    {
      if(J==0){
        std::vector<Eigen::VectorXd> w(1);
        w[0]=data;
        return std::make_pair(data,std::vector<Eigen::VectorXd>());
      }
      else{
        // Initializations
        std::vector<Eigen::VectorXd> s(J+1); s[J]=data;
        std::vector<Eigen::VectorXd> w(J);
        Eigen::VectorXd c0(L);
        Eigen::VectorXd c1(L);
        // Get elementary wavelet transform matrices
        std::vector<Eigen::MatrixXd> elem_matrix
                = get_alpert_transform_matrices(L,quadOrder);

        for(int j=J-1; j>=0; j--){
          s[j].resize(static_cast<int>(L*pow(2,j)));
          w[j].resize(static_cast<int>(L*pow(2,j)));
          for(int k=0; k<pow(2.,j); k++){
            c0=s[j+1].segment(L*2*k,L);
            c1=s[j+1].segment(L*(2*k+1),L);
            s[j].segment(k*L,L)=elem_matrix[0]*c0+elem_matrix[1]*c1;
            w[j].segment(k*L,L)=elem_matrix[2]*c0+elem_matrix[3]*c1;

          }
        }
        return std::make_pair(s[0],w);
      }
    }

    // Computes IDWT of a Multiresolution Analysis stored in w.
    // w.first     --> components in V_0
    // w.second[j] --> Components in W_j (0 <= j <= J)
    Eigen::VectorXd
    IDWT(std::pair<Eigen::VectorXd,std::vector<Eigen::VectorXd>>& w,
         size_t L,
         size_t J,
         size_t quadOrder)
    {
      // Get elementary wavelet transform matrices
      std::vector<Eigen::MatrixXd>
        elem_matrix = get_alpert_transform_matrices(L,quadOrder);
      std::vector<Eigen::VectorXd> s(J+1);
      s[0]=w.first;
      for(size_t j=1; j<=J; j++){
        s[j].resize(static_cast<int>(L*pow(2,j)));
        Eigen::VectorXd stmp(L);
        Eigen::VectorXd wtmp(L);
        for(int k=0; k<pow(2.,j); k++){
          int kdiv2=(int) k/2;
          int kmod2=k%2;
          stmp=s[j-1].segment(L*kdiv2,L);
          wtmp=w.second[j-1].segment(L*kdiv2,L);
          if(kmod2==0){
            s[j].segment(L*k,L)=elem_matrix[0].transpose()*stmp
                                +elem_matrix[2].transpose()*wtmp;
          }
          else{
            s[j].segment(L*k,L)=elem_matrix[1].transpose()*stmp
                                +elem_matrix[3].transpose()*wtmp;
          }
        }
      }
      return s[J];
    }

    template<class Function>
    Eigen::MatrixXd waveletKernelMatrix(const Function& kernelFunction,
                                        size_t wltOrder,
                                        size_t maxLevel,
                                        size_t nQuadAngle)
      {
        using namespace Eigen;
        using namespace boost::math::constants;

        const double r = pi<double>();
        const size_t num_s = (wltOrder+1) << maxLevel;
        MatrixXd kernelMatrix(num_s, num_s);

        // Get Gauss-Legendre quadrature rule of nQuadAngle points
        // i.e. order=2*nQuadAngle-1. Interval: [-1,1]
        Eigen::VectorXd quadPos,quadWeight;
        size_t quadOrder=2*nQuadAngle-1;
        getLegendreQuad(quadPos,quadWeight,quadOrder,-1.,1.);

        // Get Alpert wavelets in [-1,1]
        std::vector<Eigen::VectorXd> alpertWlt1 =
          getAlpertWlt(wltOrder+1,-1.,1.,quadPos,quadWeight);

        // Auxiliary variables to use in loop
        auto w = quadWeight * quadWeight.transpose();

        Eigen::VectorXd one =
          Eigen::VectorXd::Constant(quadPos.size(),1.);
        Eigen::VectorXd x(quadPos.size());
        Eigen::VectorXd y(quadPos.size());

        // Matrix. Terms <k, sf * sf'>
        for(size_t ly = 0; ly < wltOrder+1; ly++) {
          Eigen::VectorXd sfy
            = getLegendrePoly(quadPos,ly);
            double ymin = -r;
            double ymax = r;
            y = 0.5*(quadPos + one)*(ymax - ymin) + ymin*one;
            auto Y = y * one.transpose();

          for(size_t lx = 0; lx < wltOrder+1; lx++) {
            Eigen::VectorXd sfx
              = getLegendrePoly(quadPos,lx);

            double xmin = -r;
            double xmax = r;
            x = 0.5*(quadPos + one)*(xmax - xmin) + xmin*one;
            auto X = one * x.transpose();

            auto diff = X - Y;
            auto Keval = diff.unaryExpr(kernelFunction);

            auto prod = r*w.cwiseProduct(
              Keval.cwiseProduct( sfy * sfx.transpose()));
            kernelMatrix(lx, ly) = prod.sum();
          }
        }

        // Matrix. Terms <k, sf * wlt'> and <k, sf' * wlt>
        for(size_t ly = 0; ly < wltOrder+1; ly++) {
          Eigen::VectorXd sfy
            = getLegendrePoly(quadPos,ly);
          double ymin = -r;
          double ymax = r;
          y = 0.5*(quadPos + one)*(ymax - ymin) + ymin*one;
          auto Y = y * one.transpose();

          for(size_t jx = 0; jx < maxLevel; jx++) {
            for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {

              double xmin = r *   kx   * std::exp2(-(double)jx+1) - r;
              double xmax = r * (kx+1) * std::exp2(-(double)jx+1) - r;
              x = 0.5*(quadPos + one)*(xmax - xmin) + xmin*one;
              auto X = one * x.transpose();

              auto diff = X - Y;
              auto Keval = diff.unaryExpr(kernelFunction);

              for(size_t lx = 0; lx < wltOrder+1; lx++) {
                  auto prod = r/std::pow(2,(double)jx/2)
                  * w.cwiseProduct( Keval.cwiseProduct(
                  sfy * alpertWlt1[lx].transpose()));
                const size_t i = ((1 << jx) + kx)*(wltOrder+1)+lx;
                const size_t j = ly;
                kernelMatrix(i, j) = prod.sum();
                kernelMatrix(j, i) = kernelMatrix(i, j);
              }
            }
          }
        }

        // Matrix. Terms <k, wlt * wlt'>
        for(size_t jy = 0; jy < maxLevel; jy++) {
          for(size_t ky = 0, ky_max = 1 << jy; ky < ky_max; ky++) {

            double ymin = r *   ky   * std::exp2(-(double)jy+1) - r;
            double ymax = r * (ky+1) * std::exp2(-(double)jy+1) - r;
            y = 0.5*(quadPos + one)*(ymax - ymin) + ymin*one;
            auto Y = y * one.transpose();

            for(size_t jx = 0; jx < maxLevel; jx++) {
              for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {

                double xmin = r *   kx   * std::exp2(-(double)jx+1) - r;
                double xmax = r * (kx+1) * std::exp2(-(double)jx+1) - r;
                x = 0.5*(quadPos + one)*(xmax - xmin) + xmin*one;
                auto X = one * x.transpose();

                auto diff = X - Y;
                auto Keval = diff.unaryExpr(kernelFunction);

                for(size_t ly = 0; ly < wltOrder+1; ly++) {
                  for(size_t lx = 0; lx < wltOrder+1; lx++) {
                    auto prod = r/std::pow(2,(double)(jx+jy)/2)
                       * w.cwiseProduct( Keval.cwiseProduct(
                        alpertWlt1[ly] * alpertWlt1[lx].transpose()));
                    const size_t i = ((1 << jx) + kx)*(wltOrder+1)+lx;
                    const size_t j = ((1 << jy) + ky)*(wltOrder+1)+ly;
                    kernelMatrix(i, j) = prod.sum();
                  }
                }
              }
            }
          }
        }

        return kernelMatrix;
      }

    template<unsigned int wltOrder>
    class SVD {

      public:
        enum : unsigned int { dim = 2 };
        using Direction = FieldVector<double, dim>;

        // This is the number of quadrature points per interval.
        const unsigned int numSperInterval = wltOrder+1;

        SVD() = delete;
        SVD(const SVD&) = delete;

        template<class Function>
        SVD(const Function& kernelFunction,
            double accuracyKernel)
          : maxLevel(setLevel(accuracyKernel,100)),
            num_s((wltOrder+1) << maxLevel),
            kernelSVD(num_s, num_s, Eigen::ComputeThinU | Eigen::ComputeThinV),
            level(maxLevel),
            rank(num_s),
            nQuadAngle(30) {
          // if((1u << maxLevel) != num_s)
          //   DUNE_THROW(MathError, "You are using " << num_s
          //       << " directions, but only powers of 2 are supported.");
          std::cout << "Accuracy kernel " << accuracyKernel
            << " with wlt order " << wltOrder
            << " requires level J = " << maxLevel
            << " and " << num_s << " directions." << std::endl;

          Eigen::MatrixXd
            kernelMatrix(waveletKernelMatrix(kernelFunction,
              wltOrder, maxLevel, nQuadAngle));
          /* initialize SVD of kernel (using Eigen) */
          kernelSVD.compute(kernelMatrix);

        }

        Eigen::VectorXd
        PairToXd(const std::pair<Eigen::VectorXd,
                std::vector<Eigen::VectorXd>>& wPair) const {

          const Eigen::VectorXd w0 = wPair.first;
          const std::vector<Eigen::VectorXd> w1 = wPair.second;

          const size_t L = w0.size();
          const size_t J = w1.size();

          Eigen::VectorXd w(L << J);
          w.segment(0, L) = w0;

          size_t pos0 = L;
          size_t pos1 = L;
          for(size_t j=0; j < J; j++) {
            pos0 = pos1;
            pos1 = pos0+w1[j].size();
            w.segment(pos0, w1[j].size()) = w1[j];
          }
          return w;
        }

        std::pair<Eigen::VectorXd,std::vector<Eigen::VectorXd>>
        XdToPair(Eigen::VectorXd& v) const {
          if(level==0) {
            return std::make_pair(v,std::vector<Eigen::VectorXd>());
          }
          else {
            Eigen::VectorXd sf = v.segment(0,wltOrder+1);
            std::vector<Eigen::VectorXd> wlt(level);

            size_t pos0 = wltOrder+1;
            size_t pos1 = pos0;

            for(size_t j=0; j<level; j++) {
              pos0 = pos1;
              pos1 = pos0 + ((wltOrder+1)<<j);
              wlt[j] = v.segment(pos0, (wltOrder+1)<<j);
            }
            return std::make_pair(sf,wlt);
          }
        }

        // Given a vector of u(s_i), compute (Ku)(s_i) with SVD
        void applyToVector(Eigen::VectorXd& u) const {

          using namespace boost::math::constants;
          // Express u in Legendre basis on each segment of [-pi,pi]
          const size_t uLevel = std::ilogb(u.size()/(wltOrder+1));
          Eigen::VectorXd uLegendre
            = lagrange_to_legendre_VJ(u, uLevel, wltOrder);
          // Express u in wlt basis using DWT
          const size_t quadOrder = 2*nQuadAngle-1;
          const auto uWltPair = DWT(uLegendre, wltOrder+1, uLevel, quadOrder);
          Eigen::VectorXd uWlt = PairToXd(uWltPair);
          // Compute Ku with SVD up to level given by rank.
          // The result is expressed in the wlt basis
          Eigen::VectorXd KuWlt
            = kernelSVD.matrixU().topLeftCorner(rows, rank)
            * kernelSVD.singularValues().head(rank).asDiagonal()
            * kernelSVD.matrixV().topLeftCorner(uWlt.size(), rank).adjoint()
            * uWlt;
          // Express Ku in Legendre basis on each segment of [-pi,pi]
          // via an inverse wlt transform
          std::pair<Eigen::VectorXd,std::vector<Eigen::VectorXd>>
            KuWltPair = XdToPair(KuWlt);
          Eigen::VectorXd KuLegendre
            = IDWT(KuWltPair, wltOrder+1, level, quadOrder);
          // Express Ku in Lagrange basis basis on each segment of [-pi,pi],
          // i.e., in a vector of (Ku)(s_i)
          u = legendre_to_lagrange_VJ(KuLegendre, level, wltOrder);
        }

        std::vector<Direction> setAccuracy(double accuracy) {
          using namespace Eigen;
          VectorXd singularValues = kernelSVD.singularValues();
          size_t i = singularValues.size() - 1;
          double err = 0,
                rank_err = singularValues(i) * singularValues(i);
          while (err + rank_err < (accuracy * accuracy / 4.) && i > 0) {
            err += rank_err;
            i -= 1;
            rank_err = singularValues(i) * singularValues(i);
          }
          rank = i+1;
          // TODO: If accuracy is low enough to allow rank = 0,
          //       this gives rank = 1.

          level = setLevel(accuracy/2.,maxLevel);

          std::vector<Direction> sVector((wltOrder+1) << level);
          compute_sVector(sVector);
          rows = sVector.size();
          return sVector;
        }

        template<class Direction>
        void compute_sVector(std::vector< Direction >& sVector) {

          // Get Gauss-Legendre quadrature in [0,1]
          const int quadOrder = 2*wltOrder+1;
          GeometryType type;
          type.makeLine();
          const size_t dim = 1;
          const Dune::QuadratureRule<double, dim>& quad =
            Dune::QuadratureRules<double, dim>::rule(type,
                                                    quadOrder,
                                                    QuadratureType::GaussLegendre);
          using namespace boost::math::constants;
          size_t i=0;
          for(int k = 0; k < (1<<level); ++k)
          {
            for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
              double angle = 2*pi<double>()*(k+quad[pt].position())
                             / (1<<level) - pi<double>(); // Angle in [-pi,pi]
              sVector[i] = {cos(angle),sin(angle)};
              i++;
            }
          }
        }

        size_t getLevel() {
          return level;
        }

        size_t getMaxLevel() {
          return maxLevel;
        }

        size_t getNumS() {
          return ((wltOrder+1) << level);
        }

        size_t maxNumS() const {
          return ((wltOrder+1) << maxLevel);
        }

        size_t getWltOrder() const {
          return wltOrder;
        }

        std::string typeApprox() const {
          return "Kernel approximation with: SVD";
        }

        std::string info() const {
          return "SVD rank: "+std::to_string(rank);
        }

        std::vector<double> getSingularValues() const {
          std::vector<double> vec(kernelSVD.singularValues().data(),
            kernelSVD.singularValues().data()
            +kernelSVD.singularValues().size());
          return vec;
        }

      private:

        static inline size_t setLevel(double accuracy, size_t maxLevel) {
          size_t level=0;

          for(; level < maxLevel; ++level) {
            if(1./((1.+level*level)*( 1 << ((wltOrder+1)*level))) <= accuracy)
              break;
          }

          return level;
        }

        size_t maxLevel;
        size_t num_s;
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> kernelSVD;
        size_t level;
        size_t rows;
        size_t rank;
        size_t nQuadAngle;
      };

    template<unsigned int wltOrder>
    class MatrixTH {

      public:
        enum : unsigned int { dim = 2 };
        using Direction = FieldVector<double, dim>;

        // This is the number of quadrature points per interval.
        const unsigned int numSperInterval = (wltOrder+1);


        MatrixTH() = delete;
        MatrixTH(const MatrixTH&) = delete;

        template<class Function>
        MatrixTH(const Function& kernelFunction,
            double accuracyKernel)
          : maxLevel(setLevel(accuracyKernel,100)),
            num_s((wltOrder+1) << maxLevel),
            level(maxLevel),
            rank(num_s),
            nQuadAngle(30),
            kernelMatrix(
              waveletKernelMatrix(kernelFunction,
              wltOrder, maxLevel, nQuadAngle)),
            kernelMatrixTH(kernelMatrix),
            ofsTH("output_rad_trans_TH")
            {
          // if((1u << maxLevel) != num_s)
          //   DUNE_THROW(MathError, "You are using " << num_s
          //       << " directions, but only powers of 2 are supported.");
          std::cout << "Accuracy kernel " << accuracyKernel
            << " with wlt order " << wltOrder
            << " requires level J = " << maxLevel
            << " and " << num_s << " directions." << std::endl;

          ofsTH << "Accuracy kernel " << accuracyKernel
                << " with wlt order " << wltOrder
                << " requires level J = " << maxLevel
                << " and " << num_s << " directions."
                << std::endl << std::endl;

          ofsTH << "kernelMatrix has size: "
                << kernelMatrix.rows() << " x " << kernelMatrix.cols()
                << " . It has " << kernelMatrix.size() << " elements."
                << std::endl << std::endl;
        }

        Eigen::VectorXd
        PairToXd(const std::pair<Eigen::VectorXd,
                std::vector<Eigen::VectorXd>>& wPair) const {

          const Eigen::VectorXd w0 = wPair.first;
          const std::vector<Eigen::VectorXd> w1 = wPair.second;

          const size_t L = w0.size();
          const size_t J = w1.size();

          Eigen::VectorXd w(L << J);
          w.segment(0, L) = w0;

          size_t pos0 = L;
          size_t pos1 = L;
          for(size_t j=0; j < J; j++) {
            pos0 = pos1;
            pos1 = pos0 + w1[j].size();
            w.segment(pos0, w1[j].size()) = w1[j];
          }
          return w;
        }

        std::pair<Eigen::VectorXd,std::vector<Eigen::VectorXd>>
        XdToPair(Eigen::VectorXd& v) const {
          if(level==0) {
            return std::make_pair(v,std::vector<Eigen::VectorXd>());
          }
          else {
            Eigen::VectorXd sf = v.segment(0,wltOrder+1);
            std::vector<Eigen::VectorXd> wlt(level);

            size_t pos0 = wltOrder+1;
            size_t pos1 = pos0;

            for(size_t j=0; j<level; j++) {
              pos0 = pos1;
              pos1 = pos0 + ((wltOrder+1)<<j);
              wlt[j] = v.segment(pos0, (wltOrder+1)<<j);
            }
            return std::make_pair(sf,wlt);
          }
        }

        void applyToVector(Eigen::VectorXd& u) const {
          using namespace boost::math::constants;
          // Express u in Legendre basis on each segment of [-pi,pi]
          const size_t uLevel = std::ilogb(u.size()/(wltOrder+1));
          Eigen::VectorXd uLegendre
            = lagrange_to_legendre_VJ(u, uLevel, wltOrder);
          // Express u in wlt basis using DWT
          const size_t quadOrder = 2*nQuadAngle-1;
          const auto uWltPair = DWT(uLegendre, wltOrder+1, uLevel, quadOrder);
          Eigen::VectorXd uWlt = PairToXd(uWltPair);
          // Approx with trucated kernelMatrix
          Eigen::VectorXd KuWlt
            = kernelMatrixTH.topLeftCorner(rows,uWlt.size()) * uWlt;
          // Express Ku in Legendre basis on each segment of [-pi,pi]
          // via an inverse wlt transform
          std::pair<Eigen::VectorXd,std::vector<Eigen::VectorXd>>
            KuWltPair = XdToPair(KuWlt);
          Eigen::VectorXd KuLegendre
            = IDWT(KuWltPair, wltOrder+1, level, quadOrder);
          // Express Ku in Lagrange basis basis on each segment of [-pi,pi],
          // i.e., in a vector of (Ku)(s_i)
          u = legendre_to_lagrange_VJ(KuLegendre, level, wltOrder);
        }

        std::vector<Direction> setAccuracy(double accuracy) {
          // Compute the new wavelet level
          level = setLevel(accuracy/2.,maxLevel);

          // Directions
          std::vector<Direction> sVector((wltOrder+1) << level);
          compute_sVector(sVector);

          // Compute thresholded kernel matrix
          // using matrix compression techniques
          set_kernelMatrixTH(level);
          rows = sVector.size();  // number of rows for function applyToVector

          // Return vector of directions
          return sVector;
        }

        template<class Direction>
        void compute_sVector(std::vector< Direction >& sVector) {

          // Get Gauss-Legendre quadrature in [0,1]
          const int quadOrder = 2*wltOrder+1;
          GeometryType type;
          type.makeLine();
          const size_t dim = 1;
          const Dune::QuadratureRule<double, dim>& quad =
            Dune::QuadratureRules<double, dim>::rule(type,
                                                    quadOrder,
                                                    QuadratureType::GaussLegendre);
          using namespace boost::math::constants;
          size_t i=0;
          for(int k = 0; k < (1<<level); ++k)
          {
            for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
              double angle = 2*pi<double>()*(k+quad[pt].position())
                             / (1<<level) - pi<double>(); // Angle in [-pi,pi]
              sVector[i] = {cos(angle),sin(angle)};
              i++;
            }
          }
        }

        size_t getLevel() {
          return level;
        }

        size_t getMaxLevel() {
          return maxLevel;
        }

        size_t getNumS() const {
          return ((wltOrder+1) << level);
        }

        size_t maxNumS() const {
          return ((wltOrder+1) << maxLevel);
        }

        size_t getWltOrder() const {
          return wltOrder;
        }

        std::string typeApprox() const {
          return "Kernel approximation with: Matrix compression";
        }

        std::string info() const {
          Eigen::MatrixXd zerosFilter = kernelMatrixTH.unaryExpr( [] (double x)
            { return abs(x) <= 1.e-18 ? 1. : 0.; } );

          std::string s = "Kernel matrix is of size "
                + std::to_string(kernelMatrixTH.rows())
                +"x"
                +std::to_string(kernelMatrixTH.cols())
                +". It has "+std::to_string(kernelMatrixTH.size())
                + " elements of which " + std::to_string((int)zerosFilter.sum())
                +" are zero.";

          return s;
        }

      private:

        static inline size_t setLevel(double accuracy, size_t maxLevel) {
          size_t level=0;

          for(; level < maxLevel; ++level) {
            if(1./((1.+level*level)*( 1 << ((wltOrder+1)*level))) <= accuracy)
              break;
          }

          return level;
        }

        Eigen::MatrixXd computeArcDistanceMatrix (
          const size_t targetLevel,
          std::pair<Eigen::MatrixXd,Eigen::MatrixXd> levelsMatrix) const {

          using namespace boost::math::constants;

          // Matrix with 2^{min(j,j')}
          Eigen::MatrixXd twoPowMinLevel = (levelsMatrix.first).cwiseMin(levelsMatrix.second);
          twoPowMinLevel = twoPowMinLevel.unaryExpr( [] (double x)
            {return std::pow(2.,x);} );

          // Matrix with dist(S_\lambda,S,{\lambda'})
          const double r = pi<double>();
          const size_t num_s = (wltOrder+1) << targetLevel;
          Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(num_s, num_s);
          double d;
          for(size_t jy = 0; jy < targetLevel; jy++) {
            for(size_t ky = 0, ky_max = 1 << jy; ky < ky_max; ky++) {

              double ymin = r *   ky   * std::exp2(-(double)jy+1) - r;
              double ymax = r * (ky+1) * std::exp2(-(double)jy+1) - r;

              for(size_t jx = 0; jx < targetLevel; jx++) {
                for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {

                  double xmin = r *   kx   * std::exp2(-(double)jx+1) - r;
                  double xmax = r * (kx+1) * std::exp2(-(double)jx+1) - r;

                  if(ymax < xmin)
                    d = (xmin-ymax < 2*r+ymin-xmax) ? (xmin-ymax) : (2*r+ymin-xmax);
                  else{
                    if(xmax < ymin)
                      d = (ymin-xmax < 2*r+xmin-ymax) ? (ymin-xmax) : (2*r+xmin-ymax);
                    else d=0.; //( (xmax <= ymax and xmax>=ymin) || (xmin<=ymax and xmin>=ymin) )
                  }

                  for(size_t ly = 0; ly < wltOrder+1; ly++) {
                    for(size_t lx = 0; lx < wltOrder+1; lx++) {

                      const size_t i = ((1 << jx) + kx)*(wltOrder+1)+lx;
                      const size_t j = ((1 << jy) + ky)*(wltOrder+1)+ly;
                      dist(i,j) = d;
                    }
                  }
                }
              }
            }
          }
          return twoPowMinLevel.cwiseProduct(dist);
        }

        static std::pair<Eigen::MatrixXd,Eigen::MatrixXd>
        computeLevelsMatrix(size_t targetLevel) {

          Eigen::VectorXd ones
            = Eigen::VectorXd::Ones((wltOrder+1) << targetLevel);
          Eigen::VectorXd levelsVector = ones;

          levelsVector.segment(0,wltOrder+1)
          = Eigen::VectorXd::Zero(wltOrder+1);

          size_t pos0 = wltOrder+1;
          size_t pos1 = pos0;

          for(size_t j=0; j < targetLevel; j++) {
            pos0 = pos1;
            pos1 = pos0 + ((wltOrder+1) << j);
            levelsVector.segment(pos0,(wltOrder+1)<<j) *= j;
          }

          return std::make_pair(levelsVector*ones.transpose(),ones*levelsVector.transpose());
        }


        void truncationInScale(
          std::pair<Eigen::MatrixXd,Eigen::MatrixXd> levelsMatrix,
          size_t targetLevel)
        {
          Eigen::MatrixXd filter = (levelsMatrix.first - levelsMatrix.second).cwiseAbs();

          filter = filter.unaryExpr( [targetLevel] (double x)
            { return x <= targetLevel ? 1. : 0.; } );

          ofsTH << "Filter truncation in Scale:" << std::endl
                << filter << std::endl << std::endl;

          kernelMatrixTH = kernelMatrix.cwiseProduct(filter);
        }

        void truncationInSpace(
          std::pair<Eigen::MatrixXd,Eigen::MatrixXd> levelsMatrix,
          size_t targetLevel)
        {
          Eigen::MatrixXd lhsBound
          = computeArcDistanceMatrix(maxLevel, levelsMatrix);

          Eigen::MatrixXd rhsBound = (levelsMatrix.first - levelsMatrix.second).cwiseAbs();

          rhsBound = rhsBound.unaryExpr( [targetLevel] (double x)
            { return std::pow(2.,targetLevel-x)/(1.+x*x); });

          Eigen::MatrixXd filter = (lhsBound.array() <= rhsBound.array()).cast<double>();

          ofsTH << "Filter truncation in Space:" << std::endl
                << filter << std::endl << std::endl;

          kernelMatrixTH = kernelMatrixTH.cwiseProduct(filter);
        }

        void set_kernelMatrixTH(size_t targetLevel) {

          kernelMatrixTH = kernelMatrix;

          std::pair<Eigen::MatrixXd,Eigen::MatrixXd> levelsMatrix
          = computeLevelsMatrix(maxLevel);

          truncationInScale(levelsMatrix,targetLevel);
          truncationInSpace(levelsMatrix,targetLevel);
        }

        size_t maxLevel;
        size_t num_s;
        size_t level;
        size_t rank;
        size_t nQuadAngle;
        size_t rows;
        Eigen::MatrixXd kernelMatrix;
        Eigen::MatrixXd kernelMatrixTH;
        std::ofstream ofsTH;
      };
  } //End namespace AlpertWavelet


  namespace HaarWavelet {
      // Scaling function $\phi_{0,0}$ on the interval [-r,r)
      Eigen::VectorXd sf(const Eigen::VectorXd& x, double r) {
        double normfactor = 1./sqrt(2*r);
        Eigen::VectorXd res(x.size());
        for(size_t i=0, imax=x.size(); i<imax; i++) {
          double val = x(i);
          res(i) = normfactor * (val >= -r and val < r);
        }
        return res;
      }

      // Wavelet $\psi_{n,k}$ on the interval [-r,r)
      // Its support is in [xmin,xmax] where
      //      xmin = r*k*pow(2,-j+1)-r
      //      xmax = r*(k+1)*pow(2,-j+1)-r
      Eigen::VectorXd wlt(size_t j, size_t k,
                          const Eigen::VectorXd& x, double r) {
        double xmin = r *   k   * std::exp2(-(double)j+1) - r;
        double xmax = r * (k+1) * std::exp2(-(double)j+1) - r;
        double xmiddle = (xmin + xmax)/2;
        double normfactor = 1./sqrt(xmax-xmin);

        Eigen::VectorXd res(x.size());
        for(size_t i=0, imax=x.size(); i<imax; i++) {
          double val = x(i);
          double positivePart =     (val >= xmin and val < xmiddle);
          double negativePart = -1.*(val >= xmiddle and val <= xmax);
          res(i) = normfactor * (positivePart + negativePart);
        }
        return res;
      }

      template<class KernelFunction>
      double evalKernel(
          KernelFunction& kernelFunction,
          const Eigen::VectorXd& x, const Eigen::VectorXd& xValues,
          double quadweightx,
          const Eigen::VectorXd& y, const Eigen::VectorXd& yValues,
          double quadweighty
          ) {
        enum : unsigned int { dim = 2 };

        // evaluate the kernel function
        const size_t nquadx = x.size();
        const size_t nquady = y.size();
        const double quadweight = quadweighty * quadweightx;
        double eval = 0.;
        for(size_t i = 0; i < nquadx; i++) {
          const double s_i = i*nquadx*2*boost::math::constants::pi<double>();
          for(size_t j = 0; j < nquady; j++) {
            const double s_j = j*nquady*2*boost::math::constants::pi<double>();
            // Integral over [-pi,pi]x[-pi,pi]
            eval += kernelFunction(s_i - s_j)
                  * quadweight * xValues[i] * yValues[j];
          }
        }
        return eval;
      }

      std::tuple<Eigen::VectorXd, double> computeQuadPoints(
          size_t jx, size_t kx, double r, size_t maxLevel) {
        double xmin = r *   kx   * std::exp2(-(double)jx+1) - r;
        double xmax = r * (kx+1) * std::exp2(-(double)jx+1) - r;
        size_t nquad = 1 << (maxLevel - jx);
        Eigen::VectorXd x(nquad);
        for(size_t i=0; i < nquad; i++) {
          x[i] = xmin + (2*i+1)/(2.*nquad) * (xmax - xmin);
        }
        double quadweight = (xmax - xmin) / nquad;
        return std::make_tuple(x, quadweight);
      }

      // Discrete Haar wavelet transform
      // data is used for in- and output and is assumed to be of size 2^n
      // for a natural number n.
      void DWT(Eigen::VectorXd& data) {
        Eigen::VectorXd tmp(data.size());
        for(size_t len = data.size() >> 1; len > 0; len >>=1) {
          for(size_t i = 0; i < len; i++) {
            double sum  = data[2*i] + data[2*i+1];
            double diff = data[2*i] - data[2*i+1];
            tmp[i]     = .5 * sum;
            tmp[len+i] = .5 * diff;
          }
          data.head(2*len) = tmp.head(2*len);
        }
        double scaling_factor = sqrt(2*boost::math::constants::pi<double>());
        data.segment(0, 1) *= scaling_factor;
        for(size_t len=1, max_len=data.size()>>1; len < max_len; len <<=1) {
          scaling_factor /= sqrt(2);
          data.segment(len, len) *= scaling_factor;
        }
      }

      // Inverse discrete Haar wavelet transform
      // data is used for in- and output and is assumed to be of size 2^n
      // for a natural number n.
      void IDWT(Eigen::VectorXd& data) {
        Eigen::VectorXd tmp(data.size());
        double scaling_factor = 1./sqrt(2*boost::math::constants::pi<double>());
        data.segment(0, 1) *= scaling_factor;
        for(size_t len=1, max_len=data.size()>>1; len < max_len; len <<=1) {
          scaling_factor *= sqrt(2);
          data.segment(len, len) *= scaling_factor;
        }
        size_t max_len = data.size();
        for(size_t len = 1; len < max_len; len <<=1) {
          for(size_t i = 0; i < len; i++) {
            double sum  = data[i] + data[len+i];
            double diff = data[i] - data[len+i];
            tmp[2*i]   = sum;
            tmp[2*i+1] = diff;
          }
          data.head(2*len) = tmp.head(2*len);
        }
      }

      template<class Function>
      Eigen::MatrixXd waveletKernelMatrix(const Function& kernelFunction,
                                          size_t maxLevel)
      {
        using namespace Eigen;
        using namespace boost::math::constants;

        const size_t num_s = 1ul << maxLevel;
        MatrixXd kernelMatrix(num_s, num_s);

        /* compute wavelet transform of the kernel matrix */
        const double r = pi<double>();

        // jy, ky = 0
        {
          VectorXd y; double qwy;
          std::tie(y, qwy) = computeQuadPoints(0, 0, r, maxLevel);
          VectorXd sfy = sf(y, r);
          // jx, kx = 0
          {
            VectorXd x; double qwx;
            std::tie(x, qwx) = computeQuadPoints(0, 0, r, maxLevel);
            VectorXd sfx = sf(x, r);
            kernelMatrix(0, 0)
                = evalKernel(kernelFunction, x, sfx, qwx, y, sfy, qwy);
          }
          for(size_t jx = 0; jx < maxLevel; jx++) {
            for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {
              VectorXd x; double qwx;
              std::tie(x, qwx) = computeQuadPoints(jx, kx, r, maxLevel);
              VectorXd wltx = wlt(jx, kx, x, r);
              const size_t i = (1 << jx) + kx;
              kernelMatrix(i, 0)
                  = evalKernel(kernelFunction, x, wltx, qwx, y, sfy, qwy);
            }
          }
        }
        for(size_t jy = 0; jy < maxLevel; jy++) {
          for(size_t ky = 0, ky_max = 1 << jy; ky < ky_max; ky++) {
            VectorXd y; double qwy;
            std::tie(y, qwy) = computeQuadPoints(jy, ky, r, maxLevel);
            VectorXd wlty = wlt(jy, ky, y, r);
            // jx, kx = 0
            {
              VectorXd x; double qwx;
              std::tie(x, qwx) = computeQuadPoints(0, 0, r, maxLevel);
              VectorXd sfx = sf(x, r);
              const size_t j = (1 << jy) + ky;
              kernelMatrix(0, j)
                  = evalKernel(kernelFunction, x, sfx, qwx, y, wlty, qwy);
            }
            for(size_t jx = 0; jx < maxLevel; jx++) {
              for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {
                VectorXd x; double qwx;
                std::tie(x, qwx) = computeQuadPoints(jx, kx, r, maxLevel);
                VectorXd wltx = wlt(jx, kx, x, r);
                const size_t i = (1 << jx) + kx;
                const size_t j = (1 << jy) + ky;
                kernelMatrix(i, j)
                    = evalKernel(kernelFunction, x, wltx, qwx, y, wlty, qwy);
              }
            }
          }
        }
        return kernelMatrix;
      }

      class MatrixCompression {
      public:
        enum : unsigned int { dim = 2 };
        using Direction = FieldVector<double, dim>;

        enum : unsigned int { numSperInterval = 1 };

        MatrixCompression() = delete;
        MatrixCompression(const MatrixCompression&) = delete;

        template<class Function>
        MatrixCompression(const Function& kernelFunction, double maxAccuracy)
          : maxLevel(setLevel(maxAccuracy, 1000)),
            level(maxLevel),
            rows(1u << level),
            kernelMatrix(waveletKernelMatrix(kernelFunction, maxLevel))
          {}

        void applyToVector(Eigen::VectorXd& v) const {
          DWT(v);
          v = kernelMatrix.topLeftCorner(rows, v.size()) * v;
          IDWT(v);
        }

        std::vector<Direction> setAccuracy(double accuracy) {
          level = setLevel(accuracy, maxLevel);
          rows = 1u << level;

          // compute transport directions corresponding to quadrature points
          std::vector<Direction> sVector(rows);
          for(unsigned int i = 0; i < rows; ++i)
          {
            using namespace boost::math::constants;
            sVector[i] = {cos(2*pi<double>()*i/rows),
                          sin(2*pi<double>()*i/rows)};
          }

          return sVector;
        }

        size_t maxNumS() const {
          return 1u << maxLevel;
        }

        std::string info() const {
          std::string s = "MatrixCompression approximation with level "
                          + std::to_string(level);
          return s;
        }

      private:
        static inline size_t setLevel(double accuracy, size_t maxLevel) {
          size_t level = 1;
          for(; level < maxLevel; ++level) {
            if(1./((1 << level)*(1+level*level)) < accuracy/4.)
              break;
          }
          return level;
        }

        size_t maxLevel;
        size_t level;
        size_t rows;
        Eigen::MatrixXd kernelMatrix;
      };

      class SVD {
      public:
        enum : unsigned int { dim = 2 };
        using Direction = FieldVector<double, dim>;

        enum : unsigned int { numSperInterval = 1 };

        SVD() = delete;
        SVD(const SVD&) = delete;

        template<class Function>
        SVD(const Function& kernelFunction, double maxAccuracy)
          : maxLevel(setLevel(maxAccuracy, 1000u)),
            level(maxLevel),
            kernelSVD(maxLevel, maxLevel,
                      Eigen::ComputeThinU | Eigen::ComputeThinV),
            rows(1u << level),
            rank(rows) {
          Eigen::MatrixXd kernelMatrix(waveletKernelMatrix(kernelFunction,
                                                           maxLevel));
          /* initialize SVD of kernel (using Eigen) */
          kernelSVD.compute(kernelMatrix);
        }

        void applyToVector(Eigen::VectorXd& v) const {
          DWT(v);
          v = kernelSVD.matrixU().topLeftCorner(rows, rank)
            * kernelSVD.singularValues().head(rank).asDiagonal()
            * kernelSVD.matrixV().topLeftCorner(v.size(), rank).adjoint() * v;
          IDWT(v);
        }

        std::vector<Direction> setAccuracy(double accuracy) {
          using namespace Eigen;
          VectorXd singularValues = kernelSVD.singularValues();
          size_t i = singularValues.size() - 1;
          double err = 0,
                rank_err = singularValues(i) * singularValues(i);
          const double accuracy_squared = accuracy * accuracy / 4.;
          while (err + rank_err < accuracy_squared && i > 0) {
            err += rank_err;
            i -= 1;
            rank_err = singularValues(i) * singularValues(i);
          }
          rank = i+1;
          // TODO: If accuracy is low enough to allow rank = 0,
          //       this gives rank = 1.

          // set level according to given accuracy
          level = setLevel(accuracy, maxLevel);
          rows = 1u << level;

          // compute transport directions corresponding to quadrature points
          std::vector<Direction> sVector(rows);
          for(unsigned int i = 0; i < rows; ++i)
          {
            using namespace boost::math::constants;
            sVector[i] = {cos(2*pi<double>()*i/rows),
                          sin(2*pi<double>()*i/rows)};
          }

          return sVector;
        }

        size_t maxNumS() const {
          return 1u << maxLevel;
        }

        std::string info() const {
          std::string s = "Wavelet SVD approximation with rank "
                          + std::to_string(rank)
                          + " and level "
                          + std::to_string(level);
          return s;
        }

      private:
        static inline size_t setLevel(double accuracy, size_t maxLevel) {
          size_t level = 1;
          for(; level < maxLevel; ++level){
            if(1./((1 << level)*(1+level*level)) < accuracy/4.)
              break;
          }
          return level;
        }

        size_t maxLevel;
        size_t level;
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> kernelSVD;
        size_t rows;
        size_t rank;
      };
  }
}

}

#endif // defined DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH
