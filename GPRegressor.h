#pragma once
#include "GPKernel.h"

#include "TVector.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/IFunction.h"

class GPRegressor
{
public:
  GPRegressor(GPKernel& kern, bool normalize_y=true, TVectorD* noise = 0);

  void Fit(std::vector<GPPoint> X, TVectorD y);
  void Predict(std::vector<GPPoint> X);

  void SolveHyperParameters();


  TVectorD PosteriorMean(){ return fPosteriorMean; }
  TVectorD PosteriorStd(){ return fPosteriorStd; }

protected:

  GPKernel& fKern;
  bool doNorm;
  TVectorD* fNoise;
  
  std::vector<GPPoint> fX_T;
  mutable TVectorD fY_T;
  double fY_Tm, fY_Ts;
  
  mutable TVectorD fPosteriorMean, fPosteriorStd;
  mutable TMatrixDSym fK, fKInv;
  mutable TVectorD fAlpha;

};

class MarginalLikelihood: public ROOT::Math::IGradientFunctionMultiDim
{
  public:
    MarginalLikelihood(GPKernel& kern, std::vector<GPPoint> x, TVectorD y, TVectorD* noise = 0);

    void Solve(const double* theta) const;

    double DoEval(const double* par) const;
    double DoDerivative(const double* par, unsigned int ipar) const;
    unsigned int NDim() const { return fDim; }
    
    ROOT::Math::IGradientFunctionMultiDim* Clone() const { return new MarginalLikelihood(fKern, fX, fY, fNoise); }

  protected:
 
    GPKernel& fKern; 
    int fDim; 
    const std::vector<GPPoint> fX;
    const TVectorD fY; 
    TVectorD* fNoise;
    
    mutable TMatrixDSym fK, fKInv;
    mutable TVectorD fAlpha;
    mutable double fKDet;
};
