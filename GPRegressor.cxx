#include "GPRegressor.h"

#include "Fit/Fitter.h"
#include "TDecompChol.h"

#include <math.h>

//----------------------------------------------------------------------
GPRegressor::GPRegressor(GPKernel& kern, bool normalize_y, TVectorD* noise)
  : fKern(kern), doNorm(normalize_y), fNoise(noise)
{
}
//----------------------------------------------------------------------
void GPRegressor::Fit(std::vector<GPPoint> X, TVectorD y)
{
  fY_T.ResizeTo(y.GetNoElements());
  fY_T = y;
  fX_T = X;
  if(doNorm){
    int n = fY_T.GetNoElements(); 
    fY_Tm = fY_T.Sum()/n;
    fY_T -= fY_Tm;
    
    fY_Ts = sqrt(fY_T.Norm2Sqr()/(n-1));
    fY_T *= 1./fY_Ts;
  }
  else{
    fY_Tm = 0.;
    fY_Ts = 1.;
  }

  if(fNoise)
    assert((*fNoise).GetNoElements() == y.GetNoElements());

  SolveHyperParameters();

  fK.ResizeTo(fX_T.size(), fX_T.size());
  fK = fKern(fX_T);
  if(fNoise){
    for(int i = 0; i < fK.GetNrows(); i++)
      fK(i, i) += (*fNoise)[i];
  }

  fKInv.ResizeTo(fK.GetNrows(), fK.GetNrows());
  fAlpha.ResizeTo(fX_T.size());

  TDecompChol cho(fK);
  if(!cho.Decompose()){
    std::cerr << "Cholesky decomp failed" << std::endl;
  }
  else{
    Bool_t ok;
    fAlpha = cho.Solve(fY_T, ok);
    cho.Invert(fKInv);
    if(!ok)
      std::cerr << "something went wrong" << std::endl;
  }
}
//----------------------------------------------------------------------
void GPRegressor::SolveHyperParameters()
{
  MarginalLikelihood lml(fKern, fX_T, fY_T, fNoise);
  ROOT::Fit::Fitter fitter;
  
  int n = fKern.NPar();
  double* par = new double[n];
  std::vector<double> par_start = fKern.GetParameters();
  for(int i = 0; i < n; i++)
    par[i] = log(par_start[i]);


  fitter.SetFCN(lml, par);
  for(int i = 0; i < n; i++)
    fitter.Config().ParSettings(i).SetLimits(-5, 5);
  
  fitter.Config().SetMinimizer("GSLMultiMin");
  fitter.Config().MinimizerOptions().SetMinimizerAlgorithm("BFGS2");
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(10000000);
  fitter.Config().MinimizerOptions().SetMaxIterations(50);
  fitter.Config().MinimizerOptions().SetPrintLevel(2);
  fitter.Config().MinimizerOptions().Print();
  fitter.FitFCN();

  fKern.SetThetas(fitter.Result().Parameters());
  return;
}
//----------------------------------------------------------------------

void GPRegressor::Predict(std::vector<GPPoint> X)
{
  TMatrixD PosteriorK_T = fKern.KernelBlock(X, fX_T);
  fPosteriorMean.ResizeTo(X.size());
  fPosteriorStd.ResizeTo(X.size());
  fPosteriorMean = PosteriorK_T*fAlpha;
  fPosteriorMean *= fY_Ts;
  fPosteriorMean += fY_Tm;
 
  fPosteriorStd = fKern.KernelDiag(X);

  TMatrixD temp(X.size(), fX_T.size());
  temp.Mult(PosteriorK_T, fKInv);
  for(int i = 0; i < temp.GetNrows(); i++){
    fPosteriorStd[i] -= TVectorD(temp[i])*TVectorD(PosteriorK_T[i]);
    if(fPosteriorStd[i] < 0.)
      fPosteriorStd[i] = 0.;
  }
  fPosteriorStd *= TMath::Power(fY_Ts, 2);
  fPosteriorStd = fPosteriorStd.Sqrt();
}
//----------------------------------------------------------------------


MarginalLikelihood::MarginalLikelihood(GPKernel& kern, std::vector<GPPoint> x, TVectorD y, TVectorD* noise)
  : fKern(kern), fX(x), fY(y), fNoise(noise)
{
    fDim = fKern.NPar();
}
//----------------------------------------------------------------------
void MarginalLikelihood::Solve(const double* theta) const
{
  std::vector<double> thetas;
  for(int i = 0; i < fDim; i++)
    thetas.push_back(theta[i]);
  fKern.SetThetas(thetas);
 
  fK.ResizeTo(fX.size(), fX.size()); 
  fKInv.ResizeTo(fX.size(), fX.size()); 
  fAlpha.ResizeTo(fX.size());

  fK = fKern(fX, -1);
  if(fNoise){
    assert((*fNoise).GetNoElements() == fY.GetNoElements());
    for(int i = 0; i < fK.GetNrows(); i++)
      fK(i, i) += (*fNoise)[i];
  }

  TDecompChol cho(fK);
  if(!cho.Decompose()){
    std::cerr << "Cholesky decomp failed" << std::endl;
  }
  else{
    Bool_t ok;
    fAlpha = cho.Solve(fY, ok);
    cho.Invert(fKInv);
    if(!ok)
      std::cerr << "something went wrong" << std::endl;
    TMatrixD choU = cho.GetU();
    fKDet = 0.;
    for(int i = 0; i < choU.GetNrows(); i++)
      fKDet += 2.*log(choU[i][i]);
  }
  
  return;

}
//----------------------------------------------------------------------
double MarginalLikelihood::DoEval(const double* par) const
{
  Solve(par);
  double ret = 0.;
  for(int i = 0; i < fY.GetNoElements(); i++){
    ret += 0.5*fY[i]*fAlpha[i];
  }
  ret += 0.5*fKDet;
  ret += 0.5*fY.GetNoElements()*log(2.*M_PI);

  return ret;
}
//----------------------------------------------------------------------
double MarginalLikelihood::DoDerivative(const double* par, unsigned int ipar) const
{
  TMatrixDSym dK = fKern(fX, ipar);
  double ret = 0.;
  for(int i = 0; i < (int) fX.size(); i++){
    for(int j = 0; j < (int) fX.size(); j++){
        ret -= 0.5*((fAlpha[i]*fAlpha[j]) - fKInv[i][j])*dK[j][i];
    }
  }
  return ret;
}
//----------------------------------------------------------------------
