#include "GPKernel.h"

#include <math.h>

//----------------------------------------------------------------------
TMatrixDSym GPKernel::operator()(std::vector<GPPoint> pts, int dpar_idx) const
{
  TMatrixDSym sym(pts.size());
  for(int i = 0; i < (int)pts.size(); i++){
    for(int j = 0; j <= i; j++){
      double element = Element(pts[i], pts[j], dpar_idx);
      sym(i, j) = element;
      sym(j, i) = element;
    }
  }
  return sym;
}

//----------------------------------------------------------------------
TMatrixD GPKernel::KernelBlock(std::vector<GPPoint> pts1, std::vector<GPPoint> pts2) const
{
  TMatrixD block(pts1.size(), pts2.size());
  for(int i = 0; i < (int)pts1.size(); i++){
    for(int j = 0; j < (int)pts2.size(); j++){
      block(i, j) = Element(pts1[i], pts2[j]);
    }
  }
  return block;
}

//----------------------------------------------------------------------
TVectorD GPKernel::KernelDiag(std::vector<GPPoint> pts) const
{
  TVectorD diag(pts.size());
  for(int i = 0; i < (int)pts.size(); i++)
    diag[i] = Element(pts[i], pts[i]);
  return diag;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
double RBFKernel::Element(GPPoint pt1, GPPoint pt2, int dpar_idx) const
{
  assert(dpar_idx < 2 && "RBF Kernel has only 2 hyperparameters, index is either -1, 0 or 1");
 
  double scale = fPars[0];
  double coeff = fPars[1];
  if(pt1.Mag(pt2) == 0.){
    if(dpar_idx == 0) 
      return 0.;
    else
      return coeff * (1. + fNoise);
  }

  double val = pt1.Mag(pt2)/(scale*scale);
  if(dpar_idx == 0)
    return coeff * exp(-val/2.)*val;
  else
    return coeff * exp(-val/2.);
}

//----------------------------------------------------------------------
void RBFKernel::SetThetas(std::vector<double> thetas)
{
  assert(thetas.size() == fPars.size() && "Theta vector should be of size 2 for RBF Kernel");
  fPars[0] = exp(thetas[0]);
  fPars[1] = exp(thetas[1]);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
double RationalQuadraticKernel::Element(GPPoint pt1, GPPoint pt2, int dpar_idx) const
{
  assert(dpar_idx < 3 && "Rational Quadratic Kernel has only 3 hyperparameters, index is either -1, 0, 1 or 2");
  
  double scale = fPars[0];
  double alpha = fPars[1];
  double coeff = fPars[2];
  if(pt1.Mag(pt2) == 0.){
    if(dpar_idx >= 0 && dpar_idx < 2) 
      return 0.;
    else
      return coeff * (1. + fNoise);
  }

  double dists = pt1.Mag(pt2);
  double base = dists/(2 * alpha * TMath::Power(scale, 2));
  base += 1.;

  double val = TMath::Power(base, -alpha);

  if(dpar_idx == 0)
    return coeff * dists * val/(base * TMath::Power(scale, 2));
  else if(dpar_idx == 1)
    return coeff * val * (-alpha * log(base) + dists / (2 * TMath::Power(scale, 2) * base));
  else
    return coeff * val;
}

//----------------------------------------------------------------------
void RationalQuadraticKernel::SetThetas(std::vector<double> thetas)
{
  assert(thetas.size() == fPars.size() && "Theta vector should be of size 3 for Rational Quadratic Kernel");
  fPars[0] = exp(thetas[0]);
  fPars[1] = exp(thetas[1]);
  fPars[2] = exp(thetas[2]);
}
