#pragma once

#include <math.h>
#include "TMath.h"

class GPPoint
{
public:
  GPPoint(double x=0, double y=0, double z=0)
    : fX(x), fY(y), fZ(z)
  {
  }
  
  void SetX(double x) { fX = x; }
  void SetY(double y) { fY = y; }
  void SetZ(double z) { fZ = z; }
  
  double X() const { return fX; }
  double Y() const { return fY; }
  double Z() const { return fZ; }
  
  double Mag(GPPoint pt) const { 
    return TMath::Power(fX-pt.X(), 2) + TMath::Power(fY-pt.Y(), 2) + TMath::Power(fZ-pt.Z(), 2);
  };
  
  bool operator==(GPPoint pt) const {
    return (fX == pt.X()) && (fY == pt.Y()) && (fZ == pt.Z());
  };

protected:
  double fX, fY, fZ;
};
