/*************************************************************
 *                                                           *
 *  Simulation of a cosmic telescope composed                *
       by three scintillator counters.                       *
       Instruction usage:                                    *
       To simulate 100 cosmic rays with a polar              *
       distribution proportional to Cos(theta)^2.5           *
       root prompt:                                          *
       .L CosmicMC.cc                                      *
       CosmicMC(100,3.5)                                     *
                                                             *
                                                             *
                                *************************/

                                
#include <math.h> 
#include <iostream>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
// Units: length in meters
// Units: time in nseconds
// Axis: up vertical +z, right +x, z x x = y

class CosmicRay //costruisci una classe e inzializza le variabili con double
{
  double cosTheta;
  double phi;

  double x, y;

  double exponent;

  double size_x;
  double size_y;
  double eff;
  
public:
  // Constructor: creo un oggetto
  CosmicRay( double exponent, double size_x, double size_y ):
    exponent(exponent),
    size_x(size_x),
    size_y(size_y) { Throw(); };
 
  void Throw( void ){  
    phi = gRandom->Uniform() * 2 * M_PI;
    cosTheta = pow(gRandom->Uniform() , 1/(exponent+1) );
    x = gRandom->Uniform(size_x);
    y = gRandom->Uniform(size_y);
  }

  double CosTheta(void) const { return cosTheta ;};
  double Phi(void)      const { return phi ;     };
  double X(void)        const { return x ;       };
  double Y(void)        const { return y ;       };
  double X(double z)    const { return x -
      sqrt(1- pow(CosTheta(),2))/CosTheta() * cos( Phi()) * z ;};
  double Y(double z)    const { return y -
      sqrt(1- pow(CosTheta(),2))/CosTheta() * sin( Phi()) * z ;};
};

class Detector
{
  double xmin, xmax, ymin, ymax, z;
  double eff;//aggiunta l'efficienza
public:
  Detector( double xmin, double xmax,
            double ymin, double ymax,
            double z, double eff=0.7):
    xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), z(z), eff{eff}
  {};

  bool CheckIntersect( const CosmicRay & ray ){
    double x_on_plane = ray.X(z);
    double y_on_plane = ray.Y(z);
    return x_on_plane > xmin &&
           x_on_plane < xmax &&
           y_on_plane > ymin &&
           y_on_plane < ymax;
  };

  bool CheckEff(const CosmicRay & ray){//aggiunta l'efficienza
    return CheckIntersect( ray ) &&
            gRandom->Uniform()<0.7;
  };
};

void CosmicMC(int n, double exponent)
{
  TH1F * exponentCheck = new TH1F( "exponentCheck_pippo","Check cosTheta", 100, 0, 1);
  TH2F * intercept = new TH2F( "intercept","intercept", 300, -.4, .8, 150, -.2, .4);
  CosmicRay ray( exponent, .40, .20);
  Detector down( 0., 0.4, 0., .20, -.05);
  Detector up  ( 0., .40, 0., .20, +.05);
  
  int pairs(0), triples(0);
  for(int i = 0; i < n; i++)
    {
      ray.Throw();
      exponentCheck -> Fill( ray.CosTheta() );

      pairs   += down.CheckIntersect( ray ) ? 1 : 0;
      triples += down.CheckIntersect( ray ) && up.CheckEff( ray ) ? 1 : 0;
      if( down.CheckEff( ray ) )
        intercept -> Fill( ray.X(+.05), ray.Y(+.05) );
    }

  std::cout << triples/double(pairs) << std::endl;
  TF1 f("f","[0] * pow(x,[1])", 0, 1); 
  exponentCheck->Fit("f");
  gStyle->SetOptFit(11111);

  exponentCheck->Draw();
}