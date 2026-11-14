
/**************************************************************
 *
 * COSMIC RAY TELESCOPE - MONTE CARLO
 *
 * Geometry:
 *   - Upper scintillator = generation plane
 *   - One lower planar detector
 *   - 6 parallelepipeds:
 *       2 columns along X
 *       3 consecutive levels along Z
 *
 * Output:
 *   - Triple coincidences for each box
 *   - Geometrical acceptance for each box
 *   - One TNtuple per box
 *   - Entry and exit coordinates for each crossing
 *
 * Units: meters
 * Angles: radians
 *
 * ROOT:
 *   .L CosmicMC.cc++
 *   CosmicMC(100000,2.5)
 *
 **************************************************************/

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <vector>

#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TMath.h>


// ============================================================
// POINT IN 3D
// ============================================================

struct Point3D
{
    double x,y,z;
};


// ============================================================
// COSMIC RAY
// ============================================================

class CosmicRay
{
    double cosTheta, phi;
    double x0,y0,z0;
    double exponent;

    double xmin,xmax;
    double ymin,ymax;

public:

    CosmicRay(double exponent,
              double xmin,double xmax,
              double ymin,double ymax,
              double z)
        :
        exponent(exponent),
        xmin(xmin),xmax(xmax),
        ymin(ymin),ymax(ymax),
        z0(z)
    {}


    void Throw()
    {
        phi = gRandom->Uniform(0.,2.*TMath::Pi());

        cosTheta = pow(
            gRandom->Uniform(),
            1./(exponent+1.)
        );

        x0 = gRandom->Uniform(xmin,xmax);
        y0 = gRandom->Uniform(ymin,ymax);
    }


    double X0() const { return x0; }
    double Y0() const { return y0; }
    double Z0() const { return z0; }

    double CosTheta() const { return cosTheta; }
    double Phi() const { return phi; }


    // Direction versor: downward propagation

    double Ux() const
    {
        return sqrt(std::max(0.,1.-cosTheta*cosTheta))
               *cos(phi);
    }

    double Uy() const
    {
        return sqrt(std::max(0.,1.-cosTheta*cosTheta))
               *sin(phi);
    }

    double Uz() const
    {
        return -cosTheta;
    }


    // Position at distance s along the trajectory

    Point3D PointAt(double s) const
    {
        return {
            x0+s*Ux(),
            y0+s*Uy(),
            z0+s*Uz()
        };
    }


    // Intersection with a horizontal plane

    bool PointAtZ(double z,double &x,double &y) const
    {
        if(fabs(Uz())<1.e-14)
            return false;

        double s=(z-z0)/Uz();

        if(s<0.)
            return false;

        Point3D p=PointAt(s);

        x=p.x;
        y=p.y;

        return true;
    }

};


// ============================================================
// PLANAR DETECTOR
// ============================================================

class Detector
{
    double xmin,xmax;
    double ymin,ymax;
    double z;

public:

    Detector(double xmin,double xmax,
             double ymin,double ymax,
             double z)
        :
        xmin(xmin),xmax(xmax),
        ymin(ymin),ymax(ymax),
        z(z)
    {}


    bool Hit(const CosmicRay &ray,
             double &xHit,
             double &yHit) const
    {
        if(!ray.PointAtZ(z,xHit,yHit))
            return false;

        return xHit>=xmin && xHit<=xmax &&
               yHit>=ymin && yHit<=ymax;
    }

};


// ============================================================
// PARALLELEPIPED
// ============================================================

class Parallelepiped
{
    double xmin,xmax;
    double ymin,ymax;
    double zmin,zmax;

public:

    Parallelepiped(double xmin,double xmax,
                   double ymin,double ymax,
                   double zmin,double zmax)
        :
        xmin(xmin),xmax(xmax),
        ymin(ymin),ymax(ymax),
        zmin(zmin),zmax(zmax)
    {}


    // Ray-box intersection using the slab method.
    // Returns entry and exit Cartesian coordinates.

    bool Intersect(const CosmicRay &ray,
                   Point3D &entry,
                   Point3D &exit) const
    {
        double p[3]={
            ray.X0(),ray.Y0(),ray.Z0()
        };

        double d[3]={
            ray.Ux(),ray.Uy(),ray.Uz()
        };

        double lower[3]={
            xmin,ymin,zmin
        };

        double upper[3]={
            xmax,ymax,zmax
        };

        double sEnter=0.;

        double sExit=
            std::numeric_limits<double>::infinity();

        const double EPS=1.e-14;

        for(int i=0;i<3;i++)
        {
            if(fabs(d[i])<EPS)
            {
                if(p[i]<lower[i] || p[i]>upper[i])
                    return false;
            }
            else
            {
                double s1=(lower[i]-p[i])/d[i];
                double s2=(upper[i]-p[i])/d[i];

                if(s1>s2)
                    std::swap(s1,s2);

                sEnter=std::max(sEnter,s1);
                sExit=std::min(sExit,s2);

                if(sEnter>sExit)
                    return false;
            }
        }

        // Reject backward intersections and tangencies

        if(sExit<=0. || sExit-sEnter<=1.e-12)
            return false;

        entry=ray.PointAt(sEnter);
        exit=ray.PointAt(sExit);

        return true;
    }

};


// ============================================================
// MONTE CARLO
// ============================================================

void CosmicMC(int n,double exponent)
{
    if(n<=0 || exponent<=-1.)
    {
        std::cerr<<"Invalid n or exponent"<<std::endl;
        return;
    }


    // ========================================================
    // 1. GEOMETRY CONFIGURATION
    // ========================================================


    // --------------------------------------------------------
    // COSMIC RAY GENERATION AREA
    // Upper scintillator
    // --------------------------------------------------------

    const double genXmin=0.00;
    const double genXmax=0.40;

    const double genYmin=0.00;
    const double genYmax=0.51;

    const double genZ=+0.;


    // --------------------------------------------------------
    // LOWER PLANAR DETECTOR
    // --------------------------------------------------------

    const double lowerXmin=0.00;
    const double lowerXmax=0.40;

    const double lowerYmin=0.00;
    const double lowerYmax=0.51;

    const double lowerZ=-0.425;


    // --------------------------------------------------------
    // SIX PARALLELEPIPEDS
    //
    // Two columns along X
    // Three levels along Z
    // --------------------------------------------------------


    // Position of the first column

    const double gridXmin=0.09;

    // Common Y position

    const double gridYmin=0.09;


    // Upper face of the first level

    const double gridZmax=-0.05;


    // Dimensions of each parallelepiped

    const double boxDX=0.11;
    const double boxDY=0.38;
    const double boxDZ=0.105;


    // Gaps between adjacent boxes

    const double gapX=0.00;
    const double gapZ=0.00;


    // ========================================================
    // 2. GEOMETRY INITIALIZATION
    // ========================================================

    const int nCols=2;
    const int nLevels=3;
    const int nBoxes=nCols*nLevels;


    if(genZ<=lowerZ)
    {
        std::cerr<<"ERROR: generation plane must be above "
                 <<"lower detector"<<std::endl;
        return;
    }

    if(boxDX<=0. || boxDY<=0. || boxDZ<=0. ||
       gapX<0. || gapZ<0.)
    {
        std::cerr<<"ERROR: invalid box dimensions or gaps"
                 <<std::endl;
        return;
    }


    CosmicRay ray(
        exponent,
        genXmin,genXmax,
        genYmin,genYmax,
        genZ
    );


    Detector lower(
        lowerXmin,lowerXmax,
        lowerYmin,lowerYmax,
        lowerZ
    );


    // --------------------------------------------------------
    // Build the six boxes
    //
    // Numbering:
    //
    //       BOX 1 | BOX 4
    //       BOX 2 | BOX 5
    //       BOX 3 | BOX 6
    //
    // --------------------------------------------------------

    std::vector<Parallelepiped> boxes;

    boxes.reserve(nBoxes);

    std::cout<<"\n======= BOX GEOMETRY =======\n";

    for(int col=0;col<nCols;col++)
    {
        for(int level=0;level<nLevels;level++)
        {
            int ib=col*nLevels+level;

            double xmin=gridXmin+col*(boxDX+gapX);
            double xmax=xmin+boxDX;

            double ymin=gridYmin;
            double ymax=ymin+boxDY;

            double zmax=gridZmax-level*(boxDZ+gapZ);
            double zmin=zmax-boxDZ;

            boxes.emplace_back(
                xmin,xmax,
                ymin,ymax,
                zmin,zmax
            );

            std::cout
                <<"Box "<<ib+1
                <<": x = ["<<xmin<<", "<<xmax<<"]"
                <<", y = ["<<ymin<<", "<<ymax<<"]"
                <<", z = ["<<zmin<<", "<<zmax<<"]"
                <<std::endl;
        }
    }


    // Check whether the entire box array lies between planes

    double arrayZmin=
        gridZmax-(nLevels-1)*(boxDZ+gapZ)-boxDZ;

    if(arrayZmin<lowerZ || gridZmax>genZ)
    {
        std::cerr
            <<"WARNING: some boxes are outside the two "
            <<"scintillator planes."
            <<std::endl;
    }


    // ========================================================
    // 3. ROOT OUTPUT
    // ========================================================

    TFile output("CosmicMC_6boxes.root","RECREATE");

    if(output.IsZombie())
    {
        std::cerr<<"ERROR opening ROOT file"<<std::endl;
        return;
    }


    TH1F *exponentCheck=new TH1F(
        "exponentCheck",
        "Check cos(#theta)",
        100,0.,1.
    );

    TH2F *intercept=new TH2F(
        "intercept",
        "Intercept on lower detector",
        100,lowerXmin,lowerXmax,
        100,lowerYmin,lowerYmax
    );


    // Histograms survive after closing the ROOT file

    exponentCheck->SetDirectory(nullptr);
    intercept->SetDirectory(nullptr);


    // One TNtuple per box

    std::vector<TNtuple*> tracks(nBoxes);

    for(int ib=0;ib<nBoxes;ib++)
    {
        TString name=
            TString::Format("tracks_box%d",ib+1);

        TString title=
            TString::Format("Tracks in box %d",ib+1);

        tracks[ib]=new TNtuple(
            name.Data(),
            title.Data(),
            "event:hitLower:"
            "xGen:yGen:zGen:"
            "xLower:yLower:"
            "xIn:yIn:zIn:"
            "xOut:yOut:zOut:"
            "cosTheta:phi"
        );
    }


    // ========================================================
    // 4. COUNTERS
    // ========================================================

    int doubles=0;

    std::vector<int> boxCrossings(nBoxes,0);
    std::vector<int> triples(nBoxes,0);

    int tripleAnyBox=0;


    // ========================================================
    // 5. MONTE CARLO LOOP
    // ========================================================

    for(int i=0;i<n;i++)
    {
        ray.Throw();

        exponentCheck->Fill(ray.CosTheta());


        // ----------------------------------------------------
        // Lower planar detector
        // ----------------------------------------------------

        double xLower=-999.;
        double yLower=-999.;

        bool hitLower=lower.Hit(
            ray,xLower,yLower
        );

        if(hitLower)
        {
            doubles++;

            intercept->Fill(xLower,yLower);
        }


        // ----------------------------------------------------
        // Six parallelepipeds
        // ----------------------------------------------------

        bool hitAnyBox=false;

        for(int ib=0;ib<nBoxes;ib++)
        {
            Point3D entry,exit;

            bool hitBox=boxes[ib].Intersect(
                ray,entry,exit
            );

            if(!hitBox)
                continue;

            hitAnyBox=true;

            boxCrossings[ib]++;


            // Triple coincidence for this box

            if(hitLower)
                triples[ib]++;


            // Fill the corresponding TNtuple

            tracks[ib]->Fill(
                float(i),
                float(hitLower),

                float(ray.X0()),
                float(ray.Y0()),
                float(ray.Z0()),

                float(xLower),
                float(yLower),

                float(entry.x),
                float(entry.y),
                float(entry.z),

                float(exit.x),
                float(exit.y),
                float(exit.z),

                float(ray.CosTheta()),
                float(ray.Phi())
            );
        }


        // Count the event only once

        if(hitLower && hitAnyBox)
            tripleAnyBox++;
    }


    // ========================================================
    // 6. RESULTS
    // ========================================================

    std::cout<<"\n====================================\n";

    std::cout<<"GENERATED EVENTS: "<<n<<std::endl;

    std::cout<<"DOUBLE COINCIDENCES: "
             <<doubles<<std::endl;


    std::cout<<"\nTRIPLE COINCIDENCES PER BOX\n";

    for(int ib=0;ib<nBoxes;ib++)
    {
        std::cout<<"\nBox "<<ib+1<<std::endl;

        std::cout<<"Geometrical crossings: "
                 <<boxCrossings[ib]<<std::endl;

        std::cout<<"Triple coincidences: "
                 <<triples[ib]<<std::endl;

        if(doubles>0)
        {
            double acceptance=
                double(triples[ib])/doubles;

            std::cout<<"Acceptance: "
                     <<100.*acceptance<<" %"
                     <<std::endl;
        }
        else
        {
            std::cout<<"Acceptance undefined"
                     <<std::endl;
        }
    }


    std::cout<<"\nTriple coincidences with at least "
             <<"one box: "<<tripleAnyBox<<std::endl;

    if(doubles>0)
    {
        std::cout<<"Total acceptance: "
                 <<100.*tripleAnyBox/doubles
                 <<" %"<<std::endl;
    }


    // ========================================================
    // 7. ANGULAR DISTRIBUTION FIT
    // ========================================================

    TF1 fit(
        "fitCosTheta",
        "[0]*pow(x,[1])",
        0.,1.
    );

    fit.SetParameters(
        exponentCheck->GetMaximum(),
        exponent
    );

    gStyle->SetOptFit(1111);

    exponentCheck->Fit(&fit,"R");


    // ========================================================
    // 8. SAVE ROOT OBJECTS
    // ========================================================

    output.cd();

    exponentCheck->Write();
    intercept->Write();

    for(int ib=0;ib<nBoxes;ib++)
        tracks[ib]->Write();

    output.Close();


    // Histograms are still valid in memory

    exponentCheck->Draw();

}
