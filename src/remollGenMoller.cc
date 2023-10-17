#include "remollGenMoller.hh"

#include "Randomize.hh"

#include "G4Material.hh"
#include "G4PhysicalConstants.hh"

#include "remollEvent.hh"
#include "remollVertex.hh"
#include "remolltypes.hh"

remollGenMoller::remollGenMoller()
  : remollVEventGen("moller") {
  fThCoM_min =    30.0*deg;
  fThCoM_max =   150.0*deg;

  fApplyMultScatt = true;
}

remollGenMoller::~remollGenMoller(){
}

void remollGenMoller::SamplePhysics(remollVertex *vert, remollEvent *evt){
  // Generate Moller event

  double beamE = vert->GetBeamEnergy();
  double me    = electron_mass_c2;

  double beta_com  = sqrt( (beamE - me)/(beamE + me) );
  double gamma_com = 1.0/sqrt(1.0 - beta_com*beta_com);

  double e_com = me*gamma_com;
  double thcom = acos(G4RandFlat::shoot(cos(fThCoM_max), cos(fThCoM_min)));
  double phcom = G4RandFlat::shoot(fPh_min, fPh_max);

  double sigma = alpha*alpha*pow(3.0+cos(thcom)*cos(thcom),2.0)*hbarc*hbarc/pow(sin(thcom),4.0)/(2.0*me*beamE); // units of area

  double V = (fPh_max - fPh_min)*(cos(fThCoM_min) - cos(fThCoM_max));

  double s = 2.0*me*beamE;

  double beamE_lab = (beamE/2.0)*(1+cos(thcom));

  double cos_thcom_lab = 1-((me/beamE)*((1-cos(thcom))/(1+cos(thcom))));

  double t = -2*beamE*beamE_lab*(1-cos_thcom_lab);

  double u = -2*me*beamE_lab;

  double dSigma_Born_dOmega = (alpha*alpha/(2.0*s))*((t*t+t*u+u*u)/(t*u))*((t*t+t*u+u*u)/(t*u));

  double dSigma_phi_dOmega = -(((alpha*alpha*alpha*me/(8.0*sqrt(s)))*sin(thcom)*sin(phcom)*(1.0/(t*t*u*u))*(3.0*s*((t*(u-s)*log(-t/s))-(u*(t-s)*log(-u/s)))-(2.0*(t-u)*t*u))));

  double alpha_T = (1/sin(phcom))*((dSigma_phi_dOmega)/(dSigma_Born_dOmega));
 
  double  pi  = 3.14159265358979323846;


  //  Multiply by Z because we have Z electrons
  //  here we must also divide by two because we are double covering
  //  phasespace because of identical particles

  evt->SetEffCrossSection(sigma*V*vert->GetMaterial()->GetZ()/2.0);

  if( vert->GetMaterial()->GetNumberOfElements() != 1 ){
    G4cerr << __FILE__ << " line " << __LINE__ <<
      ": Error!  Some lazy programmer didn't account for complex materials in the moller process!" << G4endl;
    exit(1);
  }

  G4double APV = electron_mass_c2*beamE*GF*4.0*sin(thcom)*sin(thcom)*(QWe+QWe_rad)/(sqrt(2.0)*pi*alpha*pow(3.0+cos(thcom)*cos(thcom),2.0));
  G4double A_TV = alpha_T*sin(phcom);
  G4double A_TH = alpha_T*sin(phcom-CLHEP::pi/2);
  evt->SetAsymmetry(APV);
  evt->SetAsymmetryTransverseVertical(A_TV);
  evt->SetAsymmetryTransverseHorizontal(A_TH);
  evt->SetThCoM(thcom);
//  G4cout << "A_TV = " << A_TV << G4endl;
//  G4cout << "A_TH = " << A_TH << G4endl;

  //evt->SetQ2( 2.0*e_com*e_com*(1.0-cos(thcom)) );
  // Q2 is not actually well defined
  evt->SetQ2( 0.0 );

  double pperp = e_com*sin(thcom);
  double ppar  = e_com*cos(thcom);

  evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0),
                           G4ThreeVector(pperp*cos(phcom), pperp*sin(phcom), gamma_com*(ppar + e_com*beta_com) ),
                           "e-");

  evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0),
                           G4ThreeVector(-pperp*cos(phcom), -pperp*sin(phcom), gamma_com*(-ppar + e_com*beta_com) ),
                           "e-");

  return;

}
