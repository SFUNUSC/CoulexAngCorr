// nasty bit of coding practice
// but at least its modular...
#include "coulex_ang_corr.h"
#include "read_config.cpp"
#include "build_matrices.cpp"

int main(int argc, char *argv[])
{
  //Creating an instance of TApplication
  //This is evidently needed for auto-loading of ROOT libraries, 
  //otherwise the program may crash upon execution depending on how ROOT
  //is set up.
  //http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=14064
  //
  // Not always necessary for some reason... 
  //
  // int ac;
  // char* av[10];
  // theApp=new TApplication("App", &ac, av);

  FILE *list,*output;
  TFile *inp;
  randGen = new TRandom3();

  newWeight=0.;

  if(argc!=2)
    {
      printf("\ncoulex_ang_dist parameter_file\n");
      printf("-----------------------\nSorts TIGRESS add-back spectra gated by TIGRESS-CsI Doppler shift groups defined via a group map using the Coulex angular distribution weight modification based on Alder, Rev. Mod. Phys 28 (1956).\n\n");
      exit(-1);
    }
  
  readConfigFile(argv[1],"coulex_ang_dist"); //grab data from the parameter file
  
  // read in tree list file
  if((list=fopen(inp_filename,"r"))==NULL)
    {
      printf("ERROR: Cannot open the input list file %s!\n",inp_filename);
      exit(-1);
    } 
  // scan the list file for ROOT files and put their data into the output hitogram

  // build matricies for angular distribution
  build_matrices();

  // some useful precomputed quantities for the weight calculation
  // global because... well, because
  A = Ap/(Ap+Ar); 
  DEprime = (1.+Ap/Ar)*E0; // Eq. II C.4 from Alder

  while(fscanf(list,"%s",str)!=EOF)
    {
      inp = new TFile(str,"read");
      if((tree = (TTree*)inp->Get(tree_name))==NULL)
	{
	  printf("The specified tree named %s doesn't exist, trying default name 'tree'.\n",tree_name);
	  if((tree = (TTree*)inp->Get("tree"))==NULL) // try the default tree name
	    {
	      printf("ERROR: The specified tree named %s (within the ROOT file) cannot be opened!\n",tree_name);
	      exit(-1);
	    }
	}
      printf("Tree in %s read out.\n",str);

      // sort leaf (energy)
      if((sortLeaf = tree->GetLeaf(sort_path))==NULL)
	if((sortBranch = tree->GetBranch(sort_path))==NULL)
	  {
	    printf("ERROR: Sort data path '%s' doesn't correspond to a branch or leaf in the tree!\n",sort_path);
	    exit(-1);
	  }
      if(sortLeaf==NULL)
	sortLeaf = (TLeaf*)sortBranch->GetListOfLeaves()->First(); // get the first leaf from the specified branch  
      
      // position leaf
      if((posLeaf = tree->GetLeaf(pos_path))==NULL)
	if((posBranch = tree->GetBranch(pos_path))==NULL)
	  {
	    printf("ERROR: Pos data path '%s' doesn't correspond to a branch or leaf in the tree!\n",pos_path);
	    exit(-1);
	  }
      if(posLeaf==NULL)
	posLeaf = (TLeaf*)posBranch->GetListOfLeaves()->First(); // get the first leaf from the specified branch  
      
      // color leaf
      if((colLeaf = tree->GetLeaf(col_path))==NULL)
	if((colBranch = tree->GetBranch(col_path))==NULL)
	  {
	    printf("ERROR: Col data path '%s' doesn't correspond to a branch or leaf in the tree!\n",col_path);
	    exit(-1);
	  }
      if(colLeaf==NULL)
            colLeaf = (TLeaf*)colBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch  
      
      // csi leaf
      if((csiLeaf = tree->GetLeaf(csi_path))==NULL)
	if((csiBranch = tree->GetBranch(csi_path))==NULL)
	  {
	    printf("ERROR: CsI data path '%s' doesn't correspond to a branch or leaf in the tree!\n",csi_path);
	    exit(-1);
	  }
      if(csiLeaf==NULL)
	csiLeaf = (TLeaf*)csiBranch->GetListOfLeaves()->First(); // get the first leaf from the specified branch  
      
      printf("Paths to sort data set.\n");
      printf("Number of tree entries: %Ld\n",tree->GetEntries());
      /* printf("%Ld\n",tree->GetEntries()); */
      
      readGroupMap();
      
      // read and store data from tree
      readBranchData();
      // calculate new weight
      for (int i=0;i<tree->GetEntries();i++)
	{
	  newWeight = calcWeight(i);
	  addTreeDataToOutHist(i,newWeight);
	}
    }

  // write the output histogram to disk
  if((output=fopen(out_filename,"w"))==NULL)
    {
      printf("ERROR: Cannot open the output .fmca file!\n");
      exit(-1);
    }
  for (int i=0;i<NSPECT;i++)
    fwrite(dOutHist[i],S32K*sizeof(float),1,output);
  fclose(output);
  
  return 0;
}

double calcWeight(int entry)
{
  // vectors for CM transformation
  TLorentzVector pIn,rIn,pOut,pDec,gamma;
  TVector3 vCM; // vCM vector for transformation
  TVector3 gammaR,pDecR,gammaDir; // gamma ray positions and decay vector
  double df; // doppler factor for decay in LAB frame
  double xi; // xi for distribution function. See Eq. II C.12 in Alder
  double Elab; // projectile kinetic energy in LAB frame for calculation of xi
  double thetaP,phiP,thetaG,phiG; // angles for particle and gamma in the CM
  double lambda=2.;
  double modw=0.;
  double wgt=0.;

  tree->GetEntry(entry);

  // in the LAB frame
  // velocity in v/c, energy in MeV, momenta in MeV/c, masses in MeV/c^2
  // projectile and recoil
  pIn.SetPxPyPzE(lprIn_px->GetValue(0),lprIn_py->GetValue(0),lprIn_pz->GetValue(0),lprIn_E->GetValue(0)+(AMU*Ap));
  Elab = lprIn_E->GetValue(0);
  wgt = lprOut_w->GetValue(0);
  pOut.SetPxPyPzE(lprOut_px->GetValue(0),lprOut_py->GetValue(0),lprOut_pz->GetValue(0),lprOut_E->GetValue(0)+(AMU*Ap));
  vCM = A*lprIn_b->GetValue(0)*pIn.Vect()*(1./pIn.Vect().Mag()); 
  
  // gamma ray
  df = lprDec_df->GetValue(0);
  pDec.SetPxPyPzE(lprDec_px->GetValue(0),lprDec_py->GetValue(0),lprDec_pz->GetValue(0),lprDec_E->GetValue(0)+(AMU*Ap));
  pDecR.SetXYZ(lprDec_x->GetValue(0),lprDec_y->GetValue(0),lprDec_z->GetValue(0));
  gammaR.SetXYZ(lgx->GetValue(0),lgy->GetValue(0),lgz->GetValue(0));
  gammaDir = gammaR-pDecR;
  gammaDir*=(1./gammaDir.Mag());
  gammaDir*=(df*E0);
  gamma.SetPxPyPzE(gammaDir.X(),gammaDir.Y(),gammaDir.Z(),df*E0);
  
  // std::cout << std::endl;
  // std::cout << "LAB KEIn        : " << Elab << std::endl;
  // std::cout << "LAB pIN.Vect()  : "; pIn.Vect().Print();
  // std::cout << "LAB pIN         : "; pIn.Print();
  // std::cout << "LAB pOut.Vect() : "; pOut.Vect().Print();
  // std::cout << "LAB pDec.Vect() : "; pDec.Vect().Print();
  // std::cout << "LAB gamma.Vect(): "; gamma.Vect().Print();
  // std::cout << "vCM vCM.Vect()  : "; vCM.Print();
  // std::cout << std::endl;

  // boost by into CM frame moving relative to the LAB
  // pIn.Boost(-vCM);
  pOut.Boost(-vCM);
  // pDec.Boost(-vCM);
  gamma.Boost(-vCM);  

  // std::cout << "CM pIN.Vect()  : "; pIn.Vect().Print();
  // std::cout << "CM pIN         : "; pIn.Print();
  // std::cout << "CM pOut.Vect() : "; pOut.Vect().Print();
  // std::cout << "CM pDec.Vect() : "; pDec.Vect().Print();
  // std::cout << "CM gamma.Vect(): "; gamma.Vect().Print();
  // std::cout << std::endl;

  // rotation to coordinate system where the incident 
  // particle momentum vector is along the +z axis
  TVector3 ez(0.,0.,1.); // +z axis
  double angle = pIn.Angle(ez);
  // std::cout << "angle " << angle*RAD2DEG << std::endl;  
  TVector3 axis = pIn.Vect().Cross(ez); // perp. to pIn and +z axis
  // pIn.Rotate(angle,axis);
  pOut.Rotate(angle,axis);
  // pDec.Rotate(angle,axis);
  gamma.Rotate(angle,axis);

  // precision problem with the rotation so check and fix here
  // can access by index number x:0, y:1, z:2 t:3
  // see https://root.cern.ch/doc/master/classTLorentzVector.html
  // double eps = 1E-9;
  // for(int i=0;i<3;i++) // only x,y,z
  //   {
  //     if(fabs(pIn(i)) < eps)
  // 	pIn(i) = 0.;
  //   }
  
  // std::cout << "CM ROT pIN.Vect()  : "; pIn.Vect().Print();
  // std::cout << "CM ROT pOut.Vect() : "; pOut.Vect().Print();
  // std::cout << "CM ROT pDec.Vect() : "; pDec.Vect().Print();
  // std::cout << "CM ROT gamma.Vect(): "; gamma.Vect().Print();
  // std::cout << std::endl;
  
  thetaP = pOut.Vect().Theta();
  phiP = pOut.Vect().Phi();
  thetaG = gamma.Vect().Theta();
  phiG = gamma.Vect().Phi();
  xi = getXi(Elab,DEprime);
  
  // std::cout << "thetaP | phiP | thetaG | phiG | xi"<< std::endl;
  // std::cout << RAD2DEG*thetaP << " " << RAD2DEG*phiP << " " << RAD2DEG*thetaG << " " << RAD2DEG*phiG << " " << xi << std::endl;
  // std::cout << std::endl;
  
  for(int k=0;k<=4;k+=2) // sum over k = 0, 2, 4
    {
      double term = 0.;
      for(int ka=0;ka<=k;ka++) // sum over kappa = 0, 1, ... , k
	{
	  double norm = 1.;
	  if(ka == 0)
	    norm = 0.5;
	  	  
	  if(ka % 2 == 0) // kappa even
	    {
	      // std::cout << "sqrt factorial " << sqrt((double)factorial(k-ka)/(double)factorial(k+ka)) << std::endl;
	      // std::cout << "P " << associatedLegendre(k,ka,cos(thetaG)) << std::endl;
	      // std::cout << "ujm " << calcUJM(k,ka,thetaP,xi) << std::endl;
	      // std::cout << "norm " << norm << std::endl;
	      // std::cout << "cosine " << cos(ka*(phiG-phiP)) << std::endl;
	      term += sqrt((double)factorial(k-ka)/(double)factorial(k+ka))*associatedLegendre(k,ka,cos(thetaG))*calcUJM(k,ka,thetaP,xi)*norm*cos(ka*(phiG-phiP));
	      // std::cout << "term " << term << std::endl;
	    }
	  else 
	    {
	      // std::cout << "sqrt factorial " << sqrt((double)factorial(k-ka)/(double)factorial(k+ka)) << std::endl;
	      // std::cout << "P " << associatedLegendre(k,ka,cos(thetaG)) << std::endl;
	      // std::cout << "wjm " << calcWJM(k,ka,thetaP,xi) << std::endl;
	      // std::cout << "sine " << sin(ka*(phiG-phiP)) << std::endl;
	      // term += -1.*sqrt((double)factorial(k-ka)/(double)factorial(k+ka))*associatedLegendre(k,ka,cos(thetaG))*calcWJM(k,ka,thetaP,xi)*sin(ka*(phiG-phiP));
	      term += sqrt((double)factorial(k-ka)/(double)factorial(k+ka))*associatedLegendre(k,ka,cos(thetaG))*calcWJM(k,ka,thetaP,xi)*sin(ka*(phiG-phiP));
	      // std::cout << "term " << term << std::endl;
	    }
	}
      // std::cout << "coeff " << getCoeff(k) << std::endl;
      modw += getCoeff(k)*term;
      // std::cout << "norm term " << getCoeff(k)*term*-1.*CONST/sqrt(2*lambda+1)/calcVJMthxi(0,0,RAD2DEG*thetaP,xi) << std::endl;
      // std::cout << " --- calculation for k = " << k << " complete." << std::endl;
      // getc(stdin);
    }
  double v00 = calcVJMthxi(0,0,RAD2DEG*thetaP,xi);
  // std::cout << "v00 " << v00 << std::endl;
  modw *= -1.*CONST/sqrt(2*lambda+1)/v00;
  // std::cout << "original weight " << wgt << std::endl;
  // std::cout << "mod factor      " << modw << std::endl;
  if(modw < 0.)
    {
      std::cout << "Entry " << entry << " modification factor " << modw << " set to 0." << std::endl;
      modw = 0.;
      // getc(stdin);
    }
  modw *= wgt;
  // std::cout << "modified weight " << modw << std::endl;
  // getc(stdin);
  return modw;
}

/********************************************************
Calculates xi using equations from Alder.
Same as function in the TIP Coulex Geant4 code.
********************************************************/
double getXi(double KE,double DEp)
{
  double zeta,sqz;
  double eti,etf,nu,xi;
  zeta=DEp/KE;                       // Eq. II C.5 dimensionless
  sqz=sqrt(1.-zeta);
  eti=0.5*Zp*Zr*sqrt(Ap/10.008/KE);  // Eq. II C.8
  etf=eti/sqz;                       // Eq. II C.9
  nu=2.*sqrt(1./eti/eti-1./etf/etf); // Eq. II C.10 top
  xi=2./nu*sqrt(zeta)*(1./sqz-1.);   // Eq. II C.12
  return xi;
}

/********************************************************
Calculates value of VJM(theta,xi) from pre-computed
matrix for a given input J,M,theta,xi using simple 
bilinear interpolation. Returns result of interpolation.
********************************************************/
double calcVJMthxi(int J, int M, double theta, double xi)
{
  vector<double>::iterator itth = thArray.begin();
  vector<double>::iterator itxi = xiArray.begin();
  double x1,x2,y1,y2,dy,f11,f12,f21,f22;
  int shiftXi,shiftTh;
  vector< vector<double> > VJM; // VJM(theta,xi)

  // which matrix to interpolate from
  int cnd = J*M+J;

  if(cnd == 0)
    VJM = V00;
  if(cnd == 2)
    VJM = V20;
  if(cnd == 6)
    VJM = V22;
  if(cnd == 4)
    VJM = V40;
  if(cnd == 12)
    VJM = V42;
  if(cnd == 20)
    VJM = V44;

  shiftTh=0;
  shiftXi=0;

  // theta is x
  double x = theta;
  for(;itth<thArray.end()-1;itth++)
    {
    if((theta>=(*itth))&&(theta<=(*(itth+1))))
      break;
    shiftTh++;
   }
  x1=*itth;
  x2=*(itth+1);

  // xi is y
  double y = xi;
  for(;itxi<xiArray.end()-1;itxi++)
    {
    if((xi>=(*itxi))&&(xi<=(*(itxi+1))))
      break;
    shiftXi++;
   }
  y1=*itxi;
  y2=*(itxi+1);
  dy=y2-y1; // depends on bin in this case!

  vector<double>::iterator it;
  vector< vector<double> >::iterator itVJM=VJM.begin()+shiftTh;
  it = (*itVJM).begin()+shiftXi;
  f11 = *it;
  f12 = *(it+1);
  itVJM++;
  it = (*itVJM).begin()+shiftXi;
  f21 = *it;
  f22 = *(it+1);

  // interpolate in x (theta) direction first
  // dx always = 10 degrees for theta
  double fx1 = ((x2-x)/10.)*f11 + ((x-x1)/10.)*f21;
  double fx2 = ((x2-x)/10.)*f12 + ((x-x1)/10.)*f22;
  // now interpolate between fx1 and fx2 to get fxy
  double fxy = ((y2-y)/dy)*fx1 + ((y-y1)/dy)*fx2;

  return fxy;
}

/********************************************************
Calculates value of UJM(theta,xi) when M even
Lookup table for VJM assumes theta in deg.
Calculation of cosine assumes angle in rad.
********************************************************/
double calcUJM(int J, int M, double thetaP, double xi)
{
  double ujm=0.;
  for(int kap=0;kap<=J;kap+=2) // sum over kappa' (even)
    {
      // std::cout << "--> k " << J << " kappa " << M << " kappa' " << kap << std::endl;
      double norm = 1.;
      if (kap == 0)
	norm = 0.5;
      ujm += calcVJMthxi(J,kap,RAD2DEG*thetaP,xi)*wignerSmallD(J,kap,M)*norm*cos(kap*(PI2+0.5*thetaP));
      // std::cout << "ujm " << J << M << kap << " " << calcVJMthxi(J,kap,RAD2DEG*thetaP,xi)*wignerSmallD(J,kap,M)*norm*cos(kap*(PI2+0.5*thetaP)) << std::endl;

    }
  // std::cout << "ujm " << ujm << std::endl;
  return ujm;
}

/********************************************************
Calculates value of WJM(theta,xi) when M odd
Lookup table for VJM assumes theta in deg.
Calculation of sine assumes angle in rad.
********************************************************/
double calcWJM(int J, int M, double thetaP, double xi)
{
  double wjm=0.;
  for(int kap=0;kap<=J;kap+=2) // sum over kappa' (even)
    {
      // std::cout << "--> k " << J << " kappa " << M << " kappa' " << kap << std::endl;
      wjm += calcVJMthxi(J,kap,RAD2DEG*thetaP,xi)*wignerSmallD(J,kap,M)*sin(kap*(PI2+0.5*thetaP));
      // std::cout << "wjm " << J << M << kap << " " << calcVJMthxi(J,kap,RAD2DEG*thetaP,xi)*wignerSmallD(J,kap,M)*sin(kap*(PI2+0.5*thetaP)) << std::endl;
    }
  // std::cout << "wjm " << wjm << std::endl;
  return wjm;
}

/********************************************************
Gets the coefficient A_J/<j1m1j2m2|J0> in the sum for 
calculating the new weights. A_J's are given in the 
Alder review, Table II.11, pg. 468. 
********************************************************/
double getCoeff(int J)
{
  double coeff=0.;
  double clebsh=0.;

  if(J == 0)
    {
      clebsh = -1.*sqrt(1./5.);
      coeff = 1./clebsh;
    }  

  if(J == 2)
    {
      clebsh = sqrt(1./14.);
      coeff = 0.3571/clebsh;
    }

  if(J == 4)
    {
      clebsh = sqrt(8./35.);
      coeff = 1.143/clebsh;
    }

  return coeff;
}

/********************************************************
Calculates values of the Wigner small d matrices
See Varshalovich QTAM Tables 4.16 and 4.20
********************************************************/
double wignerSmallD(int J,int M1,int M2)
{
  double d=0.;
  vector< vector<double> > dJ;

  if(J == 0)
    d=1.;
  else
    {
      if(J == 2)
	dJ = d2;
      if(J == 4)
	dJ = d4;      
      
      int shiftM1 = J-M1;
      int shiftM2 = J-M2;  
      vector<double>::iterator it;
      vector< vector<double> >::iterator itdJ=dJ.begin()+shiftM1;
      it = (*itdJ).begin()+shiftM2;
      d = *it;
    }
  return d; 
}

/********************************************************
Calculates associated Legendre polynomials
P^M_J(x) for given input J, M, x=cos(theta)
********************************************************/
double associatedLegendre(int J, int M, double x)
{
  double pmj=0.;
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;

  // J and M restricted
  // J = 0, 2, 4
  // M = 0, 1, ... , J
  switch(J)
    {
    case 0:
      pmj = 1.;
      break;
    case 2:
      switch(M)
	{
	case 0:
	  pmj = 0.5*(3.*x2-1.);
	  break;
	case 1:
	  pmj = -3.*x*sqrt(1.-x2);
	  break;
	case 2:
	  pmj = 3.*(1.-x2);
	  break;
	}
      break;
    case 4:
      switch(M)
	{
	case 0:
	  pmj = 0.125*(35.*x4-30.*x2+3.);
	  break;
	case 1:
	  pmj = -2.5*(7.*x3-3.*x)*sqrt(1.-x2);
	  break;
	case 2:
	  pmj = 7.5*(7.*x2-1.)*(1.-x2);
	  break;
	case 3:
	  pmj = -105.*x*(1.-x2)*sqrt(1.-x2);
	  break;
	case 4:
	  pmj = 105.*(1.-x2)*(1.-x2);
	  break;
	}
      break;
    }
  return pmj;
}

/********************************************************
Recursively calculates factorial
Be careful since this has to return int
********************************************************/
int factorial(int n)
{
  if (n == 0)
    return 1;
  else
    return n*factorial(n-1);
}

/********************************************************
Reads position and momentum data from an input branch
This is very unseemly but it works
********************************************************/
void readBranchData()
{
  char leaf[256];
  char E[] = "E";
  char px[] = "px";
  char py[] = "py";
  char pz[] = "pz";
  char x[] = "x";
  char y[] = "y";
  char z[] = "z";  
  char beta[] = "b";
  char df[] = "df";
  char wgt[] = "w";

  // get information from the tree
  //
  // projectile reaction in
  strcpy(leaf,E);
  if((lprIn_E = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }
  strcpy(leaf,px);
  if((lprIn_px = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }
  strcpy(leaf,py);
  if((lprIn_py = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }
  strcpy(leaf,pz);
  if((lprIn_pz = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }
  strcpy(leaf,x);
  if((lprIn_x = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }
  strcpy(leaf,y);
  if((lprIn_y = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      	  exit(-1);
    }
  strcpy(leaf,z);
  if((lprIn_z = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }
  strcpy(leaf,beta);
  if((lprIn_b = tree->GetLeaf(prIn_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prIn_name,leaf);
      exit(-1);
    }

  // projectile reaction out
  strcpy(leaf,E);
  if((lprOut_E = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,px);
  if((lprOut_px = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,py);
  if((lprOut_py = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,pz);
  if((lprOut_pz = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,x);
  if((lprOut_x = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,y);
  if((lprOut_y = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      	  exit(-1);
    }
  strcpy(leaf,z);
  if((lprOut_z = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,beta);
  if((lprOut_b = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }
  strcpy(leaf,wgt);
  if((lprOut_w = tree->GetLeaf(prOut_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prOut_name,leaf);
      exit(-1);
    }

  // projectile decay
  strcpy(leaf,E);
  if((lprDec_E = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,px);
  if((lprDec_px = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,py);
  if((lprDec_py = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,pz);
  if((lprDec_pz = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,x);
  if((lprDec_x = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,y);
  if((lprDec_y = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      	  exit(-1);
    }
  strcpy(leaf,z);
  if((lprDec_z = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,beta);
  if((lprDec_b = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }
  strcpy(leaf,df);
  if((lprDec_df = tree->GetLeaf(prDec_name,leaf))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",prDec_name,leaf);
      exit(-1);
    }

  // gamma detection
  if((lgx = tree->GetLeaf(gx_name))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",gx_name,leaf);
      exit(-1);
    }
    if((lgy = tree->GetLeaf(gy_name))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",gx_name,leaf);
      	  exit(-1);
    }
      if((lgz = tree->GetLeaf(gz_name))==NULL)
    {
      printf("ERROR: Branch '%s' leaf '%s' not in the tree!\n",gx_name,leaf);
      exit(-1);
    }
}
/********************************************************
Reads group map from file.
********************************************************/
void readGroupMap()
{
  FILE *inp;
  char line[132];
  int  pos,col,csi,group;
  
  if((inp=fopen(group_file,"r"))==NULL)
    {
      printf("\nI can't open file %s\n",group_file);
      exit(1);
    }
  /* printf("\nTIGRESS-CsI group map read from the file %s\n",group_file); */
  
  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	while(fscanf(inp,"%d %d %d %d",&pos,&col,&csi,&group)!=EOF)
	  if(csi>=1 && csi<=24)
	    if(pos>=1 && pos<=16)
	      if(col>=0 && col<=3)
		group_map[pos][col][csi]=group;
		
    }  
  else
    {
      printf("Wrong structure of file %s\n",group_file);
      printf("Aborting sort\n");
      exit(1);
    }
  fclose(inp);
}

/********************************************************
Gaussian FWHM random number generator
********************************************************/
double FWHM_response(double ch_in)
{
  if(ch_in==0.)
    return ch_in;
  
  double ch_out,fwhm,sigma,ch;
  
  ch=ch_in/1000.; // ch_in is keV, ch is MeV
  fwhm=sqrt(fwhmF*fwhmF + fwhmG*fwhmG*ch + fwhmH*fwhmH*ch*ch); // fwhm is keV
  sigma=fwhm/2.35482;
  if(sigma>0)
    ch_out=randGen->Gaus(ch_in,sigma); // ch_out is keV
  else
    ch_out=ch_in;

  // std::cout << ch_in << " " << ch << " " << fwhm << " " << sigma << " " << ch_out << std::endl;
  // getc(stdin);
  
  return ch_out;
}

/********************************************************
Adds data to output histogram using new weight
calculated based on the angular distribution.
********************************************************/
void addTreeDataToOutHist(int entry, double weight)
{
  Double_t sort_value;
  Int_t pos,col,csi,group;
  double histVal;
  
  tree->GetEntry(entry);
  for(int j=0;j<sortLeaf->GetLen();j++)
    {
      sort_value = sortLeaf->GetValue(j); // in keV
      // std::cout << sort_value << std::endl;
      pos = posLeaf->GetValue(j);
      col = colLeaf->GetValue(j);
      csi = csiLeaf->GetValue(0); // recoil in csi only once per event
      group = group_map[pos][col][csi];
      
      // printf("pos %2d col %1d csi %2d group %d E %.3f w %.3f\n",pos,col,csi,group,sort_value,weight);
      // getc(stdin);
      
      // drop bad hpge crystals and csi
      /**************************************************************
      for 84Kr: 10.2 (38), 11.0 (40), 11.2 (42), 12.2 (46), no csi
      for 94Sr: 12.2 (46), drop corner CsI (11,15,19,23)
      **************************************************************/
      int hpge=(pos-1)*4+col; //0-3 pos1, 4-7 pos2, etc.
      
      if(hpge != 38)
      if(hpge != 40)
      if(hpge != 42)
      if(hpge != 46)
      // if(csi != 11)
      // if(csi != 15)
      // if(csi != 19)
      // if(csi != 23)
	{
	  if(sort_value >= 0.0)
	    {
	      if(fwhmResponse==false)
		histVal=sort_value*sort_scaling;
	      else
		{
		  histVal=FWHM_response(sort_value);
		  histVal=histVal*sort_scaling;
		}
	      
	      if(histVal>=0.0)
		if(histVal<S32K)
		  dOutHist[group][(int)(histVal)]+=(float)weight; //fill the output histogram
	    }
	}
    }
}
