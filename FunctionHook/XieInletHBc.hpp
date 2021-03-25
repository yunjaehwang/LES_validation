#ifndef XIEINLETHBC_HPP
#define XIEINLETHBC_HPP

#include "oneDProfile.hpp"

class XieInletHBc : public InletHBc {

private:

    //ibf data
    double (*Rd)[3];
    double (*Rod)[3];
    double (*up_it)[3];

    //uniform grid turbulence

    double XP;
    double LY;
    double LZ;
    double y_min, y_max, z_min, z_max;
    int NY;
    int NZ;
 
    // uniform grid resolution
    double dy;
    double dz;
   
    double Uref;
    double refD, wallD;

    // time scales  
    double lagT_u;
    double lagT_v;
    double lagT_w; 
  
    // length-scales 
    double LX_u, LX_v, LX_w;

    // vertical filter length-scale (yLu = 0.2*xLu  windTunnel slides)
    double LY_u, LY_v, LY_w;
    int NLY_u;
    int NLY_v;
    int NLY_w;

    int NLYW_u;
    int NLYW_v;
    int NLYW_w;

    // spanwise filter length-scale (zLu = 0.3*xLu  windTunnel slides)
    double LZ_u, LZ_v, LZ_w;
    int NLZ_u;
    int NLZ_v;
    int NLZ_w;

    int NLZW_u;
    int NLZW_v;
    int NLZW_w;

    //N=2n and summation over -N ~ N 
    int NLY2P1_u; 
    int NLZ2P1_u;
    int NLY2P1_v;
    int NLZ2P1_v; 
    int NLY2P1_w; 
    int NLZ2P1_w; 

    int NLYW2P1_u;
    int NLZW2P1_u;
    int NLYW2P1_v;
    int NLZW2P1_v;
    int NLYW2P1_w;
    int NLZW2P1_w;
  
    //filter coeffs
    double ** byz_u;
    double ** byz_v; 
    double ** byz_w;

    //filter coeffs near the wall
    double ** byzw_u; 
    double ** byzw_v;
    double ** byzw_w;

    double (*u1f)[3];

public:
  
  double kappa;

  XieInletHBc(BfZone* p, HelmholtzSolver* s) : InletHBc(p,s), kappa(0.41) {
    u_bc = NULL; zone_ptr->registerBfData(u_bc, "u_bc", CAN_WRITE_DATA);
    up_it = NULL; zone_ptr->registerBfData(up_it, "up_it", CAN_WRITE_DATA);
    Rd = NULL; zone_ptr->registerBfData(Rd, "Rd", CAN_WRITE_DATA);
    Rod = NULL; zone_ptr->registerBfData(Rod, "Rod", CAN_WRITE_DATA);
  }
  XieInletHBc(BfZone* p, HelmholtzSolver* s,const int icg) : InletHBc(p,s,icg), kappa(0.41) {
    u_bc = NULL; 
    up_it = NULL; 
    Rd = NULL; 
    Rod = NULL; 
  }

  ~XieInletHBc() {
    DELETE(up_it);
    DELETE(Rd);
    DELETE(Rod);

    for (int i = 0; i<NLY2P1_u; ++i){
      DELETE(byz_u[i]); 
    }
    DELETE(byz_u); 
    for (int i = 0; i<NLY2P1_v; ++i){
      DELETE(byz_v[i]); 
    }
    DELETE(byz_v); 
    for (int i = 0; i<NLY2P1_w; ++i){
      DELETE(byz_w[i]); 
    }
    DELETE(byz_w); 
    for (int i = 0; i<NLYW2P1_u; ++i){
      DELETE(byzw_u[i]);
    }
    DELETE(byzw_u);
    for (int i = 0; i<NLYW2P1_v; ++i){
      DELETE(byzw_v[i]);
    }
    DELETE(byzw_v);
    for (int i = 0; i<NLYW2P1_w; ++i){
      DELETE(byzw_w[i]);
    }
    DELETE(byzw_w);

    DELETE(u1f);
  }

  void initData() {
    initializeInflowTurbStructures();
    assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
    assert(mf   == NULL);  mf   = new double[zone_ptr->nbf];

    //pertubation stuff
    assert(up_it == NULL); up_it   = new double[zone_ptr->nbf][3];
    assert(Rd == NULL); Rd   = new double[zone_ptr->nbf][3];
    assert(Rod == NULL); Rod   = new double[zone_ptr->nbf][3];
  
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      mf[ibf] = 0.0;
      FOR_I3 {
        u_bc[ibf][i] = 0.0;
        up_it[ibf][i] = 0.0;
        Rd[ibf][i] = 0.0;
        Rod[ibf][i] = 0.0;
      }
    }

    allocateInflowTurbStructures();
  }
  
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {
    y_min = 1;
    y_max = 1;
    z_min = 1;
    z_max = 1;
    LY = 1; //LY = 4.0;
    LZ = 1; //LZ = 3.0;
    up_it = NULL;
    Rd = NULL;
    Rod = NULL;
    NY = 1;
    NZ = 1;
    dy=LY/NY;    // grid size in the y-dir of the virtual uniform mesh 
    dz=LZ/NZ;    // grid size in the z-dir of the virtual uniform mesh 

    // initialize 9 length scales
    LX_u = 0.0;
    LX_v = 0.0;
    LX_w = 0.0;

    // vertical filter legth-scale (yLu = 0.2*xLu  windTunnel slides)
    LY_u = 0.0;
    LY_v = 0.0;
    LY_w = 0.0;
    NLY_u = 0;
    NLY_v = NLY_u;
    NLY_w = NLY_u;

    NLYW_u = NLY_u;
    NLYW_v = NLY_u;
    NLYW_w = NLY_u;

    // spanwise filter legth-scale (zLu = 0.3*xLu  windTunnel slides)
    LZ_u = 0.0;
    LZ_v = 0.0;
    LZ_w = 0.0;
    NLZ_u = 0;
    NLZ_v = NLZ_u;
    NLZ_w = NLZ_u;

    NLZW_u = NLZ_u;
    NLZW_v = NLZ_u;
    NLZW_w = NLZ_u;

    //N=2n and summation over -N ~ N 
    NLY2P1_u=0;     //filter size in y-dir for u NLY_u 
    NLZ2P1_u=0;     //filter size in z-dir for u NLZ_u 
    NLY2P1_v=0;     //filter size in y-dir for v NLY_v 
    NLZ2P1_v=0;     //filter size in z-dir for v NLZ_v 
    NLY2P1_w=0;     //filter size in y-dir for w NLY_w 
    NLZ2P1_w=0;     //filter size in z-dir for w NLZ_w 

    NLYW2P1_u=0;   //filter size in y-dir near the wall NLYW_u 
    NLZW2P1_u=0;   //filter size in z-dir near the wall NLZW_u
    NLYW2P1_v=0;   //filter size in y-dir near the wall NLYW_v
    NLZW2P1_v=0;   //filter size in z-dir near the wall NLZW_v 
    NLYW2P1_w=0;   //filter size in y-dir near the wall NLYW_w 
    NLZW2P1_w=0;   //filter size in z-dir near the wall NLZW_w 
  
    //filter coeffs
    byz_u = NULL; 
    byz_v = NULL;
    byz_w = NULL; 

    //filter coeffs near the wall
    byzw_u = NULL; 
    byzw_v = NULL;
    byzw_w = NULL; 
  
    u1f = NULL;

    assert(mf   == NULL);  mf   = new double[cg->ncb];
    for (int ibf = 0; ibf < cg->ncb; ++ibf){
      mf[ibf]     = 0.0;
    }
  }
  
  void parseParam(){ 
    Param* param = getParam(getName());
    int iarg = 1;
  
    XP = 0.55;
 
    NY  = 100;
    NZ = 100;

    refD = 1.0;  
    wallD = -0.1;
    Uref = 10.0;

//    lagT_u = 0.12;
//    lagT_v = 0.033;
//    lagT_w = 0.044; 
   
    LX_u = 1.0;
    LY_u = -1.0;
    LZ_u = -1.0;
   
    LX_v = -1.0;
    LY_v = -1.0;
    LZ_v = -1.0;

    LX_w = -1.0;
    LY_w = -1.0;
    LZ_w = -1.0;


    while ( iarg < param->size()) {
      string token = param->getString(iarg++);

      if ( token == "XP") {
        XP = param->getDouble(iarg++);
      } else if ( token == "NY" ) {
        NY = param->getInt(iarg++);
      } else if ( token == "NZ") {
        NZ = param->getInt(iarg++);
      } else if ( token == "UREF"){
        Uref = param->getDouble(iarg++);
      } else if ( token == "REFD"){
        refD = param->getDouble(iarg++);
      } else if ( token == "WALLD"){
        wallD = param->getDouble(iarg++);
      } else if ( token == "xLu"){
        LX_u = param->getDouble(iarg++);
      } else if ( token == "yLu"){
        LY_u = param->getDouble(iarg++);
      } else if ( token == "zLu"){
        LZ_u = param->getDouble(iarg++);
      } else if ( token == "xLv"){
        LX_v = param->getDouble(iarg++);
      } else if ( token == "yLv"){
        LY_v = param->getDouble(iarg++);
      } else if ( token == "zLv"){
        LZ_v = param->getDouble(iarg++);
      } else if ( token == "xLw"){
        LX_w = param->getDouble(iarg++);
      } else if ( token == "yLw"){
        LY_w = param->getDouble(iarg++);
      } else if ( token == "zLw"){
        LZ_w = param->getDouble(iarg++);
      } else{ 
        CERR("Unrecognized " << getName() << " HOOK parameter: " << token);
      }
    }
    if(LY_u < 0.0){LY_u = 0.2 * LX_u; }
    if(LZ_u < 0.0){LZ_u = 0.3 * LX_u; }
   
    if(LX_v < 0.0){LX_v = 0.2 * LX_u; }
    if(LY_v < 0.0){LY_v = LY_u; }
    if(LZ_v < 0.0){LZ_v = LZ_u; }

    if(LX_w < 0.0){LX_w = 0.3 * LX_u; }
    if(LY_w < 0.0){LY_w = LY_u; }
    if(LZ_w < 0.0){LZ_w = LZ_u; }

    lagT_u = LX_u / Uref;
    lagT_v = LX_v / Uref;     
    lagT_w = LX_w / Uref;     

    if(wallD < 0.0){wallD = 0.1*refD; }

    COUT1(" > " << getName() << " Inflow Turbulence HOOK: XP " << XP << " NY " << NY << " NZ " << NZ<< " UREF " << Uref << " REF D " << refD << " WALL D " << wallD);
    COUT1(" > " << getName() << " Inflow Turbulence HOOK: LAGT_U " << lagT_u << " LAGT_V " << lagT_v << " LAGT_W " << lagT_w);
    COUT1(" > " << getName() << " Inflow Turbulence HOOK: xLu " << LX_u << " yLu " << LY_u << " zLu " << LZ_u );
    COUT1(" > " << getName() << " Inflow Turbulence HOOK: xLv " << LX_v << " yLv " << LY_v << " zLv " << LZ_v );
    COUT1(" > " << getName() << " Inflow Turbulence HOOK: xLw " << LX_w << " yLw " << LY_w << " zLw " << LZ_w );
  }
  void initializeInflowTurbStructures(){

    parseParam();
    //get inlet bounding box
    double my_bbox[4] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      for (int nob = solver->noobf_i[zone_ptr->ibf_f+ibf]; nob != solver->noobf_i[zone_ptr->ibf_f+ibf+1]; ++nob) {
        const int ino = solver->noobf_v[nob]; assert((ino >= 0)&&(ino < solver->nno));
        if (solver->x_no[ino][1] < my_bbox[0]){
          my_bbox[0] = solver->x_no[ino][1];
        }
        if (-solver->x_no[ino][1] < my_bbox[1]){
          my_bbox[1] = -solver->x_no[ino][1];
        }
        if (solver->x_no[ino][2] < my_bbox[2]){
          my_bbox[2] = solver->x_no[ino][2];
        }
        if (-solver->x_no[ino][2] < my_bbox[3]){
          my_bbox[3] = -solver->x_no[ino][2];
        }
      }
    }
    double bbox[4];
    MPI_Allreduce(my_bbox, bbox, 4, MPI_DOUBLE, MPI_MIN, mpi_comm); 
    y_min = bbox[0];
    y_max = -bbox[1];
    z_min = bbox[2];
    z_max = -bbox[3];
    LY = y_max - y_min; //LY = 4.0;
    LZ = z_max - z_min; //LZ = 3.0;
    COUT1(" > " << getName() << " Inflow Turbulence HOOK: Bounding Box y_min " << y_min << " y_max " << y_max << " LY " << LY << " z_min " << z_min << " z_max " << z_max << " LZ " << LZ);
    assert(LY>0.0 && LZ>0.0);
    up_it = NULL;
    Rd = NULL;
    Rod = NULL;

    // uniform grid resolution
    dy=LY/NY;    // grid size in the y-dir of the virtual uniform mesh 
    dz=LZ/NZ;    // grid size in the z-dir of the virtual uniform mesh 

    // vertical filter legth-scale (yLu = 0.2*xLu  windTunnel slides)
    // LY_u = 0.2*Uref*lagT_u;
    NLY_u = ceil(LY_u/dy);
    NLY_v = ceil(LY_v/dy);
    NLY_w = ceil(LY_w/dy);

    NLYW_u = NLY_u;
    NLYW_v = NLY_v;
    NLYW_w = NLY_w;

    // spanwise filter legth-scale (zLu = 0.3*xLu  windTunnel slides)
    // LZ_u = 0.3*Uref*lagT_u;
    NLZ_u = ceil(LZ_u/dz);
    NLZ_v = ceil(LZ_v/dz);
    NLZ_w = ceil(LZ_w/dz);

    NLZW_u = NLZ_u;
    NLZW_v = NLZ_v;
    NLZW_w = NLZ_w;

    //N=2n and summation over -N ~ N 
    NLY2P1_u=1+ NLY_u*4;     //filter size in y-dir for u NLY_u 
    NLZ2P1_u=1+ NLZ_u*4;     //filter size in z-dir for u NLZ_u 
    NLY2P1_v=1+ NLY_v*4;     //filter size in y-dir for v NLY_v 
    NLZ2P1_v=1+ NLZ_v*4;     //filter size in z-dir for v NLZ_v 
    NLY2P1_w=1+ NLY_w*4;     //filter size in y-dir for w NLY_w 
    NLZ2P1_w=1+ NLZ_w*4;     //filter size in z-dir for w NLZ_w 

    NLYW2P1_u=1+ NLYW_u*4;   //filter size in y-dir near the wall NLYW_u 
    NLZW2P1_u=1+ NLZW_u*4;   //filter size in z-dir near the wall NLZW_u
    NLYW2P1_v=1+ NLYW_v*4;   //filter size in y-dir near the wall NLYW_v
    NLZW2P1_v=1+ NLZW_v*4;   //filter size in z-dir near the wall NLZW_v 
    NLYW2P1_w=1+ NLYW_w*4;   //filter size in y-dir near the wall NLYW_w 
    NLZW2P1_w=1+ NLZW_w*4;   //filter size in z-dir near the wall NLZW_w 
  
    //filter coeffs
    byz_u = NULL; 
    byz_v = NULL;
    byz_w = NULL; 

    //filter coeffs near the wall
    byzw_u = NULL; 
    byzw_v = NULL;
    byzw_w = NULL; 
  
    u1f = NULL;
  }

  void allocateInflowTurbStructures(){
    //filter coeffs
    assert(byz_u==NULL); byz_u = new double*[NLY2P1_u]; 
    for (int i = 0; i<NLY2P1_u; ++i){
      byz_u[i] = new double[NLZ2P1_u]; 
      for (int j = 0; j<NLZ2P1_u; ++j){
        byz_u[i][j] = 0.0;
      }
    }
    assert(byz_v==NULL); byz_v = new double*[NLY2P1_v]; 
    for (int i = 0; i<NLY2P1_v; ++i){
      byz_v[i] = new double[NLZ2P1_v]; 
      for (int j = 0; j<NLZ2P1_v; ++j){
        byz_v[i][j] = 0.0;
      }
    }
    assert(byz_w==NULL); byz_w = new double*[NLY2P1_w]; 
    for (int i = 0; i<NLY2P1_w; ++i){
      byz_w[i] = new double[NLZ2P1_w]; 
      for (int j = 0; j<NLZ2P1_w; ++j){
        byz_w[i][j] = 0.0;
      }
    }

    //filter coeffs near the wall
    assert(byzw_u==NULL); byzw_u = new double*[NLYW2P1_u]; 
    for (int i = 0; i<NLYW2P1_u; ++i){
      byzw_u[i] = new double[NLZW2P1_u];
      for (int j = 0; j<NLZW2P1_u; ++j){
        byzw_u[i][j] = 0.0;
      }
    }
    assert(byzw_v==NULL); byzw_v = new double*[NLYW2P1_v];
    for (int i = 0; i<NLYW2P1_v; ++i){
      byzw_v[i] = new double[NLZW2P1_v];
      for (int j = 0; j<NLZW2P1_v; ++j){
        byzw_v[i][j] = 0.0;
      }
    }
    assert(byzw_w==NULL); byzw_w = new double*[NLYW2P1_w]; 
    for (int i = 0; i<NLYW2P1_w; ++i){
      byzw_w[i] = new double[NLZW2P1_w];
      for (int j = 0; j<NLZW2P1_w; ++j){
        byzw_w[i][j] = 0.0;
      }
    }

    assert(u1f==NULL); u1f = new double[NY*NZ][3];
    for(int k=0; k < NZ; k++){
      for(int j=0;j<NY;j++){
        int jk = k*NY + j;
        FOR_I3 u1f[jk][i] = 0.0;
      }
    }
  }


  void initialHook(){

    if (mpi_rank==0){
      cout << "XieInlet hook is calling initialHook()..." << endl;
    }

    oneDProfile * umProf  = new oneDProfile;
    umProf->read("inflow_data/UInlet");
    oneDProfile * uuProf  = new oneDProfile;
    uuProf->read("inflow_data/uuBarInlet");
    oneDProfile * vvProf  = new oneDProfile;
    vvProf->read("inflow_data/vvBarInlet");
    oneDProfile * wwProf  = new oneDProfile;
    wwProf->read("inflow_data/wwBarInlet");
    oneDProfile * uvProf  = new oneDProfile;
    uvProf->read("inflow_data/uvBarInlet");

    const double rho_bc = getDoubleParam("RHO"); 
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      u_bc[ibf][0]  = umProf->getData(zone_ptr->x_bf[ibf][1]);
      u_bc[ibf][1]  = 0.0;
      u_bc[ibf][2]  = 0.0;
      // Apply constant mass flow at inlet
      mf[ibf] = rho_bc*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
   
      Rd[ibf][0]  = uuProf->getData(zone_ptr->x_bf[ibf][1]);
      Rd[ibf][1]  = vvProf->getData(zone_ptr->x_bf[ibf][1]);
      Rd[ibf][2]  = wwProf->getData(zone_ptr->x_bf[ibf][1]);
      Rod[ibf][0] = 0.0;
      Rod[ibf][1] = 0.0;
      Rod[ibf][2] = uvProf->getData(zone_ptr->x_bf[ibf][1]);
    }

    delete umProf;
    delete uuProf;
    delete vvProf;
    delete wwProf;
    delete uvProf;

    // Look for an ascii file with a previous flow field on the inflow
    // turbulence grid
  
    int it_error = 0;
    if ( mpi_rank == 0){// && solver->step>0)
      stringstream filename;
      filename << "u1f." << setfill('0') << setw(8) << solver->step << ".dat";
      cout << "reading u1f inflow turbulence velocity history from file: " << filename.str().c_str() << endl;
      // first read the file to get the size
      ifstream fp;
      fp.open(filename.str().c_str());
      if ( fp.fail()){
        cout << "Warning: could not find file " << filename.str().c_str() << ", resetting inflow turbulence time history"  << endl;
      }
      else{
        string line;
        int linec = -1;
        while(!fp.eof()){
          getline(fp,line);
          stringstream ss(line);
          ++linec;
          if (linec<NY*NZ)
            ss >> u1f[linec][0] >> u1f[linec][1] >> u1f[linec][2];
        }
        fp.close();
        if (linec!=NY*NZ){
          it_error = 1;
          cout << "Error: found " << linec << " entries in " <<  filename.str().c_str() << " but was expecting " << NY*NZ << endl;
        }
      }
      cout << "  first u1f[0000] = " << u1f[0][0] << " " << u1f[0][1] << " " << u1f[0][2] << endl;
      cout << "  last  u1f[" << NY*NZ-1 << "] = " << u1f[NY*NZ-1][0] << " " << u1f[NY*NZ-1][1] << " " << u1f[NY*NZ-1][2] << endl;
    }
    MPI_Bcast(&it_error,1,MPI_INT,0,mpi_comm);
    if (it_error==1)
      throw(0);

    //send velocity history to all ranks from rank 0
    MPI_Bcast(&(u1f[0][0]),NY*NZ*3,MPI_DOUBLE,0,mpi_comm); 
    // get coefficients, by and bz : these are generated only once whole time

    double sumy_u=0.0;
    double sumy_v=0.0;
    double sumy_w=0.0;

    double sumz_u=0.0;
    double sumz_v=0.0;
    double sumz_w=0.0;

    double sumyw_u=0.0;
    double sumyw_v=0.0;
    double sumyw_w=0.0;

    double sumzw_u=0.0;
    double sumzw_v=0.0;
    double sumzw_w=0.0;

    //sumy_i
    for (int j=0; j<NLY2P1_u;j++){
      sumy_u += exp(-0.5*M_PI*fabs(float(j-2*NLY_u)/float(NLY_u))) * exp(-0.5*M_PI*fabs(float(j-2*NLY_u)/float(NLY_u)) );
    }
    for (int j=0; j<NLY2P1_v;j++){
      sumy_v += exp(-0.5*M_PI*fabs(float(j-2*NLY_v)/float(NLY_v))) * exp(-0.5*M_PI*fabs(float(j-2*NLY_v)/float(NLY_v)) );
    }
    for (int j=0; j<NLY2P1_w;j++){
      sumy_w += exp(-0.5*M_PI*fabs(float(j-2*NLY_w)/float(NLY_w))) * exp(-0.5*M_PI*fabs(float(j-2*NLY_w)/float(NLY_w)) );
    }

    //sumz_i
    for (int k=0; k<NLZ2P1_u;k++){
      sumz_u += exp(-0.5*M_PI*fabs(float(k-2*NLZ_u)/float(NLZ_u))) * exp(-0.5*M_PI*fabs(float(k-2*NLZ_u)/float(NLZ_u))); 
    }
    for (int k=0; k<NLZ2P1_v;k++){
      sumz_v += exp(-0.5*M_PI*fabs(float(k-2*NLZ_v)/float(NLZ_v))) * exp(-0.5*M_PI*fabs(float(k-2*NLZ_v)/float(NLZ_v))); 
    }
    for (int k=0; k<NLZ2P1_w;k++){
      sumz_w += exp(-0.5*M_PI*fabs(float(k-2*NLZ_w)/float(NLZ_w))) * exp(-0.5*M_PI*fabs(float(k-2*NLZ_w)/float(NLZ_w))); 
    }

    //sumyw_i
    for (int j=0; j<NLYW2P1_u;j++){
      sumyw_u +=exp(-0.5*M_PI*fabs(float(j-2*NLYW_u)/float(NLYW_u))) * exp(-0.5*M_PI*fabs(float(j-2*NLYW_u)/float(NLYW_u)) );
    }
    for (int j=0; j<NLYW2P1_v;j++){
      sumyw_v +=exp(-0.5*M_PI*fabs(float(j-2*NLYW_v)/float(NLYW_v))) * exp(-0.5*M_PI*fabs(float(j-2*NLYW_v)/float(NLYW_v)) );
    }
    for (int j=0; j<NLYW2P1_w;j++){
      sumyw_w += exp(-0.5*M_PI*fabs(float(j-2*NLYW_w)/float(NLYW_w))) * exp(-0.5*M_PI*fabs(float(j-2*NLYW_w)/float(NLYW_w)) );
    }

    //sumzw_i
    for (int k=0; k<NLZW2P1_u;k++){
      sumzw_u += exp(-0.5*M_PI*fabs(float(k-2*NLZW_u)/float(NLZW_u))) * exp(-0.5*M_PI*fabs(float(k-2*NLZW_u)/float(NLZW_u)) );
    }
    for (int k=0; k<NLZW2P1_v;k++){
      sumzw_v += exp(-0.5*M_PI*fabs(float(k-2*NLZW_v)/float(NLZW_v))) * exp(-0.5*M_PI*fabs(float(k-2*NLZW_v)/float(NLZW_v)) );
    }
    for (int k=0; k<NLZW2P1_w;k++){
      sumzw_w += exp(-0.5*M_PI*fabs(float(k-2*NLZW_w)/float(NLZW_w))) * exp(-0.5*M_PI*fabs(float(k-2*NLZW_w)/float(NLZW_w)) );
    }

    
    sumy_u=sqrt(sumy_u);
    sumy_v=sqrt(sumy_v);
    sumy_w=sqrt(sumy_w);

    sumz_u=sqrt(sumz_u);
    sumz_v=sqrt(sumz_v);
    sumz_w=sqrt(sumz_w);

    sumyw_u=sqrt(sumyw_u);
    sumyw_v=sqrt(sumyw_v);
    sumyw_w=sqrt(sumyw_w);

    sumzw_u=sqrt(sumzw_u);
    sumzw_v=sqrt(sumzw_v);
    sumzw_w=sqrt(sumzw_w);

    // by_i
    double * by_u = new double[NLY2P1_u];
    for (int j=0;j<NLY2P1_u;j++){
      by_u[j] = exp(-0.5*M_PI*fabs(float(j-2*NLY_u)/float(NLY_u)))/sumy_u;
    }

    double * by_v = new double[NLY2P1_v];
    for (int j=0;j<NLY2P1_v;j++){
      by_v[j] = exp(-0.5*M_PI*fabs(float(j-2*NLY_v)/float(NLY_v)))/sumy_v;
    }
    double * by_w = new double[NLY2P1_w];
    for (int j=0;j<NLY2P1_w;j++){
      by_w[j] = exp(-0.5*M_PI*fabs(float(j-2*NLY_w)/float(NLY_w)))/sumy_w;
    }

    //bz_i
    double * bz_u = new double[NLZ2P1_u];
    for (int k=0; k<NLZ2P1_u;k++){
      bz_u[k]=exp(-0.5*M_PI*fabs(float(k-2*NLZ_u)/float(NLZ_u)))/sumz_u;  
    }
    double * bz_v = new double[NLZ2P1_v];
    for (int k=0; k<NLZ2P1_v;k++){
      bz_v[k]= exp(-0.5*M_PI*fabs(float(k-2*NLZ_v)/float(NLZ_v)))/sumz_v;  
    }
    double * bz_w = new double[NLZ2P1_w];
    for (int k=0; k<NLZ2P1_w;k++){
      bz_w[k]= exp(-0.5*M_PI*fabs(float(k-2*NLZ_w)/float(NLZ_w)))/sumz_w;  
    }

    //byz_i
    for (int k=0; k<NLZ2P1_u;k++){
      for (int j=0;j<NLY2P1_u;j++){
        byz_u[j][k]=by_u[j]*bz_u[k];
      }
    }
    delete[] by_u;
    delete[] bz_u;
    for (int k=0; k<NLZ2P1_v;k++){
      for (int j=0;j<NLY2P1_v;j++){
        byz_v[j][k]=by_v[j]*bz_v[k];
      }
    }
    delete[] by_v;
    delete[] bz_v;
    for (int k=0; k<NLZ2P1_w;k++){
      for (int j=0;j<NLY2P1_w;j++){
        byz_w[j][k]=by_w[j]*bz_w[k];
      }
    }
    delete[] by_w;
    delete[] bz_w;

    //byw_i
    double * byw_u = new double[NLYW2P1_u];
    for (int j=0;j<NLYW2P1_u;j++){
      byw_u[j] = exp(-0.5*M_PI*fabs(float(j-2*NLYW_u)/float(NLYW_u)))/sumyw_u;
    }
    double * byw_v = new double[NLYW2P1_v];
    for (int j=0;j<NLYW2P1_v;j++){
      byw_v[j] = exp(-0.5*M_PI*fabs(float(j-2*NLYW_v)/float(NLYW_v)))/sumyw_v;
    }
    double * byw_w = new double[NLYW2P1_w];
    for (int j=0;j<NLYW2P1_w;j++){
      byw_w[j] = exp(-0.5*M_PI*fabs(float(j-2*NLYW_w)/float(NLYW_w)))/sumyw_w;
    }

    //bzw_i
    double * bzw_u = new double[NLZW2P1_u];
    for (int k=0;k<NLZW2P1_u;k++){
      bzw_u[k] = exp(-0.5*M_PI*fabs(float(k-2*NLZW_u)/float(NLZW_u)))/sumzw_u;
    }
    double * bzw_v = new double[NLZW2P1_v];
    for (int k=0;k<NLZW2P1_v;k++){
      bzw_v[k] = exp(-0.5*M_PI*fabs(float(k-2*NLZW_v)/float(NLZW_v)))/sumzw_v;
    }
    double * bzw_w = new double[NLZW2P1_w];
    for (int k=0;k<NLZW2P1_w;k++){
      bzw_w[k] = exp(-0.5*M_PI*fabs(float(k-2*NLZW_w)/float(NLZW_w)))/sumzw_w;
    }


    //byzw_i
    for (int k=0; k<NLZW2P1_u;k++){
      for (int j=0;j<NLYW2P1_u;j++){
        byzw_u[j][k]=byw_u[j]*bzw_u[k];
      }
    }
    delete[] byw_u;
    delete[] bzw_u;
    for (int k=0; k<NLZW2P1_v;k++){
      for (int j=0;j<NLYW2P1_v;j++){
        byzw_v[j][k]=byw_v[j]*bzw_v[k];
      }
    }
    delete[] byw_v;
    delete[] bzw_v;
    for (int k=0; k<NLZW2P1_w;k++){
      for (int j=0;j<NLYW2P1_w;j++){
        byzw_w[j][k]=byw_w[j]*bzw_w[k];
      }
    }  
    delete[] byw_w;
    delete[] bzw_w;
  
    // filter coefficient END
    
  }

  void setBc() {

    if (mpi_rank==0&&solver->step%solver->check_interval==0){
      cout << "XieInlet hook is calling setBc()..." << endl;
    }

    // read time step
    double timeStep= solver->dt;
    double ran,uni; 

    // streamwise filter legth-scale
    double NLX_u = lagT_u/timeStep;
    double NLX_v = lagT_v/timeStep;
    double NLX_w = lagT_w/timeStep;
  
    int NRY=3*NY;     //random array size in the y-dir NY 
    int NRZ=3*NZ;     //random array size in the z-dir NZ 
  
    double  cw_ =M_PI/4.0;

    // random numbers  
    double * rndmx1 = new double[NRY*NRZ];
    double * rndmy1 = new double[NRY*NRZ];
    double * rndmz1 = new double[NRY*NRZ];

    // correaltion functions  
    double ** pSi_u = new double*[NY];
    for (int i=0;i<NY;++i)
      pSi_u[i] = new double[NZ];
    double ** pSi_v = new double*[NY];
    for (int i=0;i<NY;++i)
      pSi_v[i] = new double[NZ];
    double ** pSi_w = new double*[NY];
    for (int i=0;i<NY;++i)
      pSi_w[i] = new double[NZ];

    // compute random numbers on rank 0 and broadcast to all ranks
    if (mpi_rank==0){
      for(int k=0; k < NRZ; k++){
        for(int j=0;j<NRY;j++){
          const int jk = k*NRY + j;  
          rndmx1[jk] = 0.0;
          rndmy1[jk] = 0.0;
          rndmz1[jk] = 0.0;

          ran = 0.0;
          for(int count=1; count<=12; ++count){
            uni = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
            ran += uni;
          }
          rndmx1[jk]=ran -6.0;

          ran = 0.0;
          for(int count=1; count<=12; ++count){
            uni = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
            ran += uni;
          }
          rndmy1[jk]=ran -6.0;

          ran = 0.0;
          for(int count=1; count<=12; ++count){
            uni = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
            ran += uni;
          }
          rndmz1[jk]=ran -6.0;
        }
      }
    }
    MPI_Bcast(rndmx1,NRY*NRZ,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast(rndmy1,NRY*NRZ,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast(rndmz1,NRY*NRZ,MPI_DOUBLE,0,mpi_comm);
    // Random Field generation End

    // cout<<"\nSpanwise length-scales 12  \n"<< "Rsuu "<< 12345 <<"\n"<< endl;
    
    // Spatial correlation (y, z-dir) implementation
    for(int kk = 0; kk < NZ; kk++){
      for(int jj = 0; jj < NY; jj++){ 
  
        pSi_u[jj][kk] =0.0;
        pSi_v[jj][kk] =0.0;
        pSi_w[jj][kk] =0.0;

        if(jj*dy < wallD ){ //near wall region
          for(int k = 0 ; k < NLZW2P1_u; k++){
            for (int j = 0; j < NLYW2P1_u; j++) {
              const int jkrw = (k+kk+NZ-2*NLZW_u)*NRY+(j+jj+NY-2*NLYW_u);
              pSi_u[jj][kk] += byzw_u[j][k]*rndmx1[jkrw];
            }
          }
          for(int k = 0 ; k < NLZW2P1_v; k++){
            for (int j = 0; j < NLYW2P1_v; j++) {
              const int jkrw = (k+kk+NZ-2*NLZW_v)*NRY+(j+jj+NY-2*NLYW_v);
              pSi_v[jj][kk] += byzw_v[j][k]*rndmy1[jkrw];
            }
          }
          for(int k = 0 ; k < NLZW2P1_w; k++){
            for (int j = 0; j < NLYW2P1_w; j++) {
              const int jkrw = (k+kk+NZ-2*NLZW_w)*NRY+(j+jj+NY-2*NLYW_w);
              pSi_w[jj][kk] += byzw_w[j][k]*rndmz1[jkrw];
            }
          }
        }
        else{
          for(int k = 0 ; k < NLZ2P1_u; k++){
            for (int j = 0; j < NLY2P1_u; j++) {
              const int jkr = (k+kk+NZ-2*NLZ_u)*NRY+(j+jj+NY-2*NLY_u);
              pSi_u[jj][kk] += byz_u[j][k]*rndmx1[jkr];
            }
          }
          for(int k = 0 ; k < NLZ2P1_v; k++){
            for (int j = 0; j < NLY2P1_v; j++) {
              const int jkr = (k+kk+NZ-2*NLZ_v)*NRY+(j+jj+NY-2*NLY_v);
              pSi_v[jj][kk] += byz_v[j][k]*rndmy1[jkr];
            }
          }
          for(int k = 0 ; k < NLZ2P1_w; k++){
            for (int j = 0; j < NLY2P1_w; j++) {
              const int jkr = (k+kk+NZ-2*NLZ_w)*NRY+(j+jj+NY-2*NLY_w);
              pSi_w[jj][kk] += byz_w[j][k]*rndmz1[jkr];
            }
          }
        }
      }
    }
    delete[] rndmx1;
    delete[] rndmy1;
    delete[] rndmz1;

    //cout<<"\nSpanwise length-scales 12 Rsuu \n"<< "RSuu  = "<< 1234 <<"\n"<< endl;

    double (*u2f)[3] = new double[NY*NZ][3];

    // correlation (x-dir) implementation
    for(int k=0; k < NZ; k++){
      for(int j=0;j<NY;j++){
        int jk = k*NY + j;
        u2f[jk][0] = u1f[jk][0]* exp(-cw_/NLX_u ) + pSi_u[j][k]*sqrt(1.0-exp(-2.0*cw_/NLX_u ) );  
        u2f[jk][1] = u1f[jk][1]* exp(-cw_/NLX_v ) + pSi_v[j][k]*sqrt(1.0-exp(-2.0*cw_/NLX_v ) ); 
        u2f[jk][2] = u1f[jk][2]* exp(-cw_/NLX_w ) + pSi_w[j][k]*sqrt(1.0-exp(-2.0*cw_/NLX_w ) );
      }
    }
    // correlation (x-dir) implementation END
   
    for (int i=0;i<NY;++i)
      delete[] pSi_u[i];
    delete[] pSi_u;
    for (int i=0;i<NY;++i)
      delete[] pSi_v[i];
    delete[] pSi_v;
    for (int i=0;i<NY;++i)
      delete[] pSi_w[i];
    delete[] pSi_w;

    // cyclic boundary manipulation (not perfect..)
    ///************************************8
  
    for(int j=0;j<NY;j++){
      const int jk0 = 0*NY + j;
      const int jk1 = 1*NY + j;
      const int jkN0 = (NZ-1)*NY + j;
      const int jkN1 = (NZ-2)*NY + j;
      FOR_I3  u2f[jk0][i] = u2f[jkN0][i];       
      FOR_I3  u2f[jk1][i] = u2f[jkN1][i];
    }
  
    // Mapping from 'UNIFORM' to 'NON-UNIFORM' mesh
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int j = min(max(int(NY*(zone_ptr->x_bf[ibf][1]-y_min)/(y_max-y_min)),0),NY-1);
      const int k = min(max(int(NZ*(zone_ptr->x_bf[ibf][2]-z_min)/(z_max-z_min)),0),NZ-1);
      int jk=k*NY + j;

      assert(Rd[ibf][0]>=0.0);
      double a11 =  sqrt(Rd[ibf][0]);
      double a21 = 0.0;
      if ( a11 > 0.0 ) a21 = Rod[ibf][2] / a11;
      double a22 = sqrt( max(Rd[ibf][1]-a21*a21,0.0) );
      double a31 = 0.0;
      if ( a11> 0.0 ) a31 = Rod[ibf][1] / a11;
      double a32 = 0.0;
      if ( a22 > 0.0 ) a32 =  (Rod[ibf][0] - a21*a31) / a22;
      double a33 = sqrt(  max(Rd[ibf][2]-a31*a31-a32*a32,0.0) );

      //set cell fluctuating velocity
      up_it[ibf][0] = a11*u2f[jk][0];
      up_it[ibf][1] = a21*u2f[jk][0] + a22*u2f[jk][1];
      up_it[ibf][2] = a31*u2f[jk][0] + a32*u2f[jk][1] + a33*u2f[jk][2];
    }

    for(int k=0; k < NZ; k++){
      for(int j=0;j<NY;j++){
        int jk = k*NY + j;
        FOR_I3   u1f[jk][i] = u2f[jk][i];
      }
    }

    delete[] u2f;


    //if the user has requested writing result files at a given interval,
    //also write the fluctuating velocity field needed for restarting.
    FOR_PARAM_MATCHING("WRITE_RESULT") {
      int interval       = -1;     
      if ( param->size() == 0 ) {
        continue;
      } else if ( param->size() == 1) {
        interval = param->getInt(0);
      } else {
        int iarg = 0;
        while ( iarg < param->size()) {
          string token = param->getString(iarg++);
          if ( token == "NAME") {
            iarg++;
            //don't use the name here
          } else if (token == "INTERVAL") {
            interval = param->getInt(iarg++);
            break;
          } else {
            //unrecognized syntax
            break;
          }
        }
      }
      if ( interval > 0 && solver->step%interval==0 ) {
        writeVelocityHistory();
      }
    }
  }

  //modify the predicted velocity field from solveU with the turbulent field
  void modifyU(){
    if (solver->step%solver->check_interval==0){
      COUT1(" > " << getName() << ":modifyU()");
      dumpRange(up_it,zone_ptr->nbf,"up_it");
    }
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      FOR_I3 solver->u[icv0][i] = u_bc[ibf][i] + up_it[ibf][i];
    }
  }

  void writeVelocityHistory(){
    if (mpi_rank==0){
      stringstream filename;
      filename << "u1f." << setfill('0') << setw(8) << solver->step << ".dat";
      cout << "writing u1f inflow turbulence velocity history to file: " << filename.str().c_str() << " ... ";
      ofstream ofile;
      ofile.open(filename.str().c_str());
      assert(ofile.is_open());
      for (int jk = 0; jk<NY*NZ; ++jk){
        ofile << u1f[jk][0] << " " << u1f[jk][1] << " " << u1f[jk][2] << endl;
      }
      ofile.close();
      cout << "complete." << endl;
    }
  }

};

class XieHelmholtzSolver : public HelmholtzSolver { 
public: 

  XieHelmholtzSolver() : HelmholtzSolver() {}
 
  ~XieHelmholtzSolver() {}

  HelmholtzBc* initHookBc(BfZone* p) {
  
    //
    // assumes that the input file required as bc specification 
    // of the form.. 
    // 
    //   x0 = HOOK 
    // 

    assert( p->getName() == "x0_turb");
    return new XieInletHBc(p,this);

  }

  //Required to support multigrid solver... 
  XieHelmholtzSolver(const int icg) : HelmholtzSolver(icg) {}

  HelmholtzBc* initHookBcMg(BfZone* p, const int icg) {
    assert( p->getName() == "x0_turb");
    return new XieInletHBc(p,this,icg);
  }

  HelmholtzSolver* getHelmholtzSolverMg(const int icg){
    return (new XieHelmholtzSolver(icg));
  }
  //End Required to support multigrid solver

  // Add calls to initialHook, temporalHook, etc...

  void finalHook() {
  
    for(vector<HelmholtzBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it) {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "HOOK") {
        XieInletHBc * xieInlet = dynamic_cast<XieInletHBc*>(*it);
        xieInlet->writeVelocityHistory();
      }
    }
  
  }

};



#endif
