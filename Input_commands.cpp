#include"Input_commands.h"
Input::Input(){cout<<"Not allowed default constructor"<<endl;}
Input::Input(string rho_in)
{
 punctualr=false;punctualp=false;scanr=false;scanp=false;scanelf=false;cps=false;indicators=false;
 integrals=false;log=false;spin_calcs=false;cas=false;cuba=false;debug=false;
 int_file=false;dmn=false;dmnp=false;dmn_integrals=false;scanindic=false;
 dmn_thresh=false;esi_int=false;dmn_indicators=false;nopath=true;dim2=false;dim3=false;quadrature=false;gnuplot=false;
 cubature=false;dmn_plots=false;rotate_grid=false;Beta_MOs=false;wfx_print=false;wfx_print_dmn=false;store_dmn=false;
 print_dm1_fchk=false;cubature2=false;cube=false;tps=false;intracule=false;Vr=false;scan_localhybs=false;dens_sim=false;
 int_pol_hyperpol=false;extracule=false;r1_moment=false;symrotdens=false;symgrad=false;intra_1rdm=false;
 multiplicity=0;ncores=1;extra_lines=0;
 string name=rho_in;
 ifstream rho_input_file;
 rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
 rho_input_file.open((rho_in).c_str());
 if(rho_input_file.good()) //Check existence of file
 {
  while(getline(rho_input_file,rho_in))
  {
   rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
   lowercase(rho_in);
   if(rho_in=="$name")
   {
    getline(rho_input_file,rho_in);
    rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    name_fchk_wfn=rho_in;
   }
   else if(rho_in=="$log")
   {
    log=true;
    getline(rho_input_file,rho_in);
    name_log=rho_in;
   }
   else if(rho_in=="$multiplicity"){rho_input_file>>multiplicity;}
   else if(rho_in=="$spin_calcs"){spin_calcs=true;}
   else if(rho_in=="$cas"){cas=true;}
   else if(rho_in=="$intra_1rdm")
   {
    intra_1rdm=true;
    rho_input_file>>order_grid_r>>order_grid_ang;
   }
   else if(rho_in=="$punctual_r")
   {
    punctualr=true;
    rho_input_file>>punctuals_r;
    coordinates_r=new double*[punctuals_r];
    for(i=0;i<punctuals_r;i++)
    {coordinates_r[i]=new double[3];}
    for(i=0;i<punctuals_r;i++)
    {rho_input_file>>coordinates_r[i][0]>>coordinates_r[i][1]>>coordinates_r[i][2];}
   }
   else if(rho_in=="$punctual_p")
   {
    punctualp=true;
    rho_input_file>>punctuals_p;
    coordinates_p=new double*[punctuals_p];
    for(i=0;i<punctuals_p;i++)
    {coordinates_p[i]=new double[3];}
    for(i=0;i<punctuals_p;i++)
    {rho_input_file>>coordinates_p[i][0]>>coordinates_p[i][1]>>coordinates_p[i][2];}
   }
   else if(rho_in=="$scan_r")
   {
    scanr=true;
    rho_input_file>>init_coord_r[0]>>init_coord_r[1]>>init_coord_r[2];
    rho_input_file>>points_scan_r;
    rho_input_file>>step_r;
    rho_input_file>>dir_r;
   }
   else if(rho_in=="$scan_p")
   {
    scanp=true;
    rho_input_file>>init_coord_p[0]>>init_coord_p[1]>>init_coord_p[2];
    rho_input_file>>points_scan_p;
    rho_input_file>>step_p;
    rho_input_file>>dir_p;
   }
   else if(rho_in=="$scan_elf")
   {
    scanelf=true;
    rho_input_file>>init_coord_elf[0]>>init_coord_elf[1]>>init_coord_elf[2];
    rho_input_file>>points_scan_elf;
    rho_input_file>>step_elf;
    rho_input_file>>dir_elf;
   }
   else if(rho_in=="$scan_indicators")
   {
    scanindic=true;
    rho_input_file>>init_coord_indic[0]>>init_coord_indic[1]>>init_coord_indic[2];
    rho_input_file>>points_scan_indic;
    rho_input_file>>step_indic;
    rho_input_file>>dir_indic;
   }
   else if(rho_in=="$scan_lhybrids")
   {
    scan_localhybs=true;
    rho_input_file>>init_coord_lh[0]>>init_coord_lh[1]>>init_coord_lh[2];
    rho_input_file>>points_scan_lh;
    rho_input_file>>step_lh;
    rho_input_file>>dir_lh;
   }
   else if(rho_in=="$vr")
   {
    Vr=true;
    rho_input_file>>Point_Vr[0]>>Point_Vr[1]>>Point_Vr[2];
    rho_input_file>>ncores;
    rho_input_file>>error_abs>>error_rel;
    rho_input_file>>minevals>>maxevals;
    if(ncores<0){ncores=0;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="vegas"){method_cuba="Vegas";}
    else if(rho_in=="suave"){method_cuba="Suave";}
    else if(rho_in=="divonne"){method_cuba="Divonne";}
    else{method_cuba="Cuhre";}
   }
   else if(rho_in=="$cps"){cps=true;}
   else if(rho_in=="$indicators"){indicators=true;}
   else if(rho_in=="$beta_mos"){Beta_MOs=true;}
   else if(rho_in=="$dmnincore"){store_dmn=true;rho_input_file>>mem;}
   else if(rho_in=="$wfx_print"){wfx_print=true;}
   else if(rho_in=="$cube")
   {
    cube=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    opcube=rho_in;
    rho_input_file>>cubex>>cubey>>cubez;
    rho_input_file>>step_cube_x>>step_cube_y>>step_cube_z;
    if(opcube=="mo" || opcube=="no")
    {rho_input_file>>mo_no_cube;}
   }
   else if(rho_in=="$intracule")
   {
    intracule=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm2=rho_in;
    rho_input_file>>Point_intra[0]>>Point_intra[1]>>Point_intra[2];
    rho_input_file>>ncores;
    rho_input_file>>error_abs>>error_rel;
    rho_input_file>>minevals>>maxevals;
    if(ncores<0){ncores=0;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="vegas"){method_cuba="Vegas";}
    else if(rho_in=="suave"){method_cuba="Suave";}
    else if(rho_in=="divonne"){method_cuba="Divonne";}
    else{method_cuba="Cuhre";}
   }
   else if(rho_in=="$extracule")
   {
    extracule=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm2=rho_in;
    rho_input_file>>Point_extra[0]>>Point_extra[1]>>Point_extra[2];
    rho_input_file>>ncores;
    rho_input_file>>error_abs>>error_rel;
    rho_input_file>>minevals>>maxevals;
    if(ncores<0){ncores=0;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="vegas"){method_cuba="Vegas";}
    else if(rho_in=="suave"){method_cuba="Suave";}
    else if(rho_in=="divonne"){method_cuba="Divonne";}
    else{method_cuba="Cuhre";}
   }
   else if(rho_in=="$tps")
   {
    tps=true;
    rho_input_file>>ncores;
    rho_input_file>>error_abs>>error_rel;
    rho_input_file>>minevals>>maxevals;
    if(ncores<0){ncores=0;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="vegas"){method_cuba="Vegas";}
    else if(rho_in=="suave"){method_cuba="Suave";}
    else if(rho_in=="divonne"){method_cuba="Divonne";}
    else{method_cuba="Cuhre";}
    for(i=0;i<6;i++)
    {
     second_moments_tps[i]=ZERO;
     rho_input_file>>second_moments_tps[i];
    }
   }
   else if(rho_in=="$dens_sim")
   {
    dens_sim=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    second_fchk_wfn=rho_in;
    rho_input_file>>ncores;
    rho_input_file>>error_abs>>error_rel;
    rho_input_file>>minevals>>maxevals;
    if(ncores<0){ncores=0;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="vegas"){method_cuba="Vegas";}
    else if(rho_in=="suave"){method_cuba="Suave";}
    else if(rho_in=="divonne"){method_cuba="Divonne";}
    else{method_cuba="Cuhre";}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="density"){symrotdens=true;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="gradient"){symgrad=true;}
   }
   else if(rho_in=="$int_pol_hyperpol")
   {
    int_pol_hyperpol=true; 
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    second_fchk_wfn=rho_in;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    third_fchk_wfn=rho_in;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    fourth_fchk_wfn=rho_in;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    fifth_fchk_wfn=rho_in;
    rho_input_file>>ncores;
    rho_input_file>>error_abs>>error_rel;
    rho_input_file>>minevals>>maxevals;
    if(ncores<0){ncores=0;}
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    if(rho_in=="vegas"){method_cuba="Vegas";}
    else if(rho_in=="suave"){method_cuba="Suave";}
    else if(rho_in=="divonne"){method_cuba="Divonne";}
    else{method_cuba="Cuhre";}
    rho_input_file>>field_F;
    rho_input_file>>dir_pol_hyper;
   }
   else if(rho_in=="$wfx_print_dmn")
   {
    wfx_print_dmn=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm1=rho_in;
   }
   else if(rho_in=="$fchk_dm1")
   {
    print_dm1_fchk=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm1=rho_in;
   }
   else if(rho_in=="$integrals")
   {
    integrals=true;
    getline(rho_input_file,cuba_cubature);
    cuba_cubature.erase(std::remove_if(cuba_cubature.begin(),cuba_cubature.end(),::isspace),cuba_cubature.end());
    lowercase(cuba_cubature);
    if(cuba_cubature=="cuba")
    {
     cuba=true;
     rho_input_file>>ncores;
     rho_input_file>>error_abs>>error_rel;
     rho_input_file>>minevals>>maxevals;
     if(ncores<0){ncores=0;}
     //Number of integral operations ops
     rho_input_file>>ops;
     integral_ops=new string[ops];
     for(i=0;i<ops;i++)
     {
      //Get the operation to compute
      do
      {
       getline(rho_input_file,rho_in);
      }while(rho_in=="");
      rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
      lowercase(rho_in);
      integral_ops[i]=rho_in;
      if(integral_ops[i]=="sij")
      {
       int_file=true;
       do
       {
        getline(rho_input_file,rho_in);
       }while(rho_in=="");
       rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
       lowercase(rho_in);
       if(rho_in=="no"){MOorNO="no";}
       rho_input_file>>nregions;
      }
     }
     //Get the method to use (Cuhre, Vegas, etc)
     do
     {
      getline(rho_input_file,rho_in);
     }while(rho_in=="");
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
     lowercase(rho_in);
     if(rho_in=="vegas"){method_cuba="Vegas";}
     else if(rho_in=="suave"){method_cuba="Suave";}
     else if(rho_in=="divonne"){method_cuba="Divonne";}
     else{method_cuba="Cuhre";}
     //Read interval
     rho_input_file>>interval_integralsCUB[0]>>interval_integralsCUB[1];
     rho_input_file>>interval_integralsCUB[2]>>interval_integralsCUB[3];
     rho_input_file>>interval_integralsCUB[4]>>interval_integralsCUB[5];
     if(int_file)
     {
      interval_integrals=new double*[nregions];
      for(i=0;i<nregions;i++)
      {interval_integrals[i]=new double[6];}
      //The first int gets the first intervale from the file!
      for(i=0;i<6;i++)
      {interval_integrals[0][i]=interval_integralsCUB[i];}
      for(i=1;i<nregions;i++)
      {
       //Read the rest of the intervals
       rho_input_file>>interval_integrals[i][0]>>interval_integrals[i][1];
       rho_input_file>>interval_integrals[i][2]>>interval_integrals[i][3];
       rho_input_file>>interval_integrals[i][4]>>interval_integrals[i][5];
      }
     }
    }
    else if(cuba_cubature=="cubature")
    {
     cubature=true;
     rho_input_file>>error_abs;
     //Number of integral operations ops
     rho_input_file>>ops;
     integral_ops=new string[ops];
     for(i=0;i<ops;i++)
     {
      //Get the operation to compute
      do
      {
       getline(rho_input_file,rho_in);
      rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
      }while(rho_in=="");
      lowercase(rho_in);
      integral_ops[i]=rho_in;
      if(integral_ops[i]=="sij")
      {
       //Get info of the orbitals to use MOs or NOs
       do
       {
        getline(rho_input_file,rho_in);
        rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
       }while(rho_in=="");
       lowercase(rho_in);
       if(rho_in=="no"){MOorNO="no";}
       rho_input_file>>orb1>>orb2;
      }
     }
     //Read interval
     rho_input_file>>interval_integralsCUB[0]>>interval_integralsCUB[1];
     rho_input_file>>interval_integralsCUB[2]>>interval_integralsCUB[3];
     rho_input_file>>interval_integralsCUB[4]>>interval_integralsCUB[5];
    }
    else if(cuba_cubature=="cubature2")
    {
     cubature2=true;
     rho_input_file>>error_abs;
     //Number of integral operations ops
     rho_input_file>>ops;
     integral_ops=new string[ops];
     for(i=0;i<ops;i++)
     {
      //Get the operation to compute
      do
      {
       getline(rho_input_file,rho_in);
      rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
      }while(rho_in=="");
      lowercase(rho_in);
      integral_ops[i]=rho_in;
      if(integral_ops[i]=="sij")
      {
       cout<<"Warning! Cubature2 does not have implemented the overlaps!"<<endl;
      }
     }
     //Read interval
     rho_input_file>>interval_integralsCUB[0]>>interval_integralsCUB[1];
     rho_input_file>>interval_integralsCUB[2]>>interval_integralsCUB[3];
     rho_input_file>>interval_integralsCUB[4]>>interval_integralsCUB[5];
    }
    else  //Quadrature
    {
     quadrature=true;
     rho_input_file>>order_grid_r>>order_grid_ang;
     if(order_grid_ang<=30){order_grid_ang=30;}
     if(order_grid_r<=30){order_grid_r=30;}
     //Number of integral operations ops
     rho_input_file>>ops;
     integral_ops=new string[ops];
     interval_integrals=new double*[ops];
     for(i=0;i<ops;i++)
     {interval_integrals[i]=new double[6];}
     for(i=0;i<ops;i++)
     {
      rho_in="";
      //Get the operation to compute
      do
      {
       getline(rho_input_file,rho_in);
       rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
      }while(rho_in=="");
      lowercase(rho_in);
      integral_ops[i]=rho_in;
      //Read interval
      rho_input_file>>interval_integrals[i][0]>>interval_integrals[i][1];
      rho_input_file>>interval_integrals[i][2]>>interval_integrals[i][3];
      rho_input_file>>interval_integrals[i][4]>>interval_integrals[i][5];
      if(integral_ops[i]=="sij" || integral_ops[i]=="r1_moment")
      {
       do
       {
        getline(rho_input_file,rho_in);
        rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
       }while(rho_in=="");
       lowercase(rho_in);
       if(rho_in=="no"){MOorNO="no";}
      }
     }
    }
   }
   else if(rho_in=="$dmn_r")
   {
    rho_input_file>>dmns;
    dmn=true;
    dmn_order=new int[dmns];
    name_dmn=new string [dmns];
    for(i=0;i<dmns;i++)
    {
     do
     {
      getline(rho_input_file,rho_in);
      rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
     }while(rho_in=="");
     name_dmn[i]=rho_in;
     rho_input_file>>dmn_order[i];
     if(dmn_order[i]==1)
     {
      rho_input_file>>pointDM1[0]>>pointDM1[1]>>pointDM1[2];
      rho_input_file>>pointprimeDM1[0]>>pointprimeDM1[1]>>pointprimeDM1[2];
     }
     else
     {}
    }
   }
   else if(rho_in=="$dmn_p")
   {
    rho_input_file>>dmnsp;
    dmnp=true;
    dmn_orderp=new int[dmnsp];
    name_dmnp=new string [dmnsp];
    for(i=0;i<dmnsp;i++)
    {
     do
     {
      getline(rho_input_file,rho_in);
      rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
     }while(rho_in=="");
     name_dmnp[i]=rho_in;
     rho_input_file>>dmn_orderp[i];
     if(dmn_orderp[i]==1)
     {
      rho_input_file>>pointDM1p[0]>>pointDM1p[1]>>pointDM1p[2];
      rho_input_file>>pointprimeDM1p[0]>>pointprimeDM1p[1]>>pointprimeDM1p[2];
     }
     else
     {}
    }
   }
   else if(rho_in=="$integrals_dmn")
   {
    dmn_integrals=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm1=rho_in;
   }
   else if(rho_in=="$integrals_rotate_grid")
   {
    rotate_grid=true;
    rho_input_file>>rotate_angle;
   }
   else if(rho_in=="$indicators_dmn")
   {
    dmn_indicators=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm1=rho_in;
   }
   else if(rho_in=="$dmn_threshold"){dmn_thresh=true;rho_input_file>>dmn_threshold;}
   else if(rho_in=="$plot")
   {
    string aux;
    gnuplot=true;
    do
    {
     getline(rho_input_file,path);
     path.erase(std::remove_if(path.begin(),path.end(),::isspace),path.end());
    }while(path=="");
    aux=path;
    lowercase(aux);
    if(!(aux=="nopath")){nopath=false;} //A path is required!
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    lowercase(rho_in);
    if(rho_in=="3d")
    {
     dim3=true;
     rho_input_file>>scan_1>>scan_2;
     rho_input_file>>points_scan_1[0]>>points_scan_1[1];
     rho_input_file>>points_scan_2[0]>>points_scan_2[1];
     rho_input_file>>grid_1>>grid_2;
     rho_input_file>>not_used_3d>>not_used_3d_doub;
     rho_input_file>>num_plot_ops;
     plot_ops=new string[num_plot_ops];
     for(i=0;i<num_plot_ops;i++)
     {
      do
      {
       getline(rho_input_file,rho_in);
       rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
      }while(rho_in=="");
      lowercase(rho_in);
      plot_ops[i]=rho_in;
     }
     rho_input_file>>extra_lines;
     if(extra_lines!=0)
     {
      extra_lines_plot=new string[extra_lines];
      for(i=0;i<extra_lines;i++)
      {
       do
       {
        getline(rho_input_file,rho_in);
        rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
       }while(rho_in=="");
       lowercase(rho_in);
       extra_lines_plot[i]=rho_in;
      }
     }
    }
    else if(rho_in=="2d")
    {
     dim2=true;
     rho_input_file>>scan_1;
     rho_input_file>>points_scan_1[0]>>points_scan_1[1];
     rho_input_file>>grid_1;
     rho_input_file>>num_plot_ops;
     plot_ops=new string[num_plot_ops];
     for(i=0;i<num_plot_ops;i++)
     {
      do
      {
       getline(rho_input_file,rho_in);
       rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
      }while(rho_in=="");
      lowercase(rho_in);
      plot_ops[i]=rho_in;
     }
     rho_input_file>>extra_lines;
     if(extra_lines!=0)
     {
      for(i=0;i<extra_lines;i++)
      {
       do
       {
        getline(rho_input_file,rho_in);
        rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
       }while(rho_in=="");
       lowercase(rho_in);
       extra_lines_plot[i]=rho_in;
      }
     }
    }
    else{}
   }
   else if(rho_in=="$plots_dmn")
   {
    dmn_plots=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    name_dm1=rho_in;
   }
   else if(rho_in=="$esi_int")
   {
    esi_int=true;
    do
    {
     getline(rho_input_file,rho_in);
     rho_in.erase(std::remove_if(rho_in.begin(),rho_in.end(),::isspace),rho_in.end());
    }while(rho_in=="");
    Sij_region=rho_in;
   }
   else if(rho_in=="$debug"){debug=true;}
   else{}
  }
  rho_input_file.close();
 }
 else
 {cout<<"File "<<name<<" not found!"<<endl;}
}
Input::~Input()
{
 if(punctualr)
 {
  for(i=0;i<punctuals_r;i++)
  {delete[] coordinates_r[i];}
  delete[] coordinates_r;
 }
 if(punctualp)
 {
  for(i=0;i<punctuals_p;i++)
  {delete[] coordinates_p[i];}
  delete[] coordinates_p;
 }
 if(cuba)
 {
  if(int_file)
  {
   for(i=0;i<nregions;i++)
   {delete[] interval_integrals[i];}
   delete[] interval_integrals;
  }
  delete[] integral_ops;
 }
 if(cubature)
 {
  delete[] integral_ops;
 }
 if(quadrature)
 {
  for(i=0;i<ops;i++)
  {delete[] interval_integrals[i];}
  delete[] interval_integrals;
  delete[] integral_ops;
 }
 if(dmn)
 {
  delete[] dmn_order;
  delete[] name_dmn;
 }
 if(dmnp)
 {
  delete[] dmn_orderp;
  delete[] name_dmnp;
 }
 if(gnuplot)
 {
  if(dim2 || dim3)
  {
   delete[] plot_ops;
   if(extra_lines!=0)
   {
    delete[] extra_lines_plot;
   }
  }
 }
}
