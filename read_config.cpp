#include <string.h>

void readConfigFile(const char * fileName,const char *configType) 
{
  
  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the parameter file %s!\n",fileName);
      exit(-1);
    }
  
  while(!(feof(config)))//go until the end of file is reached
    {
      if(fgets(cfgstr,256,config)!=NULL)
        {
          if(sscanf(cfgstr,"%s %s",str1,str2)==2) //single parameter data
	    {
	      // file list and tree name
	      if(strcmp(str1,"INPUT_FILE_LIST")==0)
		strcpy(inp_filename,str2);
	      if(strcmp(str1,"TREE_NAME")==0)
		strcpy(tree_name,str2);
	      if(strcmp(str1,"PROJ_A")==0)
		Ap=atof(str2);
	      if(strcmp(str1,"PROJ_Z")==0)
		Zp=atof(str2);
	      if(strcmp(str1,"RECL_A")==0)
		Ar=atof(str2);
	      if(strcmp(str1,"RECL_Z")==0)
		Zr=atof(str2);
	      if(strcmp(str1,"EGAMMA")==0)
		E0=atof(str2);
	      if(strcmp(str1,"PROJ_REAC_IN_PATH")==0)
		strcpy(prIn_name,str2);
	      if(strcmp(str1,"PROJ_REAC_OUT_PATH")==0)
		strcpy(prOut_name,str2);
	      if(strcmp(str1,"PROJ_DECAY_PATH")==0)
		strcpy(prDec_name,str2);
	      if(strcmp(str1,"GAMMA_X")==0)
		strcpy(gx_name,str2);
	      if(strcmp(str1,"GAMMA_Y")==0)
		strcpy(gy_name,str2);
	      if(strcmp(str1,"GAMMA_Z")==0)
		strcpy(gz_name,str2);
	      if(strcmp(str1,"OUTPUT_FILE")==0)
		strcpy(out_filename,str2);
               
	      if(strcmp(str1,"SORT_PATH")==0)
                strcpy(sort_path,str2);
	      if(strcmp(str1,"SORT_DATA_SCALING_FACTOR")==0)
                sort_scaling=atof(str2);

	      // group map stuff
	      if( (strcmp(str1,"GROUP_MAP_PATH")==0))
                strcpy(group_file,str2);
	      if(strcmp(str1,"POS_PATH")==0)
                strcpy(pos_path,str2);
	      if(strcmp(str1,"COL_PATH")==0)
                strcpy(col_path,str2);
	      if(strcmp(str1,"CSI_PATH")==0)
                strcpy(csi_path,str2);

	      // FWHM response
	      if(strcmp(str1,"SORT_DATA_FWHM_RESPONSE")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    fwhmResponse=true;
                  else
                    fwhmResponse=false;
                }
              if(strcmp(str1,"FWHM_F")==0)
                fwhmF=atof(str2);
              if(strcmp(str1,"FWHM_G")==0)
                fwhmG=atof(str2);
              if(strcmp(str1,"FWHM_H")==0)
                fwhmH=atof(str2);

	    }
	  if(sscanf(cfgstr,"%s %s",str1,str2)==1) //only one item on line
            if(strcmp(str1,"<---END_OF_PARAMETERS--->")==0)
              break;  
	}
    }
  fclose(config);
    
  //Report parameters based on the config file type used
  if(strcmp(configType,"coulex_ang_dist")==0)
    {
      printf("Input list file: %s\n",inp_filename);
      printf("Sorting from tree name: %s\n",tree_name);
      printf("Projectile Z,A: %d,%d\n",(int)Zp,(int)Ap);
      printf("Recoil Z,A: %d,%d\n",(int)Zr,(int)Ar);
      printf("Projectile Egamma: %f\n",E0);
      printf("Projectile reaction in branch : %s\n",prIn_name);
      printf("Projectile reaction out branch: %s\n",prOut_name);
      printf("Projectile decay branch: %s\n",prDec_name);
      printf("Gamma (x,y,z) branches: (%s,%s,%s)\n",gx_name,gy_name,gz_name);
      printf("Sorting from leaf with path: %s in tree: %s\n",sort_path,tree_name);
      printf("Using groups defined in %s\n",group_file);
      if(fwhmResponse==true)
        {
          printf("Will apply FWHM response function to sorted data.\n");
          printf("FWHM response function paremeters: F=%f, G=%f, H=%f.\n",fwhmF,fwhmG,fwhmH);
        }
      if(sort_scaling!=1.0)
        printf("Will scale sorted data by a factor of %f\n",sort_scaling);
        printf("Will save output data to file: %s\n",out_filename);
    }
  
}
