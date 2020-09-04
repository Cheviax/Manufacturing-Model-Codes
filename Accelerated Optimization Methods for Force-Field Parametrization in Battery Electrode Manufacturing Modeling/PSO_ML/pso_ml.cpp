#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <cassert>

#include <unistd.h>

#include <fstream>
#include <list>

#include <cstdlib>
#include <string>


#define WITH_THREADS
#ifdef WITH_THREADS
#include <thread>
#endif //WITH_THREADS


//#define USE_GPU_NODES 


#define MODE_FOR_COMPUTING_VISCOSITY 1
#define SIZE_FOR_VISCO_AVG 100

//#define ETA_EXP1 1092
//#define ETA_EXP2 293
#define DENSITY_EXP 2.14
#define VISCOSITY1_EXP 4100
#define VISCOSITY2_EXP 1550


double to_failure()
{
  return 100.0;
}


bool has_failed(double v)
{
  return v>=to_failure();
}


bool has_failed(const vector<double>& v_res)
{
  return has_failed(v_res[v_res.size()-1]);
}


bool list_isin(const list<char>& l, char x)
{
  list<char>::const_iterator i_l = l.begin();
  
  while(i_l!=l.end())
    {
      if(*i_l==x)
	return true;
      
      i_l++;
    }
  
  return false;
}


double list_average(const list<double>& l)
{
  if(l.empty()) return 0;
  
  double avg = 0;
  list<double>::const_iterator i_l = l.begin();
  while(i_l!=l.end())
    {
      avg+=*i_l;
      i_l++;
    }
  
  return avg/(double)l.size();
}

double list_standardDeviation(const list<double>& l)
{
  if(l.empty()) return 0;
  double avg = list_average(l);
  
  double s = 0;
  list<double>::const_iterator i_l = l.begin();
  while(i_l!=l.end())
    {
      double diff = ((double)*i_l)-avg;
      s+=(diff*diff);
      i_l++;
    }
  
  s = s/((double)(l.size()-1)); //l.size()-1 ???
  
  return sqrt(s);
}


bool filetools_exist(const string filename)
{
  ifstream flux(filename);
  if(not flux.is_open())
    return false;

  flux.close();
  return true;
}

bool filetools_tostring(const string& filename, string& s)
{
  ifstream flux(filename);
  if(not flux.is_open())
    {
      //cerr<<"*** ERROR *** "<<filename<<" not found"<<endl;
      return false;
    }
      
  s = std::string((std::istreambuf_iterator<char>(flux)),
		  std::istreambuf_iterator<char>());
  flux.close();
  return true;
}


bool stringtools_toDouble(const char *s, double & v) 
{
  char *p;
  v=strtod(s, &p);
  return p!=s;
}

bool stringtools_toDouble(const string s, double& v)
{
  return stringtools_toDouble(s.c_str(), v);
}


int stringtools_replace(string&  s, const string& s1, const string& s2)
{
  //cerr<<s<<" "<<s1<<" "<<s2<<endl; sleep(1);
  
  int ind = 0;
  while(true)
    {
      ind = s.find(s1, ind);
      
      if(ind==(int)std::string::npos) break;
      if(ind>=(int)s.length()) break;

      //cerr<<ind<<" - ";
      
      s.replace(ind, s1.length(), s2);
      ind = ind+s2.length();
      
      //cerr<<ind<<" "<<s.length()<<endl;
    }
  
  //cerr<<"== "<<s1<<" "<<s2<<" "<<s<<endl;
  
  return s.size();
}

int stringtools_split(const string& s, 
		      string& res1,
		      string& res2,
		      char sep = ' ')
{
  char c_inhib = ' ';
  char c_desinhib = ' ';
  int i = 0;
  bool to_stop = false;
  
  list<char> l_sep; l_sep.push_back(sep);
  
  res1 = "";
  res2 = "";
  sep = ' ';
  
  int nb = 0;
  for(; i<(int)s.size(); i++)
    {
      if(nb==0)
	{
	  if(list_isin(l_sep, s[i]))
	    {
	      sep = s[i];
	      break;
	    }
	}
      
      if(c_desinhib!=c_inhib)
	{
	  if(s[i]==c_inhib) nb++;
	  if(s[i]==c_desinhib) nb--;
	}
      
      res1+=s[i];
    }
  
  i++;
  
  if(to_stop) return i;
  
  for(; i<(int)s.size(); i++)
    {
      res2+=s[i];
    } 
  
  //cerr<<":: "<<res1<<" "<<sep<<" "<<res2<<endl;
  
  return i;
}


string stringtools_upper(const string& s1)
{
  string s = s1;
  for(int i = 0; i<(int)s.size(); i++)
    {
      if((s[i]>='a')&&(s[i]<='z'))
	{
	  s[i] = s[i]-'a'+'A';
	}
    }
  
  return s;
}


bool stringtools_extract(const string& line, string& s_param, char tag = '@')
{
  bool to_keep = false;
  for(size_t i = 0; i<line.size(); i++)
    { 
      if(line[i]==tag)
	{
	  if(to_keep)
	    {
	      s_param = s_param+tag;
	      return true;
	    }
	  else
	    to_keep = true;
	}

      if(to_keep)
	s_param = s_param+line[i]; 
    }
  
  return false;
}

double LRCS_readDensity(string& filename)
{
  ifstream flux(filename);
  if(not flux.is_open())
    {
      cerr<<"*** ERROR *** "<<filename<<" not found"<<endl;
      //assert(0);
      return 100;
    }
 
  string last_line;
  while(not flux.eof())
    {
      string line;
      std::getline(flux, line);
      if(line!="")
	last_line = line;
      //cerr<<"I read <"<<last_line<<">"<<endl;
    }
  
  //cerr<<"I read "<<last_line<<endl;
  
  string f1, f2;
  stringtools_split(last_line, f1, f2, ' ');
  
  
  double val2 = -1; 
  if(not stringtools_toDouble(f2, val2))
    {
      cerr<<"failed"<<endl;
      return 0;
    }
  
  
  //cerr<<"density = "<<val2<<endl;
  
  return val2;

}

pair<double, double> LRCS_readViscosity(string& filename, size_t mode = 1, size_t p1 = 100)
{
  ifstream flux(filename);
  if(not flux.is_open())
    { 
      cerr<<"*** ERROR *** "<<filename<<" not found"<<endl;
      //assert(0);
      return make_pair(100000, 0);
    }
  
  list<double> l_visco;
  double v = 0;
  double dev = 0; 
  while(not flux.eof())
    {
      string line;
      std::getline(flux, line);
      
      //cerr<<"I read "<<line<<endl;
      
      string f1, f2;
      stringtools_split(line, f1, f2, ' ');
      
      //cerr<<f2<<endl;
      
      double val2 = -1; 
      if(not stringtools_toDouble(f2, val2))
	{
	  //cerr<<"failed"<<endl;
	  continue;
	}
      
      
      if(mode==1)
	{
	  l_visco.push_back(val2);
	  if(l_visco.size()>p1)
	    l_visco.pop_front();
	}	
    }
  
  if(l_visco.empty())
    {
      cerr<<"*** ERROR *** "<<filename<<" no found value"<<endl;
      //assert(0);
      return make_pair(100000, 0);
    }

  if(mode==1)
    {
      list<double>::const_iterator i_visco = l_visco.begin();
      while(i_visco!=l_visco.end())
	{
	  v = v+*i_visco;
	  i_visco++;
	}
      
      v = v/(double)(l_visco.size());
      
      dev = list_standardDeviation(l_visco);
      double avg = list_average(l_visco);
      
      if(avg!=v)
	{
	  cerr<<"*** WARNING *** different averages : "<<avg<<" <> "<<v<<endl;
	}
    }
  
  //cerr<<l_visco.size()<<endl;
  
  flux.close();
  
  return make_pair(v, dev);
}


double LRCS_readTemperature(string& filename)
{
  ifstream flux(filename);
  if(not flux.is_open())
    {
      cerr<<"*** ERROR *** [ TEMPERATURE ] "<<filename<<" not found"<<endl;
      return 8000;
    }
 
  string last_line;
  while(not flux.eof())
    {
      string line;
      std::getline(flux, line);
      if(line!="")
	last_line = line;
      //cerr<<"I read <"<<last_line<<">"<<endl;
    }
  
  //cerr<<"I read "<<last_line<<endl;
  
  string f1, f2;
  stringtools_split(last_line, f1, f2, ' ');
  
  
  double val2 = -1; 
  if(not stringtools_toDouble(f2, val2))
    {
      cerr<<"Temperature : failed"<<f2<<" "<<val2<<endl;
      return 0;
    }
  
  return val2;
}



size_t LRCS_parseMetaInputFile(const string& filename, list<pair<string, string> >& l_param,
			       char tag = '@')
{
  ifstream flux(filename);
  if(not flux.is_open())
    {
      cerr<<filename<<" not found"<<endl;
      return 0;
    }
  
  while(not flux.eof())
    {
      string line;
      std::getline(flux, line);
      
      //cerr<<"I read "<<line<<endl;

      string s_param;
      if(not stringtools_extract(line, s_param, tag))
	continue;

      cerr<<"parameter : "<<s_param<<endl;
      cerr<<"line : "<<line<<endl;

      string stag = ""; stag = stag+tag;
      size_t i = line.find(stag);
      assert(i>0);
      i = i-1;
      while(true)
	{
	  if(line[i]!=' ')
	    break;
	  i--;
	  assert(i>0);
	}

      while(true)
	{
	  if(line[i]==' ')
	    break;
	  i--;
	  assert(i>0);
	}

      while(true)
	{
	  if(line[i]!=' ')
	    break;
	  i--;
	  assert(i>0);
	}

      //example
      //line : variable    epsi_am equal @P1@
      //<m> correspond au m de epsi_am

      cerr<<"<"<<line[i]<<">"<<endl;

      string s_name = "";
      while(true)
	{
	  if(line[i]==' ')
	    break;
	  
	  if(((line[i]>='0')&&(line[i]<='9'))
	     ||((line[i]>='A')&&(line[i]<='Z'))
	     ||((line[i]>='a')&&(line[i]<='z'))
	     ||((line[i]>='_')))
	    {
	      s_name = line[i]+s_name;
	    }
	  else
	    {
	      break;
	    }

	  if(i==0) break;
	  
	  i = i-1;
	}

      cerr<<"==> "<<s_name<<endl;

      l_param.push_back(make_pair(s_param, s_name));
      
      
      /*
      string f1, f2;
      stringtools_split(line, f1, f2, ' ');
      
      //cerr<<f2<<endl;
      
      double val2 = -1; 
      if(not stringtools_toDouble(f2, val2))
	{
	  //cerr<<"failed"<<endl;
	  continue;
	}
      */
    }

  
  return l_param.size();
}


void LRCS_generateInputFile(const string& metainputfilename,
			    const string& inputfilename,
			    const vector<pair<string, double> > v_values)
{
  //std::copy(metainputfilename, inputfilename);
  
  ifstream metainputflux(metainputfilename);
  string scontent( (std::istreambuf_iterator<char>(metainputflux) ),
		  (std::istreambuf_iterator<char>()    ) );
  metainputflux.close();
  
  for(size_t i = 0; i<v_values.size(); i++)
    {
      stringtools_replace(scontent, v_values[i].first, std::to_string(v_values[i].second));
    }

  ofstream inputflux(inputfilename);
  inputflux<<scontent;

  inputflux.close();
}



double simulate_with_lammps(const vector<double>& v_input, size_t num_particle,
			    const list<string>* filenamesList, vector<double>* v_res)
{

  double error_eval = to_failure();

  cerr<<"[ EVALUATION ] - SIMULATE WITH LAMMPS - PARTICLE "<<num_particle<<endl;
  
  assert(filenamesList);
  assert(not filenamesList->empty());
  
  list<string> tempFilenamesList = *filenamesList; 
  //cerr<<num_particle<<" "<<filenamesList->front()<<endl;


  string metainput_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  string metainput1_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  string metainput2_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  
  string my_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  stringtools_replace(my_filename, "@NUMPARTICLE@", std::to_string(num_particle));

  string my_filename1 = tempFilenamesList.front(); tempFilenamesList.pop_front();
  stringtools_replace(my_filename1, "@NUMPARTICLE@", std::to_string(num_particle));

  string my_filename2 = tempFilenamesList.front(); tempFilenamesList.pop_front();
  stringtools_replace(my_filename2, "@NUMPARTICLE@", std::to_string(num_particle));

  /*
  assert(not filenamesList->empty());
  string visco1_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  stringtools_replace(visco1_filename, "@NUMPARTICLE@", std::to_string(num_particle));

  assert(not filenamesList->empty());
  string density1_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  //stringtools_replace(density1_filename, "@NUMPARTICLE@", std::to_string(num_particle));
  
  //assert(not filenamesList->empty());
  //string visco2_filename = tempFilenamesList.front(); tempFilenamesList.pop_front();
  //stringtools_replace(visco2_filename, "@NUMPARTICLE@", std::to_string(num_particle));
  

  //Creer le fichier d'input

  cerr<<"visco1 filename : "<<visco1_filename<<endl;
  cerr<<"density1 filename : "<<density1_filename<<endl;
  */  


  //Question : les v_input contient les valeurs des paramètres
  vector<pair<string,double> > v_values;
  for(size_t i = 0; i<(size_t)v_input.size(); i++)
    {
      //cerr<<v_input[i]<<endl; sleep(1);
      v_values.push_back(make_pair("@P"+std::to_string(i+1)+"@",v_input[i]));
    }

  //cerr<<"==> "<<my_filename<<" "<<visco1_filename<<" "<<density1_filename<<endl;

  
  LRCS_generateInputFile(metainput_filename, my_filename, v_values);
  LRCS_generateInputFile(metainput1_filename, my_filename1, v_values);
  LRCS_generateInputFile(metainput2_filename, my_filename2, v_values);
  
  //TODO Iteration?
  //Ajouter coordin.data dans filenamesList

  string short_my_filename, res2;
  stringtools_split(my_filename, short_my_filename, res2, '_');  
  
  string repo_rootname = "XPLAMMPS";
  string repo_rootname_full = repo_rootname+std::to_string(num_particle);


  /*
  //Test de l'inexistence des outputs (density)
  cerr<<"test de l'inexistence des outputs (density)"<<endl;
  if(filetools_exist(repo_rootname_full+"/"+density1_filename))
  cerr<<repo_rootname_full+"/"+density1_filename<<" exists"<<endl;
  else cerr<<repo_rootname_full+"/"+density1_filename<<" does not exist"<<endl;
  if(filetools_exist(repo_rootname_full+"/coord_out_slurry.data"))
  cerr<<repo_rootname_full+"/coord_out_slurry.data"<<" exists"<<endl;
  else cerr<<repo_rootname_full+"/coord_out_slurry.data"<<" does not exist"<<endl;
  */
  
  size_t time_sleep = 30;

   
  string short_my_filename1, poub1;
  stringtools_split(my_filename1, short_my_filename1, poub1, '_');  
  
  string short_my_filename2, poub2;
  stringtools_split(my_filename2, short_my_filename2, poub2, '_');  


  const string short_visco1_LowerLimit_filename = "visco1_LowerLimit";
  const string short_visco2_LowerLimit_filename = "visco2_LowerLimit";
  
  //string visco1_LowerLimit_filename = repo_rootname_full+"/"+short_visco1_LowerLimit_filename;
  //string visco2_LowerLimit_filename = repo_rootname_full+"/"+short_visco2_LowerLimit_filename;
  
#ifdef USE_GPU_NODES
  cerr<<"use gpunodes"<<endl; 

  string command_lammps_slurry = "./launch_gpu.sh "+std::to_string(num_particle)+" "+short_my_filename+" coord_in.data "+repo_rootname+" "+"script_gpu.template script_gpu.sh";
  
  string command_lammps_visco1 = "./launch_gpu.sh "+std::to_string(num_particle)+" "+short_my_filename1+" coord_out_slurry.data "+repo_rootname+" "+"script1_gpu.template script1_gpu.sh";
  
  string command_lammps_visco2 = "./launch_gpu.sh "+std::to_string(num_particle)+" "+short_my_filename2+" coord_out_slurry.data "+repo_rootname+" "+"script2_gpu.template script2_gpu.sh";

  string s_clear4 = "rm -f "+repo_rootname_full+"/"+"script_gpu.sh";
  
#else
  //cerr<<"use midnodes"<<endl;
  
  string command_lammps_slurry = "./launch.sh "+std::to_string(num_particle)+" "+short_my_filename+" coord_in.data "+repo_rootname+" "+"script.template script.sh";
  
  string command_lammps_visco1 = "./launch.sh "+std::to_string(num_particle)+" "+short_my_filename1+" coord_out_slurry.data "+repo_rootname+" "+"script1.template script1.sh";
  
  string command_lammps_visco2 = "./launch.sh "+std::to_string(num_particle)+" "+short_my_filename2+" coord_out_slurry.data "+repo_rootname+" "+"script2.template script2.sh";

  string s_clear4 = "rm -f "+repo_rootname_full+"/"+"script.sh";
#endif
  
 

  
  //cerr<<command_lammps_slurry<<endl;
  cerr<<"step 1 - liquid - "<<"particle "<<num_particle<<endl;		
  system(command_lammps_slurry.c_str());
  //assert(0);

  string flagfile = "xpfinish";
  //Penser a changer dans script.template (touch xpfinish) et un peu plus loin
  //dans ce fichier (voir rm -f ==> nettoyage)
   
  string flagpath = repo_rootname_full+"/"+flagfile;

  size_t num_trial = 1;
  while(not filetools_exist(flagpath))
    {
      //cerr<<"waiting end of experiment from "<<num_trial*time_sleep<<" seconds"<<endl;
      num_trial++;
      sleep(time_sleep);
    }


  /*
  //cerr<<"je passe o la suite "<<endl;
  //Test de l'existence des outputs (density)
  cerr<<"test de l'existence des outputs (density)"<<endl;
  if(filetools_exist(repo_rootname_full+"/"+density1_filename))
  cerr<<repo_rootname_full+"/"+density1_filename<<" exists";
  else cerr<<repo_rootname_full+"/"+density1_filename<<" does not exist";
  if(filetools_exist(repo_rootname_full+"/coord_out_slurry.data"))
  cerr<<repo_rootname_full+"/coord_out_slurry.data"<<" exists"<<endl;
  else cerr<<repo_rootname_full+"/coord_out_slurry.data"<<" does not exist"<<endl;
  */

  string temperature_filename = repo_rootname_full+"/Temperature";
  double temperature = LRCS_readTemperature(temperature_filename);

  string s_coordout = "";
  bool is_ok = true;
  
  if(not filetools_tostring(repo_rootname_full+"/coord_out_slurry.data", s_coordout))
    {
      cerr<<"*** ERROR *** "<<repo_rootname_full+"/coord_out_slurry.data"<<" not found"<<endl;
      is_ok = false;
    }

   //Modifier Coord_out
  //ajouter une ligne dans le fichier de sortie coord_out
  //"0.0 0.0 0.0 xy xz yz"
  

  if(is_ok)
    {
      string s_to_find = " zhi";
      string s_new = s_to_find+"\n"+"0.0 0.0 0.0 xy xz yz";
      size_t old_size = s_coordout.size();
      stringtools_replace(s_coordout, s_to_find, s_new);
      
      if(s_coordout.size()<=old_size)
	{
	  cerr<<"*** ERROR *** "<<"the string <"<<s_to_find<<"> not found in "+repo_rootname_full+"/coord_out_slurry.data (triclinic)";
	  is_ok = false;
	}
    }
  
  
  
  double temperature_LowerLimit = 200;
  double temperature_UpperLimit = 400;
  if((temperature<temperature_LowerLimit)||(temperature>temperature_UpperLimit)||(not is_ok))
    {
      if(temperature<temperature_LowerLimit)
	{
	  cerr<<"! "<<"particle "<<num_particle<<" : "<<temperature<<"<"<<temperature_LowerLimit<<" ==> "<<"Lower bound of the temperature not respected"<<endl;
	}
      
      if(temperature>temperature_UpperLimit)
	{
	  cerr<<"! "<<"particle "<<num_particle<<" : "<<temperature<<">"<<temperature_UpperLimit<<" ==> "<<"Upper bound of the temperature not respected"<<endl;
	}
      
      
      string s_clear_all1 = "rm -Rf "+repo_rootname_full;
      string s_clear_all2 = "rm -f *.run_*";
      system(s_clear_all1.c_str());
      system(s_clear_all2.c_str());
      
      //double error_eval = 100;
      if(v_res)
	{
	  if(v_res->size()==0)
	    {
	      //A checker! Susceptible de generer des Segmentation Fault
	      v_res->resize(1+1+v_input.size()+1+3+2+1);
	    }
	  
	  (*v_res)[0] = 0; //Le numero d'iteration sera mise dans la fonction PS0...
	  (*v_res)[1] = num_particle;
	  
	  for(size_t i = 0; i<v_input.size(); i++)
	    {
	      assert(i+2<v_res->size());
	      (*v_res)[i+2] = v_input[i];
	    }
	  
	  if(2+v_input.size()+7!=v_res->size())
	    {
	      cerr<<"*** ACHECKER *** [ RESULTS ] "<<2+v_input.size()+7<<" <> "<<v_res->size()<<endl;
	    }
	  
	  assert(2+v_input.size()+6<v_res->size());
		
	  (*v_res)[2+v_input.size()] = temperature;
	  (*v_res)[2+v_input.size()+1] = 100;
	  (*v_res)[2+v_input.size()+2] = 1000000;
	  (*v_res)[2+v_input.size()+3] = 1000000;
	  (*v_res)[2+v_input.size()+4] = 0;
	  (*v_res)[2+v_input.size()+5] = 0;
	  (*v_res)[2+v_input.size()+6] = error_eval;
	}
      
      return error_eval;      
    }
  
 
  
  //s_coordout = s_coordout.replace(s_coordout.find(s_to_find),s_to_find.length(),s_new);
  
  //cerr<<s_coordout<<endl;
  
  ofstream coordoutfile(repo_rootname_full+"/coord_out_slurry.data");
  coordoutfile<<s_coordout;
  coordoutfile.close();
  sleep(1);
  
  
  //TODO Sauvegarder parametres avec density et viscosity
  
  //nettoyage
  string s_clear1 = "rm -f "+repo_rootname_full+"/"+"coord_in.data";
  string s_clear2 = "rm -f "+repo_rootname_full+"/"+"dump.com*";
  string s_clear3 = "rm -f "+repo_rootname_full+"/"+"xpfinish";
  
  string s_clear5 = "rm -f "+repo_rootname_full+"/"+"inliqreal.run";
  system(s_clear1.c_str());
  system(s_clear2.c_str());
  system(s_clear3.c_str());
  system(s_clear4.c_str());
  //system(s_clear5.c_str());
  sleep(1);
  
  
  cerr<<"step 2 - first viscosity - "<<"particle "<<num_particle<<endl;
  system(command_lammps_visco1.c_str());
  
  cerr<<"step 3 - second viscosity - "<<"particle "<<num_particle<<endl;							
  system(command_lammps_visco2.c_str());
  
  string flagfile1 = "xpfinish1";
  //Penser a changer dans script.template (touch xpfinish1) et un peu plus loin
  //dans ce fichier (voir rm -f ==> nettoyage)
  
  
  string flagpath1 = repo_rootname_full+"/"+flagfile1;
  
  //Attente du resultat viscosite 1
  num_trial = 1;
  while(not filetools_exist(flagpath1))
    {
      //cerr<<"waiting end of experiment 1 from "<<num_trial*time_sleep<<" seconds"<<endl;
      num_trial++;
      sleep(time_sleep);
    }
  
  
  string flagfile2 = "xpfinish2";
  string flagpath2 = repo_rootname_full+"/"+flagfile2;
  
  
  //Attente du resultat viscosite 2
  num_trial = 1;
  while(not filetools_exist(flagpath2))
    {
      //cerr<<"waiting end of experiment 2 from "<<num_trial*time_sleep<<" seconds"<<endl;
      num_trial++;
      sleep(time_sleep);
    }
  
  
  
  /*
    LAMMPS data file via write_data, version 31 Mar 2017, timestep = 0
    
    4077 atoms
    6 atom types
    
    -1.0000000000000000e+02 1.0000000000000000e+02 xlo xhi
    -1.0000000000000000e+02 1.0000000000000000e+02 ylo yhi
    -3.0000000000000000e+02 3.0000000000000000e+02 zlo zhi
    HERE
  */
  
  //Creer script batch pour viscosite
  cerr<<"Individu "<<num_particle<<" - input : ";
  for(size_t i = 0; i<v_input.size(); i++)
    cerr<<v_input[i]<<" ";
  cerr<<endl;
  
  //system("sbatch script");
  //sleep(1);
  
  
  string visco1_filename = repo_rootname_full+"/"+"ave_visco1";
  string visco2_filename = repo_rootname_full+"/"+"ave_visco2";
  string density1_filename = repo_rootname_full+"/"+"density_slurry";
  
  
  
  pair<double, double> viscosity1 = 
    LRCS_readViscosity(visco1_filename, MODE_FOR_COMPUTING_VISCOSITY, SIZE_FOR_VISCO_AVG);
  
  pair<double, double> viscosity2 =
    LRCS_readViscosity(visco2_filename, MODE_FOR_COMPUTING_VISCOSITY, SIZE_FOR_VISCO_AVG);
  //pair<double, double> out2 =
  //LRCS_readViscosity(visco2_filename, MODE_FOR_COMPUTING_VISCOSITY, SIZE_FOR_VISCO_AVG);
  
  double density1 =
    LRCS_readDensity(density1_filename);

  
  //double temperature = LRCS_readTemperature(temperature_filename);

  double visco1LowerLimit = 0; //default value
  ifstream visco1_LowerLimit_flux(short_visco1_LowerLimit_filename); //visco1_LowerLimit_filename);
  if(not visco1_LowerLimit_flux.is_open())
    {
      cerr<<"*** ERROR *** "<<short_visco1_LowerLimit_filename<<" not found"<<endl;
    }
  else
    {
      visco1_LowerLimit_flux>>visco1LowerLimit;
      visco1_LowerLimit_flux.close();
    }
  
  cerr<<"the lower limit of viscosity 1 is "<<visco1LowerLimit<<endl;
  
  if(viscosity1.first<visco1LowerLimit)
    {
      cerr<<"! "<<"particle "<<num_particle<<" : "<<viscosity1.first<<"<"<<visco1LowerLimit<<" ==> "<<"Lower bound of the viscosity 1 not respected"<<endl;
      
      string s_clear_all1 = "rm -Rf "+repo_rootname_full;
      string s_clear_all2 = "rm -f *.run_*";
      system(s_clear_all1.c_str());
      system(s_clear_all2.c_str());
      
      //double error_eval = 100;
      if(v_res)
	{
	  if(v_res->size()==0)
	    {
	      //A checker! Susceptible de generer des Segmentation Fault
	      v_res->resize(1+1+v_input.size()+1+3+2+1);
	    }
	  
	  (*v_res)[0] = 0; //Le numero d'iteration sera mise dans la fonction PS0...
	  (*v_res)[1] = num_particle;
	  
	  for(size_t i = 0; i<v_input.size(); i++)
	    {
	      assert(i+2<v_res->size());
	      (*v_res)[i+2] = v_input[i];
	    }

	  if(2+v_input.size()+7!=v_res->size())
	    {
	      cerr<<"*** ACHECKER1 *** [ RESULTS ] "<<2+v_input.size()+7<<" <> "<<v_res->size()<<endl;
	    }
	  
	  assert(2+v_input.size()+6<v_res->size());
	  
	  (*v_res)[2+v_input.size()] = temperature;
	  (*v_res)[2+v_input.size()+1] = 100;
	  (*v_res)[2+v_input.size()+2] = 1000001;
	  (*v_res)[2+v_input.size()+3] = 1000001;
	  (*v_res)[2+v_input.size()+4] = 0;
	  (*v_res)[2+v_input.size()+5] = 0;
	  (*v_res)[2+v_input.size()+6] = error_eval;
	}
      
      return error_eval;      
    }
    

  

  double visco2LowerLimit = 0; //default value
  ifstream visco2_LowerLimit_flux(short_visco2_LowerLimit_filename); 
  if(not visco2_LowerLimit_flux.is_open())
    {
      cerr<<"*** ERROR *** "<<short_visco2_LowerLimit_filename<<" not found"<<endl;
    }
  else
    {
      visco2_LowerLimit_flux>>visco2LowerLimit;
      visco2_LowerLimit_flux.close();
    }
  
  cerr<<"the lower limit of viscosity 2 is "<<visco2LowerLimit<<endl;
  
  
  if(viscosity2.first<visco2LowerLimit)
    {
      cerr<<"! "<<"particle "<<num_particle<<" : "<<viscosity2.first<<"<"<<visco2LowerLimit<<" ==> "<<"Lower bound of the viscosity 2 not respected"<<endl;
      
      string s_clear_all1 = "rm -Rf "+repo_rootname_full;
      string s_clear_all2 = "rm -f *.run_*";
      system(s_clear_all1.c_str());
      system(s_clear_all2.c_str());
      
      //double error_eval = 100;
      if(v_res)
	{
	  if(v_res->size()==0)
	    {
	      //A checker! Susceptible de generer des Segmentation Fault
	      v_res->resize(1+1+v_input.size()+1+3+2+1);
	    }
	  
	  (*v_res)[0] = 0; //Le numero d'iteration sera mise dans la fonction PS0...
	  (*v_res)[1] = num_particle;
	  
	  for(size_t i = 0; i<v_input.size(); i++)
	    {
	      assert(i+2<v_res->size());
	      (*v_res)[i+2] = v_input[i];
	    }
	  
	  if(2+v_input.size()+7!=v_res->size())
	    {
	      cerr<<"*** ACHECKER2 *** [ RESULTS ] "<<2+v_input.size()+7<<" <> "<<v_res->size()<<endl;
	    }
	  
	  assert(2+v_input.size()+6<v_res->size());
	  
	  (*v_res)[2+v_input.size()] = temperature;
	  (*v_res)[2+v_input.size()+1] = 100;
	  (*v_res)[2+v_input.size()+2] = 1000002;
	  (*v_res)[2+v_input.size()+3] = 1000002;
	  (*v_res)[2+v_input.size()+4] = 0;
	  (*v_res)[2+v_input.size()+5] = 0;
	  (*v_res)[2+v_input.size()+6] = error_eval;
	}
      
      return error_eval;      
    }
    




  
  
  cerr<<"particule "<<num_particle<< " - ";
  cerr<<"density "<<density1<<" - ";
  cerr<<"viscosity1 "<<viscosity1.first<<"+/-"<<viscosity1.second<<" - ";
  cerr<<"viscosity2 "<<viscosity2.first<<"+/-"<<viscosity2.second<<endl;
  
  double density_exp1 = DENSITY_EXP;
  double viscosity1_exp = VISCOSITY1_EXP;
  double viscosity2_exp = VISCOSITY2_EXP;
  
  //double ev = pow(density_exp1-density1, 2)+pow(viscosity1_exp-viscosity1.first, 2)
  //+pow(viscosity2_exp-viscosity2.first, 2);

  double density_norm = density1/density_exp1;
  double viscosity1_norm = viscosity1.first/viscosity1_exp;
  double viscosity2_norm = viscosity2.first/viscosity2_exp;

  
  //double ev = 0.9*pow(1-density_norm, 2)+0.05*pow(1-viscosity1_norm, 2)
  //+0.05*pow(1-viscosity2_norm, 2);

  double ev = 0.3333*pow(1-density_norm, 2)+0.3333*pow(1-viscosity1_norm, 2)
    +0.3333*pow(1-viscosity2_norm, 2);

  
  cerr<<"evaluation = "<<ev<<endl;


   //(1) Stocker les resultats dans un fichier (Attention besoin du numero d'iteration)
  //(2) Virer les repertoires XPLAMMPS
  string s_clear_all1 = "rm -Rf "+repo_rootname_full;
  string s_clear_all2 = "rm -f *.run_*";
  system(s_clear_all1.c_str());
  system(s_clear_all2.c_str());
  

  if(v_res)
    {
      if(v_res->size()==0)
	{
	  //A checker! Susceptible de generer des Segmentation Fault
	  v_res->resize(1+1+v_input.size()+1+3+2+1);
	}

      (*v_res)[0] = 0; //Le numero d'iteration sera mise dans la fonction PS0...
      (*v_res)[1] = num_particle;
      
      for(size_t i = 0; i<v_input.size(); i++)
	{
	  assert(i+2<v_res->size());
	  (*v_res)[i+2] = v_input[i];
	}
      
      if(2+v_input.size()+7!=v_res->size())
	{
	  cerr<<"*** ACHECKER3 *** [ RESULTS ] "<<2+v_input.size()+7<<" <> "<<v_res->size()<<endl;
	}
      
      assert(2+v_input.size()+6<v_res->size());

      (*v_res)[2+v_input.size()] = temperature;
      (*v_res)[2+v_input.size()+1] = density1;
      (*v_res)[2+v_input.size()+2] = viscosity1.first;
      (*v_res)[2+v_input.size()+3] = viscosity2.first;
      (*v_res)[2+v_input.size()+4] = viscosity1.second;
      (*v_res)[2+v_input.size()+5] = viscosity2.second;
      (*v_res)[2+v_input.size()+6] = ev;
    }
 
  return ev;
}




double sphere(const vector<double>& v_input)
{
  double s = 0;
  for(size_t i = 0; i<v_input.size(); i++)
    s = s+(v_input[i]*v_input[i]);
    //s = s+(v_input[i]*v_input[i]-2);

  return s;
}

double sphere(const vector<double>& v_input, size_t numParticle, const list<string>* l, vector<double>* v_res)
{
  double s = 0;
  for(size_t i = 0; i<v_input.size(); i++)
    s = s+(v_input[i]*v_input[i]);
  //s = s+(v_input[i]*v_input[i]-2);
  
  return s;
}


double cost_function(const vector<double>& v_input)
{
  return sphere(v_input);
}

double random_probability()
{
  //const unsigned int p = 1000000007;
  //return (rand()%(p+1))/((double)p);
  
  double ub = 1;
  double lb = 0;
  
  return (lb+((rand()%RAND_MAX)/(long double)(RAND_MAX))*(ub-lb));
}


class Particle
{
public:
  vector<double> position;
  vector<double> velocity;
  //vector<double> cost;
  double cost;
  vector<double> best_position;
  double best_cost;
  //vector<double> best_cost;
  
  Particle()
  {
    //cost = 0;
  }
  
  // void

  void init_position(unsigned int dim)
  {
    position.resize(dim, 0);
  }

  
  void init_random_position(unsigned int dim, const vector<double>& v_lb, const vector<double>& v_ub)
  {
    
    position.resize(dim);
    for(size_t i = 0; i<position.size(); i++)
      {
	double lb = v_lb[i];
	double ub = v_ub[i];
	position[i] = lb+(random_probability())*(ub-lb);
	cerr<<"position : random - "<<lb<<" "<<ub<<" "<<position[i]<<endl; //sleep(1);
      }
  }

  void init_human_position(const vector<double>& pos)
  {
    position = pos;
  }
  

  void init_velocity(unsigned int dim)
  {
    velocity.resize(dim, 0);
  }

  void init_random_velocity(unsigned int dim, const vector<double>& v_lb, const vector<double>& v_ub)
  {
    velocity.resize(dim);
    for(size_t i = 0; i<velocity.size(); i++)
      {
	double lb = v_lb[i];
	double ub = v_ub[i];
	velocity[i] = lb+(random_probability())*(ub-lb);
	cerr<<"velocity : random - "<<lb<<" "<<ub<<" "<<velocity[i]<<endl; //sleep(1);
      }
  }
  

  void eval(double (& f) (const vector<double>&))
  {
    cost = f(position); //cost_function(position);
    //cerr<<cost<<endl; sleep(1);
  }

  void eval(double (& f) (const vector<double>&, size_t, const list<string>*,  vector<double>*),
	    size_t num_particle, const list<string>* l_param, vector<double>* v_res)
  {
    cost = f(position, num_particle, l_param, v_res); //cost_function(position);
    //cerr<<cost<<endl; sleep(1);
  }

  void print(ostream& os)
  {
    os<<"position = ";
    for(size_t i = 0; i<position.size(); i++)
      {
	os<<position[i]<<" "; 
      }
    os<<"\n";
  }
};


//TODO CostFunction = problem.CostFunction;
struct PSO_Problem
{
  unsigned int nbVar; //Number of Unknown (Decision) Variables
  
  vector<double> varMin; //Lower Bound of Decision Variables
  vector<double> varMax; //Upper Bound of Decision Variables

  //double (* costFunction)(const vector<double>&);
  double (* costFunction)(const vector<double>&, size_t, const list<string>*, vector<double>*);

  list<string> filenamesList;
  vector<string> varNames;
  
  
  PSO_Problem(unsigned int nbVar, double varMin, double varMax,
	      double (* f)(const vector<double>&, size_t, const list<string>*,  vector<double>*),
	      const list<string>* p_filenamesList = 0, const vector<string>* p_varNames = 0)
  {
    this->costFunction = f;
    this->nbVar = nbVar; //Number of Unknown (Decision) Variables
    
    this->varMin.resize(nbVar, varMin); //Lower Bound of Decision Variables
    this->varMax.resize(nbVar, varMax); //Upper Bound of Decision Variables

    
    if(p_varNames)
      {
	this->varNames = *p_varNames;
      }
    
    if(p_filenamesList)
      {
	filenamesList = *p_filenamesList;
      }
  }

  PSO_Problem(const vector<double>& varMin, const vector<double>& varMax,
	      double (* f)(const vector<double>&, size_t, const list<string>*, vector<double>*),
	      const list<string>* p_filenamesList = 0, const vector<string>* p_varNames = 0)
  {
    this->costFunction = f;
    this->nbVar = varMin.size(); //Number of Unknown (Decision) Variables
    
    this->varMin = varMin; //Lower Bound of Decision Variables
    this->varMax = varMax; //Upper Bound of Decision Variables
    if(p_varNames)
      {
	this->varNames = *p_varNames;
      }
    
    if(p_filenamesList)
      {
	filenamesList = *p_filenamesList;
      }
  }

  
};

struct PSO_Parameters
{
  size_t maxIterations; //Maximum Number of Iterations
  size_t nbPop; //Population Size (Swarm Size)
  double w; //Inertia Coefficient
  double wdamp; //Damping Ration of Inertia Coefficient
  double c1; //Personal Acceleration Coefficient
  double c2; //Social Acceleration Coefficient
  
  bool showIterInfo;
  
  PSO_Parameters(size_t maxIterations, size_t nbPop, double w, double wdamp, double c1, double c2,
		 bool showIterInfo = false)
  {
    this->maxIterations = maxIterations; //Maximum Number of Iterations
    this->nbPop = nbPop; //Population Size (Swarm Size)
    this->w = w; //Inertia Coefficient
    this->wdamp = wdamp; //Damping Ration of Inertia Coefficient
    this->c1 = c1; //Personal Acceleration Coefficient
    this->c2 = c2; //Social Acceleration Coefficient

    this->showIterInfo = showIterInfo;
  }
};


struct PSO_Out
{
  vector<double> global_best_solution;
  vector<double> best_costs;
  vector<Particle> pop;

  PSO_Out() {}
  
  PSO_Out(const vector<double>& gbp, const vector<double>& gbc)
  {
    global_best_solution = gbp;
    best_costs = gbc;
  }

  PSO_Out(const vector<double>& gbp, const vector<double>& gbc, const vector<Particle>& v_particles)
  {
    global_best_solution = gbp;
    best_costs = gbc;
    pop = v_particles;
  }
  
};


void compute_one_initialization_of_pso(size_t j, unsigned int nbVar,
				       vector<Particle>& v_particles,
				       const vector<double>& varMin,
				       const vector<double>& varMax,
				       const vector<double>& minVelocity,
				       const vector<double>& maxVelocity,
				       const list<string>* filenamesList,
				       double (& costFunction) (const vector<double>&, size_t, const list<string>*,  vector<double>*),
				       const vector<string>& varNames,
				       vector<double>& v_res)
{
  assert(0); //Actually with no use when ml
  
  bool is_init = false;
  if(j==0)
    {
      //Generate an individual that we think close to the best position
      vector<double> v_human; v_human.resize(nbVar);
      assert(v_human.size()==varMin.size());
      
      //TODO filename a mettre dans filenamesList
      string particle0_filename = "particle0_parameter.txt";
      ifstream flux(particle0_filename);
      
      if((flux.is_open())&&(varNames.size()==varMin.size()))
	{
	  while(not flux.eof())
	    {
	      string line;
	      std::getline(flux, line);
	      
	      stringtools_replace(line, " ", "");
	      
	      string res1, val; 
	      stringtools_split(line, res1, val, '=');
	      
	      //stringtools_upper();
	      
	      if(res1=="")
		{ //BREAK : Attention! Un saut de ligne dans le fichier de paramètres inhibe l'initialisation
		  //des paramètres qui seraient dans les lignes suivantes
		  break;
		}
	      
	      string varname = res1;
	      //cerr<<varname<<endl;
	      
	      //cerr<<">> "<<line<<endl;
	      
	      size_t index = 0;
	      string varnameup = stringtools_upper(varname);
	      
	      for(size_t i = 0; i<varNames.size(); i++)
		{
		  string temp_varname =  varNames[i];
		  temp_varname = stringtools_upper(temp_varname);
		  
		  if(temp_varname==varnameup)
		    {
		      index = i;
		      break;
		    }
		}
	      
	      //cerr<<"@@ "<<varnameup<<" - "<<index<<endl;
	      
	      assert(index<v_human.size());
	      
	      double v; stringtools_toDouble(val, v);
	      
	      v_human[index] = v;
	    }
	  
	  v_particles[j].init_human_position(v_human);
	  is_init = true;
	  
	  
	  flux.close();
	}
      else
	{
	  if(flux.is_open())
	    {
	      flux.close();
	    }
	  else
	    {
	      cerr<<"*** ERROR *** "<<particle0_filename<<" not found"<<endl;
	    }
	  
	  cerr<<"*** ERROR *** "<<"names of variables are not correctly transmitted"<<endl; 
										      
	  //assert(0);
	}
    }
      
  if(not is_init)
    {
      //Generate Random Solution
      v_particles[j].init_random_position(nbVar,varMin, varMax);
    }
  
  //Initialize Velocity
  //v_particles[j].init_velocity(nbVar);
  v_particles[j].init_random_velocity(nbVar, minVelocity, maxVelocity); 
  
  //Evaluation
  v_particles[j].eval(costFunction, j, filenamesList, &v_res);
  
  //Update the Personal Best
  v_particles[j].best_position = v_particles[j].position;
  v_particles[j].best_cost = v_particles[j].cost;
}



void compute_one_individual_of_pso_with_ml(size_t i, size_t j, unsigned int nbVar, double c1, double c2, double w,
					   const string& mlpredict_filename, vector<Particle>& v_particles,
					   const vector<double>& global_best_position, const vector<double>& minVelocity,
					   const vector<double>& maxVelocity, const vector<double>& varMin, const vector<double>& varMax,
					   const list<string>* filenamesList, double (& costFunction) (const vector<double>&, size_t, const list<string>*, vector<double>*),
					   vector<double>& v_res)
{
  if(j==0)
    {
      v_particles[j].init_position(nbVar);

      
      size_t last_found = 0;
      size_t found = 0;
      
      //-0.04694597890955371,0.8521626808962169,-0.08025520954363234,-0.0060275342811020125,5.593586573754619,33.909558629418996,0.013135351048673126
      ifstream predict_flux(mlpredict_filename);
      if(not predict_flux.is_open())
	{
	  cerr<<"*** ERROR *** prediction file not found"<<endl;
	  assert(0);
	}
      
      string s_predict = std::string((std::istreambuf_iterator<char>(predict_flux)),
				     std::istreambuf_iterator<char>());
      predict_flux.close();

      
      cerr<<"the prediction is : "<<s_predict<<endl;
      
      size_t num_data = 1;
      while(found!=std::string::npos)
	{
	  //cerr<<"@0 "<<last_found<<" "<<s_predict<<" "<<s_predict.size()<<endl;
	  
	  
	  found = s_predict.find_first_of(",", last_found);
	  
	  
	  string s_data = s_predict.substr(last_found, found-last_found); 
	  double d_data = std::stod(s_data);
	  
			     
	  //cerr<<"==> "<<j<<" "<<(v_particles[j].position).size()<<endl;
	  assert(j<v_particles.size());
	  assert(7<(v_particles[j].position).size());
	  
	  
	  if(num_data==1) v_particles[j].position[0] = d_data;
	  if(num_data==2) v_particles[j].position[1] = d_data;
	  if(num_data==3) v_particles[j].position[2] = d_data;
	  if(num_data==4) v_particles[j].position[4] = d_data;
	  if(num_data==5) v_particles[j].position[5] = d_data;
	  if(num_data==6) v_particles[j].position[6] = d_data;
	  if(num_data==7)
	    {
	      v_particles[j].position[7] = d_data;
	      break;
	    }
	  
	  num_data++;
	  last_found = found+1;
	}
      
      v_particles[j].position[3] = 6; //varMin[3];
      
      cerr<<"[ ML ] - ITERATION "<<i+1<<" INDIVIDUAL "<<j<<" ";
      for(size_t kk = 0; kk<v_particles[j].position.size(); kk++)
	{
	  cerr<<v_particles[j].position[kk]<<" ";
	}
      cerr<<endl;
      
      v_particles[j].init_velocity(nbVar);
    }
  else
    { //j>0
      if(i==0)
	{
	  v_particles[j].init_random_position(nbVar,varMin, varMax);
	  v_particles[j].init_random_velocity(nbVar,minVelocity, maxVelocity);
	}
      else
	{ //i>0 et j>0
	  assert(v_particles[j].velocity.size()>0);
	  
	  //Update velocity
	  vector<double> v_r1; v_r1.resize(nbVar, 0);
	  vector<double> v_r2; v_r2.resize(nbVar, 0);
	  for(size_t h = 0; h<v_r1.size(); h++)
	    {
	      v_r1[h] = random_probability();
	      v_r2[h] = random_probability();
	    }
	  
	  vector<double> v_cognitive; v_cognitive.resize(nbVar, 0);
	  for(size_t h = 0; h<v_r1.size(); h++)
	    {
	      assert(j<v_particles.size());
	      assert(h<v_particles[j].best_position.size());
	      assert(h<v_particles[j].position.size());
	      
	      v_cognitive[h] =
		c1*v_r1[h]*(v_particles[j].best_position[h]-v_particles[j].position[h]);
	    }
	  
	  vector<double> v_social; v_social.resize(nbVar, 0);
	  for(size_t h = 0; h<v_r1.size(); h++)
	    {
	      v_social[h] =
		c2*v_r2[h]*(global_best_position[h]-v_particles[j].position[h]);
	    }
	  
	  
	  for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
	    {
	      assert(v_social.size()==v_particles[j].velocity.size());
	      v_particles[j].velocity[k] = w*v_particles[j].velocity[k]+
		v_cognitive[k]+v_social[k];
	    }
	  
	  //Apply velocity Limits
	  for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
	    {
	      if(v_particles[j].velocity[k]<minVelocity[k])
		v_particles[j].velocity[k] = minVelocity[k];
	      
	      if(v_particles[j].velocity[k]>maxVelocity[k])
		v_particles[j].velocity[k] = maxVelocity[k]; 
	    }

	  
	  assert(j<v_particles.size());
	  //Update position
	  for(size_t k = 0; k<v_particles[j].position.size(); k++)
	    {
	      //cerr<<j<<" "<<k<<" - velocitysize : "<<v_particles[j].velocity.size()<<endl;
	      assert(k<v_particles[j].position.size());
	      assert(k<v_particles[j].velocity.size());
	      
	      v_particles[j].position[k] =
		v_particles[j].position[k]+v_particles[j].velocity[k];
	    }
	}
    }
  
  
  /*
    if(i>0)
    {
      assert(v_particles[j].velocity.size()>0);
      
      //Update velocity
      vector<double> v_r1; v_r1.resize(nbVar, 0);
      vector<double> v_r2; v_r2.resize(nbVar, 0);
      for(size_t h = 0; h<v_r1.size(); h++)
	{
	  v_r1[h] = random_probability();
	  v_r2[h] = random_probability();
	}
      
      vector<double> v_cognitive; v_cognitive.resize(nbVar, 0);
      for(size_t h = 0; h<v_r1.size(); h++)
	{
	  assert(j<v_particles.size());
	  assert(h<v_particles[j].best_position.size());
	  assert(h<v_particles[j].position.size());
	  
	  v_cognitive[h] =
	    c1*v_r1[h]*(v_particles[j].best_position[h]-v_particles[j].position[h]);
	}
      
      vector<double> v_social; v_social.resize(nbVar, 0);
      for(size_t h = 0; h<v_r1.size(); h++)
	{
	  v_social[h] =
	    c2*v_r2[h]*(global_best_position[h]-v_particles[j].position[h]);
	}
      
      
      for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
	{
	  assert(v_social.size()==v_particles[j].velocity.size());
	  v_particles[j].velocity[k] = w*v_particles[j].velocity[k]+
	    v_cognitive[k]+v_social[k];
	}
      
      //Apply velocity Limits
      for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
	{
	  if(v_particles[j].velocity[k]<minVelocity[k])
	    v_particles[j].velocity[k] = minVelocity[k];
	  
	  if(v_particles[j].velocity[k]>maxVelocity[k])
	    v_particles[j].velocity[k] = maxVelocity[k]; 
	}
    }
  else //i==0
    {
      v_particles[j].init_random_position(nbVar,varMin, varMax);
      v_particles[j].init_random_velocity(nbVar,varMinVelocity, varMaxVelocity); //TODO
    }

      
  
  if(j==0)
    {
      
      size_t last_found = 0;
      size_t found = 0;
      
      //-0.04694597890955371,0.8521626808962169,-0.08025520954363234,-0.0060275342811020125,5.593586573754619,33.909558629418996,0.013135351048673126
      ifstream predict_flux(mlpredict_filename);
      if(not predict_flux.is_open())
	{
	  cerr<<"*** ERROR *** prediction file not found"<<endl;
	  assert(0);
	}
      
      string s_predict = std::string((std::istreambuf_iterator<char>(predict_flux)),
				     std::istreambuf_iterator<char>());
      predict_flux.close();


      cerr<<"the prediction is : "<<s_predict<<endl;
      
      size_t num_data = 1;
      while(found!=std::string::npos)
	{
	  //cerr<<"@0 "<<last_found<<" "<<s_predict<<" "<<s_predict.size()<<endl;

	  
	  found = s_predict.find_first_of(",", last_found);

	  
	  string s_data = s_predict.substr(last_found, found-last_found); 
	  double d_data = std::stod(s_data);

			     
	  
	  assert(j<v_particles.size());
	  assert(7<(v_particles[j].position).size());
										       
	   
	  if(num_data==1) v_particles[j].position[0] = d_data;
	  if(num_data==2) v_particles[j].position[1] = d_data;
	  if(num_data==3) v_particles[j].position[2] = d_data;
	  if(num_data==4) v_particles[j].position[4] = d_data;
	  if(num_data==5) v_particles[j].position[5] = d_data;
	  if(num_data==6) v_particles[j].position[6] = d_data;
	  if(num_data==7)
	    {
	      v_particles[j].position[7] = d_data;
	      break;
	    }

	  num_data++;
	  last_found = found+1;
	}
      
      v_particles[j].position[3] = 6; //varMin[3];

      cerr<<"[ ML ] - ITERATION "<<i+1<<" INDIVIDUAL "<<j<<" ";
      for(size_t kk = 0; kk<v_particles[j].position.size(); kk++)
	{
	  cerr<<v_particles[j].position[kk]<<" ";
	}
      cerr<<endl;

    }
  else
    { 
      
      if(i==0)
	{
	  assert(j<v_particles.size());
	  
	  //Generate Random Solution
	  v_particles[j].init_random_position(nbVar,varMin, varMax);
	}
      else
	{
	  //cerr<<"step 8"<<endl;

	  assert(j<v_particles.size());
	  //Update position
	  for(size_t k = 0; k<v_particles[j].position.size(); k++)
	    {
	      cerr<<j<<" "<<k<<" - velocitysize : "<<v_particles[j].velocity.size()<<endl;
	      assert(k<v_particles[j].position.size());
	      assert(k<v_particles[j].velocity.size());
	      
	      v_particles[j].position[k] =
		v_particles[j].position[k]+v_particles[j].velocity[k];
	    }
	}
    }
  */
  
  assert(j<v_particles.size());
  //Apply Lower and Upper Bound Limits
  for(size_t k = 0; k<v_particles[j].position.size(); k++)
    {
      assert(k<varMin.size());
      assert(k<varMax.size());
      
      if(v_particles[j].position[k]<varMin[k])
	v_particles[j].position[k] = varMin[k];
      
      if(v_particles[j].position[k]>varMax[k])
	v_particles[j].position[k] = varMax[k]; 
    }

  cerr<<"evaluation "<<endl;
  //Evaluation
  //v_particles[j].eval();
  v_particles[j].eval(costFunction, j, filenamesList, &v_res);
  //cerr<<v_particles[j].cost<<endl;
  
  if(i==0)
    {
      //Update the Personal Best
      v_particles[j].best_position = v_particles[j].position;
      v_particles[j].best_cost = v_particles[j].cost;
    }
}


/*
void compute_one_individual_of_pso(size_t j, unsigned int nbVar, double c1, double c2, double w,
				   vector<Particle>& v_particles,
				   const vector<double>& global_best_position, const vector<double>& minVelocity,
				   const vector<double>& maxVelocity, const vector<double>& varMin, const vector<double>& varMax,
				   const list<string>* filenamesList, double (& costFunction) (const vector<double>&, size_t, const list<string>*, vector<double>*),
				   vector<double>& v_res)
{
  assert(0);
  
  //Update velocity
  vector<double> v_r1; v_r1.resize(nbVar, 0);
  vector<double> v_r2; v_r2.resize(nbVar, 0);
  for(size_t h = 0; h<v_r1.size(); h++)
    {
      v_r1[h] = random_probability();
      v_r2[h] = random_probability();
    }
 
  vector<double> v_cognitive; v_cognitive.resize(nbVar, 0);
  for(size_t h = 0; h<v_r1.size(); h++)
    {
      v_cognitive[h] =
	c1*v_r1[h]*(v_particles[j].best_position[h]-v_particles[j].position[h]);
    }

  vector<double> v_social; v_social.resize(nbVar, 0);
  for(size_t h = 0; h<v_r1.size(); h++)
    {
      v_social[h] =
	c2*v_r2[h]*(global_best_position[h]-v_particles[j].position[h]);
    }
  
  for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
    {
      assert(v_social.size()==v_particles[j].velocity.size());
      v_particles[j].velocity[k] = w*v_particles[j].velocity[k]+
	v_cognitive[k]+v_social[k];
    }

  
  //Apply velocity Limits
  for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
    {
      if(v_particles[j].velocity[k]<minVelocity[k])
	v_particles[j].velocity[k] = minVelocity[k];
      
      if(v_particles[j].velocity[k]>maxVelocity[k])
	v_particles[j].velocity[k] = maxVelocity[k]; 
    }

  //Update position
  for(size_t k = 0; k<v_particles[j].position.size(); k++)
    v_particles[j].position[k] =
      v_particles[j].position[k]+v_particles[j].velocity[k];
  
 
  //Apply Lower and Upper Bound Limits
  for(size_t k = 0; k<v_particles[j].position.size(); k++)
    {
      if(v_particles[j].position[k]<varMin[k])
	v_particles[j].position[k] = varMin[k];
      
      if(v_particles[j].position[k]>varMax[k])
	v_particles[j].position[k] = varMax[k]; 
    }
  
  //Evaluation
  //v_particles[j].eval();
  v_particles[j].eval(costFunction, j, filenamesList, &v_res);
  //cerr<<v_particles[j].cost<<endl;
}
*/


void multithread_pso_with_ml(const PSO_Problem& problem, const PSO_Parameters& params, PSO_Out& out)
{
  const string flagfile = "mlfinish";
  const string results_filename = "results.csv";
  const string mldatasets_filename = "ML_data_sets_no_out_norm.csv";
  const string mlpredict_filename = "predicted_out.txt";
  string s_launchml = "sbatch scriptml.sh";
  
  //Lancement du reseau de neurones
  cerr<<"[ MULTITHREADED PSO ][ MACHINE LEARNING ] in progress"<<endl; 
  system(s_launchml.c_str());

  
  if(true) 
    { //Machine-learning

      //Penser a changer dans scriptml (touch mlfinish)
      
      string flagpath = flagfile;
      
      size_t time_sleep = 30;
      size_t num_trial = 1;
      while(not filetools_exist(flagpath))
	{
	  //cerr<<"waiting end of experiment from "<<num_trial*time_sleep<<" seconds"<<endl;
	  num_trial++;
	  sleep(time_sleep);
	}
      
      string rm_flag = "rm -f "+flagpath;
      system(rm_flag.c_str());
      sleep(1);

      cerr<<"Machine Learning finished"<<endl;
    }
  
  
  
   //Problem definition
  assert(problem.costFunction);
  double (& costFunction) (const vector<double>&, size_t, const list<string>*, vector<double>*) = *(problem.costFunction);
  
  unsigned int nbVar = problem.nbVar; //Number of Unknown (Decision) Variables
  vector<unsigned int> varSize; //Matrix Size of Decision Variables
  varSize.resize(2); varSize[0] = 1; varSize[1] = problem.nbVar;
  
  vector<double> varMin = problem.varMin; //Lower Bound of Decision Variables
  vector<double> varMax = problem.varMax; //Upper Bound of Decision Variables
  vector<string> varNames = problem.varNames;
  
  const list<string>* filenamesList = &(problem.filenamesList);
  
  //Parameters of PSO
  size_t maxIterations = params.maxIterations; //Maximum Number of Iterations
  size_t nbPop = params.nbPop; //Population Size (Swarm Size)
  double w = params.w; //Inertia Coefficient
  double wdamp = params.wdamp; //Damping Ration of Inertia Coefficient
  double c1 = params.c1; //Personal Acceleration Coefficient
  double c2 = params.c2; //Social Acceleration Coefficient

  vector<double> maxVelocity; //= 0.2*(varMax-varMin);
  vector<double> minVelocity; //= -maxVelocity;
  maxVelocity.resize(varMax.size(), 0);
  for(size_t i = 0; i<maxVelocity.size(); i++)
    {
      maxVelocity[i] = 0.2*(varMax[i]-varMin[i]);
    }
  
  minVelocity = maxVelocity;
  for(size_t i = 0; i<minVelocity.size(); i++)
    {
      minVelocity[i] = -minVelocity[i];
    }
  
  
  bool showIterInfo = params.showIterInfo; //The flag for Showing Iteration Information
  

  //Initialization

  // Create Population Array
  vector<Particle> v_particles; v_particles.resize(nbPop, Particle());

  // Initialize Global Best
  vector<double> global_best_position;
  double global_best_cost = -1; //Inf

  vector<double> v_poub;
  vector<vector<double> > v_results; v_results.resize(v_particles.size(), v_poub);
  for(size_t i = 0; i<v_results.size(); i++)
    {
      v_results[i].resize(2+varMin.size()+1+3+2+1, 0);
    }
  
  ofstream flux_init(results_filename, ios::app);
  if(not flux_init.is_open())
    {
      cerr<<"*** ERROR ***"<<results_filename<<" not found"<<endl;
      assert(0);
    }
  
  flux_init<<"Iteration,"<<"#Pop,";
  for(size_t i = 0; i<varNames.size(); i++)
    {
      flux_init<<varNames[i]<<",";		   
    }
  
  flux_init<<"temperature,"<<"density,"<<"viscosity1,"<<"viscosity2,"<<"std1,"<<"std2,"<<"eval"<<","<<"bestCost"<<"\n";
  flux_init.close();
  

  /*
  // Initialize Population Members
  {
    cerr<<"[ MULTITHREADED PSO ] - INITIALIZATION - #POP = "<<v_particles.size()<<endl; 
    
    vector<std::thread> v_threads0;
    for(size_t j = 0; j<v_particles.size(); j++)
      {
	//cerr<<"initialization - individu "<<j<<endl;
	v_threads0.push_back(std::thread(compute_one_initialization_of_pso, j, nbVar, std::ref(v_particles),
	varMin, varMax, minVelocity, maxVelocity,
	filenamesList,
	std::ref(costFunction),
	varNames, std::ref(v_results[j])));
      }
    
    
    //Synchronisation : Attente que tous les threads aient termine leur job
    for(size_t j = 0; j<v_threads0.size(); j++)
      {
	v_threads0[j].join();
      }
    
    cerr<<"iteration 0 : "<<endl;
   
   
    
    
    for(size_t j = 0; j<v_particles.size(); j++)
      {
	//Update Global Best
	if((global_best_position.size()==0)||(v_particles[j].best_cost<global_best_cost))
	  {
	    global_best_cost = v_particles[j].best_cost;
	    global_best_position = v_particles[j].best_position;
	  }
	v_particles[j].print(cerr);
      } 



    cerr<<"add new ml data (0)"<<endl;
    //Iteration,#Pop,Epsiam,Sigma,Epsicb,-dCBD-,RhoCBD,Kn,Gamma,Xu,temperature,density,viscosity1,viscosity2,std1,std2,eval,bestCost
    //dCBD a ete fixe comme une constante et n'est donc pas integre dans les datasets pour le Machine-Learning
    //2 3 4 6 7 8 9 11 12 13
    ofstream mldata_flux0(mldatasets_filename, ios::app);
    if(not mldata_flux0.is_open())
      {
	cerr<<"*** ERROR ***"<<mldatasets_filename<<" not found"<<endl;
	assert(0);
      }
    
    for(size_t j = 0; j<v_results.size(); j++)
      {
	
	for(size_t k = 0; k<v_results[j].size(); k++)
	  {
	    if(k==0) continue;
	    if(k==1) continue;
	    if(k==5) continue;
	    if(k==10) continue;
	    if(k>=14) break;
	    
	    cerr<<v_results[j][k]<<" ";
	    
	    mldata_flux0<<v_results[j][k];
	    mldata_flux0<<",";
	    //if(k<v_results[j].size()-1)
	    //flux<<",";
	    //else flux<<"\n";
	  }
      }
    
    mldata_flux0.close();

    //Lancement du reseau de neurones
    
    
    ofstream flux0(results_filename, ios::app);
    if(not flux0.is_open())
      {
	cerr<<"*** ERROR ***"<<results_filename<<" not found"<<endl;
	assert(0);
      }
    
    for(size_t j = 0; j<v_results.size(); j++)
      {
	for(size_t k = 0; k<v_results[j].size(); k++)
	  {
	    cerr<<v_results[j][k]<<" ";
	    
	    flux0<<v_results[j][k];
	    flux0<<",";
	    //if(k<v_results[j].size()-1)
	    //flux0<<",";
	    //else flux0<<"\n";
	  }

	flux0<<global_best_cost<<"\n";
	cerr<<endl;
	
	for(size_t k = 0; k<v_results[j].size(); k++)
	  v_results[j][k] = 0;
      }
    
    flux0.close();
  }

  cerr<<"best position = ";
  for(size_t i = 0; i<global_best_position.size(); i++)
    {
      cerr<<global_best_position[i]<<" "; 
    }
  cerr<<"\n";
  
  //Penser a changer dans scriptml (touch mlfinish)
  if(true)
    {
      string flagpath = flagfile;
      
      size_t time_sleep = 30;
      size_t num_trial = 1;
      while(not filetools_exist(flagpath))
	{
	  //cerr<<"waiting end of experiment from "<<num_trial*time_sleep<<" seconds"<<endl;
	  num_trial++;
	  sleep(time_sleep);
	}
      
      string rm_flag = "rm -f "+flagpath;
      system(rm_flag.c_str());
      sleep(1);
    }
  
  cerr<<"best cost = "<<global_best_cost<<endl;
  */  

  vector<double> best_costs; best_costs.resize(maxIterations, 0);  

  
  //Main Loop of PSO
  for(size_t i = 0; i<maxIterations; i++)
    {
      //nbPop threads
      
      cerr<<"[ MULTITHREADED PSO ] - ITERATION "<<i+1<<"/"<<maxIterations<<" - #POP = "<<nbPop<<endl; 
      
      vector<std::thread> v_threads;
      for(size_t j = 0; j<nbPop; j++) 
	{
   
	  cerr<<"iteration "<<i+1<<" - "<<"individu "<<j<<endl;
	  v_threads.push_back(std::thread(compute_one_individual_of_pso_with_ml, i, j, nbVar, c1, c2,  w, mlpredict_filename,std::ref(v_particles),
					  global_best_position, minVelocity,
					  maxVelocity, varMin,  varMax,
					  filenamesList, std::ref(costFunction),
					  std::ref(v_results[j])));
	  
	  //sleep(10);
	}

      for(int j = 0; j<(int)v_threads.size(); j++)
	{
	  v_threads[j].join();
	}
      
      for(size_t j = 0; j<nbPop; j++)
	{
	  //Update the Personal Best
	  if((i==0)||(v_particles[j].cost<v_particles[j].best_cost))
	    {
	      v_particles[j].best_position = v_particles[j].position;
	      v_particles[j].best_cost = v_particles[j].cost;
	      
	      //Update Global Best
	      if((global_best_position.size()==0)||(v_particles[j].best_cost<global_best_cost))
	      {
		global_best_cost = v_particles[j].best_cost;
		global_best_position = v_particles[j].best_position;
	      }
	    }
	}

      
      cerr<<"iteration : "<<i+1<<endl;

      cerr<<"add new ml data"<<endl;
      //Iteration,#Pop,Epsiam,Sigma,Epsicb,-dCBD-,RhoCBD,Kn,Gamma,Xu,temperature,density,viscosity1,viscosity2,std1,std2,eval,bestCost
      //dCBD a ete fixe comme une constante et n'est donc pas integre dans les datasets pour le Machine-Learning
      //2 3 4 6 7 8 9 11 12 13
      ofstream mldata_flux(mldatasets_filename, ios::app);
      if(not mldata_flux.is_open())
	{
	  cerr<<"*** ERROR ***"<<mldatasets_filename<<" not found"<<endl;
	  assert(0);
	}

      for(size_t j = 0; j<v_results.size(); j++)
	{

	  if(has_failed(v_results[j]))
	    {
	      continue;
	    }

	  
	  for(size_t k = 0; k<v_results[j].size(); k++)
	    {
	      if(k==0) continue;
	      if(k==1) continue;
	      if(k==5) continue;
	      if(k==10) continue;
	      if(k>=14) break;

	      //k==11 k==12 k==13

	      double v_to_add = v_results[j][k];
	      
	      if(k==11)
		{
		  v_to_add = v_results[j][k]/DENSITY_EXP;
		}
	      
	      if(k==12)
		{
		  v_to_add = v_results[j][k]/VISCOSITY1_EXP;
		}
	      
	      if(k==13)
		{
		  v_to_add = v_results[j][k]/VISCOSITY2_EXP;
		}
	      
	      
	      //cerr<<v_results[j][k]<<" ";
	      cerr<<v_to_add;
	      if(v_to_add!=v_results[j][k])
		cerr<<"("<<v_results[j][k]<<")";
	      cerr<<" ";
	      
	      
	      mldata_flux<<v_to_add; 
	      mldata_flux<<",";
	      //if(k<v_results[j].size()-1)
	      //flux<<",";
	      //else flux<<"\n";
	    }

	  mldata_flux<<"\n";
	}

      mldata_flux.close();
      
      //Lancement du reseau de neurones
      system(s_launchml.c_str());

      
      ofstream flux(results_filename, ios::app);
      if(not flux.is_open())
	{
	  cerr<<"*** ERROR ***"<<results_filename<<" not found"<<endl;
	  assert(0);
	}
      
      for(size_t j = 0; j<v_results.size(); j++)
	{
	  v_results[j][0] = i+1;
	  
	  for(size_t k = 0; k<v_results[j].size(); k++)
	    {
	      cerr<<v_results[j][k]<<" ";
	      
	      flux<<v_results[j][k];
	      flux<<",";
	      //if(k<v_results[j].size()-1)
	      //flux<<",";
	      //else flux<<"\n";
	    }


	  flux<<global_best_cost<<"\n";
	  cerr<<endl;
	  
	  for(size_t k = 0; k<v_results[j].size(); k++)
	    v_results[j][k] = 0;
	}

      flux.close();
      
      
      //Store the Best Cost Value
      best_costs[i] = global_best_cost;
      

      //Display Iteration information
      if(showIterInfo)
	cerr<<"iteration = "<<i+1<<" - best cost = "<<best_costs[i]<<endl;
      
      //Damping Inertia Coefficient
      w = w*wdamp;

      
     
      //Penser a changer dans scriptml (touch mlfinish)
   
      string flagpath = flagfile;

      size_t time_sleep = 30;
      size_t num_trial = 1;
      while(not filetools_exist(flagpath))
	{
	  //cerr<<"waiting end of experiment from "<<num_trial*time_sleep<<" seconds"<<endl;
	  num_trial++;
	  sleep(time_sleep);
	}

      string rm_flag = "rm -f "+flagpath;
      system(rm_flag.c_str());
      sleep(1);
      
    }

  out = PSO_Out(global_best_position, best_costs, v_particles);
  //out.global_best_position = global_best_position;
  //out.global_best_cost = global_best_cost;
}


/*
void multithread_pso(const PSO_Problem& problem, const PSO_Parameters& params, PSO_Out& out)
{
  assert(0);
  
   //Problem definition
  assert(problem.costFunction);
  double (& costFunction) (const vector<double>&, size_t, const list<string>*, vector<double>*) = *(problem.costFunction);
  
  unsigned int nbVar = problem.nbVar; //Number of Unknown (Decision) Variables
  vector<unsigned int> varSize; //Matrix Size of Decision Variables
  varSize.resize(2); varSize[0] = 1; varSize[1] = problem.nbVar;
  
  vector<double> varMin = problem.varMin; //Lower Bound of Decision Variables
  vector<double> varMax = problem.varMax; //Upper Bound of Decision Variables
  vector<string> varNames = problem.varNames;
  
  const list<string>* filenamesList = &(problem.filenamesList);
  
  //Parameters of PSO
  size_t maxIterations = params.maxIterations; //Maximum Number of Iterations
  size_t nbPop = params.nbPop; //Population Size (Swarm Size)
  double w = params.w; //Inertia Coefficient
  double wdamp = params.wdamp; //Damping Ration of Inertia Coefficient
  double c1 = params.c1; //Personal Acceleration Coefficient
  double c2 = params.c2; //Social Acceleration Coefficient

  vector<double> maxVelocity; //= 0.2*(varMax-varMin);
  vector<double> minVelocity; //= -maxVelocity;
  maxVelocity.resize(varMax.size(), 0);
  for(size_t i = 0; i<maxVelocity.size(); i++)
    {
      maxVelocity[i] = 0.2*(varMax[i]-varMin[i]);
    }
  
  minVelocity = maxVelocity;
  for(size_t i = 0; i<minVelocity.size(); i++)
    {
      minVelocity[i] = -minVelocity[i];
    }
  
  
  bool showIterInfo = params.showIterInfo; //The flag for Showing Iteration Information
  

  //Initialization

  // Create Population Array
  vector<Particle> v_particles; v_particles.resize(nbPop, Particle());

  // Initialize Global Best
  vector<double> global_best_position;
  double global_best_cost = -1; //Inf

  vector<double> v_poub;
  vector<vector<double> > v_results; v_results.resize(v_particles.size(), v_poub);
  for(size_t i = 0; i<v_results.size(); i++)
    {
      v_results[i].resize(2+varMin.size()+1+3+2+1, 0);
    }

  string results_filename = "results.csv";
  ofstream flux_init(results_filename, ios::app);
  if(not flux_init.is_open())
    {
      cerr<<"*** ERROR ***"<<results_filename<<" not found"<<endl;
      assert(0);
    }
  
  flux_init<<"Iteration,"<<"#Pop,";
  for(size_t i = 0; i<varNames.size(); i++)
    {
      flux_init<<varNames[i]<<",";		   
    }
  
  flux_init<<"temperature,"<<"density,"<<"viscosity1,"<<"viscosity2,"<<"std1,"<<"std2,"<<"eval"<<","<<"bestCost"<<"\n";
  flux_init.close();
  
  
  // Initialize Population Members
  {
    cerr<<"[ MULTITHREADED PSO ] - INITIALIZATION - #POP = "<<v_particles.size()<<endl; 
    
    vector<std::thread> v_threads0;
    for(size_t j = 0; j<v_particles.size(); j++)
      {
	//cerr<<"initialization - individu "<<j<<endl;
	v_threads0.push_back(std::thread(compute_one_initialization_of_pso, j, nbVar, std::ref(v_particles),
					 varMin, varMax, minVelocity, maxVelocity,
				       filenamesList,
					 std::ref(costFunction),
					 varNames, std::ref(v_results[j])));
      }
    
    
    //Synchronisation : Attente que tous les threads aient termine leur job
    for(size_t j = 0; j<v_threads0.size(); j++)
      {
	v_threads0[j].join();
      }
    
    cerr<<"iteration 0 : "<<endl;
   
   
    
    
    for(size_t j = 0; j<v_particles.size(); j++)
      {
	//Update Global Best
	if((global_best_position.size()==0)||(v_particles[j].best_cost<global_best_cost))
	  {
	    global_best_cost = v_particles[j].best_cost;
	    global_best_position = v_particles[j].best_position;
	  }
	v_particles[j].print(cerr);
      } 

     ofstream flux0(results_filename, ios::app);
    if(not flux0.is_open())
      {
	cerr<<"*** ERROR ***"<<results_filename<<" not found"<<endl;
	assert(0);
      }
    
    for(size_t j = 0; j<v_results.size(); j++)
      {
	for(size_t k = 0; k<v_results[j].size(); k++)
	  {
	    cerr<<v_results[j][k]<<" ";
	    
	    flux0<<v_results[j][k];
	    flux0<<",";
	    //if(k<v_results[j].size()-1)
	    //flux0<<",";
	    //else flux0<<"\n";
	  }

	flux0<<global_best_cost<<"\n";
	cerr<<endl;
	
	for(size_t k = 0; k<v_results[j].size(); k++)
	  v_results[j][k] = 0;
      }
    flux0.close();
    
    }

  cerr<<"best position = ";
  for(size_t i = 0; i<global_best_position.size(); i++)
    {
      cerr<<global_best_position[i]<<" "; 
    }
  cerr<<"\n";

  cerr<<"best cost = "<<global_best_cost<<endl;

  vector<double> best_costs; best_costs.resize(maxIterations, 0);  

  
  //Main Loop of PSO
  for(size_t i = 0; i<maxIterations; i++)
    {
      //nbPop threads
      
      cerr<<"[ MULTITHREADED PSO ] - ITERATION "<<i+1<<"/"<<maxIterations<<" - #POP = "<<nbPop<<endl; 
      
      vector<std::thread> v_threads;
      for(size_t j = 0; j<nbPop; j++) 
	{
	  
	  //cerr<<"iteration "<<i+1<<" - "<<"individu "<<j<<endl;
	  v_threads.push_back(std::thread(compute_one_individual_of_pso, j, nbVar, c1, c2,  w, std::ref(v_particles),
					  global_best_position, minVelocity,
					  maxVelocity, varMin,  varMax,
					  filenamesList, std::ref(costFunction),
					  std::ref(v_results[j])));

	}

      for(int j = 0; j<(int)v_threads.size(); j++)
	{
	  v_threads[j].join();
	}
      
      for(size_t j = 0; j<nbPop; j++)
	{
	  //Update the Personal Best
	  if(v_particles[j].cost<v_particles[j].best_cost)
	    {
	      v_particles[j].best_position = v_particles[j].position;
	      v_particles[j].best_cost = v_particles[j].cost;
	      
	      //Update Global Best
	      if((global_best_position.size()==0)||(v_particles[j].best_cost<global_best_cost))
	      {
		global_best_cost = v_particles[j].best_cost;
		global_best_position = v_particles[j].best_position;
	      }
	    }
	}

      cerr<<"iteration : "<<i+1<<endl;
      ofstream flux(results_filename, ios::app);
      if(not flux.is_open())
	{
	  cerr<<"*** ERROR ***"<<results_filename<<" not found"<<endl;
	  assert(0);
	}
      
      for(size_t j = 0; j<v_results.size(); j++)
	{
	  v_results[j][0] = i+1;
	  
	  for(size_t k = 0; k<v_results[j].size(); k++)
	    {
	      cerr<<v_results[j][k]<<" ";
	      
	      flux<<v_results[j][k];
	      flux<<",";
	      //if(k<v_results[j].size()-1)
	      //flux<<",";
	      //else flux<<"\n";
	    }


	  flux<<global_best_cost<<"\n";
	  cerr<<endl;
	  
	  for(size_t k = 0; k<v_results[j].size(); k++)
	    v_results[j][k] = 0;
	}

      flux.close();
      
      //Store the Best Cost Value
      best_costs[i] = global_best_cost;
      

      //Display Iteration information
      if(showIterInfo)
	cerr<<"iteration = "<<i+1<<" - best cost = "<<best_costs[i]<<endl;
      
      //Damping Inertia Coefficient
      w = w*wdamp;
      
    }

  out = PSO_Out(global_best_position, best_costs, v_particles);
  //out.global_best_position = global_best_position;
  //out.global_best_cost = global_best_cost;
}
*/

 /*
void pso(const PSO_Problem& problem, const PSO_Parameters& params, PSO_Out& out)
{
   //Problem definition
  assert(problem.costFunction);
  double (& costFunction) (const vector<double>&, size_t, const list<string>*, vector<double>*) = *(problem.costFunction);
  
  unsigned int nbVar = problem.nbVar; //Number of Unknown (Decision) Variables
  vector<unsigned int> varSize; //Matrix Size of Decision Variables
  varSize.resize(2); varSize[0] = 1; varSize[1] = problem.nbVar;
  
  vector<double> varMin = problem.varMin; //Lower Bound of Decision Variables
  vector<double> varMax = problem.varMax; //Upper Bound of Decision Variables

  const list<string>* filenamesList = &(problem.filenamesList);
  
  //Parameters of PSO
  size_t maxIterations = params.maxIterations; //Maximum Number of Iterations
  size_t nbPop = params.nbPop; //Population Size (Swarm Size)
  double w = params.w; //Inertia Coefficient
  double wdamp = params.wdamp; //Damping Ration of Inertia Coefficient
  double c1 = params.c1; //Personal Acceleration Coefficient
  double c2 = params.c2; //Social Acceleration Coefficient

  vector<double> maxVelocity; //= 0.2*(varMax-varMin);
  vector<double> minVelocity; //= -maxVelocity;
  maxVelocity.resize(varMax.size(), 0);
  for(size_t i = 0; i<maxVelocity.size(); i++)
    {
      maxVelocity[i] = 0.2*(varMax[i]-varMin[i]);
    }
  
  minVelocity = maxVelocity;
  for(size_t i = 0; i<minVelocity.size(); i++)
    {
      minVelocity[i] = -minVelocity[i];
    }
  
  
  bool showIterInfo = params.showIterInfo; //The flag for Showing Iteration Information
  


  
  //Initialization

  // Create Population Array
  vector<Particle> v_particles; v_particles.resize(nbPop, Particle());

  // Initialize Global Best
  vector<double> global_best_position;
  double global_best_cost = -1; //Inf
  
  // Initialize Population Members
  for(size_t i = 0; i<v_particles.size(); i++)
    {
      //Generate Random Solution
      v_particles[i].init_random_position(nbVar,varMin, varMax);

      //Initialize Velocity
      v_particles[i].init_velocity(nbVar);
      
      //Evaluation
      vector<double> v_res;
      v_particles[i].eval(costFunction, i, filenamesList, &v_res);
      
      //Update the Personal Best
      v_particles[i].best_position = v_particles[i].position;
      v_particles[i].best_cost = v_particles[i].cost;

      //Update Global Best
      if((global_best_position.size()==0)||(v_particles[i].best_cost<global_best_cost))
	{
	  global_best_cost = v_particles[i].best_cost;
	  global_best_position = v_particles[i].best_position;
	}
      v_particles[i].print(cerr);
    }

  cerr<<"best position = ";
  for(size_t i = 0; i<global_best_position.size(); i++)
    {
      cerr<<global_best_position[i]<<" "; 
    }
  cerr<<"\n";

  cerr<<"best cost = "<<global_best_cost<<endl;

  vector<double> best_costs; best_costs.resize(maxIterations, 0);
  
  //Main Loop of PSO
  
  for(size_t i = 0; i<maxIterations; i++)
    {

      for(size_t j = 0; j<nbPop; j++) 
	{
	  //Update velocity
	  vector<double> v_r1; v_r1.resize(nbVar, 0);
	  vector<double> v_r2; v_r2.resize(nbVar, 0);
	  for(size_t h = 0; h<v_r1.size(); h++)
	    {
	      v_r1[h] = random_probability();
	      v_r2[h] = random_probability();
	    }


	  
	  vector<double> v_cognitive; v_cognitive.resize(nbVar, 0);
	  for(size_t h = 0; h<v_r1.size(); h++)
	    {
	      v_cognitive[h] =
		c1*v_r1[h]*(v_particles[j].best_position[h]-v_particles[j].position[h]);
	    }
	  
	  vector<double> v_social; v_social.resize(nbVar, 0);
	  for(size_t h = 0; h<v_r1.size(); h++)
	    {
	      v_social[h] =
		c2*v_r2[h]*(global_best_position[h]-v_particles[j].position[h]);
	    }

	  for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
	    {
	      assert(v_social.size()==v_particles[j].velocity.size());
	      v_particles[j].velocity[k] = w*v_particles[j].velocity[k]+
		v_cognitive[k]+v_social[k];
	    }
	  
	  //Apply velocity Limits
	  for(size_t k = 0; k<v_particles[j].velocity.size(); k++)
	    {
	      if(v_particles[j].velocity[k]<minVelocity[k])
		v_particles[j].velocity[k] = minVelocity[k];
	      
	      if(v_particles[j].velocity[k]>maxVelocity[k])
		v_particles[j].velocity[k] = maxVelocity[k]; 
	    }
	  
	  //Update position
	  for(size_t k = 0; k<v_particles[j].position.size(); k++)
	    v_particles[j].position[k] =
	      v_particles[j].position[k]+v_particles[j].velocity[k];


	  
	  //Apply Lower and Upper Bound Limits
	  for(size_t k = 0; k<v_particles[j].position.size(); k++)
	    {
	      if(v_particles[j].position[k]<varMin[k])
		v_particles[j].position[k] = varMin[k];
	      
	      if(v_particles[j].position[k]>varMax[k])
		v_particles[j].position[k] = varMax[k]; 
	    }
	 
	  
	  //Evaluation
	  //v_particles[j].eval();
	  vector<double> v_res;
	  v_particles[j].eval(costFunction, j, filenamesList, &v_res);
	  //cerr<<v_particles[j].cost<<endl;

	  //Update the Personal Best
	  if(v_particles[j].cost<v_particles[j].best_cost)
	    {
	      v_particles[j].best_position = v_particles[j].position;
	      v_particles[j].best_cost = v_particles[j].cost;
	   
	      
	      //Update Global Best
	      if((global_best_position.size()==0)||(v_particles[j].best_cost<global_best_cost))
		{
		  global_best_cost = v_particles[j].best_cost;
		  global_best_position = v_particles[j].best_position;
		}
	    }
	}


      //Store the Best Cost Value
      best_costs[i] = global_best_cost;


      //Display Iteration information
      if(showIterInfo)
	cerr<<"iteration = "<<i<<" - best cost = "<<best_costs[i]<<endl;

      //Damping Inertia Coefficient
      w = w*wdamp;
      
    }

  out = PSO_Out(global_best_position, best_costs, v_particles);
  //out.global_best_position = global_best_position;
  //out.global_best_cost = global_best_cost;
}
 */


struct Variable
{
  double upperBound;
  double lowerBound;
  string name;
  
  Variable(double lb = 0, double ub = 1, string n = "")
  {
    upperBound = ub;
    lowerBound = lb;
    name = n;
  }

  void setLowerBound(double lb)
  {
    lowerBound = lb;
  }
  
  void setUpperBound(double ub)
  {
    upperBound = ub;
  }

  string getName() const { return name; }

  string to_string() const
  {
    return name+" belongs to [ "+std::to_string(lowerBound)+" ; "+std::to_string(upperBound)+" ]";
  }
};


struct LRCSVariables
{
  //Variable dCBD;
  //Variable rhoCBD;
  //Variable kn;
  //Variable alpha;
  //Variable xu;
  //Variable gamma;
  //Variable beta;

  vector<Variable> variables;
  
  LRCSVariables()
  {
    /*
    variables.push_back(Variable(0, 1, "dCBD"));
    variables.push_back(Variable(0, 1, "RhoCBD"));
    variables.push_back(Variable(0, 1, "Kn"));
    variables.push_back(Variable(0, 1, "alpha"));
    variables.push_back(Variable(0, 1, "Xu"));
    variables.push_back(Variable(0, 1, "Gamma"));
    variables.push_back(Variable(0, 1, "Beta"));
    */
    variables.push_back(Variable(0, 1, "Epsiam"));
    variables.push_back(Variable(0, 1, "Sigma"));
    variables.push_back(Variable(0, 1, "Epsicb"));
    variables.push_back(Variable(0, 1, "dCBD"));
    variables.push_back(Variable(0, 1, "RhoCBD"));
    variables.push_back(Variable(0, 1, "Kn"));
    variables.push_back(Variable(0, 1, "Gamma"));
    variables.push_back(Variable(0, 1, "Xu"));
  }

  double getLowerBound(size_t iVar)
  {
    assert(iVar<variables.size());
    return variables[iVar].lowerBound;
  }

  double getUpperBound(size_t iVar)
  {
    assert(iVar<variables.size());
    return variables[iVar].upperBound;
  }

  string getName(size_t iVar)
  {
    assert(iVar<variables.size());
    return variables[iVar].getName();
  }
  

  void init(string& filename)
  { 
    ifstream flux(filename);
    if(not flux.is_open())
      {
	cerr<<"*** ERROR *** file "<<filename<<" not found"<<endl;
	return;
      }
    
    while(not flux.eof())
      {
	string line;
	std::getline(flux, line);

	stringtools_replace(line, " ", "");

	string res1, val; 
	stringtools_split(line, res1, val, '=');
	string varname, minormax;
	stringtools_split(res1, varname, minormax, '_');

	//stringtools_upper();
	
	if(varname=="")
	  { //BREAK : Attention! Un saut de ligne dans le fichier de paramètres inhibe l'initialisation
	    //des paramètres qui seraient dans les lignes suivantes
	    break;
	  }
	
	
	//cerr<<varname<<" "<<minormax<<" "<<val<<endl;
	
	//cerr<<">> "<<line<<endl;


	size_t index = 0;
	string varnameup = stringtools_upper(varname);
	
	for(size_t i = 0; i<variables.size(); i++)
	  {
	    string temp_varname =  variables[i].getName();
	    temp_varname = stringtools_upper(temp_varname);
	    
	    if(temp_varname==varnameup)
	      {
		index = i;
		break;
	      }
	  }

	//cerr<<"@@ "<<varnameup<<" - "<<index<<endl;
	
	assert(index<variables.size());
	
	double v; stringtools_toDouble(val, v);
	
	if(minormax=="MAX")
	  variables[index].setUpperBound(v); 
	
	if(minormax=="MIN")
	variables[index].setLowerBound(v);
	
	//variables.
      }
    
    flux.close();
  }

  
  void print(ostream& os)
  {
    for(size_t i = 0; i<variables.size(); i++)
      {
	os<<variables[i].to_string()<<"\n";
      }
    
    os.flush();
  }
  
};





//#define DCBD_INDEX 0
//#define RHOCBD_INDEX 1
//#define KN_INDEX 2
//#define ALPHA_INDEX 3
//#define XU_INDEX 4
//#define GAMMA_INDEX 5
//#define BETA_INDEX 6

#define EPSIAM_INDEX 0
#define SIGMA_INDEX 1
#define EPSICB_INDEX 2
#define DCBD_INDEX 3
#define RHOCBD_INDEX 4
#define KN_INDEX 5
#define GAMMA_INDEX 6
#define XU_INDEX 7

int main(int argc, char* argv[])
{
  //srand(time(0));
  
  if(false)
    {
      
      string metainput_filename = "inliqreal.run.meta";

      char tag = '@';
      list<pair<string, string> > l_param;
      LRCS_parseMetaInputFile(metainput_filename, l_param, tag);
      
      string visco_filename = "ave_visco";
      pair<double, double> visco = LRCS_readViscosity(visco_filename, 1, 2);
      cerr<<"viscosity = "<<visco.first<<" - std = "<<visco.second<<endl;

      string density_filename = "density_slurry";
      double dens = LRCS_readDensity(density_filename);
      cerr<<"density = "<<dens<<endl;


      vector<pair<string,double> > v_values;
      v_values.push_back(make_pair("@P1@", 1.3));
      v_values.push_back(make_pair("@P2@", 882));
      v_values.push_back(make_pair("@P3@", 272));
      v_values.push_back(make_pair("@P4@", 1863688766278876767));
      v_values.push_back(make_pair("@P5@", 1.6277662));
      LRCS_generateInputFile(metainput_filename,
			     "titi.txt",
			     v_values);

      //lancer lammps
      
      
      return 0;
    }
  
  string filename = "paramlrcs.txt";
  
  list<string> filenamesList; //TODO Le placer dans Problem
  //filenamesList.push_back("inliq@NUMPARTICLE@.run");

  string metainput_filename = "inliqreal.run.meta"; filenamesList.push_back(metainput_filename);
  string metainput1_filename = "viscosityShearRate1.run.meta"; filenamesList.push_back(metainput1_filename);
  string metainput2_filename = "viscosityShearRate2.run.meta"; filenamesList.push_back(metainput2_filename);
  
  //ATTENTION MODIFIER AUSSI DANS SCRIPT BATCH
  filenamesList.push_back("inliqreal.run_@NUMPARTICLE@"); // ATTENTION MODIFIER '-' AUSSI DANS FICHIER DE LANCEMENT
  filenamesList.push_back("viscosityShearRate1.run_@NUMPARTICLE@");
  filenamesList.push_back("viscosityShearRate2.run_@NUMPARTICLE@");

  /*
  filenamesList.push_back("ave_visco@NUMPARTICLE@.1");
  //filenamesList.push_back("density_slurry@NUMPARTICLE@");
  filenamesList.push_back("density_slurry");
  //filenamesList.push_back("ave_visco@NUMPARTICLE@.2");
  
  filenamesList.push_back("script@NUMPARTICLE@.sh");
  filenamesList.push_back("coord_in.dat");
  */
  
  //LRCS
  LRCSVariables variables;
  variables.init(filename);
  variables.print(std::cerr);

  //Problem definition
  //double (*costFunction) (const vector<double>&) = &sphere;
  //double (*costFunction) (const vector<double>&, size_t, const list<string>*) = &sphere;
  double (*costFunction) (const vector<double>&, size_t, const list<string>*, vector<double>*) = &simulate_with_lammps;
  //unsigned int nbVar = 7; //Number of Unknown (Decision) Variables
  //double varMin = -10; //Lower Bound of Decision Variables
  //double varMax = 10; //Upper Bound of Decision Variables

  vector<double> varMin;
  vector<double> varMax;
  vector<string> varNames;
 
  //varMin.resize(nbVar, 0);
  //varMax.resize(varMin.size(), 1);
  varMin.push_back(variables.getLowerBound(EPSIAM_INDEX));
  varMax.push_back(variables.getUpperBound(EPSIAM_INDEX));
  varNames.push_back(variables.getName(EPSIAM_INDEX));
  
  varMin.push_back(variables.getLowerBound(SIGMA_INDEX));
  varMax.push_back(variables.getUpperBound(SIGMA_INDEX));
  varNames.push_back(variables.getName(SIGMA_INDEX));

  varMin.push_back(variables.getLowerBound(EPSICB_INDEX));
  varMax.push_back(variables.getUpperBound(EPSICB_INDEX));
  varNames.push_back(variables.getName(EPSICB_INDEX));
  
  varMin.push_back(variables.getLowerBound(DCBD_INDEX));
  varMax.push_back(variables.getUpperBound(DCBD_INDEX));
  varNames.push_back(variables.getName(DCBD_INDEX));
  
  varMin.push_back(variables.getLowerBound(RHOCBD_INDEX));
  varMax.push_back(variables.getUpperBound(RHOCBD_INDEX));
  varNames.push_back(variables.getName(RHOCBD_INDEX));
  
  varMin.push_back(variables.getLowerBound(KN_INDEX));
  varMax.push_back(variables.getUpperBound(KN_INDEX));
  varNames.push_back(variables.getName(KN_INDEX));
  
  varMin.push_back(variables.getLowerBound(GAMMA_INDEX));
  varMax.push_back(variables.getUpperBound(GAMMA_INDEX));
  varNames.push_back(variables.getName(GAMMA_INDEX));

  varMin.push_back(variables.getLowerBound(XU_INDEX));
  varMax.push_back(variables.getUpperBound(XU_INDEX));
  varNames.push_back(variables.getName(XU_INDEX));
  
  for(size_t i = 0; i<varMin.size(); i++)
    {
      cerr<<"parameter "<<i+1<<" between "<<varMin[i]<<" and "<<varMax[i]<<endl; 
    }

  /*
  double vvv = 0;
  double yyy = 0;
  for(size_t k = 0; k<100000; k++)
    {
      for(size_t i = 0; i<varMin.size(); i++)
	{
	  double p = random_probability();
	  double lb = varMin[i];
	  double ub = varMax[i];
	  double aaa = lb+(p*(ub-lb));
	  cerr<<"random - "<<p<<" # "<<lb<<" < "<<aaa<<" < "<<ub<<" "<<endl; //sleep(1);

	  if(p<0.5)
	    vvv++;
	  yyy++;
	}
    }

  cerr<<vvv<<"/"<<yyy<<endl;
  */
  
  //nbVar = varMax.size();
 
  //Parameters of PSO

  double kappa = 1;
  double phi1 = 2.05;
  double phi2 = 2.05;
  double phi = phi1+phi2;
  double temp = 2-phi-sqrt((phi*phi)-(4*phi)); if(temp<0) temp = -temp;  //abs
  double chi = 2*kappa/temp;
  
  size_t maxIterations = 100; //30; //Maximum Number of Iterations
  size_t nbPop = 8; //Population Size (Swarm Size)
  double w = chi; //Inertia Coefficient
  double wdamp = 1; //Damping Ration of Inertia Coefficient
  double c1 = chi*phi1; //Personal Acceleration Coefficient
  double c2 = chi*phi2; //Social Acceleration Coefficient


  
  bool showIterInfo = true;

  PSO_Problem problem(varMin, varMax, costFunction, &filenamesList, &varNames);
  PSO_Parameters params(maxIterations, nbPop, w, wdamp, c1, c2, showIterInfo);
  PSO_Out out;
  //pso(problem, params, out);
  //multithread_pso(problem, params, out);
  multithread_pso_with_ml(problem, params, out);

  
  //Results section


  cerr<<"global best position = ";
  for(size_t i = 0; i<out.global_best_solution.size(); i++)
    {
      cerr<<out.global_best_solution[i]<<" "; 
    }
  cerr<<endl;
  
  return 0;
}
