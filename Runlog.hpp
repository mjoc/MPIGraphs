#ifndef __MOC_RUNLOG_
#define __MOC_RUNLOG_

#include <time.h>
#include <cstdarg>
#include <sstream>
#include <string>
#include <iostream>

class Runlog{
  int _processRank;
  std::ofstream _logfile;
  std::string _filepath;
  std::string _filename;
public:
  Runlog(std::string filepath, std::string filename){
    _processRank = -1;
    
    std::string fullFilePath = filepath+filename;
    char filePathAndName[100];
		snprintf(filePathAndName,100,"%s", fullFilePath.c_str());
    
    _logfile.open(filePathAndName);
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      
      _logfile << formattedDate << " Opening log file" << std::endl;
    }else{
      std::cerr << "Couldn't open logfile" <<  std::endl;
    }
  }
  
  Runlog(std::string filepath, std::string filename, int rank){
    _processRank = rank;

     char filePathAndName[200];
		snprintf(filePathAndName,200,"%sr%d_%s",filepath.c_str(), rank, filename.c_str());
    
    _logfile.open(filePathAndName);
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);

      _logfile << _processRank << " " << formattedDate << " Opening log file" << std::endl;
    }else{
      std::cerr << "######Couldn't open logfile for writing on rank " << rank << "(" << filePathAndName << ") ######" <<std::endl;
    }
  }
  ~Runlog(){
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time(&rawtime);
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " Closing log file"<< std::endl;
      }else{
        _logfile << formattedDate << " Closing log file"<< std::endl;
      }
      _logfile.close();
    }
  }
  void writelog(char *message){
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << message << std::endl;
      }else{
        _logfile << formattedDate << " " << message << std::endl;
      }

    }else{
      std::cout << message << std::endl;
    }
    
  }
  void writelog(std::string message){
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << message << std::endl;
      }else{
        _logfile << formattedDate << " " << message << std::endl;
      }
      
    }else{
      std::cout << message << std::endl;
    }

    
  }
  
  void writelog(std::stringstream message){
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << message.rdbuf() << std::endl;
      }else{
        _logfile << formattedDate << " " << message.rdbuf() << std::endl;
      }
      
    }else{
      std::cout << message.rdbuf() << std::endl;
    }
    
  }

  /*
  void printVector(char* description, const std::vector<int>& vect){
    
    
    for(std::vector<int>::const_iterator it = vect.begin(); it != vect.end(); ++it){
      snprintf (description,200, "%s %d", description, *it);
    }
    
    
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << description << std::endl;
      }else{
        _logfile << formattedDate << " " << description << std::endl;
      }
      
    }else{
      std::cout << description << std::endl;
    }
    
  }

  
 
  void printVector(char* description, const std::vector<double>& vect){
    for(std::vector<double>::const_iterator it = vect.begin(); it != vect.end(); ++it){
      snprintf (description,200, "%s %f", description, *it);
    }
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << description << std::endl;
      }else{
        _logfile << formattedDate << " " << description << std::endl;
      }
      
    }else{
      std::cout << description << std::endl;
    }

  }
  */
  /*
  template<typename DataType>
  void printSet(char* description, const std::set<DataType>& setToprint){
    for(std::set<DataType>::const_iterator it = setToprint.begin(); it != setToprint.end(); ++it){
      sprintf (description, "%s %f", description, (double)(*it));
    }
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << description << std::endl;
      }else{
        _logfile << formattedDate << " " << description << std::endl;
      }
      
    }else{
      std::cout << description << std::endl;
    }
  }
  */

  
  void printSet(char* description, const std::set<int>& setToprint){
    
    for (std::set<int>::iterator it = setToprint.begin(); it != setToprint.end(); ++it)
    {
      snprintf (description,200, "%s %d", description, (*it));
    }
    
    
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << description << std::endl;
      }else{
        _logfile << formattedDate << " " << description << std::endl;
      }
      
    }else{
      std::cout << description << std::endl;
    }
    
  }
  

  
  
  void printSet(char* description, const std::set<int, int>& setToprint){
 
    for (std::set<int, int>::iterator it = setToprint.begin(); it != setToprint.end(); ++it)
      {
        snprintf (description,200, "%s %d", description, (*it));
      }
    
    
    if(_logfile.is_open()){
      time_t rawtime;
      struct tm *info;
      char formattedDate[21];
      time( &rawtime );
      info = localtime( &rawtime );
      strftime(formattedDate,21,"%F %T", info);
      if(_processRank > -1){
        _logfile << _processRank << " " << formattedDate << " " << description << std::endl;
      }else{
        _logfile << formattedDate << " " << description << std::endl;
      }
      
    }else{
      std::cout << description << std::endl;
    }

  }
  
  
  template <typename T> std::string to_string(const T& n){
    std::ostringstream stm;
    stm << n;
    return stm.str();
  }

  };

#endif
