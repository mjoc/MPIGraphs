#ifndef __MOC_GslWrapper_
#define __MOC_GslWrapper_

#include "gsl/gsl_rng.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_randist.h"
#include <vector>

//This is the maximum seed that should be generated when  gsl_rng_uniform_int
// generates seeds
#define SEEDMAX 987654321

class RngWrapper{
  const gsl_rng_type *_rngType;
  gsl_rng *_rng;
  
public:
  RngWrapper(){
    gsl_rng_env_setup();
    _rngType = gsl_rng_default;
    _rng = gsl_rng_alloc (_rngType);
  }
  
  ~RngWrapper(){
    gsl_rng_free (_rng);
  }
  
  std::vector<long> getSeeds(int n){
    std::vector<long> seeds;
    
    const gsl_rng_type * tempRngType;
    gsl_rng * tempRng;
    
    gsl_rng_env_setup();
    
    tempRngType = gsl_rng_default;
    tempRng = gsl_rng_alloc (tempRngType);
    
    for (int i = 0; i < n; i++)
    {
      seeds.push_back((gsl_rng_uniform_int(tempRng,SEEDMAX)));
    }
    
    gsl_rng_free (tempRng);
    return seeds;
  }

  void setSeed(long seed){
    gsl_rng_set (_rng, seed);
  }
  
  double getUniform01Random(){
     return gsl_rng_uniform(_rng);
  }
  
  std::vector<double> getUniform01Randoms(size_t n){
    std::vector<double> randoms(n);
    for(std::vector<double>::iterator it = randoms.begin(); it < randoms.end(); ++it)
      *it = gsl_rng_uniform(_rng);
  
    return randoms;
  }

  std::vector<int> getRandomPermutation(size_t n){
    std::vector<int> randperm(n);
    for (size_t i = 0; i < n; i++)
    {
      randperm[i] = (int)i;
    }
    gsl_ran_shuffle(_rng, &randperm[0], n, sizeof(int));
    
    return randperm;
  }
};

#endif

