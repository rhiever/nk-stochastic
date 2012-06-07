#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

#define randDouble ((double)rand()/(double)RAND_MAX)
#define discreteMode
#define ED 5

vector<double> fitnessesForUIntGenotype;
map<vector<double>,double> fitnessForGenotypes;
map<vector<double>,int> genotypeCount;

double **NKTable;
unsigned int bitMask=0;
int N=0,K=0,A=0;
int update=0,updates=0;
double replacementRate=0;
double mutationRate=0;
unsigned int nrOfGenotypes=0;
string out_file_name = "";
double dist_var=0;
double trans_error=0;

class tMutation
{
public:
    int from,to,where;
    void setup(int f, int w, int t){from=f; to=t; where=w;}
};

class tAgent
{
public:
#ifdef discreteMode
  vector<int> genome;
  unsigned int genome_uint;
#else
  vector<double> genome;
#endif
  double fitness;
  double real_fitness;
  int born;
  int nrPointingAtMe;
  tAgent *ancestor;
  void inherit(tAgent *from);
  tAgent();
  ~tAgent();
  void setupDiscrete(unsigned int newGenotype);
  void setupProbabilistic(void);
  void saveNeighborhood(FILE* lod);
  void saveLOD(FILE* lod, bool genotypic);
  void show(void);
  void computeOwnFitness(void);
};

vector<tAgent*> lodBuffer;

double genotypeToFitness(unsigned int genotype);
void showNKTable(void);
double normalDistribution(double mu, double sigma); 
void doEpistasis(FILE *f);
uint vectorToUInt(vector<int> genome);
double calc_pop_mean(vector<tAgent*> pop);
double calc_pop_mean_sq(vector<tAgent*> pop);
double calc_pop_stdv(vector<tAgent*> pop);
double calc_pop_variance(vector<tAgent*> pop);

int main (int argc, char * const argv[])
{
  int i=0,j=0;
  double maxFitness=0;
  double maxFitnessGenotypic=0;
  double lowest=1.0;
  unsigned int IDofLowest=0;
  double highest=0.0;
  unsigned int IDofHighest=0;
  bool max_found = false;
  vector<tAgent*> population;
  char outdir[50]="";
  
  strcpy(outdir, argv[1]);
  N=atoi(argv[2]);				//N
  K=atoi(argv[3]);				//K
  A=atoi(argv[4]);				//popsize;
  nrOfGenotypes=(unsigned int)1<<(unsigned int)N;
  updates=atoi(argv[5]);			//updates
  replacementRate=atof(argv[6]);	//replacement rate
  mutationRate=atof(argv[7]);		//mutation rate
  srand(atoi(argv[9]));           //set random number generator seed
  //srand((unsigned int)time(NULL));
  dist_var = .5;//atof(argv[10]);       //set variance for normal distribution
  trans_error = atof(argv[10]);       //set transcription error rate %
  
  stringstream temp_fn(stringstream::in | stringstream::out);
  temp_fn << outdir << N << "-" << K << "-" << A << "-" << updates << "-" << (replacementRate * 100) << "-" << (mutationRate * 100) << "-" << atoi(argv[9]) << "-" << (trans_error * 100);// << "-" << (dist_var * 100);
    
  out_file_name = temp_fn.str();
  
#ifdef discreteMode
  if(trans_error > -1)
    {
      out_file_name += "-nondet";//"-det";
    }
  else
    {
      out_file_name += "-det";
    }
#else
  out_file_name += "-nondet";
#endif
  
  bitMask=(1<<K)-1;
  cout<<N<<" "<<K<<" "<<bitMask<<" "<<nrOfGenotypes<<endl;
  cout<<A<<" "<<RAND_MAX<<endl;
  
  //generate new random NK table
  if(strcmp(argv[8],"RAND")==0)
    {
      //make the NKTable
      cout<<"make the table"<<endl;
      NKTable=NULL;
      NKTable=(double**)malloc(sizeof(double*)*N);
      for(i=0;i<N;i++)
        {
	  NKTable[i]=(double*)malloc(sizeof(double)*(1<<K));
	  for(j=0;j<(1<<K);j++)
            {
	      //don't produce zeros as entries
	      do
                {
		  NKTable[i][j]=randDouble;
                } while(NKTable[i][j] == 0.0);
            }
        }
      
      //showNKTable();
      
      //make all fitnesses for genotype table
      cout<<"make all fitnesses"<<endl;
      fitnessesForUIntGenotype.resize(nrOfGenotypes);
      cout<<fitnessesForUIntGenotype.size()<<endl;
      
      for(unsigned int u=0;u<nrOfGenotypes;u++)
        {
	  fitnessesForUIntGenotype[u]=genotypeToFitness(u);
          
	  // cout<<genotypeToFitness(u)<<endl;
	  if(fitnessesForUIntGenotype[u]<lowest)
            {
	      lowest=fitnessesForUIntGenotype[u];
	      IDofLowest=u;
            }
	  
	  if(fitnessesForUIntGenotype[u]>highest)
            {
	      highest=fitnessesForUIntGenotype[u];
	      IDofHighest=u;
            }
        }
      
      /*FILE *f;=fopen("fls.txt","w+t");
	for(i=0;i<nrOfGenotypes;i++)
	fprintf(f, "%f\n",fitnessesForUIntGenotype[i]);
	fclose(f);*/
    }
  
  //read in a file
  else
    {
      fitnessesForUIntGenotype.clear();
      FILE *f=fopen(out_file_name.c_str(),"r+t");
      for(unsigned int u=0;u<nrOfGenotypes;u++)
        {
	  float d;
	  fscanf(f,"%f\n",&d);
	  fitnessesForUIntGenotype.push_back((double)d);
	  if((double)d<lowest)
            {
	      lowest=(double)d;
	      IDofLowest=u;
            }
        }
      fclose(f);
    }
  
  // initialize the population
  population.clear();
  
  // temporarily turn off transcription error
  double tmp_te = trans_error;
  trans_error = 0.0;

  vector<tAgent*>::iterator popiter;

  for(i = 0; i < A; i++)
    {
      tAgent *dA=new tAgent;
#ifdef discreteMode
      //dA->setupDiscrete(IDofLowest);
      dA->setupDiscrete(rand() % nrOfGenotypes);
#else
      dA->setupProbabilistic();
#endif

      dA->computeOwnFitness();
      
      if(population.size() == 0)
	{
	  population.push_back(dA);
	}
      else if(population[population.size() - 1]->fitness <= dA->fitness)
	{
	  population.push_back(dA);
	}
      else
	{
	  for(popiter = population.begin(); popiter < population.end(); popiter++)
	    {
	      if((*popiter)->fitness > dA->fitness)
		{
		  population.insert(popiter, dA);
		  break;
		}
	    }
	}
    }

  // take the bottom 50% (fitness) of population and make clones. remove the top 50% fitness organisms.
  for(i = 0; i < A / 2; i++)
    {
      tAgent *dA = new tAgent;

      dA->setupDiscrete(population[i]->genome_uint);
      dA->computeOwnFitness();

      delete population[(A / 2) + i];
      population[(A / 2) + i] = dA;
    }

  // turn transcription error back on
  trans_error = tmp_te;
  
  cout<<lowest<<" "<<population[0]->fitness<<" "<<IDofLowest<<" "<<IDofHighest<<endl;
  cout<<"lowest fitness: "<<fitnessesForUIntGenotype[IDofLowest]<<" highest fitness: "<<fitnessesForUIntGenotype[IDofHighest]<<endl;
  cout<<"pop size: "<<population.size()<<endl;
  
  // add 1 to updates count so the < can be used instead of <= in the for loop
  updates += 1;

  string ftg_filename(out_file_name.c_str());
  ftg_filename += "_fitness_trajectory_genotypic.txt";
  FILE *ftg_data=fopen(ftg_filename.c_str(),"w+t");

  string ftp_filename(out_file_name.c_str());
  ftp_filename += "_fitness_trajectory_phenotypic.txt";
  FILE *ftp_data=fopen(ftp_filename.c_str(),"w+t");

  // main loop
  for(update = 1; update < updates; update++)
    {
      // find current highest fitness genotype in population
      maxFitness = 0.0;
      maxFitnessGenotypic = 0.0;
      for(i = 0; i < A; i++)
	{
	  // if highest peak found, end simulation
	  if(population[i]->real_fitness == highest && !max_found)
	    {
	      cout << "max found" << endl;
	      max_found = true;
	      update--; // decrement counter since highest was found in last update
	    }

	  // if higher peak is in population, store it as the current max
	  if(population[i]->fitness > maxFitness)
	    {
	      maxFitness = population[i]->fitness;
	    }

	  if(population[i]->real_fitness > maxFitnessGenotypic)
            {
	      maxFitnessGenotypic = population[i]->real_fitness;
	    }
	}

      // log trajectory information on regular intervals
      if((update - 1) % 10 == 0)
	{
	  double pop_mean = calc_pop_mean(population);
	  pop_mean *= pop_mean;

	  fprintf(ftg_data, "%d %f %f %f %f %f\n", (update - 1), highest, maxFitnessGenotypic, calc_pop_variance(population), pop_mean, calc_pop_mean_sq(population));
	  fprintf(ftp_data, "%d %f %f\n", (update - 1), highest, maxFitness);
	}

      for(i = 0; i < A; i++)
	{
	  if(randDouble < replacementRate)
	    {
	      population[i]->nrPointingAtMe--;
	      if(population[i]->nrPointingAtMe==0)
		delete population[i];
	      
	      tAgent *dA = new tAgent;
	      do
		{
		  j=rand()%A;
		} while((i == j) || (population[j]->born == update) || (randDouble > population[j]->fitness / maxFitness));
	      dA->inherit(population[j]);
	      population[i] = dA;
	    }
	}
      
      //cout << update << " " << maxFitness << " " << genotypeCount.size() << endl;
    }

  fclose(ftg_data);
  fclose(ftp_data);
    
  /*string NK_table_filename(out_file_name.c_str());
    NK_table_filename += "_NK_table.txt";
    FILE *f=fopen(NK_table_filename.c_str(),"w+t");
    
    fprintf(f,"%f   ", highest);
    
    for(j=0;j<N;j++)
    {
        fprintf(f, "%i ", (IDofHighest>>j)&1);
    }
    
    fprintf(f, "\n");
    
    for(i=0;i<N;i++)
    {
        for(j=0;j<(1<<K);j++)
            fprintf(f,"%f   ",NKTable[i][j]);
        fprintf(f,"\n");
    }
    fclose(f);*/
    
  tAgent *best_genome = population[0];
  
  // find the highest fitness genome - genotypic
  for(int i = 0; i < population.size(); i++)
    {
      tAgent *cur_ancestor = population[i];
      
      while(cur_ancestor != NULL)
	{
	  if(cur_ancestor->real_fitness > best_genome->real_fitness || (cur_ancestor->real_fitness == best_genome->real_fitness && cur_ancestor->born < best_genome->born))
	    {
	      best_genome = cur_ancestor;
	    }
	  
	  cur_ancestor = cur_ancestor->ancestor;
	}
    }
  
  tAgent *best_genome_phenotypic = population[0];
  
  // find the highest fitness genome - phenotypic
  for(int i = 0; i < population.size(); i++)
    {
      tAgent *cur_ancestor = population[i];
      
      while(cur_ancestor != NULL)
	{
	  if(cur_ancestor->fitness > best_genome_phenotypic->fitness || (cur_ancestor->fitness == best_genome_phenotypic->fitness && cur_ancestor->born < best_genome_phenotypic->born))
	    {
	      best_genome_phenotypic = cur_ancestor;
	    }
	  
	  cur_ancestor = cur_ancestor->ancestor;
	}
    }
  
  //save LOD for last most recent common ancestor
  /*string lmrca_filename(out_file_name.c_str());
    lmrca_filename += ".txt";
    FILE *lod=fopen(lmrca_filename.c_str(),"w+t");
    tAgent *lmrca = NULL;
    tAgent *dummy_var = NULL;
    
    dummy_var = population[1];
    
    while(dummy_var->ancestor != NULL)
    {
    if(dummy_var->ancestor->nrPointingAtMe != 1)
    {
    lmrca = dummy_var;
    }
    
    dummy_var = dummy_var->ancestor;
    }
    
    
    lmrca->saveLOD(lod);
    fclose(lod);*/
  
  //save LOD for best current genome - phenotypic fitnesses
  string bgp_filename(out_file_name.c_str());
  bgp_filename += "_best_genome_phenotypic.txt";
  FILE *bgp_lod=fopen(bgp_filename.c_str(),"w+t");
  
  best_genome->saveLOD(bgp_lod, false);
  fclose(bgp_lod);
  
  //save LOD for best current genome - genotypic fitnesses
  string bgg_filename(out_file_name.c_str());
  bgg_filename += "_best_genome_genotypic.txt";
  FILE *bgg_lod=fopen(bgg_filename.c_str(),"w+t");
  
  best_genome->saveLOD(bgg_lod, true);
  fclose(bgg_lod);
  
  //save genomic neighborhood for best current genome - genotypic fitnesses
  string bgn_filename(out_file_name.c_str());
  bgn_filename += "_best_genome_neighborhood.txt";
  FILE *bgn_data=fopen(bgn_filename.c_str(),"w+t");
  
  best_genome->saveNeighborhood(bgn_data);
  
  //save data for graphing
  string gd_filename(out_file_name.c_str());
  gd_filename += "_graph_data.txt";
  FILE *graph_data=fopen(gd_filename.c_str(),"w+t");
  
  fprintf(graph_data, "%f\n", highest);
  //fprintf(graph_data, "%d %f\n", lmrca->born, lmrca->real_fitness);
  fprintf(graph_data, "%d %f\n", best_genome->born, best_genome->real_fitness);
  fprintf(graph_data, "%d %f\n", best_genome_phenotypic->born, best_genome_phenotypic->fitness);
  
  fclose(graph_data);
  
  //do Epistasis
  //doEpistasis(lod);
  
  cout << "done" << endl;
  return 0;
}

//regular functions
double genotypeToFitness(unsigned int genotype)
{
	int i,j,n;
	double fitness=0.0;
	n=0;
	for(j=0;j<K;j++)
	  n=(n<<1)+((genotype>>j)&1);
	
	for(i=0;i<N;i++)
	  {
	    fitness+=log(NKTable[i][n]);
	    n=((n<<1)+((genotype>>((i+K)%N))&1))&bitMask;
	  }
	fitness=exp(fitness/(double)N);
    /* //testfitness maker
     fitness=0.0;
     for(i=0;i<N;i++)
     fitness+=1.0/(double)N*(double)((genotype>>i)&1);
     */
	return fitness;
}

void showNKTable(void)
{
	int i,j;
	for(j=0;j<(1<<K);j++)
    {
		cout<<j<<": ";
		for(i=0;i<N;i++)
			cout<<NKTable[i][j]<<" ";
		cout<<endl;
	}
}

double normalDistribution(double mu, double sigma)
{    
    double u1 = randDouble;
    double u2 = randDouble * 2;
    
    double z1 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    //double z2 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    
    double x1 = mu + z1 * sigma;
    //double x2 = mu + z2 * sigma;
    
    return x1;
}

double calc_pop_mean(vector<tAgent*> pop)
{
  double sum = 0.0;
  for(int i = 0; i < pop.size(); i++)
    {
      sum += pop[i]->real_fitness;
    }
  
  return sum / (double)pop.size();
}

double calc_pop_mean_sq(vector<tAgent*> pop)
{
    double sum2 = 0.0;
    for(int i = 0; i < pop.size(); i++)
    {
	sum2 += pow(pop[i]->real_fitness, 2);
    }

    return sum2 / (double)pop.size();
}

double calc_pop_stdv(vector<tAgent*> pop)
{
  double mean = calc_pop_mean(pop);
  double stdv = 0.0;

  for(int i = 0; i < pop.size(); i++)
    {
      stdv += pow(pop[i]->real_fitness - mean, 2);
    }

  return sqrt(stdv / (double)(pop.size() - 1));
}

double calc_pop_variance(vector<tAgent*> pop)
{
    double mean2 = calc_pop_mean_sq(pop);

    double mean = calc_pop_mean(pop);

    return mean2 - pow(mean, 2);

  //double stdv = calc_pop_stdv(pop);

  //return stdv * stdv;
}

/*
void doEpistasis(FILE *f)
{
    int i,j,n,k;
    map<string,tMutation> mutationsMap;
    vector<string> mutations;
    set<string> setM;
    set<int> setI;
    i=0;
    while(i<(lodBuffer.size()-1)){
        k=0;
        for(j=0;j<N;j++){
            if(lodBuffer[i]->genome[j]!=lodBuffer[i+1]->genome[j])
                k++;
        }
        if(k==0){
            cout<<"--"<<endl;
            lodBuffer.erase(lodBuffer.begin()+i);
        }
        else
            i++;
    }
    for(i=0;i<lodBuffer.size()-ED;i++)
    {
        mutations.clear();
        mutationsMap.clear();
        setM.clear();
        setI.clear();
        for(j=0;j<ED;j++){
            for(n=0;n<N;n++)
                if(lodBuffer[i+j]->genome[n]!=lodBuffer[i+j+1]->genome[n]){
                    string M;
                    tMutation tM;
                    char c[1000];
                    sprintf(c,"%i-%i-%i",lodBuffer[i+j]->genome[n],n,lodBuffer[i+j+1]->genome[n]);
                    tM.setup(lodBuffer[i+j]->genome[n],n,lodBuffer[i+j+1]->genome[n]);
                    M.assign(c);
                    mutations.push_back(M);
                    setM.insert(M);
                    setI.insert(n);
                    mutationsMap[M]=tM;
                }
        }
        fprintf(f,"%i   %f",lodBuffer[i]->born,lodBuffer[i]->fitness);
        for(j=0;j<N;j++)
            fprintf(f," %i",lodBuffer[i]->genome[j]);
        if((setI.size()!=ED)||(setM.size()!=mutations.size())||(mutationsMap.size()!=ED))
        {
            for(j=0;j<1<<ED;j++)
                fprintf(f,"    0.0");
            fprintf(f,"\n");
        }
        else
        {
            
            for(j=0;j<(1<<ED);j++){
                tAgent *A=new tAgent;
                A->genome=lodBuffer[i]->genome;
                for(k=0;k<ED;k++)
                    if((j&(1<<k))!=0)
                    {
                        //apply mutation
                        A->genome[mutationsMap[mutations[k]].where]=mutationsMap[mutations[k]].to;
                    }
                A->computeOwnFitness();
                fprintf(f,"    %f",A->fitness);
            }
            fprintf(f,"\n");
        }
    }
}
*/

uint vectorToUInt(vector<int> genome)
{
  // convert the genome into an unsigned int                                                                                                                                                                                   
  uint IDofGenome = 0;

  for(int j = 0; j < genome.size(); j++)
    {
      if(genome[j] == 1)
	{
	  IDofGenome |= (unsigned int)1 << (unsigned int)j;
	}
    }

  return IDofGenome;
}

// tAgent definitions
void tAgent::inherit(tAgent *from)
{
    genome.resize(from->genome.size());
    from->nrPointingAtMe++;
    ancestor=from;
    bool mutated=false;
    for(int i=0;i<from->genome.size();i++)
      {
        if(randDouble < mutationRate)
        {
	  mutated=true;

#ifdef discreteMode
	  genome[i] = rand()&1;
	  // flip the bit if mutated
	  /*if(genome[i] == 0)
	    {
	      genome[i] = 1;
	    }
	  else
	    {
	      genome[i] = 0;
	    }*/
#else
	  //randomly add or subtract a random value
	  /*double randVal = normalDistribution(dist_around, dist_var);
	    
          if(rand() % 2 == 0)
	  {
	  genome[i] += randVal;
	  }
	  else
	  {
	  genome[i] -= randVal;
	  }*/
          
	  //set the locus to a new random value
	  genome[i] = randDouble;//normalDistribution(genome[i], dist_var);
          
	  //cap the genome values in the range [0,1]
	  /*if(genome[i] > 1)
            {
	    genome[i] = 1;
            }
            else if(genome[i] < 0)
            {
	    genome[i] = 0;
            }*/
#endif
        }
        else
	  {
            genome[i]=from->genome[i];
	  }
      }
#ifdef discreteMode
    computeOwnFitness();
#else
    if(mutated)
      {
	computeOwnFitness();
      }
    else
      {
	fitness = from->fitness;
      }
#endif
}

tAgent::tAgent()
{
    nrPointingAtMe=1;
    ancestor=0;
    born=update;
    fitness=0;
    real_fitness=0;
    genome_uint=0;
}

tAgent::~tAgent()
{
    if(ancestor!=NULL)
    {
        ancestor->nrPointingAtMe--;
        if(ancestor->nrPointingAtMe==0)
            delete ancestor;
    }
    /*
     if(--genotypeCount[genome]==0)
     {
     genotypeCount.erase(genotypeCount.find(genome));
     fitnessForGenotypes.erase(fitnessForGenotypes.find(genome));
     }*/
}

void tAgent::setupDiscrete(unsigned int newGenotype)
{
  genome_uint = newGenotype;

  genome.resize(N);
  for(int i=0;i<N;i++)
    {
      if(((newGenotype>>i)&1)==1)
	genome[i]=1;
      else
	genome[i]=0;
    }

  fitness=fitnessesForUIntGenotype[newGenotype];
}

void tAgent::setupProbabilistic(void)
{
    
    genome.resize(N);
    for(int i=0;i<N;i++)
        genome[i]=0.5;
}

void tAgent::saveNeighborhood(FILE* lod)
{
  if(ancestor!=NULL)
    {
      ancestor->saveNeighborhood(lod);
      
      if(ancestor->nrPointingAtMe!=1)
        {
	  //fclose(lod);
	  //exit(0);
        }
    }
  
  lodBuffer.push_back(this);
  
  if(lod != NULL)
    {
      bool b = false;
        
      if(ancestor != NULL)
        {
	  for(int i = 0; i < genome.size(); i++)
            {
	      if(genome[i] != ancestor->genome[i])
                {
		  b = true;
		  break;
                }
            }
        }
        
      if(b || born == 0)
        {
	  fprintf(lod, "%d ", born);

	  for(int i = 0; i < genome.size(); i++)
            {
#ifdef discreteMode
	      vector<int> v_neighbor (genome);
	      if(v_neighbor[i] == 0)
		{
		  v_neighbor[i] = 1;
		}
	      else
		{
		  v_neighbor[i] = 0;
		}

	      uint neighbor = vectorToUInt(v_neighbor);
	     
	      fprintf(lod,"%f ", (fitnessesForUIntGenotype[neighbor] - fitness) / fitness);
#endif
            }
	  fprintf(lod, "\n");
        }
    }
}

void tAgent::saveLOD(FILE* lod, bool genotypic)
{
    if(ancestor!=NULL)
    {
      ancestor->saveLOD(lod, genotypic);
		
        if(ancestor->nrPointingAtMe!=1)
        {
            //fclose(lod);
            //exit(0);
        }
    }
	
    lodBuffer.push_back(this);
	
    if(lod!=NULL)
    {
        bool b = false;
        
        if(ancestor != NULL)
        {
            for(int i = 0; i < genome.size(); i++)
            {
                if(genome[i] != ancestor->genome[i])
                {
                    b = true;
		    break;
                }
            }
        }
        
        if(b || born == 0 || fitness != ancestor->fitness)
        {
	  if(genotypic)
	    {
	      fprintf(lod,"%i %f  ",born,real_fitness);
	    }
	  else
	    {
	      fprintf(lod,"%i %f  ",born,fitness);
	    }

	    for(int i=0;i<genome.size();i++)
            {
#ifdef discreteMode
	      fprintf(lod,"%i ", genome[i]);
#else
	      fprintf(lod,"%f ", genome[i]);
#endif
            }
            fprintf(lod,"\n");
        }
    }
}

void tAgent::show(void)
{
    printf("%i   %f ",born,fitness);
    for(int i=0;i<N;i++)
    {
#ifdef discreteMode
        printf("%i",(int)genome[i]);
#else
        printf("    %f",genome[i]);
#endif
    }
    printf("\n");
}
void tAgent::computeOwnFitness(void)
{
#ifdef discreteMode
 
  // compute lookahead fitness
  unsigned int G=0;
  unsigned int G_real=0;
  for(int i=0;i<genome.size();i++)
    {
      // chance of error in transcription at each locus
      if(randDouble < trans_error)
        {
	  // if error, flip the bit
	  if(genome[i]!=1)
            {
	      G|=(unsigned int)1<<(unsigned int)i;
            }
	  
	  // if genome[i]==1, then the 0 is already there since all bits are initialized to 0
        }
      // no error in transcription at this locus
      else
        {
	  // no error, so copy as-is
	  if(genome[i]==1)
            {
	      G|=(unsigned int)1<<(unsigned int)i;
            }
        }

      // convert genome for real fitness calculation
      if(genome[i]==1)
        {
          G_real|=(unsigned int)1<<(unsigned int)i;
        }
    }
  
  genome_uint = G_real;
  fitness=fitnessesForUIntGenotype[G];
  real_fitness=fitnessesForUIntGenotype[G_real];
  
#else
  fitness=0.0;
  real_fitness = 0.0;
  
  for(unsigned int u=0;u<nrOfGenotypes;u++)
    {
      double p=1.0;
      for(int i=0;i<N;i++)
	if(((u>>i)&1)==1)
	  p*=genome[i];
	else
	  p*=(1.0-genome[i]);
      
      fitness+=fitnessesForUIntGenotype[u]*p;
    }
  real_fitness = fitness;
  fitnessForGenotypes[genome]=fitness;
#endif
}
