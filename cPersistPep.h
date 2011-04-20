#ifndef _CPERSISTPEP_H
#define _CPERSISTPEP_H

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

//const int PPM=10;
//const int CONTAM=150;
//const int MAXMASS=8000;
//const int MINMASS=700;
//const int MD=12;

typedef struct sPep{
  int charge;
  double monoMass;
	double basePeak;
  float intensity;
	char formula[64];
	char mods[256];
	double xCorr;
} sPep;

typedef struct sScan {
  vector<sPep> *vPep;
	int scanNum;
	char file[256];
	float rTime;

	//Constructors & otherwise
	sScan(){vPep = new vector<sPep>;};
	sScan(const sScan& s){
		vPep = new vector<sPep>;
		for(int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
		scanNum = s.scanNum;
		strcpy(file,s.file);
		rTime=s.rTime;
	};
	~sScan(){delete vPep;};


	sScan& operator=(const sScan& s){
		if(this!=&s){
			delete vPep;
			vPep = new vector<sPep>;
			for(int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
			scanNum = s.scanNum;
			strcpy(file,s.file);
			rTime=s.rTime;
		};
		return *this;
	};

	void clear(){
		delete vPep;
		vPep = new vector<sPep>;
	};
} sScan;


typedef struct MA {
  int scanNum;
  float intensity;
//  double monoMass;
} MA;

typedef struct DBMatch {
  int scanNum;
  char sequence[64];
} DBMatch;

typedef struct realPep {
  char file[256];
  int lowScan;
  int highScan;
  int scanCount;
  int charge;
  double monoMass;
  double basePeak;
  float intensity;
  float sumIntensity;
  float rTime;
  float firstRTime;
  float lastRTime;
  double xCorr;
  char mods[256];
  char sequence[256];
  int bestScan;
  double firstMonoMass;
  double lastMonoMass;
  double avgMonoMass;
  bool MS2;
  float firstIntensity;
  float lastIntensity;
  //int massArraySize;
  //MA massArray[200];

  //for BLinks
  vector<DBMatch>* vDBM;
  vector<DBMatch>* vAM;
  vector<MA>* vMA;
  int MS2Events;

  //Constructors & otherwise
	realPep(){
    vDBM = new vector<DBMatch>;
    vAM = new vector<DBMatch>;
    vMA = new vector<MA>;
  }
	realPep(const realPep& r){
    int i;
		vDBM = new vector<DBMatch>;
    vAM = new vector<DBMatch>;
    vMA = new vector<MA>;
		for(i=0;i<r.vDBM->size();i++) vDBM->push_back(r.vDBM->at(i));
    for(i=0;i<r.vAM->size();i++) vAM->push_back(r.vAM->at(i));
    for(i=0;i<r.vMA->size();i++) vMA->push_back(r.vMA->at(i));
		strcpy(file,r.file);
    lowScan=r.lowScan;
    highScan=r.highScan;
    scanCount=r.scanCount;
    charge=r.charge;
    monoMass=r.monoMass;
    basePeak=r.basePeak;
    intensity=r.intensity;
    sumIntensity=r.sumIntensity;
    rTime=r.rTime;
    firstRTime=r.firstRTime;
    lastRTime=r.lastRTime;
    xCorr=r.xCorr;
    strcpy(mods,r.mods);
    strcpy(sequence,r.sequence);
    bestScan=r.bestScan;
    firstMonoMass=r.firstMonoMass;
    lastMonoMass=r.lastMonoMass;
    avgMonoMass=r.avgMonoMass;
    MS2=r.MS2;
    firstIntensity=r.firstIntensity;
    lastIntensity=r.lastIntensity;
    MS2Events=r.MS2Events;
	}
	~realPep(){
    delete vDBM;
    delete vAM;
    delete vMA;
  }

	realPep& operator=(const realPep& r){
		if(this!=&r){
      int i;
      delete vDBM;
      delete vAM;
      delete vMA;
			vDBM = new vector<DBMatch>;
		  for(i=0;i<r.vDBM->size();i++) vDBM->push_back(r.vDBM->at(i));
      vAM = new vector<DBMatch>;
      for(i=0;i<r.vAM->size();i++) vAM->push_back(r.vAM->at(i));
      vMA = new vector<MA>;
      for(i=0;i<r.vMA->size();i++) vMA->push_back(r.vMA->at(i));
		  strcpy(file,r.file);
      lowScan=r.lowScan;
      highScan=r.highScan;
      scanCount=r.scanCount;
      charge=r.charge;
      monoMass=r.monoMass;
      basePeak=r.basePeak;
      intensity=r.intensity;
      sumIntensity=r.sumIntensity;
      rTime=r.rTime;
      firstRTime=r.firstRTime;
      lastRTime=r.lastRTime;
      xCorr=r.xCorr;
      strcpy(mods,r.mods);
      strcpy(sequence,r.sequence);
      bestScan=r.bestScan;
      firstMonoMass=r.firstMonoMass;
      lastMonoMass=r.lastMonoMass;
      avgMonoMass=r.avgMonoMass;
      MS2=r.MS2;
      firstIntensity=r.firstIntensity;
      lastIntensity=r.lastIntensity;
      MS2Events=r.MS2Events;
		}
		return *this;
	}

} realPep;

typedef struct iTwo {
	int scanID;
	int pepID;
} iTwo;

class cPersistPep {
private:
	vector<sScan>* allScans;
	vector<realPep>* rPeps;
  vector<realPep>* contam;

	//User modifiable flags
	double PPM;	//mass accuracy threshold
  double CONTAM;  //contamination threshold
  double MAXMASS;  //maximum mass
  double MINMASS;  //minimum mass
  double PERSIST;  //persistence
  double GAP;      //maximum gap in persistence

protected:
public:

	cPersistPep();
	cPersistPep(const cPersistPep& p);
	~cPersistPep();

	cPersistPep& operator=(const cPersistPep& p);

	realPep& at(unsigned int i);
	void readHK(char *fn);
  void findPeps(char *fn="\0");
	int size();

  void add(realPep &rp);
  void clear();
	int getPPM();
	double getCONTAM();
	int getMAXMASS();
	int getMINMASS();
	void setPPM(double d);
	void setCONTAM(double d);
	void setMAXMASS(double d);
	void setMINMASS(double d);
  void setGAP(int i);
  void setPERSIST(int i);
  void sortBasePeak();
  void sortMonoMass();

  static int compareBP(const void *p1, const void *p2);
  static int compareMM(const void *p1, const void *p2);


};

#endif
