#include "cPersistPep.h"

using namespace std;

cPersistPep::cPersistPep(){
  PPM=5.0;
  CONTAM=30.0;
  MAXMASS=8000.0;
  MINMASS=600.0;
  PERSIST=3;
  GAP=2;
	allScans = new vector<sScan>;
	rPeps = new vector<realPep>;
	contam = new vector<realPep>;
};

cPersistPep::cPersistPep(const cPersistPep& p){
	int i;
  PPM=p.PPM;
  CONTAM=p.CONTAM;
  MAXMASS=p.MAXMASS;
  MINMASS=p.MINMASS;
  PERSIST=p.PERSIST;
  GAP=p.GAP;
	allScans = new vector<sScan>;
	rPeps = new vector<realPep>;
	contam = new vector<realPep>;
	for(i=0;i<p.allScans->size();i++) allScans->push_back(p.allScans->at(i));
	for(i=0;i<p.rPeps->size();i++) rPeps->push_back(p.rPeps->at(i));
	for(i=0;i<p.contam->size();i++) contam->push_back(p.contam->at(i));
};
	
cPersistPep& cPersistPep::operator=(const cPersistPep& p){
	int i;

	if(this!=&p){
    PPM=p.PPM;
    CONTAM=p.CONTAM;
    MAXMASS=p.MAXMASS;
    MINMASS=p.MINMASS;
    PERSIST=p.PERSIST;
    GAP=p.GAP;

		delete allScans;
		delete rPeps;
		delete contam;

		allScans = new vector<sScan>;
		rPeps = new vector<realPep>;
		contam = new vector<realPep>;
		for(i=0;i<p.allScans->size();i++) allScans->push_back(p.allScans->at(i));
		for(i=0;i<p.rPeps->size();i++) rPeps->push_back(p.rPeps->at(i));
		for(i=0;i<p.contam->size();i++) contam->push_back(p.contam->at(i));
	};
	return *this;
};

cPersistPep::~cPersistPep(){
	delete allScans;
	delete rPeps;
	delete contam;
}

void cPersistPep::add(realPep &rp){
	rPeps->push_back(rp);
}

realPep& cPersistPep::at(unsigned int i){
	return rPeps->at(i);
}

void cPersistPep::clear(){
	allScans->clear();
	rPeps->clear();
}

void cPersistPep::setCONTAM(double d){
  CONTAM=d;
}
void cPersistPep::setGAP(int i){
  GAP=i;
}
void cPersistPep::setMAXMASS(double d){
  MAXMASS=d;
}
void cPersistPep::setMINMASS(double d){
  MINMASS=d;
}
void cPersistPep::setPERSIST(int i){
  PERSIST=i;
}
void cPersistPep::setPPM(double d){
  PPM=d;
}
void cPersistPep::readHK(char *fn) {
	FILE *hkr;
	sScan scan;
	sPep pep;
	double td;
	char tag;
	bool firstScan;

	int pepCount=0;

	delete allScans;
	allScans=new vector<sScan>;

	firstScan=true;

	hkr = fopen(fn,"rt");
	if(hkr==NULL) {
		cout << "Problem reading file." << endl;
		return;
	};

	while(!feof(hkr)){
    
		tag=fgetc(hkr);

    if(tag=='S') {

			if(firstScan) {
				firstScan=false;
			}	else {
				allScans->push_back(scan);
				//cout << scan.vPep->size() << " vs " << allScans.at(allScans.size()-1).vPep->size() << endl;
			};
		
			scan.clear();
      fscanf(hkr,"\t%d\t%f%s\n",&scan.scanNum,&scan.rTime,scan.file);

		} else {

			pepCount++;
			fscanf(hkr,"\t%lf\t%d\t%f\t%lf\t%lf-%lf\t%lf\t%s\t%lf\n", &pep.monoMass,&pep.charge,&pep.intensity,&pep.basePeak,&td,&td,&td,pep.mods,&pep.xCorr);
			scan.vPep->push_back(pep);

		};
    
	};

	fclose(hkr);

	//cout << "Total Peptide IDs: " << pepCount << endl;

};

void cPersistPep::findPeps(char *fn){

	int i,j,k;

	FILE *out;

	vector<iTwo> hits;
	iTwo hit;

	realPep rp;
	strcpy(rp.sequence,"NULL");
	rp.MS2=false;
  rp.MS2Events=0;

	int charge;
	int scanPoint;
	int match;
	int noMatch;
	int maxIndex;
	int totalMatch=0;

	//char mods[256];
	
	bool bMatch;
	bool bFound;

	double mMass;
	double avgMonoMass;

	float maxIntensity;
	float totalIntensity;

	delete rPeps;
	rPeps = new vector<realPep>;

	delete contam;
	contam = new vector<realPep>;

  MA ma;

        //cout << "Before fn" << endl;
	if(fn[0]!=0){
                //cout << "Bad loop" << endl;
		out = fopen(fn,"wt");
		if(out==NULL) {
			cout << "Cannot write to " << fn << endl;
			return;
		};
	//} else {
        //  cout << "Good loop" << endl;
        }

	//cout << allScans->size() << endl;
	//cout << allScans->at(3000).vPep->size() << endl;

	for(i=0;i<allScans->size();i++){

		//pick a peptide
		for(j=0;j<allScans->at(i).vPep->size();j++) {

			//cout << "Checking scan " << i << " of " << allScans->size() << "  pep: " << j << " of " << allScans->at(i).vPep->size() << endl;

			//get the monoMass
			mMass = allScans->at(i).vPep->at(j).monoMass;
			charge = allScans->at(i).vPep->at(j).charge;
			//strcpy(mods,allScans->at(i).vPep->at(j).mods);
			//cout << mMass << " +" << charge << endl;

			//Set up info
			hits.clear();
			hit.scanID=i;
			hit.pepID=j;
			hits.push_back(hit);

			scanPoint=i;
			match=0;
			noMatch=0;

			bFound=false;

			while(true){

				scanPoint++;

				//cout << scanPoint << ":" << allScans.size() << endl;
				if(scanPoint>=allScans->size()) break;

				//check next scan(s)
				bMatch=false;
				for(k=0;k<allScans->at(scanPoint).vPep->size();k++){
					if( fabs(mMass - allScans->at(scanPoint).vPep->at(k).monoMass)/mMass * 1000000 < PPM &&
						  charge == allScans->at(scanPoint).vPep->at(k).charge) {

						//cout << mMass << " - " << allScans.at(scanPoint).vPep->at(k).monoMass << endl;
						//cout << fabs(mMass - allScans.at(scanPoint).vPep->at(k).monoMass)/mMass * 1000000 << endl;
						//match
					
						hit.scanID=scanPoint;
						hit.pepID=k;
						hits.push_back(hit);
						bMatch=true;
						//cout << "Match" << endl;
						break;
					};
				};

				//if match was found;
				if(bMatch) {
					match++;
					noMatch=0;
				} else {
					//no match, check to see if enough prior matches warrant keeping this one
					noMatch++;

					if(noMatch==GAP) {

						//don't keep bad matches
						if(match<(PERSIST-1)) {
							break;
						} else {

							//cout << "match: " << hits.size() << endl;
							totalMatch++;

							//Create stats
							avgMonoMass=0;
							maxIndex=0;
							maxIntensity=0;
							totalIntensity=0;
							for(k=0;k<hits.size();k++){
								avgMonoMass+=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).monoMass;
								if(allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity > maxIntensity){
									maxIntensity=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity;
									maxIndex=k;
								};
								totalIntensity+=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity;
							};

							avgMonoMass/=hits.size();

							//cout << "AVG ok" << endl;
							//cout << "maxIndex: " << maxIndex << endl;

							//output good matches
							strcpy(rp.file,allScans->at(i).file);
							rp.lowScan = allScans->at(hits.at(0).scanID).scanNum;
							rp.highScan = allScans->at(hits.at(hits.size()-1).scanID).scanNum;
							rp.bestScan = allScans->at(hits.at(maxIndex).scanID).scanNum;
							rp.scanCount = hits.size();
							rp.charge = charge;
							rp.monoMass = allScans->at(hits.at(maxIndex).scanID).vPep->at(hits.at(maxIndex).pepID).monoMass;
							rp.firstMonoMass = allScans->at(hits.at(0).scanID).vPep->at(hits.at(0).pepID).monoMass;
                                                        rp.lastMonoMass = allScans->at(hits.at(hits.size()-1).scanID).vPep->at(hits.at(hits.size()-1).pepID).monoMass;
                                                        rp.avgMonoMass = avgMonoMass;
							rp.intensity = maxIntensity;
							rp.firstIntensity = allScans->at(hits.at(0).scanID).vPep->at(hits.at(0).pepID).intensity;
		                                        rp.lastIntensity = allScans->at(hits.at(hits.size()-1).scanID).vPep->at(hits.at(hits.size()-1).pepID).intensity;
							rp.sumIntensity = totalIntensity;
							rp.firstRTime = allScans->at(hits.at(0).scanID).rTime;
							rp.lastRTime = allScans->at(hits.at(hits.size()-1).scanID).rTime;
							rp.rTime = allScans->at(hits.at(maxIndex).scanID).rTime;
							//cout << "Rtime ok" << endl;
							rp.xCorr = allScans->at(hits.at(maxIndex).scanID).vPep->at(hits.at(maxIndex).pepID).xCorr;
							//cout << "before mods" << endl;
							strcpy(rp.mods,allScans->at(hits.at(0).scanID).vPep->at(hits.at(0).pepID).mods);
							rp.basePeak=allScans->at(hits.at(maxIndex).scanID).vPep->at(hits.at(maxIndex).pepID).basePeak;
							//if(hits.size()<200){
							//	rp.massArraySize=hits.size();
              rp.vMA->clear();
							for(k=0;k<hits.size();k++){
								ma.scanNum=allScans->at(hits.at(k).scanID).scanNum;
								ma.intensity=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity;
                rp.vMA->push_back(ma);
							}
							//} else {
							//	rp.massArraySize=0;
							//}
							rPeps->push_back(rp);
							//cout << "HIT: " << rp.monoMass << ", " << rp.charge << endl;
						
							//erase all hits so they are not re-analyzed;
							//erase from last to first otherwise indexes will be incorrect. -- Should not be a factor since only one pep in each
							//scan is removed.
							for(k=0;k<hits.size();k++){
								allScans->at(hits.at(k).scanID).vPep->erase(allScans->at(hits.at(k).scanID).vPep->begin()+hits.at(k).pepID);
							};

							//cout << "Erase ok" << endl;

							//decrement j because we erased it
							j--;
							bFound=true;
							break;

						};//end if(match<3);
						bFound=true;
						break;
					};//end if(noMatch==2);

				};//end if(bMatch);

			};//end while;

			//cout << "End while" << endl;

			if(!bFound){
				if(match<(PERSIST-1)) {
					//break;
				} else {

					totalMatch++;

					//Create stats
					avgMonoMass=0;
					maxIndex=0;
					maxIntensity=0;
					totalIntensity=0;
					for(k=0;k<hits.size();k++){
						avgMonoMass+=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).monoMass;
						if(allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity > maxIntensity){
							maxIntensity=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity;
							maxIndex=k;
						};
						totalIntensity+=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity;
					};

					avgMonoMass/=hits.size();

					//output good matches
					strcpy(rp.file,allScans->at(i).file);
					rp.lowScan = allScans->at(hits.at(0).scanID).scanNum;
					rp.highScan = allScans->at(hits.at(hits.size()-1).scanID).scanNum;
					rp.bestScan = allScans->at(hits.at(maxIndex).scanID).scanNum;
					rp.scanCount = hits.size();
					rp.charge = charge;
          rp.monoMass = allScans->at(hits.at(maxIndex).scanID).vPep->at(hits.at(maxIndex).pepID).monoMass;
          rp.firstMonoMass = allScans->at(hits.at(0).scanID).vPep->at(hits.at(0).pepID).monoMass;
          rp.lastMonoMass = allScans->at(hits.at(hits.size()-1).scanID).vPep->at(hits.at(hits.size()-1).pepID).monoMass;
          rp.avgMonoMass = avgMonoMass;
					rp.intensity = maxIntensity;
          rp.firstIntensity = allScans->at(hits.at(0).scanID).vPep->at(hits.at(0).pepID).intensity;
          rp.lastIntensity = allScans->at(hits.at(hits.size()-1).scanID).vPep->at(hits.at(hits.size()-1).pepID).intensity;
					rp.sumIntensity = totalIntensity;
					rp.firstRTime = allScans->at(hits.at(0).scanID).rTime;
					rp.lastRTime = allScans->at(hits.at(hits.size()-1).scanID).rTime;
					rp.rTime = allScans->at(hits.at(maxIndex).scanID).rTime;
					rp.xCorr = allScans->at(hits.at(maxIndex).scanID).vPep->at(hits.at(maxIndex).pepID).xCorr;
					strcpy(rp.mods,allScans->at(hits.at(0).scanID).vPep->at(hits.at(0).pepID).mods);
					rp.basePeak=allScans->at(hits.at(maxIndex).scanID).vPep->at(hits.at(maxIndex).pepID).basePeak;
					//if(hits.size()<200){
	        //  rp.massArraySize=hits.size();
          rp.vMA->clear();
        	for(k=0;k<hits.size();k++){
	          ma.scanNum=allScans->at(hits.at(k).scanID).scanNum;
            ma.intensity=allScans->at(hits.at(k).scanID).vPep->at(hits.at(k).pepID).intensity;
            rp.vMA->push_back(ma);
          }
					//} else {
					//	rp.massArraySize=0;
					//}
					rPeps->push_back(rp);
					//cout << "HIT: " << rp.monoMass << ", " << rp.charge << endl;
						
					//erase all hits so they are not re-analyzed;
					//erase from last to first otherwise indexes will be incorrect. -- Should not be a factor since only one pep in each
					//scan is removed.
					for(k=0;k<hits.size();k++){
						allScans->at(hits.at(k).scanID).vPep->erase(allScans->at(hits.at(k).scanID).vPep->begin()+hits.at(k).pepID);
					}
					

					//decrement j because we erased it
					j--;
				}
			}

			//cout << "Go to next j  " << j << endl;

		}//for j (all peptides)

	}//for i (all scans)

	cout << "Persistent PIDs: " << totalMatch << endl;

	//Remove contaminants (elements that repeat too frequently).

  //Removing rolling overlaps
  /*
  for(i=0;i<rPeps->size()-1;i++){
    for(j=i+1;j<rPeps->size();j++){
      if( fabs((rPeps->at(i).monoMass-rPeps->at(j).monoMass)/rPeps->at(i).monoMass*1000000) <= PPM &&
         ((rPeps->at(i).firstRTime>=rPeps->at(j).firstRTime && rPeps->at(i).firstRTime<=rPeps->at(j).lastRTime) || 
         (rPeps->at(i).lastRTime>=rPeps->at(j).firstRTime && rPeps->at(i).lastRTime<=rPeps->at(j).lastRTime)) ){
			//merge

        rPeps->erase(rPeps->begin()+j);
        j--;
      }
    } 
	}

  cout << "Rolling errors: " << rPeps->size() << endl;
*/
  vector<realPep> tmpPeps;

  //New method just uses retention time
  for(i=0;i<rPeps->size();i++){
    if( (rPeps->at(i).lastRTime - rPeps->at(i).firstRTime) > CONTAM){
			contam->push_back(rPeps->at(i));
    } else {
      tmpPeps.push_back(rPeps->at(i));
		}
	}
	
	cout << "After removing contaminants: " << tmpPeps.size() << endl;

  delete rPeps;
  rPeps = new vector<realPep>;
  for(i=0;i<tmpPeps.size();i++){
    if(tmpPeps.at(i).monoMass<MINMASS || tmpPeps.at(i).monoMass>MAXMASS){
			contam->push_back(tmpPeps.at(i));
    } else {
      rPeps->push_back(tmpPeps.at(i));
		}
	}

	cout << "From range " << MINMASS << " to " << MAXMASS << ": " << rPeps->size() << endl;

	//Heading line
	if(fn[0]!=0){
		fprintf(out,"File\tFirst Scan\tLast Scan\tNum of Scans\tBestScan\tCharge\tMonoisotopic Mass\tBase Isotope Peak\t");
		fprintf(out,"Best Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications\n");

		for(i=0;i<rPeps->size();i++){
			fprintf(out,"%s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%f\t%f\t%f\t%f\t%f\t%lf\t%s\n",rPeps->at(i).file,
																																			 rPeps->at(i).lowScan,
																																			 rPeps->at(i).highScan,
																																			 rPeps->at(i).scanCount,rPeps->at(i).bestScan,
																																			 rPeps->at(i).charge,
																																			 rPeps->at(i).monoMass,
																																			 rPeps->at(i).basePeak,
																																			 rPeps->at(i).intensity,
																																			 rPeps->at(i).sumIntensity,
																																			 rPeps->at(i).firstRTime,
																																			 rPeps->at(i).lastRTime,
																																			 rPeps->at(i).rTime,
																																			 rPeps->at(i).xCorr,
																																			 rPeps->at(i).mods);
		};
  
		fclose(out);

	};

};

int cPersistPep::size(){
	return rPeps->size();
};

void cPersistPep::sortBasePeak(){
	qsort(&rPeps->at(0),rPeps->size(),sizeof(realPep),compareBP);
};

void cPersistPep::sortMonoMass(){
	qsort(&rPeps->at(0),rPeps->size(),sizeof(realPep),compareMM);
};

int cPersistPep::compareBP(const void *p1, const void *p2){
  const realPep d1 = *(realPep *)p1;
  const realPep d2 = *(realPep *)p2;
  if(d1.basePeak<d2.basePeak) return -1;
  else if(d1.basePeak>d2.basePeak) return 1;
  else return 0;
};

int cPersistPep::compareMM(const void *p1, const void *p2){
  const realPep d1 = *(realPep *)p1;
  const realPep d2 = *(realPep *)p2;
  if(d1.monoMass<d2.monoMass) return -1;
  else if(d1.monoMass>d2.monoMass) return 1;
  else return 0;
};
