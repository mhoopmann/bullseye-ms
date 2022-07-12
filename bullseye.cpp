#include "CKronik2.h"
#include "Spectrum.h"
#include "MSReader.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace MSToolkit;

typedef struct sBullseyeMatch {
  int index;
  int scanNumber;
} sBullseyeMatch;

typedef struct sScanInfo {
  int msLevel;
  int scanNumber;
  double mz;
  float rTime;
  //sPrecursor[] precursor;
  //sScanInfo() {
  //  precursor = new sPrecursor[100];
  //}
} sScanInfo;

string averagine(double mass);
MSFileFormat getFileFormat(char* c);
void matchMS2(CKronik2& p, char* ms2File, char* outFile, char* outFile2);
void usage();

bool compareMatch(sBullseyeMatch& p1, sBullseyeMatch& p2){
  if(p1.index==p2.index) return p1.scanNumber<p2.scanNumber;
  else return p1.index<p2.index;
}

double mean,stD;
double ppmTolerance;
double rtTolerance;
string sSumFile;
string sID;
string sDate;
string sHKFile;
string sDataFile;
bool bMatchPrecursorOnly;


int main(int argc, char* argv[]){

  int i;
  CKronik2 p1;

  sID = "Bullseye v1.32";
  sDate = "Jul 12 2022";

  cout << "Bullseye, " << sID << ", " << sDate << endl;
  cout << "Copyright 2008-2022 Mike Hoopmann, Ed Hsieh, Mike MacCoss" << endl;
  cout << "University of Washington" << endl;

  //Set global variables
  ppmTolerance=10.0;
  rtTolerance=0.5;
  bMatchPrecursorOnly=false;

	//Set default parameters for some options
	double contam=2.0;
	double maxMass=8000.0;
	double minMass=600.0;
  sSumFile.clear();

  if(argc<5){
    usage();
		exit(0);
  }

	if(argc>5){
		for(i=1;i<argc-4;i+=2){
			if(argv[i][0]!='-') {
				cout << "Invalid flag!\n" << endl;
				usage();
				exit(0);
			}
			switch(argv[i][1]){
      case 'c':
        contam=atof(argv[i+1]);
				break;
      case 'e':
        bMatchPrecursorOnly=true;
        i--;
        break;
      case 'g':
        p1.setGapTol(atoi(argv[i+1]));
        break;
      case 'm':
        maxMass=atof(argv[i+1]);
				break;
      case 'n':
        minMass=atof(argv[i+1]);
				break;
      case 'o':
        sSumFile = argv[i+1];
        break;
			case 'p':
				ppmTolerance=atof(argv[i+1]);
				break;
      case 'r':
        p1.setPPMTol(atof(argv[i+1]));
				break;
      case 's':
        p1.setMatchTol(atoi(argv[i+1])+1);
        break;
			case 't':
				rtTolerance=atof(argv[i+1]);
				break;
			default:
				cout << "Invalid flag!\n" << endl;
				usage();
				exit(0);
				break;
			}
		}
	}

  sDataFile = argv[argc - 3];
  sHKFile = argv[argc - 4];
	p1.processHK(argv[argc-4]);

	//Remove contaminants
	CKronik2 tKro;
	for(i=0;i<(int)p1.size();i++){
		if( (p1[i].lastRTime-p1[i].firstRTime) <= contam) tKro.add(p1[i]);
	}
	p1=tKro;
	cout << "Persistent peptides after removing contaminants: " << p1.size() << endl;

	//Filter based on size
	tKro.clear();
	for(i=0;i<(int)p1.size();i++){
		if( p1[i].monoMass>=minMass && p1[i].monoMass<=maxMass) tKro.add(p1[i]);
	}
	p1=tKro;
	cout << "Persistent peptides from " << minMass << " to " << maxMass << ": " << p1.size() << endl;

  //p1.findPeps();
  matchMS2(p1,argv[argc-3],argv[argc-2],argv[argc-1]);

  return 0;
  
}

string averagine(double mass) {
  string st;
  int aa = (int)(mass / 111.2137);
  int c = (int)(aa * 4.9558 + 0.5);
  int h = (int)(aa * 7.8241 + 0.5);
  int n = (int)(aa * 1.3571 + 0.5);
  int o = (int)(aa * 1.4716 + 0.5);
  int s = (int)(aa * 0.0390 + 0.5);

  double tot = c * 12 + h * 1.007825 + n * 14.003074 + o * 15.9949141 + s * 31.97207;

  h += (int)((mass - tot) / 1.007285 + 0.5);
  st = "C" + to_string(c) + "H" + to_string(h) + "N" + to_string(n) + "O" + to_string(o);
  if (s > 0) st += "S" + to_string(s);

  return st;
}

void matchMS2(CKronik2& p, char* ms2File, char* outFile, char* outFile2){

  Spectrum s;
  MSReader r,rPos,rNeg;
  MSObject o,o2;
  int i,j;
  int fragCount=0;
  int lookup[8001];
  double lowMass, highMass,ppm;
  int x,z;
  int a,b;
  int c=0;
  int d=0;
  int iPercent=0;
  int index;
  vector<int> vI;
  vector<int> vHit;
  vector<sBullseyeMatch> vMatch;
  MSFileFormat posFF, negFF;

  int ch[10];
  for(i=0;i<10;i++) ch[i]=0;

  char str1[16];
  char str2[16];

  //Check file formats for output. Make sure the user specifies the appropriate format
  posFF=getFileFormat(outFile);
  negFF=getFileFormat(outFile2);

  switch(posFF){
    case raw:
    case mzXML:
    case zs:
    case uzs:
    case ms1:
    case bms1:
    case cms1:
    case dunno:
      cout << "Output file format (pos set) not acceptable. Choose a different format." << endl;
      exit(-2);
    default:
      break;
  }

  switch(negFF){
    case raw:
    case mzXML:
    case zs:
    case uzs:
    case ms1:
    case bms1:
    case cms1:
    case dunno:
      cout << "Output file format (neg set) not acceptable. Choose a different format." << endl;
      exit(-2);
    default:
      break;
  }

  //Hardklor results are sorted and indexed to improve speed of Bullseye
  cout << "Sorting Hardklor results...";
  p.sortBasePeak();
  cout << "Done!" << endl;

  cout << "Building lookup table...";
  j=0;
  for(i=0;i<8001;i++){
    while(p.at(j).basePeak<i){
      if(j>=p.size()-1) break;
      j++;
    }
    lookup[i]=j;
  }
  cout << "Done!" << endl;

  //Read in the data
  cout << "Matching MS/MS..." << iPercent;

  b=0;
  z=0;
  a=0;

  r.setFilter(MS2);
  r.readFile(ms2File,s);

  o.setHeader(r.getHeader());
  o2.setHeader(r.getHeader());
  string sh = "FileGenerator\t"+sID+"\n\0";
  o.addToHeader(sh.c_str());
  o2.addToHeader(sh.c_str());

  rPos.setPrecisionMZ(4);
  rNeg.setPrecisionMZ(4);
  rPos.setHighResMGF(true);
  rNeg.setHighResMGF(false);
  rPos.writeFile(outFile,posFF,o);
  rNeg.writeFile(outFile2,negFF,o2);

  while(s.getScanNumber()>0){

    j=(int)(s.getMZ()+0.5);
    x=0;
    vHit.clear();
    sScanInfo scanInfo;
    scanInfo.scanNumber=s.getScanNumber();
    scanInfo.rTime=s.getRTime();
    scanInfo.msLevel=s.getMsLevel();
    scanInfo.mz=s.getMZ();
		
    //see if we can pick it up on base peak alone
    for(i=lookup[j-1];i<=lookup[j+1];i++){
      ppm = (p.at(i).basePeak-s.getMZ())/s.getMZ()*1000000;
      if( fabs(ppm)<ppmTolerance &&
          s.getRTime() > p.at(i).firstRTime-rtTolerance &&
          s.getRTime() < p.at(i).lastRTime+rtTolerance ) {
        x++;
        index=i;
        vHit.push_back(i);
      }
    }

    //if base peak wasn't enough, perhaps a different peak was isolated
    if(!bMatchPrecursorOnly){
      for(i=0;i<p.size();i++){
        lowMass = (p.at(i).monoMass+p.at(i).charge*1.00727649)/p.at(i).charge-0.05;
        switch(p.at(i).charge){
          case 1:
            highMass = (p.at(i).monoMass+p.at(i).charge*1.00727649)/p.at(i).charge + 3.10;
            break;
          case 2:
            highMass = (p.at(i).monoMass+p.at(i).charge*1.00727649)/p.at(i).charge + 2.10;
            break;
          default:
            highMass = (p.at(i).monoMass+p.at(i).charge*1.00727649)/p.at(i).charge + 4/p.at(i).charge +0.05;
            break;
        }
        if( s.getMZ() > lowMass &&
            s.getMZ() < highMass &&
            s.getRTime() > p.at(i).firstRTime-rtTolerance &&
            s.getRTime() < p.at(i).lastRTime+rtTolerance ) {
          x++;
          index=i;
          vHit.push_back(i);
        }
      }
    }

    vI.push_back(x);
    s.setFileType(MS2);

    if(x==0) {
      z++;
      o2.add(s);
      if(o2.size()>500){
        rNeg.appendFile(outFile2,o2);
        o2.clear();
      }
      ch[0]++;
    } else if(x==1) {
      a++;
      while(s.sizeZ()>0) s.eraseZ(0);
      s.addZState(p.at(index).charge,p.at(index).monoMass+1.00727649);
      s.addEZState(p.at(index).charge,p.at(index).monoMass+1.00727649,p.at(index).rTime,p.at(index).sumIntensity);
      o.add(s);
      if(o.size()>500){
        rPos.appendFile(outFile,o);
        o.clear();
      }
      c++;
    } else {
      while(s.sizeZ()>0) s.eraseZ(0);

      //erase redundancies in multiple hit list
      for(i=0;i<vHit.size()-1;i++){
        for(j=i+1;j<vHit.size();j++){
          if(p.at(vHit[i]).charge == p.at(vHit[j]).charge) {
            sprintf(str1,"%.2f\n",p.at(vHit[i]).monoMass+1.00727649);
            sprintf(str2,"%.2f\n",p.at(vHit[j]).monoMass+1.00727649);

            if(strcmp(str1,str2)==0) {
              if(p.at(vHit[i]).intensity < p.at(vHit[j]).intensity) vHit[i]=vHit[j];
              vHit.erase(vHit.begin()+j);
              j--;
            }

          }
        }
      }

      for(i=0;i<vHit.size();i++) {
        s.addZState(p.at(vHit[i]).charge,p.at(vHit[i]).monoMass+1.00727649);
        s.addEZState(p.at(vHit[i]).charge,p.at(vHit[i]).monoMass+1.00727649,p.at(vHit[i]).rTime,p.at(vHit[i]).sumIntensity);
      }

      if(vHit.size()==1) {
        a++;
        c++;
      } else {
        b++;
        d+=vHit.size();
      }

      o.add(s);
      if(o.size()>500){
        rPos.appendFile(outFile,o);
        o.clear();
      }

    }
    r.readFile(NULL,s);

    for(i=0;i<vHit.size();i++) ch[p.at(vHit[i]).charge]++;

    //export matched features
    for (i = 0; i < vHit.size(); i++)
    {
      sBullseyeMatch bm;
      bm.index=vHit[i];
      bm.scanNumber=scanInfo.scanNumber;
      vMatch.push_back(bm);
    }

    //Update file position counter
    if (r.getPercent() > iPercent){
      if(iPercent<10) cout << "\b";
      else cout << "\b\b";
      cout.flush();
      iPercent=r.getPercent();
      cout << iPercent;
      cout.flush();
    }
  }

  rPos.appendFile(outFile,o);
  rNeg.appendFile(outFile2,o2);

  cout << "Done!" << endl;

  cout << z << " scans had no visible parental distribution." << endl;
  cout << a << " scans had a single parental distribution." << endl;
  cout << b << " scans had multiple possible parental distributions." << endl;

  //Some simple tables of where parental ions are found relative to scan number
  i=0;
  a=0;
  x=0;
  z=0;
  b=1;
  cout << "\nScan:\tNeg\tPos" << endl;
  while(i<vI.size()){
    if(vI[i]==0) x++;
    else z++;
    a++;

    if(a==1000){
      cout << b*a << "\t" << x << "\t" << z << endl;
      x=0;
      z=0;
      b++;
      a=0;
    }

    i++;
  }

  cout << b*1000+a << "\t" << x << "\t" << z << endl;

  cout << "\nUnknown charge: " << ch[0] << endl;
  for(i=1;i<10;i++) cout << "+" << i << ": " << ch[i] << endl;

  cout << c << " singles and " << d << " doubles total." << endl;

  if (!sSumFile.empty()) {
    sort(vMatch.begin(),vMatch.end(),compareMatch);


    FILE* f = fopen(sSumFile.c_str(),"wt");
    //Preamble
    fprintf(f,"%s\n",sID.c_str());
    fprintf(f,"Harklor input: %s\n",sHKFile.c_str());
    fprintf(f,"MS/MS input file: %s\n",sDataFile.c_str());
    fprintf(f,"MS/MS spectra with new precursors: %s\n",outFile);
    fprintf(f,"Remaining spectra: %s\n",outFile2);
    //Heading line
    fprintf(f,"MonoisotopicMass\tz\tCharge\tStartRT\tEndRT\tApexRT\tApexAbundance\tTotalAbundance\tAveragine\tMS2Scans\n");
    for (size_t a = 0; a < vMatch.size(); a++)
    {
      fprintf(f,"%.4lf\t",p[vMatch[a].index].monoMass);
      fprintf(f,"%d\t",p[vMatch[a].index].charge);
      fprintf(f,"[M+%dH]\t",p[vMatch[a].index].charge);
      fprintf(f,"%.2f\t",p[vMatch[a].index].firstRTime);
      fprintf(f,"%.2f\t",p[vMatch[a].index].lastRTime);
      fprintf(f,"%.2f\t",p[vMatch[a].index].rTime);
      fprintf(f,"%.1f\t",p[vMatch[a].index].intensity);
      fprintf(f,"%.1f\t",p[vMatch[a].index].sumIntensity);
      fprintf(f,"%s\t",averagine(p[vMatch[a].index].monoMass).c_str());
      fprintf(f,"%d",vMatch[a].scanNumber);
      while (a < (vMatch.size() - 1) && vMatch[a + 1].index == vMatch[a].index)
      {
        a++;
        fprintf(f,";%d",vMatch[a].scanNumber);
      }
      fprintf(f,"\n");
    }
    fclose(f);
  }

}

MSFileFormat getFileFormat(char* c){

	char file[256];
	char ext[256];
	char *tok;

	strcpy(file,c);
	tok=strtok(file,".\n");
	while(tok!=NULL){
		strcpy(ext,tok);
		tok=strtok(NULL,".\n");
	}

	if(strcmp(ext,"ms1")==0 || strcmp(ext,"MS1")==0) return ms1;
	if(strcmp(ext,"ms2")==0 || strcmp(ext,"MS2")==0) return ms2;
	if(strcmp(ext,"bms1")==0 || strcmp(ext,"BMS1")==0) return bms1;
	if(strcmp(ext,"bms2")==0 || strcmp(ext,"BMS2")==0) return bms2;
	if(strcmp(ext,"cms1")==0 || strcmp(ext,"CMS1")==0) return cms1;
	if(strcmp(ext,"cms2")==0 || strcmp(ext,"CMS2")==0) return cms2;
	if(strcmp(ext,"zs")==0 || strcmp(ext,"ZS")==0) return zs;
	if(strcmp(ext,"uzs")==0 || strcmp(ext,"UZS")==0) return uzs;
  if(strcmp(ext,"mgf")==0 || strcmp(ext,"MGF")==0) return mgf;
	if(strcmp(ext,"mzXML")==0 || strcmp(ext,"MZXML")==0) return mzXML;
	if(strcmp(ext,"mzML")==0 || strcmp(ext,"MZML")==0) return mzML;
  if(strcmp(ext,"raw")==0 || strcmp(ext,"RAW")==0) return raw;
	return dunno;

}

void usage(){
	cout << "Usage: bullseye [flags] <HK file> <Data file> <Pos file> <Neg file>" << endl;
	cout << "\n  HK files are Hardklör generated results files." << endl;
  cout << "  http://proteome.gs.washington.edu/software/hardklor" << endl;
	cout << "\n  Data files contain the MS/MS data to be used." << endl;
  cout << "  Acceptable formats:" << endl;
  cout << "    Thermo RAW (Windows only), .mzXML" << endl;
  cout << "    .ms2, .bms2, .cms2 (Extracted with MakeMS2)" << endl;
  cout << "  http://proteome.gs.washington.edu/software/makems2/MakeMS2.zip" << endl;
  cout << "\n  Pos and Neg files are output files for matches and non-matches, respectively." << endl;
  cout << "  Acceptable formats (specify with file extension):" << endl;
  cout << "    .ms2, .bms2, .cms2, .mgf" << endl;
  cout << "\nExample: bullseye -p 5 Peptides.hk RawData.RAW Matches.ms2 NoMatch.ms2" << endl;
	cout << "\nFlags:" << endl;
  cout << "  -c <num>  Ignore peptides that persist for this length in time.\n"
       << "            These peptides are considered contaminants.\n"
       << "            The unit of time is whatever unit is used in your data file.\n"
       << "            Default value: 2\n" << endl;
  cout << "  -e        Use exact match to precursor ion. Rather than use wide\n"
       << "            precursor boundaries, this flag forces Bullseye to match\n"
       << "            precursors to the base isotope peak identified in Hardklor.\n"
       << "            The tolerance is set with the -p flag.\n" << endl;
  cout << "  -g <num>  Gap size tolerance when checking for peptides across consecutive\n"
       << "            scans.\n"
       << "            Default value: 1\n" << endl;
  cout << "  -m <num>  Only consider peptides below this maximum mass in daltons.\n"
       << "            Default value: 8000\n" << endl;
	cout << "  -n <num>  Only consider peptides above this minimum mass in daltons.\n"
       << "            Default value: 600\n" << endl;
  cout << "  -o <file> Output tab-delimited text summary of Peptide to MS2 matches.\n" << endl;
	cout << "  -p <num>  Sets the tolerance (+/- ppm) for exact match searches.\n"
       << "            Default value: 10\n" << endl;
  cout << "  -r <num>  Sets the tolerance (+/- ppm) for finding persistent peptides.\n"
       << "            Default value: 5\n" << endl;
  cout << "  -s <num>  Number of consecutive scans over which a peptide must be\n"
       << "            observed to be considered real. Gaps in persistence are allowed\n"
       << "            when setting the -g flag.\n"
       << "            Default value: 3\n" << endl;
	cout << "  -t <num>  Sets the tolerance (+/- minutes) around the retention\n"
		   << "            time over which a peptide can be matched to the MS/MS\n"
			 << "            spectrum.\n"
       << "            Default value: 0.5\n" << endl;
	cout << "\nPlease read the README.txt file for more information on Bullseye." << endl;

}