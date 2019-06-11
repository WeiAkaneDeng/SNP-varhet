//
//  main.cpp
//  METALEVENE
//
//  Created by Reedik Mägi on 17/10/2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//
//
//  1) no duplicate check!
//
//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>

#include "CmdLine.h"
#include "tools.h"
#include "statistics.h"

using namespace TCLAP;
using namespace std;

string switchStrand(string a)
{
    if (a=="A"){return "T";}
    if (a=="C"){return "G";}
    if (a=="G"){return "C";}
    if (a=="T"){return "A";}
    return "X";
}

class Marker
{
public:
    string name, ea, nea, presence;
    vector<double> n1, n2, n3;
    vector<double> z1, z2, z3;
    vector<double> var1, var2, var3;
    vector<double> n;           // N
    

    double metaLevene, metaStatsN, metaStatsD,p;          //Meta Levene's ML
   

bool calcMetaStatsN()           //Meta Levene's ML equasion
{
    
    metaStatsN=0;          //numerator in equasion
    metaStatsD=0;          //denominator in equasion
    
    n.push_back(0);n.push_back(0);n.push_back(0);n.push_back(0);        
    vector<double> mZn;                // mZis nis 
    mZn.push_back(0);mZn.push_back(0);mZn.push_back(0);mZn.push_back(0);
    vector<double> varZ;                //  
    varZ.push_back(0);varZ.push_back(0);varZ.push_back(0);varZ.push_back(0);
    vector<double> Zs;                //  
    Zs.push_back(0);Zs.push_back(0);Zs.push_back(0);Zs.push_back(0);
    vector<double> Z;                //  
    Z.push_back(0);Z.push_back(0);Z.push_back(0);Z.push_back(0);

    for (int i =0; i<n1.size(); i++)    //kohortid  
    {
        //for numerator
        
        n[1]+=n1[i];            //AA
        n[2]+=n2[i];            //AB
        n[3]+=n3[i];            //BB
        Z[1]+=z1[i];            //AA
        Z[2]+=z2[i];            //AB
        Z[3]+=z3[i];            //BB

    }
    n[0] = n[1]+n[2]+n[3];
    for (int i=0; i<n1.size(); i++)    //kohortid  
    {
        //for numerator
        
        Zs[1]+=(z1[i]*n1[i])/n[1];            //AA
        Zs[2]+=(z2[i]*n2[i])/n[2];            //AB
        Zs[3]+=(z3[i]*n3[i])/n[3];            //BB
    }  
    for (int j=1; j<=3; j++)
    {
        Zs[0]+=(Zs[j]*n[j])/n[0];
    }
    double zzz = n[0];
    
    double _sum1=0;
    for (int j=1; j<=3; j++)
    {
        _sum1+=n[j]*pow(Zs[j]-Zs[0],2);
    }
    double _sum2;
    
    for (int j=1; j<=3; j++)
    {
        _sum2+=((pow(Zs[j],2)*n[j]));
    }    
    double _sum3=0;

//        _sum3+=(varZ[j]*n[j])+(pow(Z[j],2)*n[j])-_sum2;
        for (int i=0; i<n1.size(); i++) 
        {
            _sum3+=((var1[i]*(n1[i]-1))+(pow(z1[i],2)*n1[i]));            //AA
            double yyy2 = ((var1[i]*(n1[i]-1))+(pow(z1[i],2)*n1[i]));
            _sum3+=((var2[i]*(n2[i]-1))+(pow(z2[i],2)*n2[i]));            //AB
            _sum3+=((var3[i]*(n3[i]-1))+(pow(z3[i],2)*n3[i]));
            double yyy = _sum3;
            
        }
    _sum3+=-_sum2;//
    //_sum3=(varZ[1]-_sum2[1])+(varZ[2]-_sum2[2])+(varZ[3]-_sum2[3]);      
    
    metaStatsN = (n[0]-3)*_sum1;
    metaStatsD = 2 * _sum3;
    metaLevene = metaStatsN / metaStatsD;
    p = 1-alglib::fdistribution(2,n[0]-3,metaLevene);

    return true;
}
};

bool readFile(string fileName, map <string, int> & markerPosition,vector <Marker> & markerList, int fileNr)
{
    map <string, int> markerDuplicate;
    ifstream F (fileName.c_str());
    if (F.is_open())
    {
        int _lineNr = 0;

        while (! F.eof() )
        {
            string line;
        	vector<string> tokens;
        	getline (F,line);
            int n = Tokenize(string(line), tokens, " ");
                
            if (_lineNr==0) //reading header
            {
            
            }
            else
            {

                if (n>2)
                {
                    string _name = tokens[0];
                    string _strand = tokens[1];
                    string _minor = tokens[5];
                    string _major = tokens[6];
                    double _n1 = atof(tokens[8].c_str());
                    double _n2 = atof(tokens[9].c_str());
                    double _n3 = atof(tokens[10].c_str());
                    double _z1 = atof(tokens[18].c_str())*atof(tokens[15].c_str());
                    double _z2 = atof(tokens[19].c_str())*atof(tokens[16].c_str());
                    double _z3 = atof(tokens[20].c_str())*atof(tokens[17].c_str());
                    double _varz1 = pow(atof(tokens[21].c_str())*atof(tokens[15].c_str()),2);
                    double _varz2 = pow(atof(tokens[22].c_str())*atof(tokens[16].c_str()),2);
                    double _varz3 = pow(atof(tokens[23].c_str())*atof(tokens[17].c_str()),2); 
 //                  double _z1 = atof(tokens[18].c_str());
 //                  double _z2 = atof(tokens[19].c_str());
 //                  double _z3 = atof(tokens[20].c_str());
 //                  double _varz1 = pow(atof(tokens[21].c_str()),2);
 //                  double _varz2 = pow(atof(tokens[22].c_str()),2);
 //                  double _varz3 = pow(atof(tokens[23].c_str()),2); 
                    
                    if (_strand == "-")
                    {
                        _minor = switchStrand(_minor);
                        _major = switchStrand(_major);
                    }
                                     
                    if (markerPosition[_name]==0 && markerDuplicate[_name]==0)
                    {
                        Marker myMarker;
                        myMarker.name=_name;
                        myMarker.ea=_minor;
                        myMarker.nea=_major;
                        myMarker.n1.push_back(_n1);
                        myMarker.n2.push_back(_n2);
                        myMarker.n3.push_back(_n3);
                        myMarker.z1.push_back(_z1);
                        myMarker.z2.push_back(_z2);
                        myMarker.z3.push_back(_z3);
                        myMarker.var1.push_back(_varz1);
                        myMarker.var2.push_back(_varz2);
                        myMarker.var3.push_back(_varz3);
                        myMarker.presence = "";
                        for (int i=0;i<fileNr; i++)myMarker.presence.append("0");
                        myMarker.presence.append("1");
                        
                        markerList.push_back(myMarker);
                        markerPosition[_name]=(int) markerList.size();
                    }
                    else if (markerDuplicate[_name]==0)
                    {
                        int pos = markerPosition[_name]-1;
                        
                        if (_minor == markerList[pos].ea && _major == markerList[pos].nea)
                        {
                            markerList[pos].n1.push_back(_n1);
                            markerList[pos].n2.push_back(_n2);
                            markerList[pos].n3.push_back(_n3);
                            markerList[pos].z1.push_back(_z1);
                            markerList[pos].z2.push_back(_z2);
                            markerList[pos].z3.push_back(_z3);
                            markerList[pos].var1.push_back(_varz1);
                            markerList[pos].var2.push_back(_varz2);
                            markerList[pos].var3.push_back(_varz3);
                            markerList[pos].presence.append("1");
                        }
                        else if (_major == markerList[pos].ea && _minor == markerList[pos].nea)
                        {
                            markerList[pos].n3.push_back(_n1);
                            markerList[pos].n2.push_back(_n2);
                            markerList[pos].n1.push_back(_n3);
                            markerList[pos].z3.push_back(_z1);
                            markerList[pos].z2.push_back(_z2);
                            markerList[pos].z1.push_back(_z3);
                            markerList[pos].var3.push_back(_varz1);
                            markerList[pos].var2.push_back(_varz2);
                            markerList[pos].var1.push_back(_varz3);
                            markerList[pos].presence.append("1");
                        }
                        else
                        {
                            cerr << "Cannot resolve marker " << _name << " - alleles dont match\n";
                        }
                    }
                    else
                    {
                        cerr << "Found duplicated " << _name << " from " << fileName << endl;
                    }
                    markerDuplicate[_name]=1;
                }
            }
            _lineNr++;        
        }
    }
    else {
        cerr << "Cannot access file " << fileName << ". Exit script" << endl; 
        return false;
    }
    return true;
}



int main (int argc,  char * argv[])
{

    vector<string> fileNames;
    map <string, int> markerPosition;
    vector <Marker> markerList;
    string outputFile;
	//resolveCommandLine(argc, argv, GLOBAL);
    try 
    {  
        CmdLine cmd("For more info: http://www.geenivaramu.ee/MetaLevene", ' ', "v.1.0.1");
        ValueArg<string> InputArg("i","input","Inputfile listing all cohort files",true,"","string", cmd);
        ValueArg<string> OutputArg("o","output","Output file name",true,"","string", cmd);
        cmd.parse(argc,argv);
        string inputFile = InputArg.getValue();
        outputFile = OutputArg.getValue();
        ifstream F (inputFile.c_str());
        if (F.is_open())
        {
            while (! F.eof() )
            {
                string line;
                vector <string> tokens;
                getline (F,line);
                int n = Tokenize(string(line), tokens, " ");
                if (n>0)
                {
                    fileNames.push_back(string(tokens[0]));
                }
            }
        }
        else
        {
            cerr << "Cannot open " << inputFile.c_str() << ". Exit program." << endl; 
        }
        
    }
    catch(ArgException &e) 
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }
    for (int i = 0; i < fileNames.size(); i++)
    {
        cout << fileNames[i]<<endl;
        readFile(fileNames[i], markerPosition, markerList, i);
        for (int j=0; j<markerList.size(); j++)
        {
            if (markerList[j].presence.length()<i+1) markerList[j].presence.append("0");
        }    
    }
    ofstream OUT (outputFile.c_str());
    cout << "Total marker count " << markerList.size() << endl;
    OUT << "Name\tEA\tNEA\tN\tN_0\tN_1\tN_2\tCohortCount\tMetaStatN\tMetaStatD\tMetaLevene\tP\tPresence" << endl;
    for (int i = 0; i < markerList.size(); i++)
    {
        markerList[i].calcMetaStatsN();
        OUT << markerList[i].name << "\t" << markerList[i].ea << "\t" << markerList[i].nea << "\t" 
        << markerList[i].n[0]  << "\t" << markerList[i].n[1] << "\t" << markerList[i].n[2]  << "\t" 
        << markerList[i].n[3] << "\t" << markerList[i].n1.size() << "\t" << markerList[i].metaStatsN << "\t" << markerList[i].metaStatsD << "\t" << markerList[i].metaLevene << "\t" << markerList[i].p << "\t" << markerList[i].presence <<endl;
        
    }
    return 0;
}

