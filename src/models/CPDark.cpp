#include <ext/alloc_traits.h> 
#include <stddef.h> 
#include <algorithm> 
#include <iostream> 
#include <memory> 
#include <iomanip>
#include <BSMPT/models/SMparam.h> 
#include "Eigen/Dense" 
#include "Eigen/Eigenvalues" 
#include "Eigen/IterativeLinearSolvers" 

#include <BSMPT/models/CPDark.h> 
#include <BSMPT/models/IncludeAllModels.h> 
#include <BSMPT/utility.h> 
using namespace Eigen;
namespace BSMPT{ 
namespace Models{ 
Class_CPDark::Class_CPDark () 
{ 
    Model = ModelID::ModelIDs::CPDARK; 
    NNeutralHiggs =5; 
    NChargedHiggs = 4; 

    NLepton = 9;
    NQuarks = 12;
    NGauge = 4;
    NHiggs = 9;
    nPar = 13; 
    nParCT = 18; 

    nVEV = 5; 
    VevOrder.resize(5); 
    VevOrder[0] = 6; 
    VevOrder[1] = 7; 
    VevOrder[2] = 8; 
    VevOrder[3] = 2; 
    VevOrder[4] = 5; 
    UseVTreeSimplified = false; 
    UseVCounterSimplified = false; 
} 

Class_CPDark::~Class_CPDark (){} 
std::vector<std::string> Class_CPDark::addLegendCT() const{
    std::vector<std::string> labels;
    labels.push_back("ms11CT"); 
    labels.push_back("ms22CT"); 
    labels.push_back("mssCT"); 
    labels.push_back("L1CT"); 
    labels.push_back("L2CT"); 
    labels.push_back("L3CT"); 
    labels.push_back("L4CT"); 
    labels.push_back("L5CT"); 
    labels.push_back("L6CT"); 
    labels.push_back("L7CT"); 
    labels.push_back("L8CT"); 
    labels.push_back("ArealCT"); 
    labels.push_back("AimagCT"); 
    labels.push_back("TCBCT"); 
    labels.push_back("TCPCT"); 
    labels.push_back("T1CT"); 
    labels.push_back("T2CT"); 
    labels.push_back("TSCT"); 
    return labels;
} 

std::vector<std::string> Class_CPDark::addLegendTemp() const{
    std::vector<std::string> labels; 
    labels.push_back("T_c"); 
    labels.push_back("v_c"); 
    labels.push_back("v_c/T_c"); 
    labels.push_back("omega1(T_c)"); 
    labels.push_back("omega2(T_c)"); 
    labels.push_back("omegaS(T_c)"); 
    labels.push_back("omegaCB(T_c)"); 
    labels.push_back("omegaCP(T_c)"); 
    return labels; 
}

std::vector<std::string> Class_CPDark::addLegendTripleCouplings() const{
    std::vector<std::string> labels; 
    std::vector<std::string> particles; 
    particles.resize(NHiggs); 
    particles.push_back("G^+");
    particles.push_back("G^-");
    particles.push_back("H^+");
    particles.push_back("H^-");
    particles.push_back("G^0");
    particles.push_back("A");
    particles.push_back("h_SM");
    particles.push_back("h_1");
    particles.push_back("h_H");
    std::string out ="Tree_";
    for(std::size_t i=0; i<NHiggs;i++)
    {
        for(std::size_t j=i;j<NHiggs;j++)
        {
            for(std::size_t k=j; k<NHiggs;k++) 
            {
                labels.push_back("Tree_"+particles.at(i)+particles.at(j)+particles.at(k));
                labels.push_back("CT_"+particles.at(i)+particles.at(j)+particles.at(k));
                labels.push_back("CW_"+particles.at(i)+particles.at(j)+particles.at(k));
            }
        }
    }
    return labels;
}

std::vector<std::string> Class_CPDark::addLegendVEV() const{
    std::vector<std::string> labels; 
    labels.push_back("omega1");
    labels.push_back("omega2");
    labels.push_back("omegaS");
    labels.push_back("omegaCB");
    labels.push_back("omegaCP");
    return labels; 
}

void Class_CPDark::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
    std::stringstream ss(linestr); 
    double tmp;
    if (UseIndexCol){
        ss >> tmp;    }
    for(int k=1; k<=13;k++)
    {
        ss >> tmp;
        if(k==1) par[0] = tmp; 
        else if(k==2) par[1] = tmp; 
        else if(k==3) par[2] = tmp; 
        else if(k==4) par[3] = tmp; 
        else if(k==5) par[4] = tmp; 
        else if(k==6) par[5] = tmp; 
        else if(k==7) par[6] = tmp; 
        else if(k==8) par[7] = tmp; 
        else if(k==9) par[8] = tmp; 
        else if(k==10) par[9] = tmp; 
        else if(k==11) par[10] = tmp; 
        else if(k==12) par[11] = tmp; 
        else if(k==13) par[12] = tmp;  
    }
    set_gen(par);
    return;
}

void Class_CPDark::set_gen(const std::vector<double>& par) {
    ms11 = par[0];
    ms22 = par[1];
    mss = par[2];
    L1 = par[3];
    L2 = par[4];
    L3 = par[5];
    L4 = par[6];
    L5 = par[7];
    L6 = par[8];
    L7 = par[9];
    L8 = par[10];
    Areal = par[11];
    Aimag = par[12];

    vh = C_vev0;
    scale = C_vev0;
  // L1 = C_MassSMHiggs*C_MassSMHiggs/(C_vev0*C_vev0);
 //   ms11 = - L1*vh*vh/0.2e1;

    vevTreeMin.resize(nVEV);
    vevTreeMin[0] =  vh;
    vevTreeMin[1] =  0;
    vevTreeMin[2] =  0;
    vevTreeMin[3] =  0;
    vevTreeMin[4] =  0;
    vevTree.resize(NHiggs);
    vevTree = MinimizeOrderVEV(vevTreeMin);
    if(!SetCurvatureDone) SetCurvatureArrays();
}

void Class_CPDark::set_CT_Pot_Par(const std::vector<double>& par){
    ms11CT = par[0]; 
    ms22CT = par[1]; 
    mssCT = par[2]; 
    L1CT = par[3]; 
    L2CT = par[4]; 
    L3CT = par[5]; 
    L4CT = par[6]; 
    L5CT = par[7]; 
    L6CT = par[8]; 
    L7CT = par[9]; 
    L8CT = par[10]; 
    ArealCT = par[11]; 
    AimagCT = par[12]; 
    TCBCT = par[13]; 
    TCPCT = par[14]; 
    T1CT = par[15]; 
    T2CT = par[16]; 
    TSCT = par[17]; 
    //ms11CT = -L1CT*vh*vh/2.;
    Curvature_Higgs_CT_L1[2]=TCBCT;
    Curvature_Higgs_CT_L1[5]=TCPCT;
    Curvature_Higgs_CT_L1[6]=T1CT;
    Curvature_Higgs_CT_L1[7]=T2CT;
    Curvature_Higgs_CT_L1[8]=TSCT;

    Curvature_Higgs_CT_L2[0][0]=ms11CT;
    Curvature_Higgs_CT_L2[1][1]=ms11CT;
    Curvature_Higgs_CT_L2[2][2]=ms22CT;
    Curvature_Higgs_CT_L2[3][3]=ms22CT;
    Curvature_Higgs_CT_L2[4][4]=ms11CT;
    Curvature_Higgs_CT_L2[5][5]=ms22CT;
    Curvature_Higgs_CT_L2[6][6]=ms11CT;
    Curvature_Higgs_CT_L2[7][7]=ms22CT;
    Curvature_Higgs_CT_L2[8][8]=mssCT;

/*
    Curvature_Higgs_CT_L2[0][0]=-L1CT*vh*vh/2.;
    Curvature_Higgs_CT_L2[1][1]=-L1CT*vh*vh/2.;
    Curvature_Higgs_CT_L2[2][2]=ms22CT;
    Curvature_Higgs_CT_L2[3][3]=ms22CT;
    Curvature_Higgs_CT_L2[4][4]=-L1CT*vh*vh/2.;
    Curvature_Higgs_CT_L2[5][5]=ms22CT;
    Curvature_Higgs_CT_L2[6][6]=-L1CT*vh*vh/2.;
    Curvature_Higgs_CT_L2[7][7]=ms22CT;
    Curvature_Higgs_CT_L2[8][8]=mssCT;
*/

    Curvature_Higgs_CT_L3[0][2][8] =ArealCT;
    Curvature_Higgs_CT_L3[0][3][8] =-AimagCT;
    Curvature_Higgs_CT_L3[0][8][2] =ArealCT;
    Curvature_Higgs_CT_L3[0][8][3] =-AimagCT;
    Curvature_Higgs_CT_L3[1][2][8] =AimagCT;
    Curvature_Higgs_CT_L3[1][3][8] =ArealCT;
    Curvature_Higgs_CT_L3[1][8][2] =AimagCT;
    Curvature_Higgs_CT_L3[1][8][3] =ArealCT;
    Curvature_Higgs_CT_L3[2][0][8] =ArealCT;
    Curvature_Higgs_CT_L3[2][1][8] =AimagCT;
    Curvature_Higgs_CT_L3[2][8][0] =ArealCT;
    Curvature_Higgs_CT_L3[2][8][1] =AimagCT;
    Curvature_Higgs_CT_L3[3][0][8] =-AimagCT;
    Curvature_Higgs_CT_L3[3][1][8] =ArealCT;
    Curvature_Higgs_CT_L3[3][8][0] =-AimagCT;
    Curvature_Higgs_CT_L3[3][8][1] =ArealCT;
    Curvature_Higgs_CT_L3[4][5][8] =ArealCT;
    Curvature_Higgs_CT_L3[4][7][8] =AimagCT;
    Curvature_Higgs_CT_L3[4][8][5] =ArealCT;
    Curvature_Higgs_CT_L3[4][8][7] =AimagCT;
    Curvature_Higgs_CT_L3[5][4][8] =ArealCT;
    Curvature_Higgs_CT_L3[5][6][8] =-AimagCT;
    Curvature_Higgs_CT_L3[5][8][4] =ArealCT;
    Curvature_Higgs_CT_L3[5][8][6] =-AimagCT;
    Curvature_Higgs_CT_L3[6][5][8] =-AimagCT;
    Curvature_Higgs_CT_L3[6][7][8] =ArealCT;
    Curvature_Higgs_CT_L3[6][8][5] =-AimagCT;
    Curvature_Higgs_CT_L3[6][8][7] =ArealCT;
    Curvature_Higgs_CT_L3[7][4][8] =AimagCT;
    Curvature_Higgs_CT_L3[7][6][8] =ArealCT;
    Curvature_Higgs_CT_L3[7][8][4] =AimagCT;
    Curvature_Higgs_CT_L3[7][8][6] =ArealCT;
    Curvature_Higgs_CT_L3[8][0][2] =ArealCT;
    Curvature_Higgs_CT_L3[8][0][3] =-AimagCT;
    Curvature_Higgs_CT_L3[8][1][2] =AimagCT;
    Curvature_Higgs_CT_L3[8][1][3] =ArealCT;
    Curvature_Higgs_CT_L3[8][2][0] =ArealCT;
    Curvature_Higgs_CT_L3[8][2][1] =AimagCT;
    Curvature_Higgs_CT_L3[8][3][0] =-AimagCT;
    Curvature_Higgs_CT_L3[8][3][1] =ArealCT;
    Curvature_Higgs_CT_L3[8][4][5] =ArealCT;
    Curvature_Higgs_CT_L3[8][4][7] =AimagCT;
    Curvature_Higgs_CT_L3[8][5][4] =ArealCT;
    Curvature_Higgs_CT_L3[8][5][6] =-AimagCT;
    Curvature_Higgs_CT_L3[8][6][5] =-AimagCT;
    Curvature_Higgs_CT_L3[8][6][7] =ArealCT;
    Curvature_Higgs_CT_L3[8][7][4] =AimagCT;
    Curvature_Higgs_CT_L3[8][7][6] =ArealCT;

    Curvature_Higgs_CT_L4[0][0][0][0]=3*L1CT;
    Curvature_Higgs_CT_L4[0][0][1][1]=L1CT;
    Curvature_Higgs_CT_L4[0][0][2][2]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[0][0][3][3]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[0][0][4][4]=L1CT;
    Curvature_Higgs_CT_L4[0][0][5][5]=L3CT;
    Curvature_Higgs_CT_L4[0][0][6][6]=L1CT;
    Curvature_Higgs_CT_L4[0][0][7][7]=L3CT;
    Curvature_Higgs_CT_L4[0][0][8][8]=L7CT;
    Curvature_Higgs_CT_L4[0][1][0][1]=L1CT;
    Curvature_Higgs_CT_L4[0][1][1][0]=L1CT;
    Curvature_Higgs_CT_L4[0][1][2][3]=L5CT;
    Curvature_Higgs_CT_L4[0][1][3][2]=L5CT;
    Curvature_Higgs_CT_L4[0][2][0][2]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[0][2][1][3]=L5CT;
    Curvature_Higgs_CT_L4[0][2][2][0]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[0][2][3][1]=L5CT;
    Curvature_Higgs_CT_L4[0][2][4][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][2][5][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][2][6][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][2][7][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][3][0][3]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[0][3][1][2]=L5CT;
    Curvature_Higgs_CT_L4[0][3][2][1]=L5CT;
    Curvature_Higgs_CT_L4[0][3][3][0]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[0][3][4][7]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][3][5][6]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[0][3][6][5]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[0][3][7][4]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][4][0][4]=L1CT;
    Curvature_Higgs_CT_L4[0][4][2][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][4][3][7]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][4][4][0]=L1CT;
    Curvature_Higgs_CT_L4[0][4][5][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][4][7][3]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][5][0][5]=L3CT;
    Curvature_Higgs_CT_L4[0][5][2][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][5][3][6]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[0][5][4][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][5][5][0]=L3CT;
    Curvature_Higgs_CT_L4[0][5][6][3]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[0][6][0][6]=L1CT;
    Curvature_Higgs_CT_L4[0][6][2][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][6][3][5]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[0][6][5][3]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[0][6][6][0]=L1CT;
    Curvature_Higgs_CT_L4[0][6][7][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][7][0][7]=L3CT;
    Curvature_Higgs_CT_L4[0][7][2][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][7][3][4]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][7][4][3]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][7][6][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[0][7][7][0]=L3CT;
    Curvature_Higgs_CT_L4[0][8][0][8]=L7CT;
    Curvature_Higgs_CT_L4[0][8][8][0]=L7CT;
    Curvature_Higgs_CT_L4[1][0][0][1]=L1CT;
    Curvature_Higgs_CT_L4[1][0][1][0]=L1CT;
    Curvature_Higgs_CT_L4[1][0][2][3]=L5CT;
    Curvature_Higgs_CT_L4[1][0][3][2]=L5CT;
    Curvature_Higgs_CT_L4[1][1][0][0]=L1CT;
    Curvature_Higgs_CT_L4[1][1][1][1]=3*L1CT;
    Curvature_Higgs_CT_L4[1][1][2][2]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[1][1][3][3]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[1][1][4][4]=L1CT;
    Curvature_Higgs_CT_L4[1][1][5][5]=L3CT;
    Curvature_Higgs_CT_L4[1][1][6][6]=L1CT;
    Curvature_Higgs_CT_L4[1][1][7][7]=L3CT;
    Curvature_Higgs_CT_L4[1][1][8][8]=L7CT;
    Curvature_Higgs_CT_L4[1][2][0][3]=L5CT;
    Curvature_Higgs_CT_L4[1][2][1][2]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[1][2][2][1]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[1][2][3][0]=L5CT;
    Curvature_Higgs_CT_L4[1][2][4][7]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[1][2][5][6]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][2][6][5]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][2][7][4]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[1][3][0][2]=L5CT;
    Curvature_Higgs_CT_L4[1][3][1][3]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[1][3][2][0]=L5CT;
    Curvature_Higgs_CT_L4[1][3][3][1]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[1][3][4][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][3][5][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][3][6][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][3][7][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][4][1][4]=L1CT;
    Curvature_Higgs_CT_L4[1][4][2][7]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[1][4][3][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][4][4][1]=L1CT;
    Curvature_Higgs_CT_L4[1][4][5][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][4][7][2]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[1][5][1][5]=L3CT;
    Curvature_Higgs_CT_L4[1][5][2][6]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][5][3][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][5][4][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][5][5][1]=L3CT;
    Curvature_Higgs_CT_L4[1][5][6][2]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][6][1][6]=L1CT;
    Curvature_Higgs_CT_L4[1][6][2][5]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][6][3][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][6][5][2]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][6][6][1]=L1CT;
    Curvature_Higgs_CT_L4[1][6][7][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][7][1][7]=L3CT;
    Curvature_Higgs_CT_L4[1][7][2][4]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[1][7][3][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][7][4][2]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[1][7][6][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[1][7][7][1]=L3CT;
    Curvature_Higgs_CT_L4[1][8][1][8]=L7CT;
    Curvature_Higgs_CT_L4[1][8][8][1]=L7CT;
    Curvature_Higgs_CT_L4[2][0][0][2]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[2][0][1][3]=L5CT;
    Curvature_Higgs_CT_L4[2][0][2][0]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[2][0][3][1]=L5CT;
    Curvature_Higgs_CT_L4[2][0][4][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][0][5][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][0][6][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][0][7][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][1][0][3]=L5CT;
    Curvature_Higgs_CT_L4[2][1][1][2]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[2][1][2][1]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[2][1][3][0]=L5CT;
    Curvature_Higgs_CT_L4[2][1][4][7]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[2][1][5][6]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][1][6][5]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][1][7][4]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[2][2][0][0]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[2][2][1][1]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[2][2][2][2]=3*L2CT;
    Curvature_Higgs_CT_L4[2][2][3][3]=L2CT;
    Curvature_Higgs_CT_L4[2][2][4][4]=L3CT;
    Curvature_Higgs_CT_L4[2][2][5][5]=L2CT;
    Curvature_Higgs_CT_L4[2][2][6][6]=L3CT;
    Curvature_Higgs_CT_L4[2][2][7][7]=L2CT;
    Curvature_Higgs_CT_L4[2][2][8][8]=L8CT;
    Curvature_Higgs_CT_L4[2][3][0][1]=L5CT;
    Curvature_Higgs_CT_L4[2][3][1][0]=L5CT;
    Curvature_Higgs_CT_L4[2][3][2][3]=L2CT;
    Curvature_Higgs_CT_L4[2][3][3][2]=L2CT;
    Curvature_Higgs_CT_L4[2][4][0][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][4][1][7]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[2][4][2][4]=L3CT;
    Curvature_Higgs_CT_L4[2][4][4][2]=L3CT;
    Curvature_Higgs_CT_L4[2][4][5][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][4][7][1]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[2][5][0][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][5][1][6]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][5][2][5]=L2CT;
    Curvature_Higgs_CT_L4[2][5][4][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][5][5][2]=L2CT;
    Curvature_Higgs_CT_L4[2][5][6][1]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][6][0][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][6][1][5]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][6][2][6]=L3CT;
    Curvature_Higgs_CT_L4[2][6][5][1]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][6][6][2]=L3CT;
    Curvature_Higgs_CT_L4[2][6][7][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][7][0][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][7][1][4]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[2][7][2][7]=L2CT;
    Curvature_Higgs_CT_L4[2][7][4][1]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[2][7][6][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[2][7][7][2]=L2CT;
    Curvature_Higgs_CT_L4[2][8][2][8]=L8CT;
    Curvature_Higgs_CT_L4[2][8][8][2]=L8CT;
    Curvature_Higgs_CT_L4[3][0][0][3]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[3][0][1][2]=L5CT;
    Curvature_Higgs_CT_L4[3][0][2][1]=L5CT;
    Curvature_Higgs_CT_L4[3][0][3][0]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[3][0][4][7]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][0][5][6]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[3][0][6][5]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[3][0][7][4]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][1][0][2]=L5CT;
    Curvature_Higgs_CT_L4[3][1][1][3]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[3][1][2][0]=L5CT;
    Curvature_Higgs_CT_L4[3][1][3][1]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[3][1][4][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][1][5][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][1][6][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][1][7][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][2][0][1]=L5CT;
    Curvature_Higgs_CT_L4[3][2][1][0]=L5CT;
    Curvature_Higgs_CT_L4[3][2][2][3]=L2CT;
    Curvature_Higgs_CT_L4[3][2][3][2]=L2CT;
    Curvature_Higgs_CT_L4[3][3][0][0]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[3][3][1][1]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[3][3][2][2]=L2CT;
    Curvature_Higgs_CT_L4[3][3][3][3]=3*L2CT;
    Curvature_Higgs_CT_L4[3][3][4][4]=L3CT;
    Curvature_Higgs_CT_L4[3][3][5][5]=L2CT;
    Curvature_Higgs_CT_L4[3][3][6][6]=L3CT;
    Curvature_Higgs_CT_L4[3][3][7][7]=L2CT;
    Curvature_Higgs_CT_L4[3][3][8][8]=L8CT;
    Curvature_Higgs_CT_L4[3][4][0][7]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][4][1][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][4][3][4]=L3CT;
    Curvature_Higgs_CT_L4[3][4][4][3]=L3CT;
    Curvature_Higgs_CT_L4[3][4][5][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][4][7][0]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][5][0][6]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[3][5][1][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][5][3][5]=L2CT;
    Curvature_Higgs_CT_L4[3][5][4][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][5][5][3]=L2CT;
    Curvature_Higgs_CT_L4[3][5][6][0]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[3][6][0][5]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[3][6][1][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][6][3][6]=L3CT;
    Curvature_Higgs_CT_L4[3][6][5][0]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[3][6][6][3]=L3CT;
    Curvature_Higgs_CT_L4[3][6][7][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][7][0][4]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][7][1][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][7][3][7]=L2CT;
    Curvature_Higgs_CT_L4[3][7][4][0]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][7][6][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[3][7][7][3]=L2CT;
    Curvature_Higgs_CT_L4[3][8][3][8]=L8CT;
    Curvature_Higgs_CT_L4[3][8][8][3]=L8CT;
    Curvature_Higgs_CT_L4[4][0][0][4]=L1CT;
    Curvature_Higgs_CT_L4[4][0][2][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][0][3][7]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][0][4][0]=L1CT;
    Curvature_Higgs_CT_L4[4][0][5][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][0][7][3]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][1][1][4]=L1CT;
    Curvature_Higgs_CT_L4[4][1][2][7]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[4][1][3][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][1][4][1]=L1CT;
    Curvature_Higgs_CT_L4[4][1][5][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][1][7][2]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[4][2][0][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][2][1][7]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[4][2][2][4]=L3CT;
    Curvature_Higgs_CT_L4[4][2][4][2]=L3CT;
    Curvature_Higgs_CT_L4[4][2][5][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][2][7][1]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[4][3][0][7]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][3][1][5]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][3][3][4]=L3CT;
    Curvature_Higgs_CT_L4[4][3][4][3]=L3CT;
    Curvature_Higgs_CT_L4[4][3][5][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][3][7][0]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][4][0][0]=L1CT;
    Curvature_Higgs_CT_L4[4][4][1][1]=L1CT;
    Curvature_Higgs_CT_L4[4][4][2][2]=L3CT;
    Curvature_Higgs_CT_L4[4][4][3][3]=L3CT;
    Curvature_Higgs_CT_L4[4][4][4][4]=3*L1CT;
    Curvature_Higgs_CT_L4[4][4][5][5]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[4][4][6][6]=L1CT;
    Curvature_Higgs_CT_L4[4][4][7][7]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[4][4][8][8]=L7CT;
    Curvature_Higgs_CT_L4[4][5][0][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][5][1][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][5][2][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][5][3][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][5][4][5]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[4][5][5][4]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[4][5][6][7]=L5CT;
    Curvature_Higgs_CT_L4[4][5][7][6]=L5CT;
    Curvature_Higgs_CT_L4[4][6][4][6]=L1CT;
    Curvature_Higgs_CT_L4[4][6][5][7]=L5CT;
    Curvature_Higgs_CT_L4[4][6][6][4]=L1CT;
    Curvature_Higgs_CT_L4[4][6][7][5]=L5CT;
    Curvature_Higgs_CT_L4[4][7][0][3]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][7][1][2]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[4][7][2][1]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[4][7][3][0]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[4][7][4][7]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[4][7][5][6]=L5CT;
    Curvature_Higgs_CT_L4[4][7][6][5]=L5CT;
    Curvature_Higgs_CT_L4[4][7][7][4]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[4][8][4][8]=L7CT;
    Curvature_Higgs_CT_L4[4][8][8][4]=L7CT;
    Curvature_Higgs_CT_L4[5][0][0][5]=L3CT;
    Curvature_Higgs_CT_L4[5][0][2][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][0][3][6]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[5][0][4][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][0][5][0]=L3CT;
    Curvature_Higgs_CT_L4[5][0][6][3]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[5][1][1][5]=L3CT;
    Curvature_Higgs_CT_L4[5][1][2][6]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][1][3][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][1][4][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][1][5][1]=L3CT;
    Curvature_Higgs_CT_L4[5][1][6][2]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][2][0][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][2][1][6]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][2][2][5]=L2CT;
    Curvature_Higgs_CT_L4[5][2][4][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][2][5][2]=L2CT;
    Curvature_Higgs_CT_L4[5][2][6][1]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][3][0][6]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[5][3][1][4]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][3][3][5]=L2CT;
    Curvature_Higgs_CT_L4[5][3][4][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][3][5][3]=L2CT;
    Curvature_Higgs_CT_L4[5][3][6][0]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[5][4][0][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][4][1][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][4][2][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][4][3][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][4][4][5]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[5][4][5][4]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[5][4][6][7]=L5CT;
    Curvature_Higgs_CT_L4[5][4][7][6]=L5CT;
    Curvature_Higgs_CT_L4[5][5][0][0]=L3CT;
    Curvature_Higgs_CT_L4[5][5][1][1]=L3CT;
    Curvature_Higgs_CT_L4[5][5][2][2]=L2CT;
    Curvature_Higgs_CT_L4[5][5][3][3]=L2CT;
    Curvature_Higgs_CT_L4[5][5][4][4]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[5][5][5][5]=3*L2CT;
    Curvature_Higgs_CT_L4[5][5][6][6]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[5][5][7][7]=L2CT;
    Curvature_Higgs_CT_L4[5][5][8][8]=L8CT;
    Curvature_Higgs_CT_L4[5][6][0][3]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[5][6][1][2]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][6][2][1]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[5][6][3][0]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[5][6][4][7]=L5CT;
    Curvature_Higgs_CT_L4[5][6][5][6]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[5][6][6][5]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[5][6][7][4]=L5CT;
    Curvature_Higgs_CT_L4[5][7][4][6]=L5CT;
    Curvature_Higgs_CT_L4[5][7][5][7]=L2CT;
    Curvature_Higgs_CT_L4[5][7][6][4]=L5CT;
    Curvature_Higgs_CT_L4[5][7][7][5]=L2CT;
    Curvature_Higgs_CT_L4[5][8][5][8]=L8CT;
    Curvature_Higgs_CT_L4[5][8][8][5]=L8CT;
    Curvature_Higgs_CT_L4[6][0][0][6]=L1CT;
    Curvature_Higgs_CT_L4[6][0][2][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][0][3][5]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[6][0][5][3]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[6][0][6][0]=L1CT;
    Curvature_Higgs_CT_L4[6][0][7][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][1][1][6]=L1CT;
    Curvature_Higgs_CT_L4[6][1][2][5]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][1][3][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][1][5][2]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][1][6][1]=L1CT;
    Curvature_Higgs_CT_L4[6][1][7][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][2][0][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][2][1][5]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][2][2][6]=L3CT;
    Curvature_Higgs_CT_L4[6][2][5][1]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][2][6][2]=L3CT;
    Curvature_Higgs_CT_L4[6][2][7][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][3][0][5]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[6][3][1][7]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][3][3][6]=L3CT;
    Curvature_Higgs_CT_L4[6][3][5][0]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[6][3][6][3]=L3CT;
    Curvature_Higgs_CT_L4[6][3][7][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][4][4][6]=L1CT;
    Curvature_Higgs_CT_L4[6][4][5][7]=L5CT;
    Curvature_Higgs_CT_L4[6][4][6][4]=L1CT;
    Curvature_Higgs_CT_L4[6][4][7][5]=L5CT;
    Curvature_Higgs_CT_L4[6][5][0][3]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[6][5][1][2]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][5][2][1]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][5][3][0]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[6][5][4][7]=L5CT;
    Curvature_Higgs_CT_L4[6][5][5][6]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[6][5][6][5]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[6][5][7][4]=L5CT;
    Curvature_Higgs_CT_L4[6][6][0][0]=L1CT;
    Curvature_Higgs_CT_L4[6][6][1][1]=L1CT;
    Curvature_Higgs_CT_L4[6][6][2][2]=L3CT;
    Curvature_Higgs_CT_L4[6][6][3][3]=L3CT;
    Curvature_Higgs_CT_L4[6][6][4][4]=L1CT;
    Curvature_Higgs_CT_L4[6][6][5][5]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[6][6][6][6]=3*L1CT;
    Curvature_Higgs_CT_L4[6][6][7][7]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[6][6][8][8]=L7CT;
    Curvature_Higgs_CT_L4[6][7][0][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][7][1][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][7][2][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][7][3][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[6][7][4][5]=L5CT;
    Curvature_Higgs_CT_L4[6][7][5][4]=L5CT;
    Curvature_Higgs_CT_L4[6][7][6][7]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[6][7][7][6]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[6][8][6][8]=L7CT;
    Curvature_Higgs_CT_L4[6][8][8][6]=L7CT;
    Curvature_Higgs_CT_L4[7][0][0][7]=L3CT;
    Curvature_Higgs_CT_L4[7][0][2][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][0][3][4]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][0][4][3]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][0][6][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][0][7][0]=L3CT;
    Curvature_Higgs_CT_L4[7][1][1][7]=L3CT;
    Curvature_Higgs_CT_L4[7][1][2][4]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[7][1][3][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][1][4][2]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[7][1][6][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][1][7][1]=L3CT;
    Curvature_Higgs_CT_L4[7][2][0][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][2][1][4]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[7][2][2][7]=L2CT;
    Curvature_Higgs_CT_L4[7][2][4][1]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[7][2][6][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][2][7][2]=L2CT;
    Curvature_Higgs_CT_L4[7][3][0][4]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][3][1][6]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][3][3][7]=L2CT;
    Curvature_Higgs_CT_L4[7][3][4][0]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][3][6][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][3][7][3]=L2CT;
    Curvature_Higgs_CT_L4[7][4][0][3]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][4][1][2]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[7][4][2][1]=(L4CT - L5CT)/2.;
    Curvature_Higgs_CT_L4[7][4][3][0]=(-L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][4][4][7]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[7][4][5][6]=L5CT;
    Curvature_Higgs_CT_L4[7][4][6][5]=L5CT;
    Curvature_Higgs_CT_L4[7][4][7][4]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[7][5][4][6]=L5CT;
    Curvature_Higgs_CT_L4[7][5][5][7]=L2CT;
    Curvature_Higgs_CT_L4[7][5][6][4]=L5CT;
    Curvature_Higgs_CT_L4[7][5][7][5]=L2CT;
    Curvature_Higgs_CT_L4[7][6][0][2]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][6][1][3]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][6][2][0]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][6][3][1]=(L4CT + L5CT)/2.;
    Curvature_Higgs_CT_L4[7][6][4][5]=L5CT;
    Curvature_Higgs_CT_L4[7][6][5][4]=L5CT;
    Curvature_Higgs_CT_L4[7][6][6][7]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[7][6][7][6]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[7][7][0][0]=L3CT;
    Curvature_Higgs_CT_L4[7][7][1][1]=L3CT;
    Curvature_Higgs_CT_L4[7][7][2][2]=L2CT;
    Curvature_Higgs_CT_L4[7][7][3][3]=L2CT;
    Curvature_Higgs_CT_L4[7][7][4][4]=L3CT + L4CT - L5CT;
    Curvature_Higgs_CT_L4[7][7][5][5]=L2CT;
    Curvature_Higgs_CT_L4[7][7][6][6]=L3CT + L4CT + L5CT;
    Curvature_Higgs_CT_L4[7][7][7][7]=3*L2CT;
    Curvature_Higgs_CT_L4[7][7][8][8]=L8CT;
    Curvature_Higgs_CT_L4[7][8][7][8]=L8CT;
    Curvature_Higgs_CT_L4[7][8][8][7]=L8CT;
    Curvature_Higgs_CT_L4[8][0][0][8]=L7CT;
    Curvature_Higgs_CT_L4[8][0][8][0]=L7CT;
    Curvature_Higgs_CT_L4[8][1][1][8]=L7CT;
    Curvature_Higgs_CT_L4[8][1][8][1]=L7CT;
    Curvature_Higgs_CT_L4[8][2][2][8]=L8CT;
    Curvature_Higgs_CT_L4[8][2][8][2]=L8CT;
    Curvature_Higgs_CT_L4[8][3][3][8]=L8CT;
    Curvature_Higgs_CT_L4[8][3][8][3]=L8CT;
    Curvature_Higgs_CT_L4[8][4][4][8]=L7CT;
    Curvature_Higgs_CT_L4[8][4][8][4]=L7CT;
    Curvature_Higgs_CT_L4[8][5][5][8]=L8CT;
    Curvature_Higgs_CT_L4[8][5][8][5]=L8CT;
    Curvature_Higgs_CT_L4[8][6][6][8]=L7CT;
    Curvature_Higgs_CT_L4[8][6][8][6]=L7CT;
    Curvature_Higgs_CT_L4[8][7][7][8]=L8CT;
    Curvature_Higgs_CT_L4[8][7][8][7]=L8CT;
    Curvature_Higgs_CT_L4[8][8][0][0]=L7CT;
    Curvature_Higgs_CT_L4[8][8][1][1]=L7CT;
    Curvature_Higgs_CT_L4[8][8][2][2]=L8CT;
    Curvature_Higgs_CT_L4[8][8][3][3]=L8CT;
    Curvature_Higgs_CT_L4[8][8][4][4]=L7CT;
    Curvature_Higgs_CT_L4[8][8][5][5]=L8CT;
    Curvature_Higgs_CT_L4[8][8][6][6]=L7CT;
    Curvature_Higgs_CT_L4[8][8][7][7]=L8CT;
    Curvature_Higgs_CT_L4[8][8][8][8]=6*L6CT;

    for(std::size_t k1=0;k1<NHiggs;k1++)
    {
        for(std::size_t k2=k1;k2<NHiggs;k2++)
        {
            Curvature_Higgs_CT_L2[k2][k1] = Curvature_Higgs_CT_L2[k1][k2];
            for(std::size_t k3=k2;k3<NHiggs;k3++)
            {
                Curvature_Higgs_CT_L3[k1][k3][k2] = Curvature_Higgs_CT_L3[k1][k2][k3];
                Curvature_Higgs_CT_L3[k2][k1][k3]= Curvature_Higgs_CT_L3[k1][k2][k3];
                Curvature_Higgs_CT_L3[k2][k3][k1] = Curvature_Higgs_CT_L3[k1][k2][k3];
                Curvature_Higgs_CT_L3[k3][k1][k2] = Curvature_Higgs_CT_L3[k1][k2][k3];
                Curvature_Higgs_CT_L3[k3][k2][k1] = Curvature_Higgs_CT_L3[k1][k2][k3];
                for(std::size_t k4=k3;k4<NHiggs;k4++)
                {
                    Curvature_Higgs_CT_L4[k2][k3][k4][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k3][k4][k1][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k4][k1][k2][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k2][k1][k3][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k4][k2][k1][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k3][k4][k2][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k1][k3][k4][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k3][k2][k1][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k4][k3][k2][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k1][k4][k3][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k2][k1][k4][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k4][k2][k3][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k1][k4][k2][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k3][k1][k4][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k2][k3][k1][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k1][k3][k2][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k4][k1][k3][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k2][k4][k1][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k3][k2][k4][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k1][k2][k4][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k3][k1][k2][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k4][k3][k1][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                    Curvature_Higgs_CT_L4[k2][k4][k3][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];
                }
            }
        }
    }

    return;
}

void Class_CPDark::write() const { 
    typedef std::numeric_limits<double> dbl;
    std::cout.precision(dbl::max_digits10);
    std::cout << "Model = " << Model << " \n ";
    std::cout << "The parameters are : \n";
    std::cout <<"Renorm Scale = " << scale << "\n";
    std::cout <<"ms11 = " << ms11 << std::endl;
    std::cout <<"ms22 = " << ms22 << std::endl;
    std::cout <<"mss = " << mss << std::endl;
    std::cout <<"L1 = " << L1 << std::endl;
    std::cout <<"L2 = " << L2 << std::endl;
    std::cout <<"L3 = " << L3 << std::endl;
    std::cout <<"L4 = " << L4 << std::endl;
    std::cout <<"L5 = " << L5 << std::endl;
    std::cout <<"L6 = " << L6 << std::endl;
    std::cout <<"L7 = " << L7 << std::endl;
    std::cout <<"L8 = " << L8 << std::endl;
    std::cout <<"Areal = " << Areal << std::endl;
    std::cout <<"Aimag = " << Aimag << std::endl;
    std::cout <<"vh = " << vh << std::endl;
    std::cout << "The counterterm parameters are :  \n";
    std::cout << "ms11CT = " << ms11CT << std::endl;
    std::cout << "ms22CT = " << ms22CT << std::endl;
    std::cout << "mssCT = " << mssCT << std::endl;
    std::cout << "L1CT = " << L1CT << std::endl;
    std::cout << "L2CT = " << L2CT << std::endl;
    std::cout << "L3CT = " << L3CT << std::endl;
    std::cout << "L4CT = " << L4CT << std::endl;
    std::cout << "L5CT = " << L5CT << std::endl;
    std::cout << "L6CT = " << L6CT << std::endl;
    std::cout << "L7CT = " << L7CT << std::endl;
    std::cout << "L8CT = " << L8CT << std::endl;
    std::cout << "ArealCT = " << ArealCT << std::endl;
    std::cout << "AimagCT = " << AimagCT << std::endl;
    std::cout << "TCBCT = " << TCBCT << std::endl;
    std::cout << "TCPCT = " << TCPCT << std::endl;
    std::cout << "T1CT = " << T1CT << std::endl;
    std::cout << "T2CT = " << T2CT << std::endl;
    std::cout << "TSCT = " << TSCT << std::endl;
}


std::vector<double>Class_CPDark::calc_CT() const {
    std::vector<double> parCT;
    if(!SetCurvatureDone){
        std::string retmes = __func__;
        throw std::runtime_error(retmes);
    }
    if(!CalcCouplingsdone){
        std::string retmes = __func__;
        throw std::runtime_error(retmes);
    }
    std::vector<double> WeinbergNabla, WeinbergHesse, WeinbergThird, WeinbergFourth;
    WeinbergNabla = WeinbergFirstDerivative();
    WeinbergHesse = WeinbergSecondDerivative();
    WeinbergThird = WeinbergThirdDerivative();
    WeinbergFourth = WeinbergFourthDerivative();
    VectorXd NW(NHiggs);
    MatrixXd HW(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
    
    std::vector<std::vector<std::vector<double>>> TW(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                        std::vector<double>(NHiggs)));
    std::vector<std::vector<std::vector<std::vector<double>>>> FW(NHiggs, 
		std::vector<std::vector<std::vector<double>>>(NHiggs, std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs))));
    
    for(std::size_t i=0;i<NHiggs;i++){
    	for(std::size_t j=0; j<NHiggs;j++){
    		for(std::size_t k=0; k<NHiggs; k++){
    			TW[i][j][k] = 0.0;
    			for(std::size_t l=0; l<NHiggs; l++){
    				FW[i][j][k][l] = 0.0;
    			}
    		}
    	}
    }
    
    
    for(std::size_t i=0;i<NHiggs;i++){	
        NW(i) = WeinbergNabla[i];
        for(std::size_t j=0;j<NHiggs;j++){
            HW(i,j) = WeinbergHesse.at(j*NHiggs+i);
	    for(std::size_t k=0; k<NHiggs;k++){
	        TW[i][j][k] = WeinbergThird.at(k*NHiggs*NHiggs + j*NHiggs + i);
	        for(std::size_t l=0; l<NHiggs; l++){
	            FW[i][j][k][l] = WeinbergFourth.at(l*NHiggs*NHiggs*NHiggs + k*NHiggs*NHiggs + j*NHiggs +i);
	        }
	    }
        }
    }

//    ms11CT = par[0]; 
//    ms22CT = par[1]; 
//    mssCT = par[2]; 
//    L1CT = par[3]; 
//    L2CT = par[4]; 
//    L3CT = par[5]; 
//    L4CT = par[6]; 
//    L5CT = par[7]; 
//    L6CT = par[8]; 
//    L7CT = par[9]; 
//    L8CT = par[10]; 
//    ArealCT = par[11]; 
//    AimagCT = par[12]; 
//    TCBCT = par[13]; 
//    TCPCT = par[14]; 
//    T1CT = par[15]; 
//    T2CT = par[16]; 
//    TSCT = par[17];

    /*
    normal BSMPT constraints i.e up to second derivative
    */
    /*
    double t3 = 0., t7 = 0.;
    double t2 = 0., t6 = 0., t8 = 0.;
    parCT.push_back((double)((-0.3e1*HW(0,0) + HW(6,6))/0.2e1));
    parCT.push_back((double)(-(t3*(vh*vh))/0.2e1 - HW(2,2)));
    parCT.push_back((double)(-(t7*(vh*vh))/0.2e1 - HW(8,8)));
    parCT.push_back((double)((HW(0,0) - HW(6,6))/(vh*vh)));
    parCT.push_back((double)(t2));
    parCT.push_back((double)(t3));
    parCT.push_back((double)((2*HW(2,2) - HW(5,5) - HW(7,7))/(vh*vh)));
    parCT.push_back((double)((HW(5,5) - HW(7,7))/(vh*vh)));
    parCT.push_back((double)(t6));
    parCT.push_back((double)(t7));
    parCT.push_back((double)(t8));
    parCT.push_back((double)(-HW(7,8)/vh));
    parCT.push_back((double)(HW(5,8)/vh));
    parCT.push_back((double)(-NW(2)));
    parCT.push_back((double)( -NW(5)));
    parCT.push_back((double)( vh*HW(0,0) - NW(6)));
    parCT.push_back((double)( -NW(7)));
    parCT.push_back((double)(-NW(8)));
    */
    /*
    third and fourth derivatives !
    */
    parCT.push_back((double)(-HW(6,6) + (3*vh*TW[0][0][6])/2.));
    parCT.push_back((double)(-HW(2,2) + (vh*TW[2][2][6])/2.));
    parCT.push_back((double)(-HW(8,8) + (vh*TW[6][8][8])/2.));
    parCT.push_back((double)(-(TW[0][0][6]/vh)));
    parCT.push_back((double)(-FW[2][2][2][2]/3.));
    parCT.push_back((double)(-(TW[2][2][6]/vh)));
    parCT.push_back((double)(-(HW(2,2) - HW(7,7) - vh*TW[0][3][5])/(vh*vh)));
    parCT.push_back((double)(-(HW(2,2) - HW(7,7) + vh*TW[0][3][5])/(vh*vh)));
    parCT.push_back((double)(-FW[8][8][8][8]/6.));
    parCT.push_back((double)(-(TW[6][8][8]/vh)));
    parCT.push_back((double)(-FW[2][2][8][8]));
    parCT.push_back((double)(-(HW(7,8)/vh)));
    parCT.push_back((double)(-HW(5,8)/vh));
    parCT.push_back((double)(-NW(2)));
    parCT.push_back((double)(-NW(5)));
    parCT.push_back((double)(-vh*HW(6,6) - NW(6) - (vh*vh)*TW[0][0][6]));
    parCT.push_back((double)(-NW(7)));
    parCT.push_back((double)(-NW(8)));
    
    
    //for(std::size_t i=0; i<nParCT; i++) if(std::abs(parCT[i]) <= 1e-8) parCT[i] = 0;
     
    return parCT;
}
void Class_CPDark::TripleHiggsCouplings()
{
	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();
	//std::vector<double> HiggsOrder(NHiggs);
	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)
	// example for keeping the mass order
	//for(std::size_t i=0;i<NHiggs;i++) {
	//	HiggsOrder[i]=i;
	//}
	std::vector<double> TripleDeriv;
	std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,std::vector<double>(NHiggs)));     TripleDeriv=WeinbergThirdDerivative();
	for(std::size_t i=0;i<NHiggs;i++)
	  {
		for(std::size_t j=0;j<NHiggs;j++)
		{
		  for(std::size_t k=0;k<NHiggs;k++)
			{
			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);
			}
		}
	  }
	MatrixXd HiggsRot(NHiggs,NHiggs);
	for(std::size_t i=0;i<NHiggs;i++)
	{
		for(std::size_t j=0;j<NHiggs;j++)
		{
			HiggsRot(i,j) = HiggsRotationMatrix[i][j];
		}
	}

        int posMHCS1=0, posMHCS2=0;
        int posN[3];
        int countposN=0;
        int posG1=0, posG2=0, posG0=0;
        int posA=0;
        double testsum = 0;

        for(std::size_t i=0;i<3;i++)
        {
            testsum = std::abs(HiggsRot(i,0)) + std::abs(HiggsRot(i,1));
            if(testsum != 0) posG1 = i;
            testsum = std::abs(HiggsRot(i,2)) + std::abs(HiggsRot(i,3));
            if(testsum != 0) posG2 = i;
            testsum = std::abs(HiggsRot(i,4)) + std::abs(HiggsRot(i,5));
            if(testsum != 0) posG0 = i;
        }

        for(std::size_t i=3;i<NHiggs;i++)
        {
            testsum = std::abs(HiggsRot(i,0)) + std::abs(HiggsRot(i,1));
            if(testsum != 0) posMHCS1 = i;
            testsum = std::abs(HiggsRot(i,2)) + std::abs(HiggsRot(i,3));
            if(testsum != 0) posMHCS2 = i;
            testsum = std::abs(HiggsRot(i,6)) + std::abs(HiggsRot(i,7))+std::abs(HiggsRot(i,8));
            if(testsum != 0)
            {
                posN[countposN] = i;
                countposN++;
            }
            testsum = std::abs(HiggsRot(i,4)) + std::abs(HiggsRot(i,5));
            if(testsum != 0) posA = i;
        }

        std::vector<double> HiggsMasses;
        HiggsMasses=HiggsMassesSquared(vevTree, 0);

        double NeutralHiggs[3];
        double MSMlocal=0,MhUplocal=0,MhDownlocal=0;
        for(int i=0;i<3;i++){
            NeutralHiggs[i] = HiggsMasses[posN[i]];
        }

        for(int i=0; i<3; i++)
        {
            if(std::sqrt(NeutralHiggs[i]) < 126 and std::sqrt(NeutralHiggs[i]) > 124) MSMlocal = std::sqrt(NeutralHiggs[i]);
        }

        if(std::sqrt(NeutralHiggs[0]) == MSMlocal)
        {
            MhUplocal = std::sqrt(NeutralHiggs[2]);
            MhDownlocal = std::sqrt(NeutralHiggs[1]);
        }
        else if(std::sqrt(NeutralHiggs[1]) == MSMlocal)
        {
            MhUplocal = std::sqrt(NeutralHiggs[2]);
            MhDownlocal = std::sqrt(NeutralHiggs[0]);
        }
        else{
             MhUplocal = std::sqrt(NeutralHiggs[1]);
             MhDownlocal = std::sqrt(NeutralHiggs[0]);
        }
        if(MSMlocal > MhUplocal)
        {
            double tmp = posN[1];
            posN[1] = posN[2];
            posN[2] = tmp;
        }
        if(MSMlocal > MhDownlocal)
        {
            double tmp = posN[0];
            posN[0] = posN[1];
            posN[1] = tmp;
        }

        // this is for RNH2DM when all 3 CP-even Higgs and additional S are mixed
        // In case of CP in dark, h_SM are kept seperated, A, h1, h2 and S are mixed
        // The situtation is different --> need to understand the position chosen in the code.
        std::vector<double> HiggsOrder(NHiggs);
        HiggsOrder[0] = posG1;
        HiggsOrder[1] = posG2;
        HiggsOrder[2] = posMHCS1;
        HiggsOrder[3] = posMHCS2;
        HiggsOrder[4] = posG0;
        HiggsOrder[5] = posA;
        HiggsOrder[6] = posN[0];
        HiggsOrder[7] = posN[1];
        HiggsOrder[8] = posN[2];

	MatrixXd HiggsRotSort(NHiggs,NHiggs);
	for(std::size_t i=0;i<NHiggs;i++)
	{
		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
	}
	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
	for(std::size_t i=0;i<NHiggs;i++) {
		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
		for(std::size_t j=0;j<NHiggs;j++) {
			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
		}
	}
	for(std::size_t i=0;i<NHiggs;i++)
	{
		for(std::size_t j=0;j<NHiggs;j++)
		{
			for(std::size_t k=0;k<NHiggs;k++)
			{
			    TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
			    TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
			    TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
			    for(std::size_t l=0;l<NHiggs;l++)
			    {
				    for(std::size_t m=0;m<NHiggs;m++)
				    {
					    for(std::size_t n=0;n<NHiggs;n++)
					    {
						     double RotFac = HiggsRotSort(i,l)*HiggsRotSort(j,m)*HiggsRotSort(k,n);
						     TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac*GaugeBasis[l][m][n];
						     TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac*LambdaHiggs_3[l][m][n];
						     TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac*LambdaHiggs_3_CT[l][m][n];
					    }
				    }
			    }
			}
		}
	}
}
void Class_CPDark::SetCurvatureArrays(){
    initVectors();
    SetCurvatureDone=true;
    for(std::size_t i=0; i<NHiggs;i++) {
        HiggsVev[i] = vevTree[i];
        Curvature_Higgs_L1[i] = 0;
    }

    Curvature_Higgs_L2[0][0]=ms11;
    Curvature_Higgs_L2[1][1]=ms11;
    Curvature_Higgs_L2[2][2]=ms22;
    Curvature_Higgs_L2[3][3]=ms22;
    Curvature_Higgs_L2[4][4]=ms11;
    Curvature_Higgs_L2[5][5]=ms22;
    Curvature_Higgs_L2[6][6]=ms11;
    Curvature_Higgs_L2[7][7]=ms22;
    Curvature_Higgs_L2[8][8]=mss;

    Curvature_Higgs_L3[0][2][8] =Areal;
    Curvature_Higgs_L3[0][3][8] =-Aimag;
    Curvature_Higgs_L3[0][8][2] =Areal;
    Curvature_Higgs_L3[0][8][3] =-Aimag;
    Curvature_Higgs_L3[1][2][8] =Aimag;
    Curvature_Higgs_L3[1][3][8] =Areal;
    Curvature_Higgs_L3[1][8][2] =Aimag;
    Curvature_Higgs_L3[1][8][3] =Areal;
    Curvature_Higgs_L3[2][0][8] =Areal;
    Curvature_Higgs_L3[2][1][8] =Aimag;
    Curvature_Higgs_L3[2][8][0] =Areal;
    Curvature_Higgs_L3[2][8][1] =Aimag;
    Curvature_Higgs_L3[3][0][8] =-Aimag;
    Curvature_Higgs_L3[3][1][8] =Areal;
    Curvature_Higgs_L3[3][8][0] =-Aimag;
    Curvature_Higgs_L3[3][8][1] =Areal;
    Curvature_Higgs_L3[4][5][8] =Areal;
    Curvature_Higgs_L3[4][7][8] =Aimag;
    Curvature_Higgs_L3[4][8][5] =Areal;
    Curvature_Higgs_L3[4][8][7] =Aimag;
    Curvature_Higgs_L3[5][4][8] =Areal;
    Curvature_Higgs_L3[5][6][8] =-Aimag;
    Curvature_Higgs_L3[5][8][4] =Areal;
    Curvature_Higgs_L3[5][8][6] =-Aimag;
    Curvature_Higgs_L3[6][5][8] =-Aimag;
    Curvature_Higgs_L3[6][7][8] =Areal;
    Curvature_Higgs_L3[6][8][5] =-Aimag;
    Curvature_Higgs_L3[6][8][7] =Areal;
    Curvature_Higgs_L3[7][4][8] =Aimag;
    Curvature_Higgs_L3[7][6][8] =Areal;
    Curvature_Higgs_L3[7][8][4] =Aimag;
    Curvature_Higgs_L3[7][8][6] =Areal;
    Curvature_Higgs_L3[8][0][2] =Areal;
    Curvature_Higgs_L3[8][0][3] =-Aimag;
    Curvature_Higgs_L3[8][1][2] =Aimag;
    Curvature_Higgs_L3[8][1][3] =Areal;
    Curvature_Higgs_L3[8][2][0] =Areal;
    Curvature_Higgs_L3[8][2][1] =Aimag;
    Curvature_Higgs_L3[8][3][0] =-Aimag;
    Curvature_Higgs_L3[8][3][1] =Areal;
    Curvature_Higgs_L3[8][4][5] =Areal;
    Curvature_Higgs_L3[8][4][7] =Aimag;
    Curvature_Higgs_L3[8][5][4] =Areal;
    Curvature_Higgs_L3[8][5][6] =-Aimag;
    Curvature_Higgs_L3[8][6][5] =-Aimag;
    Curvature_Higgs_L3[8][6][7] =Areal;
    Curvature_Higgs_L3[8][7][4] =Aimag;
    Curvature_Higgs_L3[8][7][6] =Areal;

    Curvature_Higgs_L4[0][0][0][0]=3*L1;
    Curvature_Higgs_L4[0][0][1][1]=L1;
    Curvature_Higgs_L4[0][0][2][2]=L3 + L4 + L5;
    Curvature_Higgs_L4[0][0][3][3]=L3 + L4 - L5;
    Curvature_Higgs_L4[0][0][4][4]=L1;
    Curvature_Higgs_L4[0][0][5][5]=L3;
    Curvature_Higgs_L4[0][0][6][6]=L1;
    Curvature_Higgs_L4[0][0][7][7]=L3;
    Curvature_Higgs_L4[0][0][8][8]=L7;
    Curvature_Higgs_L4[0][1][0][1]=L1;
    Curvature_Higgs_L4[0][1][1][0]=L1;
    Curvature_Higgs_L4[0][1][2][3]=L5;
    Curvature_Higgs_L4[0][1][3][2]=L5;
    Curvature_Higgs_L4[0][2][0][2]=L3 + L4 + L5;
    Curvature_Higgs_L4[0][2][1][3]=L5;
    Curvature_Higgs_L4[0][2][2][0]=L3 + L4 + L5;
    Curvature_Higgs_L4[0][2][3][1]=L5;
    Curvature_Higgs_L4[0][2][4][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][2][5][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][2][6][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][2][7][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][3][0][3]=L3 + L4 - L5;
    Curvature_Higgs_L4[0][3][1][2]=L5;
    Curvature_Higgs_L4[0][3][2][1]=L5;
    Curvature_Higgs_L4[0][3][3][0]=L3 + L4 - L5;
    Curvature_Higgs_L4[0][3][4][7]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[0][3][5][6]=(L4 - L5)/2.;
    Curvature_Higgs_L4[0][3][6][5]=(L4 - L5)/2.;
    Curvature_Higgs_L4[0][3][7][4]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[0][4][0][4]=L1;
    Curvature_Higgs_L4[0][4][2][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][4][3][7]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[0][4][4][0]=L1;
    Curvature_Higgs_L4[0][4][5][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][4][7][3]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[0][5][0][5]=L3;
    Curvature_Higgs_L4[0][5][2][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][5][3][6]=(L4 - L5)/2.;
    Curvature_Higgs_L4[0][5][4][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][5][5][0]=L3;
    Curvature_Higgs_L4[0][5][6][3]=(L4 - L5)/2.;
    Curvature_Higgs_L4[0][6][0][6]=L1;
    Curvature_Higgs_L4[0][6][2][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][6][3][5]=(L4 - L5)/2.;
    Curvature_Higgs_L4[0][6][5][3]=(L4 - L5)/2.;
    Curvature_Higgs_L4[0][6][6][0]=L1;
    Curvature_Higgs_L4[0][6][7][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][7][0][7]=L3;
    Curvature_Higgs_L4[0][7][2][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][7][3][4]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[0][7][4][3]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[0][7][6][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[0][7][7][0]=L3;
    Curvature_Higgs_L4[0][8][0][8]=L7;
    Curvature_Higgs_L4[0][8][8][0]=L7;
    Curvature_Higgs_L4[1][0][0][1]=L1;
    Curvature_Higgs_L4[1][0][1][0]=L1;
    Curvature_Higgs_L4[1][0][2][3]=L5;
    Curvature_Higgs_L4[1][0][3][2]=L5;
    Curvature_Higgs_L4[1][1][0][0]=L1;
    Curvature_Higgs_L4[1][1][1][1]=3*L1;
    Curvature_Higgs_L4[1][1][2][2]=L3 + L4 - L5;
    Curvature_Higgs_L4[1][1][3][3]=L3 + L4 + L5;
    Curvature_Higgs_L4[1][1][4][4]=L1;
    Curvature_Higgs_L4[1][1][5][5]=L3;
    Curvature_Higgs_L4[1][1][6][6]=L1;
    Curvature_Higgs_L4[1][1][7][7]=L3;
    Curvature_Higgs_L4[1][1][8][8]=L7;
    Curvature_Higgs_L4[1][2][0][3]=L5;
    Curvature_Higgs_L4[1][2][1][2]=L3 + L4 - L5;
    Curvature_Higgs_L4[1][2][2][1]=L3 + L4 - L5;
    Curvature_Higgs_L4[1][2][3][0]=L5;
    Curvature_Higgs_L4[1][2][4][7]=(L4 - L5)/2.;
    Curvature_Higgs_L4[1][2][5][6]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[1][2][6][5]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[1][2][7][4]=(L4 - L5)/2.;
    Curvature_Higgs_L4[1][3][0][2]=L5;
    Curvature_Higgs_L4[1][3][1][3]=L3 + L4 + L5;
    Curvature_Higgs_L4[1][3][2][0]=L5;
    Curvature_Higgs_L4[1][3][3][1]=L3 + L4 + L5;
    Curvature_Higgs_L4[1][3][4][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][3][5][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][3][6][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][3][7][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][4][1][4]=L1;
    Curvature_Higgs_L4[1][4][2][7]=(L4 - L5)/2.;
    Curvature_Higgs_L4[1][4][3][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][4][4][1]=L1;
    Curvature_Higgs_L4[1][4][5][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][4][7][2]=(L4 - L5)/2.;
    Curvature_Higgs_L4[1][5][1][5]=L3;
    Curvature_Higgs_L4[1][5][2][6]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[1][5][3][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][5][4][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][5][5][1]=L3;
    Curvature_Higgs_L4[1][5][6][2]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[1][6][1][6]=L1;
    Curvature_Higgs_L4[1][6][2][5]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[1][6][3][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][6][5][2]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[1][6][6][1]=L1;
    Curvature_Higgs_L4[1][6][7][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][7][1][7]=L3;
    Curvature_Higgs_L4[1][7][2][4]=(L4 - L5)/2.;
    Curvature_Higgs_L4[1][7][3][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][7][4][2]=(L4 - L5)/2.;
    Curvature_Higgs_L4[1][7][6][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[1][7][7][1]=L3;
    Curvature_Higgs_L4[1][8][1][8]=L7;
    Curvature_Higgs_L4[1][8][8][1]=L7;
    Curvature_Higgs_L4[2][0][0][2]=L3 + L4 + L5;
    Curvature_Higgs_L4[2][0][1][3]=L5;
    Curvature_Higgs_L4[2][0][2][0]=L3 + L4 + L5;
    Curvature_Higgs_L4[2][0][3][1]=L5;
    Curvature_Higgs_L4[2][0][4][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][0][5][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][0][6][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][0][7][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][1][0][3]=L5;
    Curvature_Higgs_L4[2][1][1][2]=L3 + L4 - L5;
    Curvature_Higgs_L4[2][1][2][1]=L3 + L4 - L5;
    Curvature_Higgs_L4[2][1][3][0]=L5;
    Curvature_Higgs_L4[2][1][4][7]=(L4 - L5)/2.;
    Curvature_Higgs_L4[2][1][5][6]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[2][1][6][5]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[2][1][7][4]=(L4 - L5)/2.;
    Curvature_Higgs_L4[2][2][0][0]=L3 + L4 + L5;
    Curvature_Higgs_L4[2][2][1][1]=L3 + L4 - L5;
    Curvature_Higgs_L4[2][2][2][2]=3*L2;
    Curvature_Higgs_L4[2][2][3][3]=L2;
    Curvature_Higgs_L4[2][2][4][4]=L3;
    Curvature_Higgs_L4[2][2][5][5]=L2;
    Curvature_Higgs_L4[2][2][6][6]=L3;
    Curvature_Higgs_L4[2][2][7][7]=L2;
    Curvature_Higgs_L4[2][2][8][8]=L8;
    Curvature_Higgs_L4[2][3][0][1]=L5;
    Curvature_Higgs_L4[2][3][1][0]=L5;
    Curvature_Higgs_L4[2][3][2][3]=L2;
    Curvature_Higgs_L4[2][3][3][2]=L2;
    Curvature_Higgs_L4[2][4][0][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][4][1][7]=(L4 - L5)/2.;
    Curvature_Higgs_L4[2][4][2][4]=L3;
    Curvature_Higgs_L4[2][4][4][2]=L3;
    Curvature_Higgs_L4[2][4][5][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][4][7][1]=(L4 - L5)/2.;
    Curvature_Higgs_L4[2][5][0][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][5][1][6]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[2][5][2][5]=L2;
    Curvature_Higgs_L4[2][5][4][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][5][5][2]=L2;
    Curvature_Higgs_L4[2][5][6][1]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[2][6][0][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][6][1][5]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[2][6][2][6]=L3;
    Curvature_Higgs_L4[2][6][5][1]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[2][6][6][2]=L3;
    Curvature_Higgs_L4[2][6][7][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][7][0][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][7][1][4]=(L4 - L5)/2.;
    Curvature_Higgs_L4[2][7][2][7]=L2;
    Curvature_Higgs_L4[2][7][4][1]=(L4 - L5)/2.;
    Curvature_Higgs_L4[2][7][6][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[2][7][7][2]=L2;
    Curvature_Higgs_L4[2][8][2][8]=L8;
    Curvature_Higgs_L4[2][8][8][2]=L8;
    Curvature_Higgs_L4[3][0][0][3]=L3 + L4 - L5;
    Curvature_Higgs_L4[3][0][1][2]=L5;
    Curvature_Higgs_L4[3][0][2][1]=L5;
    Curvature_Higgs_L4[3][0][3][0]=L3 + L4 - L5;
    Curvature_Higgs_L4[3][0][4][7]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[3][0][5][6]=(L4 - L5)/2.;
    Curvature_Higgs_L4[3][0][6][5]=(L4 - L5)/2.;
    Curvature_Higgs_L4[3][0][7][4]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[3][1][0][2]=L5;
    Curvature_Higgs_L4[3][1][1][3]=L3 + L4 + L5;
    Curvature_Higgs_L4[3][1][2][0]=L5;
    Curvature_Higgs_L4[3][1][3][1]=L3 + L4 + L5;
    Curvature_Higgs_L4[3][1][4][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][1][5][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][1][6][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][1][7][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][2][0][1]=L5;
    Curvature_Higgs_L4[3][2][1][0]=L5;
    Curvature_Higgs_L4[3][2][2][3]=L2;
    Curvature_Higgs_L4[3][2][3][2]=L2;
    Curvature_Higgs_L4[3][3][0][0]=L3 + L4 - L5;
    Curvature_Higgs_L4[3][3][1][1]=L3 + L4 + L5;
    Curvature_Higgs_L4[3][3][2][2]=L2;
    Curvature_Higgs_L4[3][3][3][3]=3*L2;
    Curvature_Higgs_L4[3][3][4][4]=L3;
    Curvature_Higgs_L4[3][3][5][5]=L2;
    Curvature_Higgs_L4[3][3][6][6]=L3;
    Curvature_Higgs_L4[3][3][7][7]=L2;
    Curvature_Higgs_L4[3][3][8][8]=L8;
    Curvature_Higgs_L4[3][4][0][7]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[3][4][1][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][4][3][4]=L3;
    Curvature_Higgs_L4[3][4][4][3]=L3;
    Curvature_Higgs_L4[3][4][5][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][4][7][0]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[3][5][0][6]=(L4 - L5)/2.;
    Curvature_Higgs_L4[3][5][1][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][5][3][5]=L2;
    Curvature_Higgs_L4[3][5][4][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][5][5][3]=L2;
    Curvature_Higgs_L4[3][5][6][0]=(L4 - L5)/2.;
    Curvature_Higgs_L4[3][6][0][5]=(L4 - L5)/2.;
    Curvature_Higgs_L4[3][6][1][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][6][3][6]=L3;
    Curvature_Higgs_L4[3][6][5][0]=(L4 - L5)/2.;
    Curvature_Higgs_L4[3][6][6][3]=L3;
    Curvature_Higgs_L4[3][6][7][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][7][0][4]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[3][7][1][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][7][3][7]=L2;
    Curvature_Higgs_L4[3][7][4][0]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[3][7][6][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[3][7][7][3]=L2;
    Curvature_Higgs_L4[3][8][3][8]=L8;
    Curvature_Higgs_L4[3][8][8][3]=L8;
    Curvature_Higgs_L4[4][0][0][4]=L1;
    Curvature_Higgs_L4[4][0][2][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][0][3][7]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[4][0][4][0]=L1;
    Curvature_Higgs_L4[4][0][5][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][0][7][3]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[4][1][1][4]=L1;
    Curvature_Higgs_L4[4][1][2][7]=(L4 - L5)/2.;
    Curvature_Higgs_L4[4][1][3][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][1][4][1]=L1;
    Curvature_Higgs_L4[4][1][5][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][1][7][2]=(L4 - L5)/2.;
    Curvature_Higgs_L4[4][2][0][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][2][1][7]=(L4 - L5)/2.;
    Curvature_Higgs_L4[4][2][2][4]=L3;
    Curvature_Higgs_L4[4][2][4][2]=L3;
    Curvature_Higgs_L4[4][2][5][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][2][7][1]=(L4 - L5)/2.;
    Curvature_Higgs_L4[4][3][0][7]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[4][3][1][5]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][3][3][4]=L3;
    Curvature_Higgs_L4[4][3][4][3]=L3;
    Curvature_Higgs_L4[4][3][5][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][3][7][0]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[4][4][0][0]=L1;
    Curvature_Higgs_L4[4][4][1][1]=L1;
    Curvature_Higgs_L4[4][4][2][2]=L3;
    Curvature_Higgs_L4[4][4][3][3]=L3;
    Curvature_Higgs_L4[4][4][4][4]=3*L1;
    Curvature_Higgs_L4[4][4][5][5]=L3 + L4 + L5;
    Curvature_Higgs_L4[4][4][6][6]=L1;
    Curvature_Higgs_L4[4][4][7][7]=L3 + L4 - L5;
    Curvature_Higgs_L4[4][4][8][8]=L7;
    Curvature_Higgs_L4[4][5][0][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][5][1][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][5][2][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][5][3][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[4][5][4][5]=L3 + L4 + L5;
    Curvature_Higgs_L4[4][5][5][4]=L3 + L4 + L5;
    Curvature_Higgs_L4[4][5][6][7]=L5;
    Curvature_Higgs_L4[4][5][7][6]=L5;
    Curvature_Higgs_L4[4][6][4][6]=L1;
    Curvature_Higgs_L4[4][6][5][7]=L5;
    Curvature_Higgs_L4[4][6][6][4]=L1;
    Curvature_Higgs_L4[4][6][7][5]=L5;
    Curvature_Higgs_L4[4][7][0][3]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[4][7][1][2]=(L4 - L5)/2.;
    Curvature_Higgs_L4[4][7][2][1]=(L4 - L5)/2.;
    Curvature_Higgs_L4[4][7][3][0]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[4][7][4][7]=L3 + L4 - L5;
    Curvature_Higgs_L4[4][7][5][6]=L5;
    Curvature_Higgs_L4[4][7][6][5]=L5;
    Curvature_Higgs_L4[4][7][7][4]=L3 + L4 - L5;
    Curvature_Higgs_L4[4][8][4][8]=L7;
    Curvature_Higgs_L4[4][8][8][4]=L7;
    Curvature_Higgs_L4[5][0][0][5]=L3;
    Curvature_Higgs_L4[5][0][2][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][0][3][6]=(L4 - L5)/2.;
    Curvature_Higgs_L4[5][0][4][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][0][5][0]=L3;
    Curvature_Higgs_L4[5][0][6][3]=(L4 - L5)/2.;
    Curvature_Higgs_L4[5][1][1][5]=L3;
    Curvature_Higgs_L4[5][1][2][6]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[5][1][3][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][1][4][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][1][5][1]=L3;
    Curvature_Higgs_L4[5][1][6][2]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[5][2][0][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][2][1][6]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[5][2][2][5]=L2;
    Curvature_Higgs_L4[5][2][4][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][2][5][2]=L2;
    Curvature_Higgs_L4[5][2][6][1]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[5][3][0][6]=(L4 - L5)/2.;
    Curvature_Higgs_L4[5][3][1][4]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][3][3][5]=L2;
    Curvature_Higgs_L4[5][3][4][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][3][5][3]=L2;
    Curvature_Higgs_L4[5][3][6][0]=(L4 - L5)/2.;
    Curvature_Higgs_L4[5][4][0][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][4][1][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][4][2][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][4][3][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[5][4][4][5]=L3 + L4 + L5;
    Curvature_Higgs_L4[5][4][5][4]=L3 + L4 + L5;
    Curvature_Higgs_L4[5][4][6][7]=L5;
    Curvature_Higgs_L4[5][4][7][6]=L5;
    Curvature_Higgs_L4[5][5][0][0]=L3;
    Curvature_Higgs_L4[5][5][1][1]=L3;
    Curvature_Higgs_L4[5][5][2][2]=L2;
    Curvature_Higgs_L4[5][5][3][3]=L2;
    Curvature_Higgs_L4[5][5][4][4]=L3 + L4 + L5;
    Curvature_Higgs_L4[5][5][5][5]=3*L2;
    Curvature_Higgs_L4[5][5][6][6]=L3 + L4 - L5;
    Curvature_Higgs_L4[5][5][7][7]=L2;
    Curvature_Higgs_L4[5][5][8][8]=L8;
    Curvature_Higgs_L4[5][6][0][3]=(L4 - L5)/2.;
    Curvature_Higgs_L4[5][6][1][2]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[5][6][2][1]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[5][6][3][0]=(L4 - L5)/2.;
    Curvature_Higgs_L4[5][6][4][7]=L5;
    Curvature_Higgs_L4[5][6][5][6]=L3 + L4 - L5;
    Curvature_Higgs_L4[5][6][6][5]=L3 + L4 - L5;
    Curvature_Higgs_L4[5][6][7][4]=L5;
    Curvature_Higgs_L4[5][7][4][6]=L5;
    Curvature_Higgs_L4[5][7][5][7]=L2;
    Curvature_Higgs_L4[5][7][6][4]=L5;
    Curvature_Higgs_L4[5][7][7][5]=L2;
    Curvature_Higgs_L4[5][8][5][8]=L8;
    Curvature_Higgs_L4[5][8][8][5]=L8;
    Curvature_Higgs_L4[6][0][0][6]=L1;
    Curvature_Higgs_L4[6][0][2][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][0][3][5]=(L4 - L5)/2.;
    Curvature_Higgs_L4[6][0][5][3]=(L4 - L5)/2.;
    Curvature_Higgs_L4[6][0][6][0]=L1;
    Curvature_Higgs_L4[6][0][7][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][1][1][6]=L1;
    Curvature_Higgs_L4[6][1][2][5]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[6][1][3][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][1][5][2]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[6][1][6][1]=L1;
    Curvature_Higgs_L4[6][1][7][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][2][0][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][2][1][5]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[6][2][2][6]=L3;
    Curvature_Higgs_L4[6][2][5][1]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[6][2][6][2]=L3;
    Curvature_Higgs_L4[6][2][7][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][3][0][5]=(L4 - L5)/2.;
    Curvature_Higgs_L4[6][3][1][7]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][3][3][6]=L3;
    Curvature_Higgs_L4[6][3][5][0]=(L4 - L5)/2.;
    Curvature_Higgs_L4[6][3][6][3]=L3;
    Curvature_Higgs_L4[6][3][7][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][4][4][6]=L1;
    Curvature_Higgs_L4[6][4][5][7]=L5;
    Curvature_Higgs_L4[6][4][6][4]=L1;
    Curvature_Higgs_L4[6][4][7][5]=L5;
    Curvature_Higgs_L4[6][5][0][3]=(L4 - L5)/2.;
    Curvature_Higgs_L4[6][5][1][2]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[6][5][2][1]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[6][5][3][0]=(L4 - L5)/2.;
    Curvature_Higgs_L4[6][5][4][7]=L5;
    Curvature_Higgs_L4[6][5][5][6]=L3 + L4 - L5;
    Curvature_Higgs_L4[6][5][6][5]=L3 + L4 - L5;
    Curvature_Higgs_L4[6][5][7][4]=L5;
    Curvature_Higgs_L4[6][6][0][0]=L1;
    Curvature_Higgs_L4[6][6][1][1]=L1;
    Curvature_Higgs_L4[6][6][2][2]=L3;
    Curvature_Higgs_L4[6][6][3][3]=L3;
    Curvature_Higgs_L4[6][6][4][4]=L1;
    Curvature_Higgs_L4[6][6][5][5]=L3 + L4 - L5;
    Curvature_Higgs_L4[6][6][6][6]=3*L1;
    Curvature_Higgs_L4[6][6][7][7]=L3 + L4 + L5;
    Curvature_Higgs_L4[6][6][8][8]=L7;
    Curvature_Higgs_L4[6][7][0][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][7][1][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][7][2][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][7][3][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[6][7][4][5]=L5;
    Curvature_Higgs_L4[6][7][5][4]=L5;
    Curvature_Higgs_L4[6][7][6][7]=L3 + L4 + L5;
    Curvature_Higgs_L4[6][7][7][6]=L3 + L4 + L5;
    Curvature_Higgs_L4[6][8][6][8]=L7;
    Curvature_Higgs_L4[6][8][8][6]=L7;
    Curvature_Higgs_L4[7][0][0][7]=L3;
    Curvature_Higgs_L4[7][0][2][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][0][3][4]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[7][0][4][3]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[7][0][6][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][0][7][0]=L3;
    Curvature_Higgs_L4[7][1][1][7]=L3;
    Curvature_Higgs_L4[7][1][2][4]=(L4 - L5)/2.;
    Curvature_Higgs_L4[7][1][3][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][1][4][2]=(L4 - L5)/2.;
    Curvature_Higgs_L4[7][1][6][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][1][7][1]=L3;
    Curvature_Higgs_L4[7][2][0][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][2][1][4]=(L4 - L5)/2.;
    Curvature_Higgs_L4[7][2][2][7]=L2;
    Curvature_Higgs_L4[7][2][4][1]=(L4 - L5)/2.;
    Curvature_Higgs_L4[7][2][6][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][2][7][2]=L2;
    Curvature_Higgs_L4[7][3][0][4]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[7][3][1][6]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][3][3][7]=L2;
    Curvature_Higgs_L4[7][3][4][0]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[7][3][6][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][3][7][3]=L2;
    Curvature_Higgs_L4[7][4][0][3]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[7][4][1][2]=(L4 - L5)/2.;
    Curvature_Higgs_L4[7][4][2][1]=(L4 - L5)/2.;
    Curvature_Higgs_L4[7][4][3][0]=(-L4 + L5)/2.;
    Curvature_Higgs_L4[7][4][4][7]=L3 + L4 - L5;
    Curvature_Higgs_L4[7][4][5][6]=L5;
    Curvature_Higgs_L4[7][4][6][5]=L5;
    Curvature_Higgs_L4[7][4][7][4]=L3 + L4 - L5;
    Curvature_Higgs_L4[7][5][4][6]=L5;
    Curvature_Higgs_L4[7][5][5][7]=L2;
    Curvature_Higgs_L4[7][5][6][4]=L5;
    Curvature_Higgs_L4[7][5][7][5]=L2;
    Curvature_Higgs_L4[7][6][0][2]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][6][1][3]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][6][2][0]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][6][3][1]=(L4 + L5)/2.;
    Curvature_Higgs_L4[7][6][4][5]=L5;
    Curvature_Higgs_L4[7][6][5][4]=L5;
    Curvature_Higgs_L4[7][6][6][7]=L3 + L4 + L5;
    Curvature_Higgs_L4[7][6][7][6]=L3 + L4 + L5;
    Curvature_Higgs_L4[7][7][0][0]=L3;
    Curvature_Higgs_L4[7][7][1][1]=L3;
    Curvature_Higgs_L4[7][7][2][2]=L2;
    Curvature_Higgs_L4[7][7][3][3]=L2;
    Curvature_Higgs_L4[7][7][4][4]=L3 + L4 - L5;
    Curvature_Higgs_L4[7][7][5][5]=L2;
    Curvature_Higgs_L4[7][7][6][6]=L3 + L4 + L5;
    Curvature_Higgs_L4[7][7][7][7]=3*L2;
    Curvature_Higgs_L4[7][7][8][8]=L8;
    Curvature_Higgs_L4[7][8][7][8]=L8;
    Curvature_Higgs_L4[7][8][8][7]=L8;
    Curvature_Higgs_L4[8][0][0][8]=L7;
    Curvature_Higgs_L4[8][0][8][0]=L7;
    Curvature_Higgs_L4[8][1][1][8]=L7;
    Curvature_Higgs_L4[8][1][8][1]=L7;
    Curvature_Higgs_L4[8][2][2][8]=L8;
    Curvature_Higgs_L4[8][2][8][2]=L8;
    Curvature_Higgs_L4[8][3][3][8]=L8;
    Curvature_Higgs_L4[8][3][8][3]=L8;
    Curvature_Higgs_L4[8][4][4][8]=L7;
    Curvature_Higgs_L4[8][4][8][4]=L7;
    Curvature_Higgs_L4[8][5][5][8]=L8;
    Curvature_Higgs_L4[8][5][8][5]=L8;
    Curvature_Higgs_L4[8][6][6][8]=L7;
    Curvature_Higgs_L4[8][6][8][6]=L7;
    Curvature_Higgs_L4[8][7][7][8]=L8;
    Curvature_Higgs_L4[8][7][8][7]=L8;
    Curvature_Higgs_L4[8][8][0][0]=L7;
    Curvature_Higgs_L4[8][8][1][1]=L7;
    Curvature_Higgs_L4[8][8][2][2]=L8;
    Curvature_Higgs_L4[8][8][3][3]=L8;
    Curvature_Higgs_L4[8][8][4][4]=L7;
    Curvature_Higgs_L4[8][8][5][5]=L8;
    Curvature_Higgs_L4[8][8][6][6]=L7;
    Curvature_Higgs_L4[8][8][7][7]=L8;
    Curvature_Higgs_L4[8][8][8][8]=6*L6;

    Curvature_Gauge_G2H2[0][0][0][0]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][1][1]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][2][2]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][3][3]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][4][4]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][5][5]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][6][6]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][0][7][7]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[0][3][0][6]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][1][4]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][2][7]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][3][5]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][4][1]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][5][3]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][6][0]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[0][3][7][2]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][1][0][0]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][1][1]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][2][2]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][3][3]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][4][4]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][5][5]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][6][6]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][1][7][7]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[1][3][0][4]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][1][6]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][2][5]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][3][7]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][4][0]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][5][2]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][6][1]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[1][3][7][3]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][2][0][0]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][1][1]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][2][2]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][3][3]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][4][4]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][5][5]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][6][6]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][2][7][7]=(0.1e1/2.) * (C_g*C_g);
    Curvature_Gauge_G2H2[2][3][0][0]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][1][1]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][2][2]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][3][3]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][4][4]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][5][5]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][6][6]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[2][3][7][7]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][0][6]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][1][4]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][2][7]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][3][5]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][4][1]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][5][3]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][6][0]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][0][7][2]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][0][4]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][1][6]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][2][5]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][3][7]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][4][0]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][5][2]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][6][1]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][1][7][3]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][0][0]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][1][1]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][2][2]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][3][3]=(0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][4][4]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][5][5]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][6][6]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][2][7][7]=(-0.1e1/2.) * ((C_g) * (C_gs));
    Curvature_Gauge_G2H2[3][3][0][0]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][1][1]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][2][2]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][3][3]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][4][4]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][5][5]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][6][6]=(0.1e1/2.) * (C_gs*C_gs);
    Curvature_Gauge_G2H2[3][3][7][7]=(0.1e1/2.) * (C_gs*C_gs);


    std::complex<double> A(0,1);
    Curvature_Lepton_F2H1[0][1][4]=(0.1e1*A) * ((C_MassElectron) * (0.1e1/vh));
    Curvature_Lepton_F2H1[0][1][6]=(C_MassElectron) * (0.1e1/vh);
    Curvature_Lepton_F2H1[1][0][4]=(0.1e1*A) * ((C_MassElectron) * (0.1e1/vh));
    Curvature_Lepton_F2H1[1][0][6]=(C_MassElectron) * (0.1e1/vh);
    Curvature_Lepton_F2H1[1][6][0]=(C_MassElectron) * (0.1e1/vh);
    Curvature_Lepton_F2H1[1][6][1]=(0.1e1*A) * ((C_MassElectron) * (0.1e1/vh));
    Curvature_Lepton_F2H1[2][3][4]=(0.1e1*A) * ((C_MassMu) * (0.1e1/vh));
    Curvature_Lepton_F2H1[2][3][6]=(C_MassMu) * (0.1e1/vh);
    Curvature_Lepton_F2H1[3][2][4]=(0.1e1*A) * ((C_MassMu) * (0.1e1/vh));
    Curvature_Lepton_F2H1[3][2][6]=(C_MassMu) * (0.1e1/vh);
    Curvature_Lepton_F2H1[3][7][0]=(C_MassMu) * (0.1e1/vh);
    Curvature_Lepton_F2H1[3][7][1]=(0.1e1*A) * ((C_MassMu) * (0.1e1/vh));
    Curvature_Lepton_F2H1[4][5][4]=(0.1e1*A) * ((C_MassTau) * (0.1e1/vh));
    Curvature_Lepton_F2H1[4][5][6]=(C_MassTau) * (0.1e1/vh);
    Curvature_Lepton_F2H1[5][4][4]=(0.1e1*A) * ((C_MassTau) * (0.1e1/vh));
    Curvature_Lepton_F2H1[5][4][6]=(C_MassTau) * (0.1e1/vh);
    Curvature_Lepton_F2H1[5][8][0]=(C_MassTau) * (0.1e1/vh);
    Curvature_Lepton_F2H1[5][8][1]=(0.1e1*A) * ((C_MassTau) * (0.1e1/vh));
    Curvature_Lepton_F2H1[6][1][0]=(C_MassElectron) * (0.1e1/vh);
    Curvature_Lepton_F2H1[6][1][1]=(0.1e1*A) * ((C_MassElectron) * (0.1e1/vh));
    Curvature_Lepton_F2H1[7][3][0]=(C_MassMu) * (0.1e1/vh);
    Curvature_Lepton_F2H1[7][3][1]=(0.1e1*A) * ((C_MassMu) * (0.1e1/vh));
    Curvature_Lepton_F2H1[8][5][0]=(C_MassTau) * (0.1e1/vh);
    Curvature_Lepton_F2H1[8][5][1]=(0.1e1*A) * ((C_MassTau) * (0.1e1/vh));


    std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
    V11 = C_Vud;
    V12 = C_Vus;
    V13 = C_Vub;
    V21 = C_Vcd;
    V22 = C_Vcs;
    V23 = C_Vcb;
    V31 = C_Vtd;
    V32 = C_Vts;
    V33 = C_Vtb;
    Curvature_Quark_F2H1[0][6][4]=(-0.1e1*A) * ((C_MassUp) * (0.1e1/vh));
    Curvature_Quark_F2H1[0][6][6]=(C_MassUp) * (0.1e1/vh);
    Curvature_Quark_F2H1[0][9][0]=(-0.1e1) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V11))));
    Curvature_Quark_F2H1[0][9][1]=(0.1e1*A) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V11))));
    Curvature_Quark_F2H1[0][10][0]=(-0.1e1) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V12))));
    Curvature_Quark_F2H1[0][10][1]=(0.1e1*A) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V12))));
    Curvature_Quark_F2H1[0][11][0]=(-0.1e1) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V13))));
    Curvature_Quark_F2H1[0][11][1]=(0.1e1*A) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V13))));
    Curvature_Quark_F2H1[1][7][4]=(-0.1e1*A) * ((C_MassCharm) * (0.1e1/vh));
    Curvature_Quark_F2H1[1][7][6]=(C_MassCharm) * (0.1e1/vh);
    Curvature_Quark_F2H1[1][9][0]=(-0.1e1) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V21))));
    Curvature_Quark_F2H1[1][9][1]=(0.1e1*A) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V21))));
    Curvature_Quark_F2H1[1][10][0]=(-0.1e1) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V22))));
    Curvature_Quark_F2H1[1][10][1]=(0.1e1*A) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V22))));
    Curvature_Quark_F2H1[1][11][0]=(-0.1e1) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V23))));
    Curvature_Quark_F2H1[1][11][1]=(0.1e1*A) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V23))));
    Curvature_Quark_F2H1[2][8][4]=(-0.1e1*A) * ((C_MassTop) * (0.1e1/vh));
    Curvature_Quark_F2H1[2][8][6]=(C_MassTop) * (0.1e1/vh);
    Curvature_Quark_F2H1[2][9][0]=(-0.1e1) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V31))));
    Curvature_Quark_F2H1[2][9][1]=(0.1e1*A) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V31))));
    Curvature_Quark_F2H1[2][10][0]=(-0.1e1) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V32))));
    Curvature_Quark_F2H1[2][10][1]=(0.1e1*A) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V32))));
    Curvature_Quark_F2H1[2][11][0]=(-0.1e1) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V33))));
    Curvature_Quark_F2H1[2][11][1]=(0.1e1*A) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V33))));
    Curvature_Quark_F2H1[3][6][0]=(C_MassDown) * ((V11) * (0.1e1/vh));
    Curvature_Quark_F2H1[3][6][1]=(0.1e1*A) * ((C_MassDown) * ((V11) * (0.1e1/vh)));
    Curvature_Quark_F2H1[3][7][0]=(C_MassDown) * ((V21) * (0.1e1/vh));
    Curvature_Quark_F2H1[3][7][1]=(0.1e1*A) * ((C_MassDown) * ((V21) * (0.1e1/vh)));
    Curvature_Quark_F2H1[3][8][0]=(C_MassDown) * ((V31) * (0.1e1/vh));
    Curvature_Quark_F2H1[3][8][1]=(0.1e1*A) * ((C_MassDown) * ((V31) * (0.1e1/vh)));
    Curvature_Quark_F2H1[3][9][4]=(0.1e1*A) * ((C_MassDown) * (0.1e1/vh));
    Curvature_Quark_F2H1[3][9][6]=(C_MassDown) * (0.1e1/vh);
    Curvature_Quark_F2H1[4][6][0]=(C_MassStrange) * ((V12) * (0.1e1/vh));
    Curvature_Quark_F2H1[4][6][1]=(0.1e1*A) * ((C_MassStrange) * ((V12) * (0.1e1/vh)));
    Curvature_Quark_F2H1[4][7][0]=(C_MassStrange) * ((V22) * (0.1e1/vh));
    Curvature_Quark_F2H1[4][7][1]=(0.1e1*A) * ((C_MassStrange) * ((V22) * (0.1e1/vh)));
    Curvature_Quark_F2H1[4][8][0]=(C_MassStrange) * ((V32) * (0.1e1/vh));
    Curvature_Quark_F2H1[4][8][1]=(0.1e1*A) * ((C_MassStrange) * ((V32) * (0.1e1/vh)));
    Curvature_Quark_F2H1[4][10][4]=(0.1e1*A) * ((C_MassStrange) * (0.1e1/vh));
    Curvature_Quark_F2H1[4][10][6]=(C_MassStrange) * (0.1e1/vh);
    Curvature_Quark_F2H1[5][6][0]=(C_MassBottom) * ((V13) * (0.1e1/vh));
    Curvature_Quark_F2H1[5][6][1]=(0.1e1*A) * ((C_MassBottom) * ((V13) * (0.1e1/vh)));
    Curvature_Quark_F2H1[5][7][0]=(C_MassBottom) * ((V23) * (0.1e1/vh));
    Curvature_Quark_F2H1[5][7][1]=(0.1e1*A) * ((C_MassBottom) * ((V23) * (0.1e1/vh)));
    Curvature_Quark_F2H1[5][8][0]=(C_MassBottom) * ((V33) * (0.1e1/vh));
    Curvature_Quark_F2H1[5][8][1]=(0.1e1*A) * ((C_MassBottom) * ((V33) * (0.1e1/vh)));
    Curvature_Quark_F2H1[5][11][4]=(0.1e1*A) * ((C_MassBottom) * (0.1e1/vh));
    Curvature_Quark_F2H1[5][11][6]=(C_MassBottom) * (0.1e1/vh);
    Curvature_Quark_F2H1[6][0][4]=(-0.1e1*A) * ((C_MassUp) * (0.1e1/vh));
    Curvature_Quark_F2H1[6][0][6]=(C_MassUp) * (0.1e1/vh);
    Curvature_Quark_F2H1[6][3][0]=(C_MassDown) * ((V11) * (0.1e1/vh));
    Curvature_Quark_F2H1[6][3][1]=(0.1e1*A) * ((C_MassDown) * ((V11) * (0.1e1/vh)));
    Curvature_Quark_F2H1[6][4][0]=(C_MassStrange) * ((V12) * (0.1e1/vh));
    Curvature_Quark_F2H1[6][4][1]=(0.1e1*A) * ((C_MassStrange) * ((V12) * (0.1e1/vh)));
    Curvature_Quark_F2H1[6][5][0]=(C_MassBottom) * ((V13) * (0.1e1/vh));
    Curvature_Quark_F2H1[6][5][1]=(0.1e1*A) * ((C_MassBottom) * ((V13) * (0.1e1/vh)));
    Curvature_Quark_F2H1[7][1][4]=(-0.1e1*A) * ((C_MassCharm) * (0.1e1/vh));
    Curvature_Quark_F2H1[7][1][6]=(C_MassCharm) * (0.1e1/vh);
    Curvature_Quark_F2H1[7][3][0]=(C_MassDown) * ((V21) * (0.1e1/vh));
    Curvature_Quark_F2H1[7][3][1]=(0.1e1*A) * ((C_MassDown) * ((V21) * (0.1e1/vh)));
    Curvature_Quark_F2H1[7][4][0]=(C_MassStrange) * ((V22) * (0.1e1/vh));
    Curvature_Quark_F2H1[7][4][1]=(0.1e1*A) * ((C_MassStrange) * ((V22) * (0.1e1/vh)));
    Curvature_Quark_F2H1[7][5][0]=(C_MassBottom) * ((V23) * (0.1e1/vh));
    Curvature_Quark_F2H1[7][5][1]=(0.1e1*A) * ((C_MassBottom) * ((V23) * (0.1e1/vh)));
    Curvature_Quark_F2H1[8][2][4]=(-0.1e1*A) * ((C_MassTop) * (0.1e1/vh));
    Curvature_Quark_F2H1[8][2][6]=(C_MassTop) * (0.1e1/vh);
    Curvature_Quark_F2H1[8][3][0]=(C_MassDown) * ((V31) * (0.1e1/vh));
    Curvature_Quark_F2H1[8][3][1]=(0.1e1*A) * ((C_MassDown) * ((V31) * (0.1e1/vh)));
    Curvature_Quark_F2H1[8][4][0]=(C_MassStrange) * ((V32) * (0.1e1/vh));
    Curvature_Quark_F2H1[8][4][1]=(0.1e1*A) * ((C_MassStrange) * ((V32) * (0.1e1/vh)));
    Curvature_Quark_F2H1[8][5][0]=(C_MassBottom) * ((V33) * (0.1e1/vh));
    Curvature_Quark_F2H1[8][5][1]=(0.1e1*A) * ((C_MassBottom) * ((V33) * (0.1e1/vh)));
    Curvature_Quark_F2H1[9][0][0]=(-0.1e1) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V11))));
    Curvature_Quark_F2H1[9][0][1]=(0.1e1*A) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V11))));
    Curvature_Quark_F2H1[9][1][0]=(-0.1e1) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V21))));
    Curvature_Quark_F2H1[9][1][1]=(0.1e1*A) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V21))));
    Curvature_Quark_F2H1[9][2][0]=(-0.1e1) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V31))));
    Curvature_Quark_F2H1[9][2][1]=(0.1e1*A) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V31))));
    Curvature_Quark_F2H1[9][3][4]=(0.1e1*A) * ((C_MassDown) * (0.1e1/vh));
    Curvature_Quark_F2H1[9][3][6]=(C_MassDown) * (0.1e1/vh);
    Curvature_Quark_F2H1[10][0][0]=(-0.1e1) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V12))));
    Curvature_Quark_F2H1[10][0][1]=(0.1e1*A) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V12))));
    Curvature_Quark_F2H1[10][1][0]=(-0.1e1) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V22))));
    Curvature_Quark_F2H1[10][1][1]=(0.1e1*A) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V22))));
    Curvature_Quark_F2H1[10][2][0]=(-0.1e1) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V32))));
    Curvature_Quark_F2H1[10][2][1]=(0.1e1*A) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V32))));
    Curvature_Quark_F2H1[10][4][4]=(0.1e1*A) * ((C_MassStrange) * (0.1e1/vh));
    Curvature_Quark_F2H1[10][4][6]=(C_MassStrange) * (0.1e1/vh);
    Curvature_Quark_F2H1[11][0][0]=(-0.1e1) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V13))));
    Curvature_Quark_F2H1[11][0][1]=(0.1e1*A) * ((C_MassUp) * ((0.1e1/vh) * (std::conj(V13))));
    Curvature_Quark_F2H1[11][1][0]=(-0.1e1) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V23))));
    Curvature_Quark_F2H1[11][1][1]=(0.1e1*A) * ((C_MassCharm) * ((0.1e1/vh) * (std::conj(V23))));
    Curvature_Quark_F2H1[11][2][0]=(-0.1e1) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V33))));
    Curvature_Quark_F2H1[11][2][1]=(0.1e1*A) * ((C_MassTop) * ((0.1e1/vh) * (std::conj(V33))));
    Curvature_Quark_F2H1[11][5][4]=(0.1e1*A) * ((C_MassBottom) * (0.1e1/vh));
    Curvature_Quark_F2H1[11][5][6]=(C_MassBottom) * (0.1e1/vh);
    SetCurvatureDone = true;
}

bool Class_CPDark::CalculateDebyeSimplified(){
    return false;
}

bool Class_CPDark::CalculateDebyeGaugeSimplified(){
    return false;
}

double Class_CPDark::VTreeSimplified(const std::vector<double>& v) const
{
    (void) v;
    return 0;
}

double Class_CPDark::VCounterSimplified(const std::vector<double>& v) const
{
    (void) v;
    return 0;
}

void Class_CPDark::Debugging(const std::vector<double>& input, std::vector<double>& output) const{
    (void) input;
    (void) output;
}
}
}
