/*
 * ClassPotentialOrigin.cpp
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */




#include <iomanip>
#include <gsl/gsl_sf_gamma.h>
#include <random>

#include "Eigen/Dense"

#include <gsl/gsl_sf_zeta.h>

#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/utility.h>
using namespace Eigen;


namespace BSMPT{
Class_Potential_Origin::Class_Potential_Origin ()
{
  // TODO Auto-generated constructor stub

}

Class_Potential_Origin::~Class_Potential_Origin ()
{
  // TODO Auto-generated destructor stub
}

/**
   * This will call set_gen(par), SetCurvatureArrays, set_CT_Pot_Par(parCT), CalculateDebye() as well as
   * CalculateDebyeGauge()
   */
void Class_Potential_Origin::set_All(const std::vector<double>& par,const std::vector<double>& parCT)
{


  set_gen(par);
  if(!SetCurvatureDone) SetCurvatureArrays();
  set_CT_Pot_Par(parCT);
  CalculateDebye();
  CalculateDebyeGauge();

}


void Class_Potential_Origin::Prepare_Triple()
{
    for(std::size_t a=0;a<NHiggs;a++)
	{
        for(std::size_t b=0;b<NHiggs;b++)
		{
            for(std::size_t i=0;i<NHiggs;i++)
			{
				LambdaHiggs_3_CT[a][b][i] = Curvature_Higgs_CT_L3[a][b][i];
                for(std::size_t j=0;j<NHiggs;j++)
				{
					LambdaHiggs_3_CT[a][b][i] += Curvature_Higgs_CT_L4[a][b][i][j]*HiggsVev[j];
				}
			}
		}
	}

}


double Class_Potential_Origin::FCW(double MassSquared) const
{
	double res=0;
	double x;
	double Boarder=std::pow(10,-200);
	if(std::isnan(MassSquared)) x = Boarder;
	else if(std::abs(MassSquared) < Boarder) x = Boarder;
	else x = std::abs(MassSquared);
	res = std::log(x) - 2*std::log(scale);
	return res;

}



double Class_Potential_Origin::CWTerm(double MassSquared, double cb, int diff) const
{

    if(std::abs(MassSquared) < C_threshold) return 0;
   double LogTerm=0, PotVal=0;
   LogTerm = FCW(MassSquared);
   if(diff == 0) PotVal = 1.0/(64*M_PI*M_PI) * MassSquared*MassSquared*(LogTerm-cb);
   else if(diff > 0){
    PotVal = 1.0/(32*M_PI*M_PI)*MassSquared*(LogTerm-cb+0.5);
   }
   return PotVal;
}





double Class_Potential_Origin::boson(double MassSquared, double Temp, double cb, int diff) const{
	double res = 0;
	if(diff >= 0) res = CWTerm(std::abs(MassSquared),cb,diff);
	if(Temp == 0) return res;
	double Ratio = MassSquared/std::pow(Temp,2);
	if(diff == 0) {
        res += std::pow(Temp,4)/(2*std::pow(M_PI,2)) * ThermalFunctions::JbosonInterpolated(Ratio);
	}
	else if(diff == 1){
        res += std::pow(Temp,2)/(2*std::pow(M_PI,2)) * ThermalFunctions::JbosonNumericalIntegration(Ratio,1);
	}
	else if(diff == -1){
        res += 1.0/(2*std::pow(M_PI,2))*(4*std::pow(Temp,3)*ThermalFunctions::JbosonNumericalIntegration(Ratio,0)
                - 2*Temp*MassSquared * ThermalFunctions::JbosonNumericalIntegration(Ratio,1) );
	}
	return res;
}





double Class_Potential_Origin::fermion(double MassSquared, double Temp, int diff) const{
	double res = 0;
	if(diff >=0) res = CWTerm(std::abs(MassSquared),C_CWcbFermion,diff);
	double Ratio = MassSquared/std::pow(Temp,2);
	if(Temp == 0) return res;
	if(diff == 0) {
        res += std::pow(Temp,4)/(2*std::pow(M_PI,2)) * ThermalFunctions::JfermionInterpolated(Ratio);
	}
	else if(diff == 1){
        res += std::pow(Temp,2)/(2*std::pow(M_PI,2)) * ThermalFunctions::JfermionNumericalIntegration(Ratio,1);
	}
	else if(diff == -1){
        res += 1.0/(2*std::pow(M_PI,2)) *(4*std::pow(Temp,3)*ThermalFunctions::JfermionNumericalIntegration(Ratio)
                - 2*Temp*MassSquared * ThermalFunctions::JfermionNumericalIntegration(Ratio,1) );

	}
	return res;
}



std::vector<double> Class_Potential_Origin::FirstDerivativeOfEigenvalues(const Ref<MatrixXcd> M, const Ref<MatrixXcd> MDiff) const
{
    std::vector<double> res;
    const std::size_t nRows = M.rows();
    const std::size_t nCols = M.cols();

    const double EVThres = std::pow(10,-6);

    if(nCols != nRows) {
        throw std::runtime_error("ERROR ! M needs to be an quadratic Matrix for calculating the derivatives !\n");
    }

    const std::size_t nSize = nRows;

    SelfAdjointEigenSolver<MatrixXcd> es;
    es.compute(M);


    std::vector<std::complex<double>> Eigenvalues(nSize);
    std::vector<std::complex<double>> Derivatives(nSize);
    std::vector<double> AlreadyCalculated(nSize);//Array to check which EVs already been calculated.
    for(std::size_t i=0;i<nSize;i++) {
        Eigenvalues[i] = es.eigenvalues()[i];
        // std::cout << Eigenvalues[i] << sep << std::abs(Eigenvalues[i]) << std::endl;
        if(std::abs(Eigenvalues[i]) < EVThres) {
            // std::cout << "hallo" << std::endl;
            Eigenvalues[i] = 0;
        }
    }

    std::vector<std::vector<double>> Mapping(nSize,std::vector<double>(nSize));

    for(std::size_t i=0;i<nSize;i++)
    {
        AlreadyCalculated[i] = -1;
        for(std::size_t j=i;j<nSize;j++)
        {
            if(std::abs(Eigenvalues[i]-Eigenvalues[j]) > EVThres)
            {
                Mapping[i][j] = 0;
            }
            else{
                Mapping[i][j] = 1;
            }
        }
    }
    for(std::size_t i=1;i<nSize;i++)
    {
        for(std::size_t j=0;j<i;j++) Mapping[i][j] = Mapping[j][i];
    }

    for(std::size_t p=0;p<nSize;p++)
    {
        if(AlreadyCalculated[p] == -1)
        {
            std::size_t NumOfReps = 0;
            for(std::size_t i=p+1;i<nSize;i++)
            {
                NumOfReps += Mapping[p][i];
            }
            if(NumOfReps == 0)
            {
                VectorXcd v(nSize);
                v = es.eigenvectors().col(p);
                // Derivatives[p] = (v.transpose()*MDiff*v).value();
                Derivatives[p] = (v.adjoint()*MDiff*v).value();
                AlreadyCalculated[p] = 1;
            }
            else{
                MatrixXcd Phi(nSize,NumOfReps+1);
                std::size_t helpCol = 0;
                MatrixXcd MXWork(NumOfReps+1,NumOfReps+1);

                for(std::size_t i=p;i<nSize;i++)
                {
                    if(Mapping[p][i] == 1)
                    {
                        Phi.col(helpCol) = es.eigenvectors().col(i);
                        helpCol++;
                    }
                }
                MXWork = Phi.adjoint()*MDiff*Phi;
                SelfAdjointEigenSolver<MatrixXcd> esWork;
                esWork.compute(MXWork);
                helpCol = 0;
                for(std::size_t i=p;i<nSize;i++)
                {
                    if(Mapping[p][i] == 1)
                    {
                        AlreadyCalculated[i] = 1;
                        Derivatives[i] = esWork.eigenvalues()[helpCol];
                        helpCol++;
                    }
                }

            }
        }

    }
    for(std::size_t i=0;i<nSize;i++)
    {
        if(std::abs(Derivatives[i]) < EVThres) Derivatives[i] = 0;
    }

    for(std::size_t i=0;i<nSize;i++)
    {
        res.push_back(Eigenvalues[i].real());
    }
    for(std::size_t i=0;i<nSize;i++) res.push_back(Derivatives[i].real());
    return res;
}

double Class_Potential_Origin::fbasequad(double MassSquaredA, double MassSquaredB, double MassSquaredC, double MassSquaredD) const
{
	double res = 0;
	double mas = MassSquaredA;
	double mbs = MassSquaredB;
	double mcs = MassSquaredC;
	double mds = MassSquaredD;
	double LogA = 0, LogB = 0, LogC =  0, LogD = 0;

	double Thres = 1e-6;
	if(std::abs(mas) < Thres) mas = 0;
	if(std::abs(mbs) < Thres) mbs = 0;
	if(std::abs(mcs) < Thres) mcs = 0;
	if(std::abs(mds) < Thres) mds = 0;
	if(std::abs(mas - mbs) < Thres) mbs = mas;
	if(std::abs(mas - mcs) < Thres) mcs = mas;
	if(std::abs(mas - mds) < Thres) mds = mas;
	if(std::abs(mbs - mcs) < Thres) mbs = mcs;
	if(std::abs(mbs - mds) < Thres) mbs = mds;
	if(std::abs(mcs - mds) < Thres) mcs = mds;

	if(mas != 0) LogA = std::log(mas) - 2*std::log(scale);
	if(mbs != 0) LogB = std::log(mbs) - 2*std::log(scale);
	if(mcs != 0) LogC = std::log(mcs) - 2*std::log(scale);
	if(mds != 0) LogD = std::log(mds) - 2*std::log(scale);
/*
	double f0000 = 0.0, fa000 = LogA/(mas*mas), f0b00 = LogB/(mbs*mbs), f00c0 = LogC/(mcs*mcs), f000d = LogD/(mds*mds);
	double fab00 = (mbs*LogA - mas*LogB)/((mas-mbs)*mas*mbs);
	double faa00 = (1.0-LogA)/(mas*mas);
	double fa0c0 = (mcs*LogA - mas*LogC)/(mas*mcs*(mas - mcs));
	double fa0a0 = (1.0-LogA)/(mas*mas);
	double fa00d = (mds*LogA - mas*LogD)/(mas*mds*(mas-mds));
	double fa00a = (1.0 - LogA)/(mas*mas);
	double f0bc0 = (mcs*LogB - mbs*LogC)/(mbs*mcs*(mbs-mcs));
	double f0bb0 = (1.0 - LogB)/(mbs*mbs);
	double f0b0d = (mds*LogB - mbs*LogD)/(mbs*mds*(mbs-mds));
	double f0b0b = (1.0 - LogB)/(mbs*mbs);
	double f00cd = (mds*LogC - mcs*LogD)/(mcs*mds*(mcs-mds));
	double f00cc = (1.0 - LogC)/(mcs*mcs);
	double fabc0 = ((mbs-mcs)*LogA + (mcs-mas)*LogB + (mas-mbs)*LogC)/((mas - mbs)*(mas - mcs)*(mbs - mcs));
	double faac0 = (mas - mcs - mas*LogA + mas*LogC)/((mas-mcs)*(mas-mcs)*mas);
	double faba0 = (mas - mbs - mas*LogA + mas*LogB)/(mas*(mas-mbs)*(mas-mbs));
	double fabb0 = (-mas + mbs + mbs*LogA - mbs*LogB)/(mbs*(mas-mbs)*(mas-mbs));
	double faaa0 = -1.0/(2.0*mas*mas);
	double fab0d = ((mbs-mds)*LogA + (mds - mas)*LogB + (mas - mbs)*LogD)/((mas - mbs)*(mas-mds)*(mbs-mds));
	double faa0d = (mas -mds -mas*LogA + mas*LogD)/(mas*(mas-mds)*(mas-mds));
	double fab0a = (mas - mbs - mas*LogA + mas*LogB)/(mas*(mas-mbs)*(mas-mbs));
	double fab0b = (-mas + mbs + mbs*LogB - mbs*LogB)/(mbs*(mas-mbs)*(mas-mbs));
	double faa0a = -1.0/(2.0*mas*mas);

	double fa0cd = ((mcs - mds)*LogA + (mds - mas)*LogC + (mas - mcs)*LogD)/((mas-mcs)*(mas-mds)*(mcs-mds));
	double fa0ad = (mas - mds - mas*LogA + mas*LogD)/(mas*(mas - mds)*(mas-mds));
	double fa0ca = (mas -mcs - mas*LogA + mas*LogC)/(mas*(mas-mcs)*(mas-mcs));
	double fa0cc = -1.0/(2.0*mcs*mcs);
	double fa0aa = -1.0/(2.0*mas*mas);
	double f0bcd = ((mcs-mds)*LogB + (mds-mbs)*LogC + (mbs-mcs)*LogD)/((mbs-mcs)*(mbs-mds)*(mcs-mds));
	double f0bbd = (mbs - mds - mbs*LogB + mbs*LogD)/(mbs*(mbs-mds)*(mbs-mds));
	double f0bcb = (mbs - mcs - mbs*LogB + mbs*LogC)/(mbs*(mbs-mcs)*(mbs-mcs));
	double f0bcc = (-mbs + mcs + mcs*LogB - mcs*LogC)/((mbs-mcs)*(mbs-mcs)*mcs);
	double f0bbb = -1.0/(2.0*mbs*mbs);

	double faacd = ((mcs-mds)*(-mas*mas + mcs*mds)*LogA + mcs*(mas-mds)*(mas-mds)*LogC + (mas-mcs)*((mas-mds)*(mcs-mds) + (mcs-mas)*mds*LogD))/((mas-mcs)*(mas-mcs)*(mas-mds)*(mas-mds)*(mcs-mds));
	double fabad = ((mbs-mds)*(-mas*mas + mbs*mds)*LogA + mbs*(mas-mds)*(mas-mds)*LogB + (mas-mbs)*((mas-mds)*(mbs-mds) + (-mas + mds)*mds*LogD))/((mas-mbs)*(mas-mbs)*(mas-mds)*(mas-mds)*(mbs-mds));
	double fabca = ((mbs-mcs)*(-mas*mas + mbs*mcs)*LogA + mbs*(mas-mcs)*(mas-mcs)*LogB + (mas-mbs)*((mas-mcs)*(mbs-mcs) + (-mas + mbs)*mcs*LogC))/((mas-mbs)*(mas-mbs)*(mas-mcs)*(mas-mcs)*(mbs-mcs));
	double fabbd = mas*LogA/((mas-mbs)*(mas-mbs)*(mas-mds)) + (-(mas-mbs)*(mbs-mds) + (-mbs*mbs + mas*mds)*LogB)/((mas-mbs)*(mas-mbs)*(mbs-mds)*(mbs-mds)) + (mds*LogD)/((mbs-mds)*(mbs-mds)*(mds-mas));
	double fabcb = mas*LogA/((mas-mbs)*(mas-mbs)*(mas-mcs)) + (-(mas-mbs)*(mbs-mcs) + (-mbs*mbs + mas*mcs)*LogB)/((mas-mbs)*(mas-mbs)*(mbs-mcs)*(mbs-mcs)) + mcs*LogC/((mbs-mcs)*(mbs-mcs)*(-mas+mcs));
	double fabcc = (mas*(mbs-mcs)*(mbs-mcs)*LogA - mbs*(mas-mcs)*(mas-mcs)*LogB + (mas-mbs)*((mas-mcs)*(mbs-mcs) + (mas*mbs - mcs*mcs)*LogC))/((mas-mbs)*(mas-mcs)*(mas-mcs)*(mbs-mcs)*(mbs-mcs));

	double faaad = (-mas*mas + mds*mds + 2.0*mas*mds*LogA - 2.0*mas*mds*LogD)/(2.0*mas*(mas-mds)*(mas-mds)*(mas-mds));
	double faaca = (-mas*mas + mcs*mcs + 2.0*mas*mcs*LogA - 2.0*mas*mcs*LogC)/(2.0*mas*(mas-mcs)*(mas-mcs)*(mas-mcs));
	double fabaa = (-mas*mas + mbs*mbs + 2.0*mas*mbs*LogA - 2.0*mas*mbs*LogB)/(2.0*mas*(mas-mbs)*(mas-mbs)*(mas-mbs));
	double fabbb = (-mas*mas + mbs*mbs + 2.0*mas*mbs*LogA - 2.0*mas*mbs*LogB)/(2.0*mbs*(mas-mbs)*(mas-mbs)*(mas-mbs));

	double faacc = (2.0*(mas-mcs) - (mas+mcs)*LogA + (mas+mcs)*LogC)/((mas-mcs)*(mas-mcs)*(mas-mcs));
	double fabab = (2.0*(mas-mbs) - (mas+mbs)*LogA + (mas+mbs)*LogB)/((mas-mbs)*(mas-mbs)*(mas-mbs));
	double fabba = (2.0*(mas-mbs) - (mas+mbs)*LogA + (mas+mbs)*LogB)/((mas-mbs)*(mas-mbs)*(mas-mbs));

	double faaaa = -1.0/(6.0*mas*mas);
	double fabcd = mas*LogA/((mas-mbs)*(mas-mcs)*(mas-mds)) + mbs*LogB/((mbs-mas)*(mbs-mcs)*(mbs-mds)) + mcs*LogC/((mcs-mas)*(mcs-mbs)*(mcs-mds)) + mds*LogD/((mds-mas)*(mds-mbs)*(mds-mcs));
*/
	std::size_t C = 1;
	if(mas == 0 and mbs == 0 and mcs ==  0 and mds == 0){
		C =1;
		//res = f0000;
		res =  0.0;
	}
	else if(mas != 0 and mbs == 0 and mcs == 0 and mds == 0){
		C = 2;
		//res = fa000;
		//res = 0.0; 
		res = LogA/(mas*mas);
	}
	else if(mas == 0 and mbs != 0 and mcs == 0 and mds == 0){
		C = 3;
		//res = f0b00;
		//res = 0.0; 
		res =  LogB/(mbs*mbs);
	}
	else if(mas == 0 and mbs == 0 and mcs != 0 and mds == 0){
		C = 4;
		//res = f00c0;
		//res = 0.0; 
		res = LogC/(mcs*mcs);
	}
	else if(mas == 0 and mbs == 0 and mcs == 0 and mds != 0){
		C = 5;
		//res = f000d;
		//res = 0.0; 
		res = LogD/(mds*mds);
	}
	else if(mas != 0 and mbs != 0 and mcs == 0 and mds == 0 and mas != mbs){
		C =  6;
		//res = fab00;
		//res = (mbs*LogA - mas*LogB)/((mas-mbs)*mas*mbs);
		res = std::pow(mas*mbs,-1) + (mbs*LogA-mas*LogB)/(mas*mbs*(mas-mbs));
	}
	else if(mas != 0 and mbs == mas and mcs == 0 and mds == 0){
		C = 7;
		//res = faa00;
		res = (2.0*mas - mas*LogA)/(mas*mas*mas); 
		//res = (1.0-LogA)/(mas*mas);
	}
	else if(mas != 0 and mbs == 0 and mcs !=0 and mds == 0 and mas != mcs){
		C =  8;
		//res = fa0c0;
		//res = (mcs*LogA - mas*LogC)/(mas*mcs*(mas - mcs));
		res = std::pow(mas*mcs,-1) + (mcs*LogA-mas*LogC)/(mas*mcs*(mas-mcs));
	}
	else if(mas != 0 and mbs == 0 and mcs == mas and mds == 0){
		C = 9;
		//res = fa0a0;
		res = (2.0*mas - mas*LogA)/(mas*mas*mas); 
		//res = (1.0-LogA)/(mas*mas);
	}
	else if(mas != 0 and mbs == 0 and mcs == 0 and mds != 0 and mds != mas){
		C = 10;
		//res = fa00d;
		//res =(mds*LogA - mas*LogD)/(mas*mds*(mas-mds));
		res = std::pow(mas*mds,-1) + (mds*LogA-mas*LogD)/(mas*mds*(mas-mds));
	}
	else if(mas != 0 and mbs == 0 and mcs == 0 and mds == mas){
		C = 11;
		//res = fa00a;
		res = (2.0*mas - mas*LogA)/(mas*mas*mas); 
		//res = (1.0 - LogA)/(mas*mas);
	}
	else if(mas == 0 and mbs != 0 and mcs != 0 and mds == 0 and mbs != mcs){
		C = 12;
		//res = f0bc0;
		res = std::pow(mbs*mcs,-1) + (mcs*LogB-mbs*LogC)/(mbs*mcs*(mbs-mcs));
		//res =(mcs*LogB - mbs*LogC)/(mbs*mcs*(mbs-mcs));
	}
	else if(mas == 0 and mbs !=  0 and mcs == mbs and mds == 0){
		C = 13;
		//res = f0bb0;
		res =  (2.0*mbs + mbs*LogB)/(mbs*mbs*mbs); 
		//res = (1.0 - LogB)/(mbs*mbs);
	}
	else if(mas == 0 and mbs != 0 and mcs == 0 and mds != 0 and mbs != mds){
		C = 14;
		//res = f0b0d;
		res = std::pow(mbs*mds,-1) + (mds*LogB-mbs*LogD)/(mbs*mds*(mbs-mds));
		//res = (mds*LogB - mbs*LogD)/(mbs*mds*(mbs-mds));
	}
	else if(mas == 0 and mbs != 0 and mcs == 0 and mds != 0 and mds == mbs){
		C = 15;
		//res = f0b0b;
		res = (2.0*mbs + mbs*LogB)/(mbs*mbs*mbs); 
		//res = (1.0 - LogB)/(mbs*mbs);
	}
	else if(mas == 0 and mbs == 0 and mcs != 0 and mds != 0 and mcs != mds){
		C = 16;
		//res = f00cd;
		res = std::pow(mcs*mds,-1) + (mds*LogC-mcs*LogD)/(mcs*mds*(mcs-mds));
		//res =(mds*LogC - mcs*LogD)/(mcs*mds*(mcs-mds));
	}
	else if(mas == 0 and mbs == 0 and mcs != 0 and mds != 0 and mds == mcs){
		C = 17;
		//res = f00cc;
		res = (2.0*mcs + mcs*LogC)/(mcs*mcs*mcs); 
		//res = (1.0 - LogC)/(mcs*mcs);
	}
	else if(mas != 0 and mbs != 0 and mcs !=0 and mds == 0 and mas != mbs and mas != mcs and mbs!= mcs){
		C = 18;
		//res = fabc0;
		res =((mbs-mcs)*LogA + (mcs-mas)*LogB + (mas-mbs)*LogC)/((mas - mbs)*(mas - mcs)*(mbs - mcs));
	}
	else if(mas != 0 and mbs == mas and mcs != 0 and mds == 0 and mcs != mas){
		C = 19;
		//res = faac0;
		res =(mas - mcs - mas*LogA + mas*LogC)/((mas-mcs)*(mas-mcs)*mas);
	}
	else if(mas != 0 and mbs != 0 and mcs == mas and mds == 0 and mas != mbs){
		C = 20;
		//res = faba0;
		res =(mas - mbs - mas*LogA + mas*LogB)/(mas*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs != 0 and mcs == mbs and mds == 0 and mbs != mas){
		C = 21;
		//res = fabb0;
		res =(-mas + mbs + mbs*LogA - mbs*LogB)/(mbs*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs == mas and mcs == mas and mds == 0){
		C = 22;
		//res = faaa0;
		res = -1.0/(2.0*mas*mas);
	}
	else if(mas != 0 and mbs != 0 and mcs == 0 and mds != 0 and mas != mbs and mas != mds and mbs != mds){
		C = 23;
		//res = fab0d;
		res =((mbs-mds)*LogA + (mds - mas)*LogB + (mas - mbs)*LogD)/((mas - mbs)*(mas-mds)*(mbs-mds));
	}
	else if(mas != 0 and mbs == mas and mcs == 0 and mds != 0 and mds != mas){
		C = 24;
		//res = faa0d;
		res = (mas -mds -mas*LogA + mas*LogD)/(mas*(mas-mds)*(mas-mds));
	}
	else if(mas != 0 and mbs != 0 and mcs == 0 and mds == mas and mbs != mas){
		C = 25;
		//res = fab0a;
		res =(mas - mbs - mas*LogA + mas*LogB)/(mas*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs != 0 and mcs == 0 and mds == mbs and mbs != mas){
		C = 26;
		//res = fab0b;
		res =(-mas + mbs + mbs*LogB - mbs*LogB)/(mbs*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs == mas and mcs == 0 and mds == mas){
		C = 27;
		//res = faa0a;
		res =-1/(2*mas*mas);
	}
	else if(mas != 0 and mbs == 0 and mcs !=  0 and mds != 0 and mas!= mcs and mas!= mds and mcs != mds){
		C = 28;
		//res = fa0cd;
		res =((mcs - mds)*LogA + (mds - mas)*LogC + (mas - mcs)*LogD)/((mas-mcs)*(mas-mds)*(mcs-mds));
	}
	else if(mas != 0 and mbs == 0 and mcs == mas and mds != 0 and mds != mas){
		C = 29;
		//res = fa0ad;
		res =(mas - mds - mas*LogA + mas*LogD)/(mas*(mas - mds)*(mas-mds));
	}
	else if(mas != 0 and mbs == 0 and mcs != 0 and mds == mas and mcs != mas){
		C = 30;
		//res = fa0ca;
		res =(mas -mcs - mas*LogA + mas*LogC)/(mas*(mas-mcs)*(mas-mcs));
	}
	else if(mas != 0 and mbs == 0 and mcs != 0 and mds == mcs and mcs != mas){
		C = 31;
		//res = fa0cc;
		res = (mcs - mas + mcs * std::log(mas/mcs)) / (mcs * std::pow((mcs - mas),2));
		//res =-1.0/(2.0*mcs*mcs);
	}
	else if(mas != 0 and mbs == 0 and mcs == mas and mds == mas){
		C = 32;
		//res = fa0aa;
		res = -1/(2*mas*mas);
	}
	else if(mas == 0 and mbs != 0 and mcs !=0 and mds != 0 and mbs != mcs and mbs != mds and mcs!= mds){
		C = 33;
		//res = f0bcd;
		res = ((mcs - mds)*LogB + (mds - mbs)*LogC + (mbs - mcs)*LogD)/((mbs-mcs)*(mbs-mds)*(mcs-mds));
	}
	else if( mas == 0 and mbs != 0 and mcs == mbs and mds !=0 and mbs != mds){
		C = 34;
		//res = f0bbd;
		res =(mbs - mds - mbs*LogB + mbs*LogD)/(mbs*(mbs-mds)*(mbs-mds));
	}
	else if(mas == 0 and mbs != 0 and mcs != 0 and mds == mbs and mbs != mcs){
		C = 35;
		//res = f0bcb;
		res =(mbs - mcs - mbs*LogB + mbs*LogC)/(mbs*(mbs-mcs)*(mbs-mcs));
	}
	else if(mas == 0 and mbs != 0 and mcs != 0 and mds == mcs and mcs != mbs){
		C = 36;
		//res = f0bcc;
		res =(-mbs + mcs + mcs*LogB - mcs*LogC)/((mbs-mcs)*(mbs-mcs)*mcs);
	}
	else if(mas == 0 and mbs != 0 and mcs == mbs and mds == mbs){
		C = 37;
		//res = f0bbb;
		res =-1.0/(2.0*mbs*mbs);
	}
	else if(mas != 0 and mbs == mas and mcs != 0 and mds != 0 and mcs != mas and mas != mds and mcs != mds){
		C = 38;
		//res = faacd;
		res =((mcs-mds)*(-mas*mas + mcs*mds)*LogA + mcs*(mas-mds)*(mas-mds)*LogC + (mas-mcs)*((mas-mds)*(mcs-mds) + (mcs-mas)*mds*LogD))/((mas-mcs)*(mas-mcs)*(mas-mds)*(mas-mds)*(mcs-mds));
	}
	else if(mas != 0 and mbs != 0 and mcs == mas and mds != 0 and mas != mbs and mas != mds and mbs != mds){
		C = 39;
		//res = fabad;
		//res =((mbs-mds)*(-mas*mas + mbs*mds)*LogA + mbs*(mas-mds)*(mas-mds)*LogB + (mas-mbs)*((mas-mds)*(mbs-mds) + (-mas + mds)*mds*LogD))/((mas-mbs)*(mas-mbs)*(mas-mds)*(mas-mds)*(mbs-mds));
		res = (mbs*std::pow((mas-mds),2)*LogB-mds*std::pow((mas-mbs),2)*LogD+(mbs-mds)*((mas-mbs)*(mas-mds)+(-std::pow(mas,2)+mbs*mds)*LogA))/(std::pow((mas-mbs),2)*std::pow((mas-mds),2)*(mbs-mds));
	}
	else if(mas != 0 and mbs != 0 and mcs != 0 and mds == mas and mas != mbs and mas != mcs and mbs != mcs){
		C = 40;
		//res = fabca;
		res =((mbs-mcs)*(-mas*mas + mbs*mcs)*LogA + mbs*(mas-mcs)*(mas-mcs)*LogB + (mas-mbs)*((mas-mcs)*(mbs-mcs) + (-mas + mbs)*mcs*LogC))/((mas-mbs)*(mas-mbs)*(mas-mcs)*(mas-mcs)*(mbs-mcs));
	}
	else if(mas != 0 and mbs != 0 and mcs == mbs and mds != 0 and mas != mbs and mbs != mds and mas != mds){
		C = 41;
		//res = fabbd;
		res = mas*LogA/((mas-mbs)*(mas-mbs)*(mas-mds)) + (-(mas-mbs)*(mbs-mds) + (-mbs*mbs + mas*mds)*LogB)/((mas-mbs)*(mas-mbs)*(mbs-mds)*(mbs-mds)) + (mds*LogD)/((mbs-mds)*(mbs-mds)*(mds-mas));
	}
	else if(mas != 0 and mbs != 0 and mcs != 0 and mds == mbs and mas != mbs and mas != mcs and mbs != mcs){
		C = 42;
		//res = fabcb;
		res = mas*LogA/((mas-mbs)*(mas-mbs)*(mas-mcs)) + (-(mas-mbs)*(mbs-mcs) + (-mbs*mbs + mas*mcs)*LogB)/((mas-mbs)*(mas-mbs)*(mbs-mcs)*(mbs-mcs)) + mcs*LogC/((mbs-mcs)*(mbs-mcs)*(-mas+mcs));
	}
	else if(mas != 0 and mbs != 0 and mcs != 0 and mds == mcs and mas != mbs and mas != mcs and mbs != mcs){
		C = 43;
		//res = fabcc;
		res =(mas*(mbs-mcs)*(mbs-mcs)*LogA - mbs*(mas-mcs)*(mas-mcs)*LogB + (mas-mbs)*((mas-mcs)*(mbs-mcs) + (mas*mbs - mcs*mcs)*LogC))/((mas-mbs)*(mas-mcs)*(mas-mcs)*(mbs-mcs)*(mbs-mcs));
	}
	else if(mas != 0 and mbs == mas and mcs == mas and mds != 0 and mds != mas){
		C = 44;
		//res = faaad;
		res = (-mas*mas + mds*mds + 2.0*mas*mds*LogA - 2.0*mas*mds*LogD)/(2.0*mas*(mas-mds)*(mas-mds)*(mas-mds));
	}
	else if(mas != 0 and mbs == mas and mcs != 0 and mds == mas and mas != mcs){
		C = 45;
		//res = faaca;
		res = (-mas*mas + mcs*mcs + 2.0*mas*mcs*LogA - 2.0*mas*mcs*LogC)/(2.0*mas*(mas-mcs)*(mas-mcs)*(mas-mcs));
	}
	else if(mas != 0 and mbs != 0 and mcs == mas and mds == mas and mas != mbs){
		C = 46;
		//res = fabaa;
		res = (-mas*mas + mbs*mbs + 2.0*mas*mbs*LogA - 2.0*mas*mbs*LogB)/(2.0*mas*(mas-mbs)*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs != 0 and mcs == mbs and mds == mbs and mas != mbs){
		C = 47;
		//res = fabbb;
		res = (-mas*mas + mbs*mbs + 2.0*mas*mbs*LogA - 2.0*mas*mbs*LogB)/(2.0*mbs*(mas-mbs)*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs == mas and mcs != 0 and mds == mcs and mas != mcs){
		C = 48;
		//res = faacc;
		res = (2.0*(mas-mcs) - (mas+mcs)*LogA + (mas+mcs)*LogC)/((mas-mcs)*(mas-mcs)*(mas-mcs));
	}
	else if(mas != 0 and mbs != 0 and mcs == mas and mds == mbs and mas != mbs){
		C = 49;
		//res = fabab;
		res = (2.0*(mas-mbs) - (mas+mbs)*LogA + (mas+mbs)*LogB)/((mas-mbs)*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs != 0 and mcs == mbs and mds == mas and mas != mbs){
		C = 50;
		//res = fabba;
		res = (2.0*(mas-mbs) - (mas+mbs)*LogA + (mas+mbs)*LogB)/((mas-mbs)*(mas-mbs)*(mas-mbs));
	}
	else if(mas != 0 and mbs == mas and mcs == mas and mds == mas){
		C = 51;
		//res = faaaa;
		res = -1.0/(6.0*mas*mas);
	}
	else if(mas != 0 and mbs != 0 and mcs != 0 and mds != 0 and mas != mbs and mas != mcs and mas != mds and mbs != mcs and mbs != mds and mcs != mds)
	{
		C = 52;
		//res = fabcd;
		res = mas*LogA/((mas-mbs)*(mas-mcs)*(mas-mds)) + mbs*LogB/((mbs-mas)*(mbs-mcs)*(mbs-mds)) + mcs*LogC/((mcs-mas)*(mcs-mbs)*(mcs-mds)) + mds*LogD/((mds-mas)*(mds-mbs)*(mds-mcs));
	}
	if(std::isnan(res) or std::isinf(res)) {
		std::cout << "mas = " << mas << std::endl;
		std::cout << "mbs = " << mbs << std::endl;
		std::cout << "mcs = " << mcs << std::endl;
		std::cout << "mds = " << mds << std::endl;
		std::cout << "C = " << C << std::endl;
		std::cout << "res = " << res << std::endl;
		std::string throwstring = "Found nan at line = ";
		throwstring += std::to_string(InputLineNumber);
		throwstring += " in function ";
		throwstring+= __func__;
		throwstring+= "\n";
		std::cerr << "Found nan at line = " << InputLineNumber << " in function " << __func__ << std::endl;
		std::cerr << mas << sep << mbs << sep << mcs << sep << res << sep << C << std::endl;
		throw std::runtime_error(throwstring.c_str());
	}
	
	return res;
}

double Class_Potential_Origin::fbaseTri(double MassSquaredA, double MassSquaredB, double MassSquaredC) const
{
    double res = 0;
    double mas = MassSquaredA;
    double mbs = MassSquaredB;
    double mcs = MassSquaredC;
    double LogA = 0, LogB = 0, LogC = 0;
    double Thres = 1e-8;
    if(std::abs(mas) < Thres) mas = 0;
    if(std::abs(mbs) < Thres) mbs = 0;
    if(std::abs(mcs) < Thres) mcs = 0;
    if(std::abs(mas-mbs) < Thres) mas = mbs;
    if(std::abs(mas-mcs) < Thres) mas = mcs;
    if(std::abs(mbs-mcs) < Thres) mbs = mcs;

    if(mas != 0) LogA = std::log(mas) - 2*std::log(scale);
    if(mbs != 0) LogB = std::log(mbs) - 2*std::log(scale);
    if(mcs != 0) LogC = std::log(mcs) - 2*std::log(scale);

    std::size_t C = 1;

    // 1) all masses are zero
    if(mas == 0 and mbs == 0 and mcs == 0) res = 0;

    // 2) one mass is non-zero
    else if(mas != 0 and mbs == 0 and mcs == 0)
    {
        C = 2;
        res = (LogA-1.0)/ mas;
        //res = LogA/mas;     
    }
    else if(mas == 0 and mbs != 0 and mcs == 0)
    {
        C = 3;
        res = (LogB-1.0)/ mbs;
        //res = LogB/mbs;
    }
    else if(mas == 0 and mbs == 0 and mcs != 0)
    {
        C = 4;
        res = (LogC-1.0)/ mcs;
        //res = LogC/mcs;
    }

    // 3) two masses are non-zero
    // 3.1) non-zero masses are equal
    else if (mas == mbs and mas != 0 and mcs == 0)
    {
        C = 5;
        res = std::pow(mas,-1);
    }
    else if (mas == mcs and mas != 0 and mbs == 0)
    {
        C = 6;
        res = std::pow(mas,-1);
    }
    else if (mbs == mcs and mbs != 0 and mas == 0)
    {
        C = 7;
        res = std::pow(mbs,-1);
    }
    // 3.2) non-zero masses are not equal
    else if (mas != 0 and mbs != 0 and mas != mbs and mcs == 0)
    {
        C = 8;
        res = (std::log(mas) - std::log(mbs)) / (mas - mbs);    //DIFFERENT!!!! (not included in fbaseTri)
    }
    else if (mas != 0 and mcs != 0 and mas != mcs and mbs == 0)
    {
        C = 9;
        res = (std::log(mas) - std::log(mcs)) / (mas - mcs);    //DIFFERENT!!!!
    }
    else if (mas == 0 and mbs != 0 and mcs != 0 and mbs != mcs)
    {
        C = 10;
        res = (std::log(mbs) - std::log(mcs)) / (mbs - mcs);    //DIFFERENT!!!!
    }
    // 4) all three masses are non-zero
    // 4.1) all masses are equal
    else if (mas == mbs and mbs == mcs and mas != 0)
    {
        C = 11;
        res = 1.0 /(2*mas);
    }
    // 4.2) only two masses are equal
    else if (mas == mbs and mas != 0 and mas != mcs and mcs != 0)
    {
        C = 12;
        res = (mas - mcs + mcs * std::log(mcs/mas)) * std::pow((mas-mcs),-2);
    }
    else if (mas == mcs and mas != 0 and mas != mbs and mbs != 0)
    {
        C = 13;
        res = (mas - mbs + mbs * std::log(mbs/mas)) * std::pow((mas-mbs),-2);
    }
    else if (mbs == mcs and mbs != 0 and mbs != mas and mas != 0)
    {
        C = 14;
        res = (mbs - mas + mas * std::log(mas/mbs)) * std::pow((mas-mbs),-2);
    }
    // 4.3) all masses are not equal
    else
    {
        C = 15;
        res = mas*LogA/((mas-mbs)*(mas-mcs)) + mbs*LogB/((mbs-mas)*(mbs-mcs)) + mcs*LogC/((mcs-mas)*(mcs-mbs));
    }

    if(std::isnan(res) or std::isinf(res)) {
    	std::string throwstring = "Found nan at line = ";
    	throwstring += std::to_string(InputLineNumber);
    	throwstring += " in function ";
    	throwstring+= __func__;
    	throwstring+= "\n";
    	std::cerr << "Found nan at line = " << InputLineNumber << " in function " << __func__ << std::endl;
        std::cerr << mas << sep << mbs << sep << mcs << sep << res << sep << C << sep << std::endl;
    	throw std::runtime_error(throwstring.c_str());
    }

    return res;
}


/*
double Class_Potential_Origin::fbaseTri(double MassSquaredA, double MassSquaredB, double MassSquaredC) const
{
    double res = 0;
    double mas = MassSquaredA;
    double mbs = MassSquaredB;
    double mcs = MassSquaredC;
    double LogA = 0, LogB = 0, LogC = 0;
    double Thres = 1e-8;
    if(std::abs(mas) < Thres) mas = 0;
    if(std::abs(mbs) < Thres) mbs = 0;
    if(std::abs(mcs) < Thres) mcs = 0;
    if(std::abs(mas-mbs) < Thres) mas = mbs;
    if(std::abs(mas-mcs) < Thres) mas = mcs;
    if(std::abs(mbs-mcs) < Thres) mbs = mcs;

    if(mas != 0) LogA = std::log(mas) - 2*std::log(scale);
    if(mbs != 0) LogB = std::log(mbs) - 2*std::log(scale);
    if(mcs != 0) LogC = std::log(mcs) - 2*std::log(scale);

    std::size_t C = 1;
    if(mas == 0 and mbs == 0 and mcs == 0) res = 0;
    else if(mas != 0 and mbs == 0 and mcs == 0)
      {
	C=2;
	res = 1.0/mas*(std::log(mas) - 1);
	res = 1.0/mas*(LogA-1);
      }
    else if(mas == 0 and mbs != 0 and mcs == 0)
      {
	C=3;
	res = 1.0/mbs*(std::log(mbs)-1);
	res = (LogB-1)/mbs;
      }
    else if(mas == 0 and mbs == 0 and mcs != 0)
      {
	C=4;
	res = 1.0/mcs*(std::log(mcs)-1);
	res = (LogC-1)/mcs;
      }
    else if(mas == mbs and mas != 0 and mas != mcs and mcs != 0)
      {
	C=6;
	res = (mbs-mcs+mcs*std::log(mcs/mbs))/std::pow(mbs-mcs,2);
      }
    else if(mas == mcs and mas != 0 and mas != mbs and mbs != 0 )
      {
	C=7;
	res = (mbs*log(mbs/mcs) - mbs + mcs)/std::pow(mbs-mcs,2);
      }
    else if(mbs == mcs and mas != 0 and mbs != mas and mbs != 0 )
      {
	C=8;
	res = (mas*std::log(mas/mcs) - mas + mcs)/std::pow(mas-mcs,2);
      }
    else if(mas == mbs and mas == mcs and mas != 0)
      {
	C=9;
	res = 1.0/(2*mcs);
      }
    else if(mas == mbs and mas != mcs and mas != 0 and mcs == 0 )
      {
	C=10;
	res = 1.0/mbs;
      }
    else if(mas == mcs and mas != mbs and mas != 0 and mbs == 0)
      {
	C=11;
	res = 1.0/mas;
      }
    else if(mbs == mcs and mbs != 0 and mbs != mas and mas == 0 )
      {
	C=12;
	res = 1.0/mbs;
      }
    else
      {
    	C=5;
    	res = mas*LogA/((mas-mbs)*(mas-mcs)) + mbs*LogB/((mbs-mas)*(mbs-mcs));
    	res += mcs*LogC/((mcs-mas)*(mcs-mbs));
      }

    if(std::isnan(res) or std::isinf(res)) {
    	std::string throwstring = "Found nan at line = ";
    	throwstring += std::to_string(InputLineNumber);
    	throwstring += " in function ";
    	throwstring+= __func__;
    	throwstring+= "\n";
    	std::cerr << "Found nan at line = " << InputLineNumber << " in function " << __func__ << std::endl;
        std::cerr << mas << sep << mbs << sep << mcs << sep << res << sep << C << std::endl;
    	throw std::runtime_error(throwstring.c_str());
    }



    return res;
}
*/

double Class_Potential_Origin::fbase(double MassSquaredA,double MassSquaredB) const
{
    double res = 0;
    double LogA = 0;
    if(MassSquaredA == 0 and MassSquaredB == 0 ) return 1;
    double ZB = std::pow(10,-5);
    if(MassSquaredA != 0) LogA = std::log(MassSquaredA) - 2*std::log(scale);
    if(std::abs(MassSquaredA-MassSquaredB) > ZB)
    {
        double LogB = 0;
        if(MassSquaredB != 0)LogB = std::log(MassSquaredB) - 2*std::log(scale);
        if(MassSquaredA == 0) res = LogB;
        else if(MassSquaredB == 0 ) res = LogA;
        else res = (LogA*MassSquaredA - LogB*MassSquaredB)/(MassSquaredA-MassSquaredB);
    }
    else{
        res = 1 + LogA;
    }
    return res;
}


std::vector<double> Class_Potential_Origin::SecondDerivativeOfEigenvaluesNonRepeated(
        const Eigen::Ref<Eigen::MatrixXd> M,
        const Eigen::Ref<Eigen::MatrixXd> MDiffX,
        const Eigen::Ref<Eigen::MatrixXd> MDiffY,
        const Eigen::Ref<Eigen::MatrixXd> MDiffXY) const
{
    std::vector<double> res;
    const std::size_t nRows = M.rows();
    const std::size_t nCols = M.cols();

    const std::size_t EVThres = std::pow(10,-6);

    if(nCols != nRows) {
        throw std::runtime_error("ERROR ! M needs to be an quadratic Matrix for calculating the derivatives !\n");
    }

    const std::size_t nSize = nRows;

    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(M);

    std::vector<double> Eigenvalues(nSize);
    for(std::size_t i=0;i<nSize;i++) Eigenvalues[i] = es.eigenvalues()[i];
    for(std::size_t i=0;i<nSize-1;i++)
    {
        if(std::abs(Eigenvalues[i]-Eigenvalues[i+1]) < EVThres) {
            std::cerr << "ERROR ! repeated eigenvalues. \n";
        }
    }

    std::vector<std::vector<double>> Deriv(nSize,std::vector<double>(4));
    VectorXd v(nSize);
    MatrixXd C(nSize,nSize),E(nSize,nSize),Identity(nSize,nSize);
    Identity = MatrixXd::Identity(nSize,nSize);
    VectorXd vDiffX(nSize),vDiffY(nSize);

    for(std::size_t i=0;i<nSize;i++)
    {
        Deriv[i][0] = Eigenvalues[i];
        v = es.eigenvectors().col(i);
        Deriv[i][1] = v.transpose()*MDiffX*v;
        Deriv[i][2] = v.transpose()*MDiffY*v;

        C = (M-Deriv[i][0]*Identity).transpose()*(M-Deriv[i][0]*Identity) + v*v.transpose();
        E = (M-Deriv[i][0]*Identity).transpose()*(MDiffX-Deriv[i][1]*Identity);

        vDiffX = C.colPivHouseholderQr().solve(-E*v);

        E = (M-Deriv[i][0]*Identity).transpose()*(MDiffY-Deriv[i][2]*Identity);
        vDiffY = C.colPivHouseholderQr().solve(-E*v);

        Deriv[i][3] = v.transpose()*MDiffXY*v;
        Deriv[i][3] += v.transpose()*(MDiffX-Deriv[i][1]*Identity)*vDiffY;
        Deriv[i][3] += v.transpose()*(MDiffY-Deriv[i][2]*Identity)*vDiffX;

    }

    for(const auto& x : Eigenvalues) res.push_back(x);

    for(std::size_t i=0;i<nSize;i++)
    {
        for(std::size_t j=0;j<4;j++) res.push_back(Deriv[i][j]);
    }
    return res;
}


void Class_Potential_Origin::CalculatePhysicalCouplings()
{
    if(!SetCurvatureDone) SetCurvatureArrays();
    const double ZeroMass = std::pow(10,-5);
    MatrixXd MassHiggs(NHiggs,NHiggs),MassGauge(NGauge,NGauge);
    MatrixXcd MassQuark(NQuarks,NQuarks),MassLepton(NLepton,NLepton);
    MassHiggs = MatrixXd::Zero(NHiggs,NHiggs);
    MassGauge = MatrixXd::Zero(NGauge,NGauge);
    MassQuark = MatrixXcd::Zero(NQuarks,NQuarks);
    MassLepton = MatrixXcd::Zero(NLepton,NLepton);

    MassSquaredGauge.resize(NGauge);
    MassSquaredHiggs.resize(NHiggs);
    MassSquaredQuark.resize(NQuarks);
    MassSquaredLepton.resize(NLepton);
    HiggsRotationMatrix.resize(NHiggs);
    for(std::size_t i=0;i<NHiggs;i++) HiggsRotationMatrix[i].resize(NHiggs);


    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            MassHiggs(i,j) += Curvature_Higgs_L2[i][j];
            for(std::size_t k=0;k<NHiggs;k++)
            {
                MassHiggs(i,j) += Curvature_Higgs_L3[i][j][k]*HiggsVev[k];
                for(std::size_t l=0;l<NHiggs;l++) MassHiggs(i,j) += 0.5*Curvature_Higgs_L4[i][j][k][l]*HiggsVev[k]*HiggsVev[l];
            }
        }
    }

    for(std::size_t a=0;a<NGauge;a++)
    {
        for(std::size_t b=0;b<NGauge;b++)
        {
            for(std::size_t i=0;i<NHiggs;i++)
            {
                for(std::size_t j=0;j<NHiggs;j++)
                {
                    MassGauge(a,b) += 0.5*Curvature_Gauge_G2H2[a][b][i][j]*HiggsVev[i]*HiggsVev[j];
                }
            }
        }
    }

    MatrixXcd MIJQuarks(NQuarks,NQuarks);
    MIJQuarks = MatrixXcd::Zero(NQuarks,NQuarks);
    for(std::size_t a=0;a<NQuarks;a++)
    {
        for(std::size_t b=0;b<NQuarks;b++)
        {
            // MIJQuarks(a,b) = 0;
            for(std::size_t k=0;k<NHiggs;k++) MIJQuarks(a,b) += Curvature_Quark_F2H1[a][b][k]*HiggsVev[k];
        }
    }

    MassQuark = MIJQuarks.conjugate()*MIJQuarks;

    MatrixXcd MIJLeptons(NLepton,NLepton);
    MIJLeptons=MatrixXcd::Zero(NLepton,NLepton);
    for(std::size_t a=0;a<NLepton;a++)
    {
        for(std::size_t b=0;b<NLepton;b++)
        {
            // MIJLeptons(a,b) = 0;
            for(std::size_t k=0;k<NHiggs;k++) MIJLeptons(a,b) += Curvature_Lepton_F2H1[a][b][k]*HiggsVev[k];
        }
    }

    MassLepton = MIJLeptons.conjugate()*MIJLeptons;

    MatrixXd HiggsRot(NHiggs,NHiggs),GaugeRot(NGauge,NGauge),QuarkRot(NQuarks,NQuarks),LepRot(NLepton,NLepton);
    HiggsRot = MatrixXd::Identity(NHiggs,NHiggs);
    GaugeRot = MatrixXd::Identity(NGauge,NGauge);
    QuarkRot = MatrixXd::Identity(NQuarks,NQuarks);
    LepRot = MatrixXd::Identity(NLepton,NLepton);

    SelfAdjointEigenSolver<MatrixXd> es;

	es.compute(MassHiggs);
	HiggsRot = es.eigenvectors().transpose();
    for(std::size_t i=0;i<NHiggs;i++)
	{
        for(std::size_t j=0;j<NHiggs;j++)
		{
			if(std::abs(HiggsRot(i,j)) < std::pow(10,-10)) HiggsRot(i,j) = 0;
		}
	}

    for(std::size_t i=0;i<NHiggs;i++)
    {
        MassSquaredHiggs[i] = es.eigenvalues()[i];
        if(std::abs(MassSquaredHiggs[i]) < ZeroMass) MassSquaredHiggs[i] = 0;
    }


	es.compute(MassGauge);
	GaugeRot = es.eigenvectors().transpose();

    for(std::size_t i=0;i<NGauge;i++)
    {
        MassSquaredGauge[i] = es.eigenvalues()[i];
        if(std::abs(MassSquaredGauge[i]) < ZeroMass) MassSquaredGauge[i] = 0;
    }


	SelfAdjointEigenSolver<MatrixXcd> esQuark(MassQuark);
	QuarkRot = esQuark.eigenvectors().transpose().real();
    for(std::size_t i=0;i<NQuarks;i++) MassSquaredQuark[i] = esQuark.eigenvalues().real()[i];

	SelfAdjointEigenSolver<MatrixXcd> esLepton(MassLepton);
	LepRot = esLepton.eigenvectors().transpose().real();
    for(std::size_t i=0;i<NLepton;i++) MassSquaredLepton[i] = esLepton.eigenvalues().real()[i];

    for(std::size_t a=0;a<NGauge;a++)
    {
        for(std::size_t b=0;b<NGauge;b++)
        {
            for(std::size_t i=0;i<NHiggs;i++)
            {
                LambdaGauge_3[a][b][i] = 0;
                for(std::size_t j=0;j<NHiggs;j++) LambdaGauge_3[a][b][i] += Curvature_Gauge_G2H2[a][b][i][j]*HiggsVev[j];
            }
        }
    }
    for(std::size_t a=0;a<NHiggs;a++)
    {
        for(std::size_t b=0;b<NHiggs;b++)
        {
            for(std::size_t i=0;i<NHiggs;i++)
            {
                LambdaHiggs_3[a][b][i] = Curvature_Higgs_L3[a][b][i];

                for(std::size_t j=0;j<NHiggs;j++) {
                	LambdaHiggs_3[a][b][i] += Curvature_Higgs_L4[a][b][i][j]*HiggsVev[j];

                }
            }
        }
    }

    for(std::size_t i=0;i<NQuarks;i++)
    {
        for(std::size_t j=0;j<NQuarks;j++)
        {
            for(std::size_t k=0;k<NHiggs;k++)
            {
                LambdaQuark_3[i][j][k] = 0;
                for(std::size_t l=0;l<NQuarks;l++){
                    LambdaQuark_3[i][j][k] += conj(Curvature_Quark_F2H1[i][l][k])*MIJQuarks(l,j);
                    LambdaQuark_3[i][j][k] += conj(MIJQuarks(i,l)) * Curvature_Quark_F2H1[l][j][k];
                }
                for(std::size_t m=0;m<NHiggs;m++)
                {
                    LambdaQuark_4[i][j][k][m] = 0;
                    for(std::size_t l=0;l<NQuarks;l++)
                    {
                        LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][k])*Curvature_Quark_F2H1[l][j][m];
                        LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][m])*Curvature_Quark_F2H1[l][j][k];
                    }
                }
            }
        }
    }


    for(std::size_t i=0;i<NLepton;i++)
    {
        for(std::size_t j=0;j<NLepton;j++)
        {
            for(std::size_t k=0;k<NHiggs;k++)
            {
                LambdaLepton_3[i][j][k] = 0;
                for(std::size_t l=0;l<NLepton;l++){
                    LambdaLepton_3[i][j][k] += conj(Curvature_Lepton_F2H1[i][l][k])*MIJLeptons(l,j);
                    LambdaLepton_3[i][j][k] += conj(MIJLeptons(i,l)) * Curvature_Lepton_F2H1[l][j][k];
                }
                for(std::size_t m=0;m<NHiggs;m++)
                {
                    LambdaLepton_4[i][j][k][m] = 0;
                    for(std::size_t l=0;l<NLepton;l++)
                    {
                        LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][k])*Curvature_Lepton_F2H1[l][j][m];
                        LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][m])*Curvature_Lepton_F2H1[l][j][k];
                    }
                }
            }
        }
    }

    // Rotate and save std::size_to corresponding vectors

    Couplings_Higgs_Quartic.resize(NHiggs);
    Couplings_Higgs_Triple.resize(NHiggs);
    for(std::size_t i=0;i<NHiggs;i++) {
        Couplings_Higgs_Quartic[i].resize(NHiggs);
        Couplings_Higgs_Triple[i].resize(NHiggs);
        for(std::size_t j=0;j<NHiggs;j++){
            Couplings_Higgs_Quartic[i][j].resize(NHiggs);
            Couplings_Higgs_Triple[i][j].resize(NHiggs);
            for(std::size_t k=0;k<NHiggs;k++) Couplings_Higgs_Quartic[i][j][k].resize(NHiggs);
        }
    }



    Couplings_Gauge_Higgs_22.resize(NGauge);
    Couplings_Gauge_Higgs_21.resize(NGauge);
    for(std::size_t a=0;a<NGauge;a++)
    {
        Couplings_Gauge_Higgs_22[a].resize(NGauge);
        Couplings_Gauge_Higgs_21[a].resize(NGauge);
        for(std::size_t b=0;b<NGauge;b++)
        {
            Couplings_Gauge_Higgs_22[a][b].resize(NHiggs);
            Couplings_Gauge_Higgs_21[a][b].resize(NHiggs);
            for(std::size_t i=0;i<NHiggs;i++)
            {
                Couplings_Gauge_Higgs_22[a][b][i].resize(NHiggs);
            }
        }
    }



    Couplings_Quark_Higgs_22.resize(NQuarks);
    Couplings_Quark_Higgs_21.resize(NQuarks);
    for(std::size_t a=0;a<NQuarks;a++)
    {
        Couplings_Quark_Higgs_22[a].resize(NQuarks);
        Couplings_Quark_Higgs_21[a].resize(NQuarks);
        for(std::size_t b=0;b<NQuarks;b++) {
            Couplings_Quark_Higgs_22[a][b].resize(NHiggs);
            Couplings_Quark_Higgs_21[a][b].resize(NHiggs);
            for(std::size_t i=0;i<NHiggs;i++) Couplings_Quark_Higgs_22[a][b][i].resize(NHiggs);
        }
    }




    Couplings_Lepton_Higgs_22.resize(NLepton);
    Couplings_Lepton_Higgs_21.resize(NLepton);
    for(std::size_t a=0;a<NLepton;a++)
    {
        Couplings_Lepton_Higgs_22[a].resize(NLepton);
        Couplings_Lepton_Higgs_21[a].resize(NLepton);
        for(std::size_t b=0;b<NLepton;b++)
        {
            Couplings_Lepton_Higgs_22[a][b].resize(NHiggs);
            Couplings_Lepton_Higgs_21[a][b].resize(NHiggs);
            for(std::size_t i=0;i<NHiggs;i++) Couplings_Lepton_Higgs_22[a][b][i].resize(NHiggs);
        }
    }


  for(std::size_t i=0;i<NHiggs;i++)
	{
        for(std::size_t j=0;j<NHiggs;j++)
		{
            for(std::size_t k=0;k<NHiggs;k++)
			{
				Couplings_Higgs_Triple[i][j][k] = 0;
                for(std::size_t is=0;is<NHiggs;is++)
				{
                    for(std::size_t js=0;js<NHiggs;js++)
					{
                        for(std::size_t ks=0;ks<NHiggs;ks++)
						{
							Couplings_Higgs_Triple[i][j][k] += HiggsRot(i,is)*HiggsRot(j,js)*HiggsRot(k,ks)*LambdaHiggs_3[is][js][ks];
						}
					}
				}
                for(std::size_t l=0;l<NHiggs;l++)
				{
					Couplings_Higgs_Quartic[i][j][k][l] = 0;
                    for(std::size_t is=0;is<NHiggs;is++)
					{
                        for(std::size_t js=0;js<NHiggs;js++)
						{
                            for(std::size_t ks=0;ks<NHiggs;ks++)
							{
                                for(std::size_t ls=0;ls<NHiggs;ls++)
								{
									Couplings_Higgs_Quartic[i][j][k][l] += HiggsRot(i,is)*HiggsRot(j,js)*HiggsRot(k,ks)*HiggsRot(l,ls)*Curvature_Higgs_L4[is][js][ks][ls];
								}
							}
						}
					}
				}
			}
		}
	}




    // Gauge Rot
    for(std::size_t a=0;a<NGauge;a++)
	{
        for(std::size_t b=0;b<NGauge;b++)
		{
            for(std::size_t i=0;i<NHiggs;i++)
			{
				Couplings_Gauge_Higgs_21[a][b][i] = 0;
                for(std::size_t as=0;as<NGauge;as++)
				{
                    for(std::size_t bs=0;bs<NGauge;bs++)
					{
                        for(std::size_t is=0;is<NHiggs;is++)
						Couplings_Gauge_Higgs_21[a][b][i] += GaugeRot(a,as)*GaugeRot(b,bs)*HiggsRot(i,is)*LambdaGauge_3[as][bs][is];
					}
				}
                for(std::size_t j=0;j<NHiggs;j++)
				{
					Couplings_Gauge_Higgs_22[a][b][i][j] = 0;
                    for(std::size_t as=0;as<NGauge;as++)
					{
                        for(std::size_t bs=0;bs<NGauge;bs++)
						{
                            for(std::size_t is=0;is<NHiggs;is++)
							{
                                for(std::size_t js=0;js<NHiggs;js++){
									double RotFac = GaugeRot(a,as)*GaugeRot(b,bs)*HiggsRot(i,is)*HiggsRot(j,js);
									Couplings_Gauge_Higgs_22[a][b][i][j] += RotFac*Curvature_Gauge_G2H2[as][bs][is][js];
								}
							}

						}
					}
				}
			}
		}
	}


    // Quark

    for(std::size_t a=0;a<NQuarks;a++)
	{
        for(std::size_t b=0;b<NQuarks;b++)
		{
            for(std::size_t i=0;i<NHiggs;i++)
			{
				Couplings_Quark_Higgs_21[a][b][i] = 0;
                for(std::size_t as = 0;as<NQuarks;as++)
				{
                    for(std::size_t bs=0;bs<NQuarks;bs++)
					{
                        for(std::size_t is=0;is<NHiggs;is++)
						{
							double RotFac = QuarkRot(a,as)*QuarkRot(b,bs)*HiggsRot(i,is);
							Couplings_Quark_Higgs_21[a][b][i] += RotFac*LambdaQuark_3[as][bs][is];
						}
					}
				}
                for(std::size_t j=0;j<NHiggs;j++)
				{
					Couplings_Quark_Higgs_22[a][b][i][j] = 0;
                    for(std::size_t as=0;as<NQuarks;as++)
					{
                        for(std::size_t bs=0;bs<NQuarks;bs++)
						{
                            for(std::size_t is=0;is<NHiggs;is++)
							{
                                for(std::size_t js=0;js<NHiggs;js++)
								{
									double RotFac = QuarkRot(a,as)*QuarkRot(b,bs)*HiggsRot(i,is)*HiggsRot(j,js);
									Couplings_Quark_Higgs_22[a][b][i][j] += RotFac*LambdaQuark_4[as][bs][is][js];
								}
							}
						}
					}
				}
			}
		}
	}


    // Lepton


    for(std::size_t a=0;a<NLepton;a++)
	{
        for(std::size_t b=0;b<NLepton;b++)
		{
            for(std::size_t i=0;i<NHiggs;i++)
			{
				Couplings_Lepton_Higgs_21[a][b][i] = 0;
                for(std::size_t as = 0;as<NLepton;as++)
				{
                    for(std::size_t bs=0;bs<NLepton;bs++)
					{
                        for(std::size_t is=0;is<NHiggs;is++)
						{
							double RotFac = LepRot(a,as)*LepRot(b,bs)*HiggsRot(i,is);
							Couplings_Lepton_Higgs_21[a][b][i] += RotFac*LambdaLepton_3[as][bs][is];
						}
					}
				}
                for(std::size_t j=0;j<NHiggs;j++)
				{
					Couplings_Lepton_Higgs_22[a][b][i][j] = 0;
                    for(std::size_t as=0;as<NLepton;as++)
					{
                        for(std::size_t bs=0;bs<NLepton;bs++)
						{
                            for(std::size_t is=0;is<NHiggs;is++)
							{
                                for(std::size_t js=0;js<NHiggs;js++)
								{
									double RotFac = LepRot(a,as)*LepRot(b,bs)*HiggsRot(i,is)*HiggsRot(j,js);
									Couplings_Lepton_Higgs_22[a][b][i][j] += RotFac*LambdaLepton_4[as][bs][is][js];
								}
							}
						}
					}
				}
			}
		}
	}

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            HiggsRotationMatrix[i][j] = HiggsRot(i,j);
        }
    }


    CalcCouplingsdone = true;

     return;


}


std::vector<double> Class_Potential_Origin::WeinbergFirstDerivative() const{
    std::vector<double> res;
    if(!CalcCouplingsdone) {
//        CalculatePhysicalCouplings();
        std::string retmes = __func__;
        retmes += " tries to use Physical couplings but they are not initialised.";
        throw std::runtime_error(retmes);
    }
    const double NumZero = std::pow(10,-10);
    VectorXd FirstDeriv(NHiggs),FirstDerivGauge(NHiggs),FirstDerivHiggs(NHiggs),FirstDerivQuark(NHiggs),FirstDerivLepton(NHiggs);
    FirstDeriv = VectorXd::Zero(NHiggs);
    FirstDerivGauge = VectorXd::Zero(NHiggs);
    FirstDerivHiggs = VectorXd::Zero(NHiggs);
    FirstDerivQuark = VectorXd::Zero(NHiggs);
    FirstDerivLepton = VectorXd::Zero(NHiggs);
    double epsilon = 1.0/(16*M_PI*M_PI);

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t a=0;a<NGauge;a++)
        {
            if(MassSquaredGauge[a] != 0)
            {
                FirstDerivGauge[i] += MassSquaredGauge[a]*Couplings_Gauge_Higgs_21[a][a][i]*(std::log(MassSquaredGauge[a]/std::pow(scale,2))-C_CWcbGB + 0.5 );
            }
        }

        for(std::size_t a = 0;a<NHiggs;a++)
        {
            if(MassSquaredHiggs[a] != 0)
            {

                FirstDerivHiggs[i] += MassSquaredHiggs[a]*Couplings_Higgs_Triple[a][a][i]*(std::log(MassSquaredHiggs[a]/std::pow(scale,2)) - C_CWcbHiggs + 0.5);
            }
        }
        for(std::size_t a=0;a<NQuarks;a++)
        {
            if(MassSquaredQuark[a] != 0)
            {
                double Coup = Couplings_Quark_Higgs_21[a][a][i].real();
                FirstDerivQuark[i] += MassSquaredQuark[a]*Coup*(std::log(MassSquaredQuark[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);
            }
        }
        for(std::size_t a=0;a<NLepton;a++)
        {
            if(MassSquaredLepton[a] != 0)
            {
                double Coup = Couplings_Lepton_Higgs_21[a][a][i].real();
                FirstDerivLepton[i] += MassSquaredLepton[a]*Coup*(std::log(MassSquaredLepton[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);
            }
        }
    }
    FirstDerivGauge *= 1.5;
    FirstDerivHiggs *= 0.5;
    FirstDerivQuark *= -3;
    FirstDerivLepton *= -1;

    MatrixXd HiggsRot(NHiggs,NHiggs);
    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++) HiggsRot(i,j) = HiggsRotationMatrix[i][j];
    }

    FirstDeriv = HiggsRot.transpose()*(FirstDerivGauge+FirstDerivHiggs+FirstDerivQuark+FirstDerivLepton);
    FirstDeriv *= epsilon;

    for(std::size_t i=0;i<NHiggs;i++)
    {
        if(std::abs(FirstDeriv[i]) < NumZero) FirstDeriv[i] = 0;
    }

    for(std::size_t i=0;i<NHiggs;i++)    res.push_back(FirstDeriv[i]);

    return  res;
}


std::vector<double> Class_Potential_Origin::WeinbergSecondDerivative() const
{
    if(!CalcCouplingsdone) {
//        CalculatePhysicalCouplings();
        std::string retmes = __func__;
        retmes += " tries to use Physical couplings but they are not initialised.";
        throw std::runtime_error(retmes);
    }
    std::vector<double> res;
    const double NumZero = std::pow(10,-10);
    MatrixXd GaugePart(NHiggs,NHiggs),HiggsPart(NHiggs,NHiggs),QuarkPart(NHiggs,NHiggs),LeptonPart(NHiggs,NHiggs);
    GaugePart = MatrixXd::Zero(NHiggs,NHiggs);
    HiggsPart = MatrixXd::Zero(NHiggs,NHiggs);
    QuarkPart = MatrixXd::Zero(NHiggs,NHiggs);
    LeptonPart = MatrixXd::Zero(NHiggs,NHiggs);

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            for(std::size_t a=0;a<NGauge;a++)
            {
                for(std::size_t b=0;b<NGauge;b++)
                {
                    double Coup1 = Couplings_Gauge_Higgs_21[a][b][i];
                    double Coup2 = Couplings_Gauge_Higgs_21[b][a][j];
                    double Br = fbase(MassSquaredGauge[a],MassSquaredGauge[b]) - C_CWcbGB + 0.5;
                    GaugePart(i,j) += Coup1*Coup2*Br;
                }
                if(MassSquaredGauge[a] != 0)
                {
                    GaugePart(i,j) += MassSquaredGauge[a]*Couplings_Gauge_Higgs_22[a][a][i][j]*(std::log(MassSquaredGauge[a]/std::pow(scale,2)) - C_CWcbGB + 0.5);
                }
            }


            for(std::size_t a=0;a<NHiggs;a++)
            {
                for(std::size_t b=0;b<NHiggs;b++)
                {
                    double Coup1 = Couplings_Higgs_Triple[a][b][i];
                    double Coup2 = Couplings_Higgs_Triple[b][a][j];
                    double Br = fbase(MassSquaredHiggs[a],MassSquaredHiggs[b]) - C_CWcbHiggs + 0.5;
                    HiggsPart(i,j) += Coup1*Coup2*Br;
                }
                if(MassSquaredHiggs[a] != 0)
                {
                    HiggsPart(i,j) += MassSquaredHiggs[a]*Couplings_Higgs_Quartic[a][a][i][j]*(std::log(MassSquaredHiggs[a]/std::pow(scale,2)) - C_CWcbHiggs + 0.5);
                }
            }

            for(std::size_t a=0;a<NQuarks;a++)
            {
                for(std::size_t b=0;b<NQuarks;b++)
                {
                    double Coup = (Couplings_Quark_Higgs_21[a][b][i] * Couplings_Quark_Higgs_21[b][a][j]).real();
                    double Br = fbase(MassSquaredQuark[a],MassSquaredQuark[b]) - C_CWcbFermion + 0.5;
                    QuarkPart(i,j) += Coup*Br;
                }
                if(MassSquaredQuark[a] != 0)
                {
                    double Coup = Couplings_Quark_Higgs_22[a][a][i][j].real();
                    QuarkPart(i,j) += MassSquaredQuark[a]*Coup*(std::log(MassSquaredQuark[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);

                }
            }

            for(std::size_t a=0;a<NLepton;a++)
            {
                for(std::size_t b=0;b<NLepton;b++)
                {
                    double Coup = (Couplings_Lepton_Higgs_21[a][b][i]*Couplings_Lepton_Higgs_21[b][a][j]).real();
                    double Br = fbase(MassSquaredLepton[a],MassSquaredLepton[b]) - C_CWcbFermion + 0.5;
                    LeptonPart(i,j) += Coup*Br;
                }
                if(MassSquaredLepton[a] != 0)
                {
                    double Coup = Couplings_Lepton_Higgs_22[a][a][i][j].real();
                    LeptonPart(i,j) += Coup*MassSquaredLepton[a]*(std::log(MassSquaredLepton[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);
                }
            }
        }
    }

    HiggsPart *= 0.5;
    GaugePart *= 1.5;
    QuarkPart*= -3;
    LeptonPart *= -1;

    MatrixXd Storage(NHiggs,NHiggs);
    Storage = HiggsPart+GaugePart+QuarkPart+LeptonPart;

    for(std::size_t i=0;i<NHiggs;i++){
      for(std::size_t j=0;j<NHiggs;j++){
        if(std::abs(Storage(i,j)) < NumZero) Storage(i,j) = 0;
      }
    }
    MatrixXd ResMatrix;
    MatrixXd HiggsRot(NHiggs,NHiggs);
    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            HiggsRot(i,j) = HiggsRotationMatrix[i][j];
        }
    }

    ResMatrix = 0.5*HiggsRot.transpose()*(Storage+Storage.transpose())*HiggsRot;
    double epsilon = 1.0/(16.0*M_PI*M_PI);
    ResMatrix *= epsilon;

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            if(std::abs(ResMatrix(i,j)) < NumZero) ResMatrix(i,j) = 0;
        }
    }

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            res.push_back(ResMatrix(j,i));
        }
    }

    return  res;
}


std::vector<double> Class_Potential_Origin::WeinbergThirdDerivative() const{

    if(not CalcCouplingsdone){
        std::string retmes = __func__;
        retmes += " tries to use Physical couplings but they are not initialised.";
        throw std::runtime_error(retmes);
    }
    const double NumZero = std::pow(10,-10);
    double epsilon = 1.0/(16.0*M_PI*M_PI);

    std::vector<double> res;

    std::vector<std::vector<std::vector< std::complex<double>>>> restmp;
    std::vector<std::vector<std::vector< std::complex<double>>>> QuarkPart;
    std::vector<std::vector<std::vector< std::complex<double>>>> LeptonPart;
    std::vector<std::vector<std::vector< std::complex<double>>>> QuarkPartSym;
    std::vector<std::vector<std::vector< std::complex<double>>>> LeptonPartSym;
    restmp.resize(NHiggs);
    QuarkPart.resize(NHiggs);
    LeptonPart.resize(NHiggs);
    QuarkPartSym.resize(NHiggs);
    LeptonPartSym.resize(NHiggs);
    for(std::size_t i=0;i<NHiggs;i++)
    {
        restmp[i].resize(NHiggs);
        QuarkPart[i].resize(NHiggs);
        LeptonPart[i].resize(NHiggs);
        QuarkPartSym[i].resize(NHiggs);
        LeptonPartSym[i].resize(NHiggs);
        for(std::size_t j=0;j<NHiggs;j++) {
            restmp[i][j].resize(NHiggs);
            QuarkPart[i][j].resize(NHiggs );
            LeptonPart[i][j].resize(NHiggs);
            QuarkPartSym[i][j].resize(NHiggs);
            LeptonPartSym[i][j].resize(NHiggs);
        }
    }

    std::vector<std::vector<std::vector<double>>> resGaugeBase(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                        std::vector<double>(NHiggs)));
    std::vector<std::vector<std::vector<double>>> Higgspart(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                     std::vector<double>(NHiggs)));


    std::vector<std::vector<std::vector<double>>> GaugePart(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                     std::vector<double>(NHiggs)));

    std::vector<std::vector<std::vector<double>>> HiggspartSym(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                        std::vector<double>(NHiggs)));

    std::vector<std::vector<std::vector<double>>> GaugePartSym(NHiggs, std::vector<std::vector<double>>(NHiggs,
                                                                                                        std::vector<double>(NHiggs)));


    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            for(std::size_t k=0;k<NHiggs;k++)
            {
                Higgspart[i][j][k] = 0;
                for(std::size_t a=0;a<NHiggs;a++)
                {
                    for(std::size_t b=0;b<NHiggs;b++)
                    {
                        for(std::size_t c=0;c<NHiggs;c++)
                        {
                            double f1 = fbaseTri(MassSquaredHiggs[a],MassSquaredHiggs[b],MassSquaredHiggs[c]);
                            double f2 = Couplings_Higgs_Triple[a][b][i];
                            double f3 = Couplings_Higgs_Triple[b][c][j];
                            double f4 = Couplings_Higgs_Triple[c][a][k];
                            Higgspart[i][j][k] += 2*f1*f2*f3*f4;
                        }
                        double f1 = Couplings_Higgs_Quartic[a][b][i][j];
                        double f2 = Couplings_Higgs_Triple[b][a][k];
                        double f3 = fbase(MassSquaredHiggs[a],MassSquaredHiggs[b]) - C_CWcbHiggs + 0.5;
                        Higgspart[i][j][k] += 3*f1*f2*f3;
                    }
                }

                GaugePart[i][j][k] = 0;
                for(std::size_t a=0;a<NGauge;a++)
                {
                    for(std::size_t b=0;b<NGauge;b++)
                    {
                        for(std::size_t c=0;c<NGauge;c++)
                        {
                            double f1 = fbaseTri(MassSquaredGauge[a],MassSquaredGauge[b],MassSquaredGauge[c]);
                            double f2 = Couplings_Gauge_Higgs_21[a][b][i];
                            double f3 = Couplings_Gauge_Higgs_21[b][c][j];
                            double f4 = Couplings_Gauge_Higgs_21[c][a][k];
                            GaugePart[i][j][k] += 2*f1*f2*f3*f4;
                        }
                        double f1 = Couplings_Gauge_Higgs_22[a][b][i][j];
                        double f2 = Couplings_Gauge_Higgs_21[b][a][k];
                        double f3 = fbase(MassSquaredGauge[a],MassSquaredGauge[b]) - C_CWcbGB + 0.5;
                        GaugePart[i][j][k] += 3*f1*f2*f3;
                    }
                }

                QuarkPart[i][j][k] = 0;
                for(std::size_t a=0;a<NQuarks;a++)
                {
                    for(std::size_t b=0;b<NQuarks;b++)
                    {
                        for(std::size_t c=0;c<NQuarks;c++)
                        {
                            std::complex<double> f1 = fbaseTri(MassSquaredQuark[a],MassSquaredQuark[b],MassSquaredQuark[c]);
                            std::complex<double> f2 = Couplings_Quark_Higgs_21[a][b][i];
                            std::complex<double> f3 = Couplings_Quark_Higgs_21[b][c][j];
                            std::complex<double> f4 = Couplings_Quark_Higgs_21[c][a][k];
                            QuarkPart[i][j][k] += 2.0*f1*f2*f3*f4;
                        }
                        std::complex<double> f1 = Couplings_Quark_Higgs_22[a][b][i][j];
                        std::complex<double> f2 = Couplings_Quark_Higgs_21[b][a][k];
                        std::complex<double> f3 = fbase(MassSquaredQuark[a],MassSquaredQuark[b]) - C_CWcbFermion + 0.5;
                        QuarkPart[i][j][k] += 3.0*f1*f2*f3;
                    }
                }
                LeptonPart[i][j][k] = 0;
                for(std::size_t a=0;a<NLepton;a++)
                {
                    for(std::size_t b=0;b<NLepton;b++)
                    {
                        for(std::size_t c=0;c<NLepton;c++)
                        {
                            std::complex<double> f1 = fbaseTri(MassSquaredLepton[a],MassSquaredLepton[b],MassSquaredLepton[c]);
                            std::complex<double> f2 = Couplings_Lepton_Higgs_21[a][b][i];
                            std::complex<double> f3 = Couplings_Lepton_Higgs_21[b][c][j];
                            std::complex<double> f4 = Couplings_Lepton_Higgs_21[c][a][k];
                            LeptonPart[i][j][k] += 2.0*f1*f2*f3*f4;
                        }
                        std::complex<double> f1 = Couplings_Lepton_Higgs_22[a][b][i][j];
                        std::complex<double> f2 = Couplings_Lepton_Higgs_21[b][a][k];
                        std::complex<double> f3 = fbase(MassSquaredLepton[a],MassSquaredLepton[b]) - C_CWcbFermion + 0.5;
                        LeptonPart[i][j][k] += 3.0*f1*f2*f3;
                    }
                }
            }
        }
    }

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            for(std::size_t k=0;k<NHiggs;k++)
            {
                HiggspartSym[i][j][k] = Higgspart[i][j][k] + Higgspart[i][k][j];
                HiggspartSym[i][j][k] += Higgspart[j][i][k] + Higgspart[j][k][i];
                HiggspartSym[i][j][k] += Higgspart[k][i][j] + Higgspart[k][j][i];
                HiggspartSym[i][j][k] *= 1.0/6.0;

                GaugePartSym[i][j][k] = GaugePart[i][j][k] + GaugePart[i][k][j];
                GaugePartSym[i][j][k] += GaugePart[j][i][k] + GaugePart[j][k][i];
                GaugePartSym[i][j][k] += GaugePart[k][i][j] + GaugePart[k][j][i];
                GaugePartSym[i][j][k] *= 1.0/6.0;

                QuarkPartSym[i][j][k] = QuarkPart[i][j][k] + QuarkPart[i][k][j];
                QuarkPartSym[i][j][k] += QuarkPart[j][i][k] + QuarkPart[j][k][i];
                QuarkPartSym[i][j][k] += QuarkPart[k][i][j] + QuarkPart[k][j][i];
                QuarkPartSym[i][j][k] *= 1.0/6.0;

                LeptonPartSym[i][j][k] = LeptonPart[i][j][k] + LeptonPart[i][k][j];
                LeptonPartSym[i][j][k] += LeptonPart[j][i][k] + LeptonPart[j][k][i];
                LeptonPartSym[i][j][k] += LeptonPart[k][i][j] + LeptonPart[k][j][i];
                LeptonPartSym[i][j][k] *= 1.0/6.0;
            }
        }
    }

    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            for(std::size_t k=0;k<NHiggs;k++)
            {
                restmp[i][j][k] = 0.5*HiggspartSym[i][j][k];
                restmp[i][j][k] += 1.5*GaugePartSym[i][j][k];
                restmp[i][j][k] += -1.0*LeptonPartSym[i][j][k];
                restmp[i][j][k] += -3.0*QuarkPartSym[i][j][k];
            }
        }
    }

    for(std::size_t l=0;l<NHiggs;l++)
    {
        for(std::size_t m=0;m<NHiggs;m++)
        {
            for(std::size_t n=0;n<NHiggs;n++)
            {
                resGaugeBase[l][m][n] = 0;
                for(std::size_t i=0;i<NHiggs;i++)
                {
                    for(std::size_t j=0;j<NHiggs;j++)
                    {
                        for(std::size_t k=0;k<NHiggs;k++)
                        {
                            double RotFac = HiggsRotationMatrix[i][l]*HiggsRotationMatrix[j][m]*HiggsRotationMatrix[k][n];
                            resGaugeBase[l][m][n] += RotFac*restmp[i][j][k].real();
                        }
                    }
                }
                resGaugeBase[l][m][n] *= epsilon;
                if(std::abs(resGaugeBase[l][m][n]) < NumZero) resGaugeBase[l][m][n] = 0;
            }
        }
    }

    for(std::size_t l=0;l<NHiggs;l++)
    {
        for(std::size_t m=0;m<NHiggs;m++)
        {
            for(std::size_t n=0;n<NHiggs;n++)
            {
                res.push_back(resGaugeBase[l][m][n]);
            }
        }
    }

    return res;


}


std::vector<double> Class_Potential_Origin::WeinbergFourthDerivative() const{
	if(not CalcCouplingsdone){
		std::string retmes = __func__;
		retmes += "tries to use Physical couplings but they are not initialised.";
		throw std::runtime_error(retmes);
	}

	const double NumZero = std::pow(10,-10);
	double epsilon = 1.0/(10.0*M_PI*M_PI);

	std::vector<double> res;

	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> restmp;
	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> QuarkPart;
	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> LeptonPart;
	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> QuarkPartSym;
	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> LeptonPartSym;

	restmp.resize(NHiggs);
	QuarkPart.resize(NHiggs);
	LeptonPart.resize(NHiggs);
	QuarkPartSym.resize(NHiggs);
	LeptonPartSym.resize(NHiggs);
	for(std::size_t i=0; i<NHiggs;i++)
	{
		restmp[i].resize(NHiggs);
		QuarkPart[i].resize(NHiggs);
		LeptonPart[i].resize(NHiggs);
		QuarkPartSym[i].resize(NHiggs);
		LeptonPartSym[i].resize(NHiggs);
		for(std::size_t j=0; j<NHiggs;j++){
			restmp[i][j].resize(NHiggs);
			QuarkPart[i][j].resize(NHiggs);
			LeptonPart[i][j].resize(NHiggs);
			QuarkPartSym[i][j].resize(NHiggs);
			LeptonPartSym[i][j].resize(NHiggs);
			for(std::size_t k=0; k<NHiggs; k++){
				restmp[i][j][k].resize(NHiggs);
				QuarkPart[i][j][k].resize(NHiggs);
				LeptonPart[i][j][k].resize(NHiggs);
				QuarkPartSym[i][j][k].resize(NHiggs);
				LeptonPartSym[i][j][k].resize(NHiggs);
			}
		}
	}
	
	std::vector<std::vector<std::vector<std::vector<double>>>> resGaugeBase(NHiggs, 
		std::vector<std::vector<std::vector<double>>>(NHiggs, std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs))));
	std::vector<std::vector<std::vector<std::vector<double>>>> HiggsPart(NHiggs, 
		std::vector<std::vector<std::vector<double>>>(NHiggs, std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs))));
	std::vector<std::vector<std::vector<std::vector<double>>>> GaugePart(NHiggs, 
		std::vector<std::vector<std::vector<double>>>(NHiggs, std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs))));
	std::vector<std::vector<std::vector<std::vector<double>>>> HiggsPartSym(NHiggs, 
		std::vector<std::vector<std::vector<double>>>(NHiggs, std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs))));
	std::vector<std::vector<std::vector<std::vector<double>>>> GaugePartSym(NHiggs, 
		std::vector<std::vector<std::vector<double>>>(NHiggs, std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs))));


	for(std::size_t i1=0;i1<NHiggs;i1++)
	{
		for(std::size_t i2=0;i2<NHiggs;i2++)
		{
			for(std::size_t i3=0;i3<NHiggs;i3++)
			{
				for(std::size_t i4=0;i4<NHiggs;i4++)
				{
					HiggsPart[i1][i2][i3][i4] = 0;
					for(std::size_t a=0;a<NHiggs;a++)
					{
						for(std::size_t b=0;b<NHiggs;b++)
						{
							for(std::size_t c=0;c<NHiggs;c++)
							{
								double f1h = fbaseTri(MassSquaredHiggs[a], MassSquaredHiggs[b], MassSquaredHiggs[c]);
								double f2h = Couplings_Higgs_Triple[a][b][i4];
								double f3h = Couplings_Higgs_Quartic[b][c][i1][i2];
								double f4h = Couplings_Higgs_Triple[c][a][i3];
								HiggsPart[i1][i2][i3][i4] += 0.2e1*f1h*f2h*f3h*f4h;
								for(std::size_t d=0;d<NHiggs;d++)
								{
									double f1h1 = fbasequad(MassSquaredHiggs[a],MassSquaredHiggs[b],MassSquaredHiggs[c],MassSquaredHiggs[d]);
									double f2h2 = Couplings_Higgs_Triple[a][b][i4];
									double f3h3 = Couplings_Higgs_Triple[b][c][i1];
									double f4h4 = Couplings_Higgs_Triple[c][d][i2];
									double f5h5 = Couplings_Higgs_Triple[d][a][i3];
									HiggsPart[i1][i2][i3][i4] += f1h1*f2h2*f3h3*f4h4*f5h5;
								}
							}
							double Br = fbase(MassSquaredHiggs[a],MassSquaredHiggs[b]) - C_CWcbHiggs + 0.5;
							double coup1 = Couplings_Higgs_Quartic[a][b][i1][i4];
							double coup2 = Couplings_Higgs_Quartic[b][a][i2][i3];
							HiggsPart[i1][i2][i3][i4] += coup1*coup2*Br;
						}
					}
					GaugePart[i1][i2][i3][i4] = 0;
					for(std::size_t a=0;a<NGauge;a++)
                                        {
                                                for(std::size_t b=0;b<NGauge;b++)
                                                {
                                                        for(std::size_t c=0;c<NGauge;c++)
                                                        {
                                                                double f1g = fbaseTri(MassSquaredGauge[a], MassSquaredGauge[b], MassSquaredGauge[c]);
                                                                double f2g = Couplings_Gauge_Higgs_21[a][b][i4];
                                                                double f3g = Couplings_Gauge_Higgs_22[b][c][i1][i2];
                                                                double f4g = Couplings_Gauge_Higgs_21[c][a][i3];
                                                                GaugePart[i1][i2][i3][i4] +=  0.2e1*f1g*f2g*f3g*f4g;
                                                                for(std::size_t d=0;d<NGauge;d++)
                                                                {
                                                                        double f1g1 = fbasequad(MassSquaredGauge[a],MassSquaredGauge[b],MassSquaredGauge[c], MassSquaredGauge[d]);
                                                                        double f2g2 = Couplings_Gauge_Higgs_21[a][b][i4];
                                                                        double f3g3 = Couplings_Gauge_Higgs_21[b][c][i1];
                                                                        double f4g4 = Couplings_Gauge_Higgs_21[c][d][i2];
                                                                        double f5g5 = Couplings_Gauge_Higgs_21[d][a][i3];
                                                                        GaugePart[i1][i2][i3][i4] += f1g1*f2g2*f3g3*f4g4*f5g5;
                                                                }
                                                        }
                                                        double Br = fbase(MassSquaredGauge[a],MassSquaredGauge[b]) - C_CWcbGB + 0.5;
                                                        double coup1 = Couplings_Gauge_Higgs_22[a][b][i1][i4];
                                                        double coup2 = Couplings_Gauge_Higgs_22[b][a][i2][i3];
                                                        GaugePart[i1][i2][i3][i4] += coup1*coup2*Br;
						}
					}
					QuarkPart[i1][i2][i3][i4] = 0;
                                        for(std::size_t a=0;a<NQuarks;a++)
					{
						for(std::size_t b=0;b<NQuarks;b++)
						{
							for(std::size_t c=0;c<NQuarks;c++)
							{
								std::complex<double> f1q = fbaseTri(MassSquaredQuark[a], MassSquaredQuark[b], MassSquaredQuark[c]);
								std::complex<double> f2q = Couplings_Quark_Higgs_21[a][b][i4];
								std::complex<double> f3q = Couplings_Quark_Higgs_22[b][c][i1][i2];
								std::complex<double> f4q = Couplings_Quark_Higgs_21[c][a][i3];
								QuarkPart[i1][i2][i3][i4] = QuarkPart[i1][i2][i3][i4] + 0.2e1*f1q*f2q*f3q*f4q;
								for(std::size_t d=0;d<NQuarks;d++)
								{
									std::complex<double> f1q1 = fbasequad(MassSquaredQuark[a],MassSquaredQuark[b],MassSquaredQuark[c], MassSquaredQuark[d]);
									std::complex<double> f2q2 = Couplings_Quark_Higgs_21[a][b][i4];
									std::complex<double> f3q3 = Couplings_Quark_Higgs_21[b][c][i1];
									std::complex<double> f4q4 = Couplings_Quark_Higgs_21[c][d][i2];
									std::complex<double> f5q5 = Couplings_Quark_Higgs_21[d][a][i3];
									QuarkPart[i1][i2][i3][i4] += f1q1*f2q2*f3q3*f4q4*f5q5;
								}
							}
							std::complex<double> Br = fbase(MassSquaredQuark[a],MassSquaredQuark[b]) - C_CWcbFermion + 0.5;
							std::complex<double> coup1 = Couplings_Quark_Higgs_22[a][b][i1][i4];
							std::complex<double> coup2 = Couplings_Quark_Higgs_22[b][a][i2][i3];
							QuarkPart[i1][i2][i3][i4] += coup1*coup2*Br;
						}
					}
					LeptonPart[i1][i2][i3][i4] = 0;
					for(std::size_t a=0;a<NLepton;a++)
					{
						for(std::size_t b=0;b<NLepton;b++)
						{
							for(std::size_t c=0;c<NLepton;c++)
							{
								std::complex<double> f1l = fbaseTri(MassSquaredLepton[a],MassSquaredLepton[b],MassSquaredLepton[c]);
								std::complex<double> f2l = Couplings_Lepton_Higgs_21[a][b][i4];
								std::complex<double> f3l = Couplings_Lepton_Higgs_22[b][c][i1][i2];
								std::complex<double> f4l = Couplings_Lepton_Higgs_21[c][a][i3];
								LeptonPart[i1][i2][i3][i4] += 0.2e1*f1l*f2l*f3l*f4l;
								for(std::size_t d=0;d<NLepton;d++)
								{
									std::complex<double> f1l1 = fbasequad(MassSquaredLepton[a],MassSquaredLepton[b],MassSquaredLepton[c], MassSquaredLepton[d]);
									std::complex<double> f2l1 = Couplings_Lepton_Higgs_21[a][b][i4];
									std::complex<double> f3l1 = Couplings_Lepton_Higgs_21[b][c][i1];
									std::complex<double> f4l1 = Couplings_Lepton_Higgs_21[c][d][i2];
									std::complex<double> f5l1 = Couplings_Lepton_Higgs_21[d][a][i3];
									LeptonPart[i1][i2][i3][i4] += f1l1*f2l1*f3l1*f4l1*f5l1;
								}
 							}
							std::complex<double> Br = fbase(MassSquaredLepton[a],MassSquaredLepton[b]) - C_CWcbFermion + 0.5;
							std::complex<double> coup1 = Couplings_Lepton_Higgs_22[a][b][i1][i4];
							std::complex<double> coup2 = Couplings_Lepton_Higgs_22[b][a][i2][i3];
							LeptonPart[i1][i2][i3][i4] += coup1*coup2*Br;
						}
					}
				}
			}
		}
	}

	for(std::size_t i1=0;i1<NHiggs;i1++)
	{
		for(std::size_t i2=0;i2<NHiggs;i2++)
		{
			for(std::size_t i3=0;i3<NHiggs;i3++)
			{
				for(std::size_t i4=0;i4<NHiggs;i4++)
				{
					HiggsPartSym[i1][i2][i3][i4] = HiggsPart[i1][i2][i3][i4] + HiggsPart[i1][i2][i4][i3] + HiggsPart[i1][i4][i2][i3] + HiggsPart[i1][i4][i3][i2] + HiggsPart[i1][i3][i4][i2] + HiggsPart[i1][i3][i2][i4];
					HiggsPartSym[i1][i2][i3][i4] += HiggsPart[i2][i1][i3][i4] + HiggsPart[i2][i1][i4][i3] + HiggsPart[i2][i3][i1][i4] + HiggsPart[i2][i3][i4][i1] + HiggsPart[i2][i4][i1][i3] + HiggsPart[i2][i4][i3][i1];
					HiggsPartSym[i1][i2][i3][i4] += HiggsPart[i3][i1][i2][i4] + HiggsPart[i3][i1][i4][i2] + HiggsPart[i3][i2][i1][i4] + HiggsPart[i3][i2][i4][i1] + HiggsPart[i3][i4][i1][i2] + HiggsPart[i3][i4][i2][i1];
					HiggsPartSym[i1][i2][i3][i4] += HiggsPart[i4][i1][i2][i3] + HiggsPart[i4][i1][i3][i2] + HiggsPart[i4][i2][i3][i1] + HiggsPart[i4][i2][i1][i3] + HiggsPart[i4][i3][i1][i2]+ HiggsPart[i4][i3][i2][i1];
					HiggsPartSym[i1][i2][i3][i4] *= 1.0/24.0;

					GaugePartSym[i1][i2][i3][i4] = GaugePart[i1][i2][i3][i4] + GaugePart[i1][i2][i4][i3] + GaugePart[i1][i4][i2][i3] + GaugePart[i1][i4][i3][i2] + GaugePart[i1][i3][i4][i2]+ GaugePart[i1][i3][i2][i4];
					GaugePartSym[i1][i2][i3][i4] += GaugePart[i2][i1][i3][i4] + GaugePart[i2][i1][i4][i3] + GaugePart[i2][i3][i1][i4] + GaugePart[i2][i3][i4][i1] + GaugePart[i2][i4][i1][i3]+ GaugePart[i2][i4][i3][i1];
					GaugePartSym[i1][i2][i3][i4] += GaugePart[i3][i1][i2][i4] + GaugePart[i3][i1][i4][i2] + GaugePart[i3][i2][i1][i4] + GaugePart[i3][i2][i4][i1] + GaugePart[i3][i4][i1][i2]+ GaugePart[i3][i4][i2][i1];
					GaugePartSym[i1][i2][i3][i4] += GaugePart[i4][i1][i2][i3] + GaugePart[i4][i1][i3][i2] + GaugePart[i4][i2][i3][i1] + GaugePart[i4][i2][i1][i3] + GaugePart[i4][i3][i1][i2]+ GaugePart[i4][i3][i2][i1];
					GaugePartSym[i1][i2][i3][i4] *= 1.0/24.0;

					QuarkPartSym[i1][i2][i3][i4] = QuarkPart[i1][i2][i3][i4] + QuarkPart[i1][i2][i4][i3] + QuarkPart[i1][i4][i2][i3] + QuarkPart[i1][i4][i3][i2] + QuarkPart[i1][i3][i4][i2] + QuarkPart[i1][i3][i2][i4];
					QuarkPartSym[i1][i2][i3][i4] += QuarkPart[i2][i1][i3][i4] + QuarkPart[i2][i1][i4][i3] + QuarkPart[i2][i3][i1][i4] + QuarkPart[i2][i3][i4][i1] + QuarkPart[i2][i4][i1][i3]+ QuarkPart[i2][i4][i3][i1];
					QuarkPartSym[i1][i2][i3][i4] += QuarkPart[i3][i1][i2][i4] + QuarkPart[i3][i1][i4][i2] + QuarkPart[i3][i2][i1][i4] + QuarkPart[i3][i2][i4][i1] + QuarkPart[i3][i4][i1][i2]+ QuarkPart[i3][i4][i2][i1];
					QuarkPartSym[i1][i2][i3][i4] += QuarkPart[i4][i1][i2][i3] + QuarkPart[i4][i1][i3][i2] + QuarkPart[i4][i2][i3][i1] + QuarkPart[i4][i2][i1][i3] + QuarkPart[i4][i3][i1][i2]+ QuarkPart[i4][i3][i2][i1];
					QuarkPartSym[i1][i2][i3][i4] *= 1.0/24.0;

					LeptonPartSym[i1][i2][i3][i4] = LeptonPart[i1][i2][i3][i4] + LeptonPart[i1][i2][i4][i3] + LeptonPart[i1][i4][i2][i3] + LeptonPart[i1][i4][i3][i2] + LeptonPart[i1][i3][i4][i2]+ LeptonPart[i1][i3][i2][i4];
					LeptonPartSym[i1][i2][i3][i4] += LeptonPart[i2][i1][i3][i4] + LeptonPart[i2][i1][i4][i3] + LeptonPart[i2][i3][i1][i4] + LeptonPart[i2][i3][i4][i1] + LeptonPart[i2][i4][i1][i3]+ LeptonPart[i2][i4][i3][i1];
					LeptonPartSym[i1][i2][i3][i4] += LeptonPart[i3][i1][i2][i4] + LeptonPart[i3][i1][i4][i2] + LeptonPart[i3][i2][i1][i4] + LeptonPart[i3][i2][i4][i1] + LeptonPart[i3][i4][i1][i2]+ LeptonPart[i3][i4][i2][i1];
					LeptonPartSym[i1][i2][i3][i4] += LeptonPart[i4][i1][i2][i3] + LeptonPart[i4][i1][i3][i2] + LeptonPart[i4][i2][i3][i1] + LeptonPart[i4][i2][i1][i3] + LeptonPart[i4][i3][i1][i2] + LeptonPart[i4][i3][i2][i1];
					LeptonPartSym[i1][i2][i3][i4] *= 1.0/24.0;
				}
			}
		}
	}



	for(std::size_t i1=0;i1<NHiggs;i1++)
	{
		for(std::size_t i2=0;i2<NHiggs;i2++)
		{
			for(std::size_t i3=0;i3<NHiggs;i3++)
			{
				for(std::size_t i4=0;i4<NHiggs;i4++)
				{
					restmp[i1][i2][i3][i4] = 0.5*HiggsPartSym[i1][i2][i3][i4];
					restmp[i1][i2][i3][i4] += 1.5*GaugePartSym[i1][i2][i3][i4];
					restmp[i1][i2][i3][i4] += -1.0*LeptonPartSym[i1][i2][i3][i4];
					restmp[i1][i2][i3][i4] += -3.0*QuarkPartSym[i1][i2][i3][i4];
				}
			}
		}
	}

	for(std::size_t j1=0;j1<NHiggs;j1++)
	{
		for(std::size_t j2=0;j2<NHiggs;j2++)
		{
			for(std::size_t j3=0;j3<NHiggs;j3++)
 			{
				for(std::size_t j4=0; j4<NHiggs;j4++)
				{
					resGaugeBase[j1][j2][j3][j4] = 0;
					for(std::size_t i1=0;i1<NHiggs;i1++)
					{
						for(std::size_t i2=0;i2<NHiggs;i2++)
						{
							for(std::size_t i3=0;i3<NHiggs;i3++)
							{
								for(std::size_t i4=0;i4<NHiggs;i4++)
								{
									double RotFac = HiggsRotationMatrix[i1][j1]*HiggsRotationMatrix[i2][j2]*HiggsRotationMatrix[i3][j3]*HiggsRotationMatrix[i4][j4];
									resGaugeBase[j1][j2][j3][j4] += 3*RotFac*restmp[i1][i2][i3][i4].real();
								}
							}
						}
					}
					resGaugeBase[j1][j2][j3][j4] *= epsilon;
					if(std::abs(resGaugeBase[j1][j2][j3][j4]) < NumZero) resGaugeBase[j1][j2][j3][j4] = 0;
				}
			}
		}
	}


	for(std::size_t j1=0;j1<NHiggs;j1++)
	{
		for(std::size_t j2=0;j2<NHiggs;j2++)
		{
			for(std::size_t j3=0;j3<NHiggs;j3++)
			{
				for(std::size_t j4=0;j4<NHiggs;j4++)
				{
					res.push_back(resGaugeBase[j1][j2][j3][j4]);
				}
			}
 		}
	}

	return res;
}


std::vector<double> Class_Potential_Origin::HiggsMassesSquared(
        const std::vector<double>& v,
        const double& Temp,
        const int& diff) const
{
    std::vector<double> res;
    if(v.size() != nVEV and v.size() != NHiggs){
      std::string ErrorString = std::string("You have called ") + std::string(__func__)
        + std::string(" with an invalid vev configuration. Your vev is of dimension ")
      + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
      throw std::runtime_error(ErrorString);
    }
    if(v.size() == nVEV and nVEV != NHiggs){
    	std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
    	    			<< nVEV << " while it should have " << NHiggs
    					<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
      std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      res = HiggsMassesSquared(Transformedv,Temp,diff);
      return res;
    }
    if(!SetCurvatureDone) {
//        SetCurvatureArrays();
        throw std::runtime_error("SetCurvatureDone is not set. The Model is not initiliased correctly");
    }
    MatrixXd MassMatrix(NHiggs,NHiggs);
    double ZeroMass = std::pow(10,-5);
    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=0;j<NHiggs;j++)
        {
            MassMatrix(i,j) = Curvature_Higgs_L2[i][j];
            for(std::size_t k=0;k<NHiggs;k++)
            {
                MassMatrix(i,j) += Curvature_Higgs_L3[i][j][k]*v[k];
                for(std::size_t l=0;l<NHiggs;l++) MassMatrix(i,j) += 0.5*Curvature_Higgs_L4[i][j][k][l] * v[k]*v[l];
            }

            if(Temp != 0)
            {
                MassMatrix(i,j) += DebyeHiggs[i][j]*std::pow(Temp,2);
            }
        }
    }

    if(diff == 0 and res.size() == 0)
    {
        SelfAdjointEigenSolver<MatrixXd> es(MassMatrix,EigenvaluesOnly);
        for(std::size_t i =0;i<NHiggs;i++)
        {
            double tmp = es.eigenvalues()[i];
            if(std::abs(tmp) < ZeroMass ) res.push_back(0);
            else res.push_back(tmp);
        }
    }
    else if(diff == 0 and res.size() == NHiggs){
    	SelfAdjointEigenSolver<MatrixXd> es(MassMatrix,EigenvaluesOnly);
        for(std::size_t i =0;i<NHiggs;i++)
  		{
  			double tmp = es.eigenvalues()[i];
  			if(std::abs(tmp) < ZeroMass ) tmp = 0;
  			res[i] = tmp;
  		}
    }
    else if(diff == 0 and res.size()!= 0 and res.size() != NHiggs){
    	std::cout << "Something went wrong in " << __func__ << std::endl;
    	std::cout << __func__ << "Is calculating the mass for " << NHiggs << "fields but the resolution vector has a size of "
    			<< res.size() << ". This should be zero or " << NHiggs << std::endl;
    }
    else if((size_t) diff <= NHiggs and diff >0)
    {
        MatrixXd Diff(NHiggs,NHiggs);
        std::size_t x0 = diff -1;
        for(std::size_t i=0;i<NHiggs;i++)
        {
            for(std::size_t j=0;j<NHiggs;j++)
            {
                Diff(i,j) = Curvature_Higgs_L3[i][j][x0];
                for(std::size_t k=0;k<NHiggs;k++)
                {
                    Diff(i,j) += Curvature_Higgs_L4[i][j][x0][k]*v[k];
                }
            }
        }
        MatrixXcd MassCast(NHiggs,NHiggs);
        MassCast = MassMatrix;
        MatrixXcd DiffCast(NHiggs,NHiggs);
        DiffCast = Diff;
        res = FirstDerivativeOfEigenvalues(MassCast,DiffCast);
    }
    else if(diff == -1){
    	MatrixXd Diff(NHiggs,NHiggs);
        for(std::size_t i=0;i<NHiggs;i++){
            for(std::size_t j=0;j<NHiggs;j++){
    			Diff(i,j) = 2* DebyeHiggs[i][j] * Temp;
    		}
    	}

    	MatrixXcd MassCast(NHiggs,NHiggs);
		MassCast = MassMatrix;
		MatrixXcd DiffCast(NHiggs,NHiggs);
		DiffCast = Diff;
        res=FirstDerivativeOfEigenvalues(MassCast,DiffCast);
    }

    return  res;


}


std::vector<double> Class_Potential_Origin::GaugeMassesSquared(
        const std::vector<double>& v,
        const double& Temp,
        const int& diff) const
{
    std::vector<double> res;
    if(v.size() != nVEV and v.size() != NHiggs){
      std::string ErrorString = std::string("You have called ") + std::string(__func__)
        + std::string(" with an invalid vev configuration. Your vev is of dimension ")
      + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
      throw std::runtime_error(ErrorString);
    }
    if(v.size() == nVEV and nVEV != NHiggs){
    	std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
    			<< nVEV << " while it should have " << NHiggs
				<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
      std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      res =GaugeMassesSquared(Transformedv,Temp,diff);
      return res ;
    }
    if(!SetCurvatureDone) {
//        SetCurvatureArrays();
        std::string retmes = __func__;
        retmes += "was called while the model was not initialised correctly.\n";
        throw std::runtime_error(retmes);
    }
    MatrixXd MassMatrix(NGauge,NGauge);
    double ZeroMass = std::pow(10,-5);
    for(std::size_t i=0;i<NGauge;i++)
    {
        for(std::size_t j=0;j<NGauge;j++)
        {
            MassMatrix(i,j) = 0;
            for(std::size_t k=0;k<NHiggs;k++)
            {
                for(std::size_t l=0;l<NHiggs;l++) MassMatrix(i,j) += 0.5*Curvature_Gauge_G2H2[i][j][k][l] * v[k]*v[l];
            }
            if(Temp != 0)
            {
                MassMatrix(i,j) += DebyeGauge[i][j]*std::pow(Temp,2);
            }
        }
    }


    if(diff == 0)
    {

        SelfAdjointEigenSolver<MatrixXd> es(MassMatrix,EigenvaluesOnly);
        for(std::size_t i =0;i<NGauge;i++)
        {
            double tmp = es.eigenvalues()[i];
            if(std::abs(tmp) < ZeroMass ) res.push_back(0);
            else res.push_back(tmp);
        }
    }
    else if(diff > 0 and static_cast<size_t>(diff) <= NHiggs)
    {
        std::size_t i = diff -1;
        MatrixXd Diff(NGauge,NGauge);
        Diff = MatrixXd::Zero(NGauge,NGauge);
        for(std::size_t a=0;a<NGauge;a++)
        {
            for(std::size_t b=0;b<NGauge;b++)
            {
                for(std::size_t j=0;j<NHiggs;j++) Diff(a,b) += Curvature_Gauge_G2H2[a][b][i][j]*v[j];
            }
        }
        MatrixXcd MassCast(NGauge,NGauge);
        MassCast = MassMatrix;
        MatrixXcd DiffCast(NGauge,NGauge);
        DiffCast = Diff;
        res=FirstDerivativeOfEigenvalues(MassCast,DiffCast);
    }
    else if(diff == -1){
		MatrixXd Diff(NGauge,NGauge);
        for(std::size_t i=0;i<NGauge;i++){
            for(std::size_t j=0;j<NGauge;j++){
				Diff(i,j) = 2* DebyeGauge[i][j] * Temp;
			}
		}

		MatrixXcd MassCast(NGauge,NGauge);
		MassCast = MassMatrix;
		MatrixXcd DiffCast(NGauge,NGauge);
		DiffCast = Diff;
        res=FirstDerivativeOfEigenvalues(MassCast,DiffCast);
	}

    return res;


}


std::vector<double> Class_Potential_Origin::QuarkMassesSquared(const std::vector<double>& v, const int& diff) const
{
    std::vector<double> res;
	if(v.size() != nVEV and v.size() != NHiggs){
	  std::string ErrorString = std::string("You have called ") + std::string(__func__)
		+ std::string(" with an invalid vev configuration. Your vev is of dimension ")
	  + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
	  throw std::runtime_error(ErrorString);
	}
	if(v.size() == nVEV and nVEV != NHiggs){
		std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
		    			<< nVEV << " while it should have " << NHiggs
						<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
	  std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      res = QuarkMassesSquared(Transformedv,diff);
      return res;
	}
    if(!SetCurvatureDone) {
//        SetCurvatureArrays(); std::cout << "Reset of Set Curvature " << std::endl;
        std::string retmes = __func__;
        retmes += "was called while the model was not initialised correctly.\n";
        throw std::runtime_error(retmes);
    };
    MatrixXcd MassMatrix(NQuarks,NQuarks),MIJ(NQuarks,NQuarks);
    MIJ = MatrixXcd::Zero(NQuarks,NQuarks);
    double ZeroMass = std::pow(10,-10);
    for(std::size_t i=0;i<NQuarks;i++)
    {
        for(std::size_t j=0;j<NQuarks;j++)
        {

            for(std::size_t k=0;k<NHiggs;k++)
            {
                MIJ(i,j) += Curvature_Quark_F2H1[i][j][k]*v[k];
            }
        }
    }

    MassMatrix = MIJ.conjugate()*MIJ;

    if(diff <= 0) // No temperature dependent part here
    {
        SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix,EigenvaluesOnly);
        for(std::size_t i =0;i<NQuarks;i++)
        {
            double tmp = es.eigenvalues().real()[i];
            if(std::abs(tmp) < ZeroMass ) res.push_back(0);
            else res.push_back(tmp);
        }
    }
    else if(diff > 0 and static_cast<size_t>(diff) <= NHiggs){
      std::size_t m = diff-1;
      MatrixXcd Diff(NQuarks,NQuarks);
      Diff = MatrixXcd::Zero(NQuarks,NQuarks);
      for(std::size_t a=0;a<NQuarks;a++){
        for(std::size_t b=0;b<NQuarks;b++){
          for(std::size_t i=0;i<NQuarks;i++){
            for(std::size_t l=0;l<NHiggs;l++){
              Diff(a,b) += std::conj(Curvature_Quark_F2H1[a][i][m]) * Curvature_Quark_F2H1[i][b][l] * v[l];
              Diff(a,b) += std::conj(Curvature_Quark_F2H1[a][i][l]) * Curvature_Quark_F2H1[i][b][m] * v[l];
            }
          }
        }
      }

      res=FirstDerivativeOfEigenvalues(MassMatrix,Diff);

      for(std::size_t j=0;j<res.size();j++){
    	  if(std::isnan(res.at(j))){
    		  std::cout << "MassMatrix = \n" << MassMatrix << "\nDiff = \n" << Diff << std::endl;
    		  std::cout << "Fermion Masses : " ;
              for(std::size_t i=0;i<NQuarks;i++) std::cout << std::sqrt(std::abs(res.at(i))) << sep;
    		  std::cout << std::endl;
    		  std::cout << "VEV fields : ";
              for(std::size_t i=0;i<v.size();i++) std::cout << v.at(i) << sep;
    		  std::cout << std::endl;




              for(std::size_t l=0;l<NHiggs;l++){

				  std::cout << "Curvature_Quark * v an Higgs  =  :" << l << "\n";
                  for(std::size_t a=0;a<NQuarks;a++){
                      for(std::size_t i=0;i<NQuarks;i++){
                          std::cout << Curvature_Quark_F2H1[a][i][l] * v[l] << sep;
					  }
					  std::cout << std::endl;
				  }
				  std::cout << "conj Curvature_Quark an Higgs = :" << l << "\n";
                  for(std::size_t a=0;a<NQuarks;a++){
                      for(std::size_t i=0;i<NQuarks;i++){
                          std::cout << std::conj(Curvature_Quark_F2H1[a][i][l]) *v[l] << sep;
					  }
					  std::cout << std::endl;
				  }
    		  }

    		  std::string retmessage = "Nan found in ";
			  retmessage+= __func__;
			  retmessage+= " at deriv number ";
			  retmessage+= std::to_string(j);
			  retmessage+= " and m = ";
			  retmessage+= std::to_string(m);
			  throw std::runtime_error(retmessage);
    	  }
      }

    }

    return  res;
}


std::vector<double> Class_Potential_Origin::LeptonMassesSquared(const std::vector<double>& v, const int& diff) const
{
    std::vector<double> res;
	if(v.size() != nVEV and v.size() != NHiggs){
	      std::string ErrorString = std::string("You have called ") + std::string(__func__)
	        + std::string(" with an invalid vev configuration. Your vev is of dimension ")
	      + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
	      throw std::runtime_error(ErrorString);
	}
	if(v.size() == nVEV and nVEV != NHiggs){
		std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
		    			<< nVEV << " while it should have " << NHiggs
						<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
	  std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      res = LeptonMassesSquared(Transformedv,diff);
      return res;
	}
    if(!SetCurvatureDone) {
//        SetCurvatureArrays();
        std::string retmes = __func__;
        retmes += "was called while the model was not initialised correctly.\n";
        throw std::runtime_error(retmes);
    }
    MatrixXcd MassMatrix(NLepton,NLepton),MIJ(NLepton,NLepton);
    MIJ = MatrixXcd::Zero(NLepton,NLepton);
    double ZeroMass = std::pow(10,-10);
    for(std::size_t i=0;i<NLepton;i++)
    {
        for(std::size_t j=0;j<NLepton;j++)
        {

            for(std::size_t k=0;k<NHiggs;k++)
            {
                MIJ(i,j) += Curvature_Lepton_F2H1[i][j][k]*v[k];
            }
        }
    }

    MassMatrix = MIJ.conjugate()*MIJ;

    if(diff <= 0) // no temperature part here
	{
		SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix,EigenvaluesOnly);
        for(std::size_t i =0;i<NLepton;i++)
		{
			double tmp = es.eigenvalues().real()[i];
			if(std::abs(tmp) < ZeroMass ) res.push_back(0);
			else res.push_back(tmp);
		}
	}
    else if(diff > 0 and static_cast<size_t>(diff) <= NHiggs){
      std::size_t m = diff-1;
	  MatrixXcd Diff(NLepton,NLepton);
	  Diff = MatrixXcd::Zero(NLepton,NLepton);
      for(std::size_t a=0;a<NLepton;a++){
        for(std::size_t b=0;b<NLepton;b++){
          for(std::size_t i=0;i<NLepton;i++){
            for(std::size_t l=0;l<NHiggs;l++){
			  Diff(a,b) += std::conj(Curvature_Lepton_F2H1[a][i][m]) * Curvature_Lepton_F2H1[i][b][l] * v[l];
			  Diff(a,b) += std::conj(Curvature_Lepton_F2H1[a][i][l]) * Curvature_Lepton_F2H1[i][b][m] * v[l];
			}
		  }
		}
	  }

      res=FirstDerivativeOfEigenvalues(MassMatrix,Diff);

      for(std::size_t j=0;j<res.size();j++){
		  if(std::isnan(res.at(j))){
			  std::cout << "MassMatrix = \n" << MassMatrix << "\nDiff = \n" << Diff << std::endl;
			  std::cout << "Fermion Masses : " ;
              for(std::size_t i=0;i<NLepton;i++) std::cout << std::sqrt(std::abs(res.at(i))) << sep;
			  std::cout << std::endl;
			  std::cout << "VEV fields : ";
              for(std::size_t i=0;i<v.size();i++) std::cout << v.at(i) << sep;
			  std::cout << std::endl;




              for(std::size_t l=0;l<NHiggs;l++){

				  std::cout << "Curvature_Lepton * v an Higgs  =  :" << l << "\n";
                  for(std::size_t a=0;a<NLepton;a++){
                      for(std::size_t i=0;i<NLepton;i++){
                          std::cout << Curvature_Quark_F2H1[a][i][l] * v[l] << sep;
					  }
					  std::cout << std::endl;
				  }
				  std::cout << "conj Curvature_Lepton an Higgs = :" << l << "\n";
                  for(std::size_t a=0;a<NLepton;a++){
                      for(std::size_t i=0;i<NLepton;i++){
                          std::cout << std::conj(Curvature_Lepton_F2H1[a][i][l]) *v[l] << sep;
					  }
					  std::cout << std::endl;
				  }
			  }

			  std::string retmessage = "Nan found in ";
			  retmessage+= __func__;
			  retmessage+= " at deriv number ";
			  retmessage+= std::to_string(j);
			  retmessage+= " and m = ";
			  retmessage+= std::to_string(m);
			  throw std::runtime_error(retmessage);
		  }
	  }

	}

    return res;


}


double Class_Potential_Origin::VTree(const std::vector<double>& v, int diff, bool ForceExplicitCalculation) const
{
    double res = 0;

    if(not ForceExplicitCalculation){
		res = VTreeSimplified(v) ;
		if(UseVTreeSimplified and diff == 0)
		{

			return res;
		}
    }
    res = 0;


    if(diff == 0)
    {
        for(std::size_t i=0;i<NHiggs;i++)
        {
        	if(v[i] != 0)
        	{
				res += Curvature_Higgs_L1[i]*v[i];
                for(std::size_t j=0;j<NHiggs;j++)
				{
					if(v[j] !=0)
					{
						res += 0.5*Curvature_Higgs_L2[i][j]*v[i]*v[j];
                        for(std::size_t k=0;k<NHiggs;k++)
						{
							res += 1.0/6.0 * Curvature_Higgs_L3[i][j][k]*v[i]*v[j]*v[k];
                            for(std::size_t l=0;l<NHiggs;l++)
							{
								res += 1.0/24.0*Curvature_Higgs_L4[i][j][k][l]*v[i]*v[j]*v[k]*v[l];
							}
						}
					}
				}
        	}
        }
    }
    else if(diff > 0 and static_cast<size_t>(diff) <= NHiggs)
    {
        std::size_t i = diff -1;
        res = Curvature_Higgs_L1[i];
        for(std::size_t j=0;j<NHiggs;j++)
        {
            res += Curvature_Higgs_L2[i][j]*v[j];
            for(std::size_t k=0;k<NHiggs;k++)
            {
                res += 0.5*Curvature_Higgs_L3[i][j][k]*v[j]*v[k];
                for(std::size_t l=0;l<NHiggs;l++)
                {
                    res += 1.0/6.0*Curvature_Higgs_L4[i][j][k][l]*v[j]*v[k]*v[l];
                }
            }
        }
    }

    return res;

}



double Class_Potential_Origin::CounterTerm(const std::vector<double>& v, int diff, bool ForceExplicitCalculation) const
{
    double res = 0;
    if(not ForceExplicitCalculation){
    	res = VCounterSimplified(v) ;
    	if(UseVCounterSimplified and diff == 0) return res;
    }




    res =0;
    if(diff == 0)
    {
        for(std::size_t i=0;i<NHiggs;i++)
        {
            res += Curvature_Higgs_CT_L1[i]*v[i];
            for(std::size_t j=0;j<NHiggs;j++)
            {
                res += 0.5*Curvature_Higgs_CT_L2[i][j]*v[i]*v[j];
                for(std::size_t k=0;k<NHiggs;k++)
                {
                    res += 1.0/6.0 * Curvature_Higgs_CT_L3[i][j][k]*v[i]*v[j]*v[k];
                    for(std::size_t l=0;l<NHiggs;l++)
                    {
                        res += 1.0/24.0*Curvature_Higgs_CT_L4[i][j][k][l]*v[i]*v[j]*v[k]*v[l];
                    }
                }
            }
        }
    }
    else if(diff > 0 and static_cast<size_t>(diff) <= NHiggs)
    {
        std::size_t i = diff -1;
        res = Curvature_Higgs_CT_L1[i];
        for(std::size_t j=0;j<NHiggs;j++)
        {
            res += Curvature_Higgs_CT_L2[i][j]*v[j];
            for(std::size_t k=0;k<NHiggs;k++)
            {
                res += 0.5*Curvature_Higgs_CT_L3[i][j][k]*v[j]*v[k];
                for(std::size_t l=0;l<NHiggs;l++)
                {
                    res += 1.0/6.0*Curvature_Higgs_CT_L4[i][j][k][l]*v[j]*v[k]*v[l];
                }
            }
        }
    }

    return res;
}



double Class_Potential_Origin::VEff(const std::vector<double>& v,
        double Temp, int diff, int Order) const {
    if(v.size() != nVEV and v.size() != NHiggs){
      std::string ErrorString = std::string("You have called ") + std::string(__func__)
        + std::string(" with an invalid vev configuration. Your vev is of dimension ")
      + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
      throw std::runtime_error(ErrorString);
    }
    if(v.size() == nVEV and nVEV != NHiggs){
    	std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
    	    			<< nVEV << " while it should have " << NHiggs
    					<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
      std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      return VEff(Transformedv,Temp,diff);
    }

	double resOut = 0;
	resOut = VTree(v,diff);
	if(Order != 0 and not UseTreeLevel){
		resOut+= CounterTerm(v,diff);
		resOut+= V1Loop(v,Temp,diff);
	}
// for(std::size_t i=0;i<NHiggs;i++) resOut += DebyeHiggs[i][i]*0.5*std::pow(v.at(i),2)*std::pow(Temp,2);
	return resOut;
}



double Class_Potential_Origin::V1Loop(const std::vector<double>& v, double Temp, int diff) const
{
	double res = 0;

	/**
	 * In case of diff != 0 the mass vectors will directly have the derivatives of the masses
	 */
    std::vector<double> HiggsMassesVec,QuarkMassesVec,GaugeMassesVec,LeptonMassesVec,HiggsMassesZeroTempVec,GaugeMassesZeroTempVec;
    HiggsMassesVec=HiggsMassesSquared(v,Temp,diff);
    GaugeMassesVec=GaugeMassesSquared(v,Temp,diff);
    GaugeMassesZeroTempVec= GaugeMassesSquared(v,0,diff);
    QuarkMassesVec=QuarkMassesSquared(v,diff);
    LeptonMassesVec=LeptonMassesSquared(v,diff);


	if(diff == 0)
	{
        if(C_UseParwani)
        {
            for(std::size_t k=0;k<NHiggs;k++) res += boson(HiggsMassesVec[k],Temp,C_CWcbHiggs,0);
            for(std::size_t k=0;k<NGauge;k++) res += boson(GaugeMassesVec[k],Temp,C_CWcbGB,0);
            for(std::size_t k=0;k<NGauge;k++) res += 2*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,0);
            for(std::size_t k=0;k<NQuarks;k++) res += -6*fermion(QuarkMassesVec[k],Temp,0);
            for(std::size_t k=0;k<NLepton;k++) res += -2*fermion(LeptonMassesVec[k],Temp,0);
        }
        else{
            HiggsMassesZeroTempVec=HiggsMassesSquared(v,0);
            for(std::size_t k=0;k<NHiggs;k++) res += boson(HiggsMassesZeroTempVec[k],Temp,C_CWcbHiggs,0);
            for(std::size_t k=0;k<NGauge;k++) res += 3*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,0);
            double AddContQuark=0;
            for(std::size_t k=0;k<NQuarks;k++) AddContQuark += -2*fermion(QuarkMassesVec[k],Temp,0);
            for(std::size_t k=0;k<NColour;k++) res += AddContQuark;
            for(std::size_t k=0;k<NLepton;k++) res += -2*fermion(LeptonMassesVec[k],Temp,0);


            double VDebay = 0;
            for(std::size_t k=0;k<NHiggs;k++)
            {
                if( HiggsMassesVec[k] > 0 ) VDebay += std::pow(HiggsMassesVec[k],1.5);
                if( HiggsMassesZeroTempVec[k] > 0) VDebay += -std::pow(HiggsMassesZeroTempVec[k],1.5);
            }
            for(std::size_t k=0;k<NGauge;k++)
            {
                if(GaugeMassesVec[k] > 0 ) VDebay += std::pow(GaugeMassesVec[k],1.5);
                if(GaugeMassesZeroTempVec[k] > 0 ) VDebay += -std::pow(GaugeMassesZeroTempVec[k],1.5);
            }

            VDebay *= -Temp/(12*M_PI);
            res += VDebay;
        }
	}
    else if(diff > 0 and static_cast<size_t>(diff) <= NHiggs){
		if(C_UseParwani){
            for(std::size_t k=0;k<NHiggs;k++) {
				res += HiggsMassesVec.at(k+NHiggs)*boson(HiggsMassesVec.at(k),Temp,C_CWcbHiggs,diff);
			}
            for(std::size_t k=0;k<NGauge;k++) {
				res += GaugeMassesVec.at(k+NGauge)*boson(GaugeMassesVec.at(k),Temp,C_CWcbHiggs,diff);
			}
            for(std::size_t k=0;k<NGauge;k++) {
				res += 2*GaugeMassesZeroTempVec.at(k+NGauge)*boson(GaugeMassesZeroTempVec.at(k),Temp,C_CWcbHiggs,diff);
			}
            for(std::size_t k=0;k<NQuarks;k++) {
				res += -6*QuarkMassesVec.at(k+NQuarks)*fermion(QuarkMassesVec[k],Temp,diff);
			}
            for(std::size_t k=0;k<NLepton;k++) {
				res += -2*LeptonMassesVec.at(k+NLepton)*fermion(LeptonMassesVec[k],Temp,diff);
			}
		}
		else{
            HiggsMassesZeroTempVec=HiggsMassesSquared(v,0,diff);
            for(std::size_t k=0;k<NHiggs;k++) {
				res += HiggsMassesZeroTempVec.at(k+NHiggs)*boson(HiggsMassesZeroTempVec[k],Temp,C_CWcbHiggs,diff);
			}
            for(std::size_t k=0;k<NGauge;k++) {
				res += 3*GaugeMassesZeroTempVec.at(k+NGauge)*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,diff);
			}
			double AddContQuark=0;
            for(std::size_t k=0;k<NQuarks;k++) {
				AddContQuark += -2*QuarkMassesVec.at(k+NQuarks)*fermion(QuarkMassesVec[k],Temp,diff);
			}
            for(std::size_t k=0;k<NColour;k++) res += AddContQuark;
            for(std::size_t k=0;k<NLepton;k++) {
				res += -2*LeptonMassesVec.at(k+NLepton)*fermion(LeptonMassesVec[k],Temp,diff);
			}


			double VDebay = 0;
            for(std::size_t k=0;k<NHiggs;k++)
			{
				if( HiggsMassesVec[k] > 0 ) {
					VDebay += 1.5 * HiggsMassesVec.at(k+NHiggs) * std::pow(HiggsMassesVec.at(k),0.5);
				}
				if( HiggsMassesZeroTempVec[k] > 0) {
					VDebay += -1.5 * HiggsMassesZeroTempVec.at(k+NHiggs)*std::pow(HiggsMassesZeroTempVec[k],0.5);
				}
			}
            for(std::size_t k=0;k<NGauge;k++)
			{
				if(GaugeMassesVec[k] > 0 ) {
					VDebay += 1.5*GaugeMassesVec.at(k+NGauge)*std::pow(GaugeMassesVec[k],0.5);
				}
				if(GaugeMassesZeroTempVec[k] > 0 ) {
					VDebay += -1.5*GaugeMassesZeroTempVec.at(k+NGauge)*std::pow(GaugeMassesZeroTempVec[k],0.5);
				}
			}

			VDebay *= -Temp/(12*M_PI);
			res += VDebay;
		}
	}
	else if(diff == -1){
		if(C_UseParwani){
            for(std::size_t k=0;k<NHiggs;k++) {
				res += HiggsMassesVec.at(k+NHiggs)*boson(HiggsMassesVec[k],Temp,C_CWcbHiggs,0);
				res += boson(HiggsMassesVec[k],Temp,C_CWcbHiggs,-1);
			}
            for(std::size_t k=0;k<NGauge;k++){
				res += GaugeMassesVec.at(k+NGauge)*boson(GaugeMassesVec[k],Temp,C_CWcbGB,0);
				res += boson(GaugeMassesVec.at(k),Temp,C_CWcbGB,-1);
			}
            for(std::size_t k=0;k<NGauge;k++){
				res += 2*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,-1);
			}
            for(std::size_t k=0;k<NQuarks;k++) {
				res += -6*fermion(QuarkMassesVec[k],Temp,-1);
			}
            for(std::size_t k=0;k<NLepton;k++){
				res += -2*fermion(LeptonMassesVec[k],Temp,-1);
			}
		}
		else{
            HiggsMassesZeroTempVec=HiggsMassesSquared(v,0,diff);
            for(std::size_t k=0;k<NHiggs;k++) {
				res += boson(HiggsMassesZeroTempVec[k],Temp,C_CWcbHiggs,-1);
			}
            for(std::size_t k=0;k<NGauge;k++) {
				res += 3*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,-1);
			}
			double AddContQuark=0;
            for(std::size_t k=0;k<NQuarks;k++) AddContQuark += -2*fermion(QuarkMassesVec[k],Temp,-1);
            for(std::size_t k=0;k<NColour;k++) res += AddContQuark;
            for(std::size_t k=0;k<NLepton;k++) res += -2*fermion(LeptonMassesVec[k],Temp,-1);


			double VDebay = 0;
            for(std::size_t k=0;k<NHiggs;k++)
			{
				if( HiggsMassesVec[k] > 0 ) VDebay += std::pow(HiggsMassesVec[k],1.5);
				if( HiggsMassesZeroTempVec[k] > 0) VDebay += -std::pow(HiggsMassesZeroTempVec[k],1.5);
			}
            for(std::size_t k=0;k<NGauge;k++)
			{
				if(GaugeMassesVec[k] > 0 ) VDebay += std::pow(GaugeMassesVec[k],1.5);
				if(GaugeMassesZeroTempVec[k] > 0 ) VDebay += -std::pow(GaugeMassesZeroTempVec[k],1.5);
			}

			VDebay *= -1.0/(12*M_PI);
			res += VDebay;
		}
	}



        return res;




}


void Class_Potential_Origin::CalculateDebye()
{
    if(!SetCurvatureDone) SetCurvatureArrays();

    bool Done = CalculateDebyeSimplified();
    if(!Done)
    {
        for(std::size_t i=0;i<NHiggs;i++)
        {
            for(std::size_t j=i;j<NHiggs;j++)
            {
                DebyeHiggs[i][j] = 0;
                for(std::size_t k=0;k<NHiggs;k++)
                {
                    DebyeHiggs[i][j] += 0.5*Curvature_Higgs_L4[i][j][k][k]/12.0;
                }
                for(std::size_t k=0;k<NGauge;k++)
                {
                    DebyeHiggs[i][j] += 3*0.5*Curvature_Gauge_G2H2[k][k][i][j]/12.0;
                }

                for(std::size_t a=0;a<NQuarks;a++)
                {
                    for(std::size_t b=0;b<NQuarks;b++)
                    {
                        double tmp = 0.5*(std::conj(Curvature_Quark_F2H1[a][b][j])*Curvature_Quark_F2H1[a][b][i] + std::conj(Curvature_Quark_F2H1[a][b][i])*Curvature_Quark_F2H1[a][b][j]).real();
                        DebyeHiggs[i][j] += 6.0/24.0*tmp;
                    }
                }

                for(std::size_t a=0;a<NLepton;a++)
                {
                    for(std::size_t b=0;b<NLepton;b++)
                    {
                        double tmp = 0.5*(std::conj(Curvature_Lepton_F2H1[a][b][j])*Curvature_Lepton_F2H1[a][b][i] + std::conj(Curvature_Lepton_F2H1[a][b][i])*Curvature_Lepton_F2H1[a][b][j]).real();
                        DebyeHiggs[i][j] += 2.0/24.0*tmp;
                    }
                }

                //	            if(i==j) DebyeHiggs[i][j] *= 0.5;
            }
        }

        for(std::size_t i=0;i<NHiggs;i++){
            for(std::size_t j=i;j<NHiggs;j++) {
                if(std::abs(DebyeHiggs[i][j])<=1e-5) DebyeHiggs[i][j] = 0;
            }
        }


        for(std::size_t i=0;i<NHiggs;i++)
        {
            for(std::size_t j=0;j<i;j++) {
                DebyeHiggs[i][j] = DebyeHiggs[j][i];
            }
        }
    }


    std::size_t nHiggsGauge = 0;
    for(std::size_t i=0;i<NHiggs;i++)
    {
        if(Curvature_Gauge_G2H2[0][0][i][i] != 0) nHiggsGauge++;
    }



}


void Class_Potential_Origin::CalculateDebyeGauge(){
    for(std::size_t i=0;i<NGauge;i++){
        for(std::size_t j=0;j<NGauge;j++) DebyeGauge[i][j] = 0;
    }

    bool Done = CalculateDebyeGaugeSimplified();
    if(Done) return ;

    std::size_t nGaugeHiggs = 0;



    for(std::size_t i=0;i<NHiggs;i++) {
	if(Curvature_Gauge_G2H2[0][0][i][i] != 0){
	    nGaugeHiggs++;
	}
    }
    for(std::size_t i=0;i<NGauge;i++)
      {
    	double GaugeFac = 0;
        for(std::size_t k=0;k<NHiggs;k++)
    	{
    		GaugeFac += Curvature_Gauge_G2H2[i][i][k][k];
    	}
    	GaugeFac*= 1.0/nGaugeHiggs;
    	DebyeGauge[i][i] = 2.0/3.0*(nGaugeHiggs/8.0 + 5)*GaugeFac;
      }

    for(std::size_t i=0;i<NGauge;i++){
        for(std::size_t j=0;j<NGauge;j++)
    	{
    		if(std::abs(DebyeGauge[i][j]) <=1e-5) DebyeGauge[i][j] = 0;
    	}
    }
}


void Class_Potential_Origin::initVectors(){
    using vec2 = std::vector<std::vector<double>>;
    using vec3 = std::vector<std::vector<std::vector<double>>>;
    using vec4 = std::vector<std::vector<std::vector<std::vector<double>>>>;

    using vec1Complex = std::vector<std::complex<double>>;
    using vec2Complex = std::vector<std::vector<std::complex<double>>>;
    using vec3Complex = std::vector<std::vector<std::vector<std::complex<double>>>>;
    using vec4Complex = std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>;

    Curvature_Higgs_L1 = std::vector<double>(NHiggs,0);
    Curvature_Higgs_L2 = vec2{NHiggs,std::vector<double>(NHiggs,0)};
    Curvature_Higgs_L3 = vec3{NHiggs,vec2{NHiggs,std::vector<double>(NHiggs,0)}};
    Curvature_Higgs_L4 = vec4{NHiggs,vec3{NHiggs,vec2{NHiggs,std::vector<double>(NHiggs,0)}}};

    Curvature_Higgs_CT_L1 = std::vector<double>(NHiggs,0);
    Curvature_Higgs_CT_L2 = vec2{NHiggs,std::vector<double>(NHiggs,0)};
    Curvature_Higgs_CT_L3 = vec3{NHiggs,vec2{NHiggs,std::vector<double>(NHiggs,0)}};
    Curvature_Higgs_CT_L4 = vec4{NHiggs,vec3{NHiggs,vec2{NHiggs,std::vector<double>(NHiggs,0)}}};

    VEVSymmetric = std::vector<double>(NHiggs,0);

    DebyeHiggs = vec2{NHiggs,std::vector<double>(NHiggs,0)};

    LambdaHiggs_3 = vec3{NHiggs,vec2{NHiggs,std::vector<double>(NHiggs,0)}};
    LambdaHiggs_3_CT = vec3{NHiggs,vec2{NHiggs,std::vector<double>(NHiggs,0)}};

    Curvature_Gauge_G2H2 = vec4{NGauge,vec3{NGauge,vec2{NHiggs,std::vector<double>(NHiggs,0)}}};
    DebyeGauge = vec2{NGauge,std::vector<double>(NGauge,0)};
    LambdaGauge_3 = vec3{NGauge,vec2{NGauge,std::vector<double>(NHiggs,0)}};

    Curvature_Lepton_F2H1 = vec3Complex{NLepton,vec2Complex{NLepton,vec1Complex(NHiggs,0)}};
    LambdaLepton_3 = vec3Complex{NLepton,vec2Complex{NLepton,vec1Complex(NHiggs,0)}};
    LambdaLepton_4 = vec4Complex{NLepton,vec3Complex{NLepton,vec2Complex{NHiggs,vec1Complex(NHiggs,0)}}};

    Curvature_Quark_F2H1 = vec3Complex{NQuarks,vec2Complex{NQuarks,vec1Complex(NHiggs,0)}};
    LambdaQuark_3 = vec3Complex{NQuarks,vec2Complex{NQuarks,vec1Complex(NHiggs,0)}};
    LambdaQuark_4 = vec4Complex{NQuarks,vec3Complex{NQuarks,vec2Complex{NHiggs,vec1Complex(NHiggs,0)}}};

    HiggsVev = std::vector<double>(NHiggs,0);
}



void Class_Potential_Origin::resetbools()
{
  SetCurvatureDone = false;
  CalcCouplingsdone=false;
  CalculatedTripleCopulings = false;
  parStored.clear();
  parCTStored.clear();

}


bool Class_Potential_Origin::CheckNLOVEV(const std::vector<double>& v) const
{
	//std::vector<double> vPotential;
	double MaxDiff = 0;
	double AllowedDifference = 1;
    for(std::size_t i=0;i<nVEV;i++)
	{
		double tmp=std::abs(std::abs(v[i])-std::abs(vevTreeMin[i]));
		if(tmp > MaxDiff) MaxDiff = tmp;
	}

	return (MaxDiff < AllowedDifference);
}



double Class_Potential_Origin::EWSBVEV(const std::vector<double>& v) const
{
  double res=0;
  for(std::size_t i=0;i<NHiggs;i++)
  {
	  double checkgauge=0;
      for(std::size_t j=0;j<NGauge;j++)
	  {
		  checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][i][i]);
	  }
	  if(checkgauge != 0) res += std::pow(v.at(i),2);
  }
  res = std::sqrt(res);

  if(res <= 0.5)
    {
      return 0;
    }
  return res;
}


void Class_Potential_Origin::setUseIndexCol(std::string legend){
    UseIndexCol = legend.rfind(sep,0) == 0;
}

void Class_Potential_Origin::CheckImplementation(
        const int& WhichMinimizer ) const{
	using std::pow;

    const std::string Pass = "Pass";
    const std::string Fail = "Fail";


    std::cout << "The tested Model is the " << Model<< std::endl;

    std::vector<std::string> TestNames, TestResults;

    TestNames.push_back("CT number/label match");
    TestResults.push_back(Pass);
    if(addLegendCT().size() != nParCT){
        std::cerr << "WARNING: The number of labels for the Counterterms does not match the number of Counterterms."
                  << " If you don't fix this, then your header will not match the numerical output in the output file."
                  << std::endl;
        TestResults.at(TestResults.size()-1) = Fail;

    }

    TestNames.push_back("VEV number/label match");
    TestResults.push_back(Pass);
    if(addLegendVEV().size() != nVEV){
        std::cerr << "WARNING: The number of labels for the VEVs does not match the number of VEVs defined in the model."
                  << " If you don't fix this, then your header will not match the numerical output in the output file."
                  << std::endl;
        TestResults.at(TestResults.size()-1) = Fail;
    }

    TestNames.push_back("addLegendTemp number/label match");
    TestResults.push_back(Pass);
    if(addLegendTemp().size() != nVEV + 3){
        std::cerr << "WARNING: The number of labels in addLegendTemp does not match the number of VEVs + 3."
                  << " If you don't fix this, then your header will not match the numerical output in the output file."
                  << " It is expected to be label for the critical temperature, the cirtical VEV, the ratio of VEV and temperature and the labels\
                     for the VEVs."
                  << std::endl;
        TestResults.at(TestResults.size()-1) = Fail;
    }

    TestNames.push_back("addLegendTripleCouplings number/label match");
    TestResults.push_back(Pass);
    {
        std::size_t ExpectedTripleHiggs = 0;
        for(std::size_t i=0;i<NHiggs;i++){
            for(std::size_t j=i;j<NHiggs;j++){
                for(std::size_t k=j;k<NHiggs;k++){
                    ExpectedTripleHiggs++;
                }
            }
        }
        ExpectedTripleHiggs*=3;
        if(addLegendTripleCouplings().size() != ExpectedTripleHiggs){
            std::cerr << "WARNING: The number of labels in addLegendTripleCouplings does not match the number of calculated Triple Higgs Couplings."
                      << " If you don't fix this, then your header will not match the numerical output in the output file."
                      << std::endl;
            TestResults.at(TestResults.size()-1) = Fail;
        }
    }

    TestNames.push_back("CKM matrix unitarity");
    TestResults.push_back(Pass);
	double ZeroMass = std::pow(10,-5);
    {
        // Checks if the implemented CKM Matrix is unitary
        MatrixXcd VCKM(3,3);
        VCKM(0,0) = C_Vud;
        VCKM(0,1) = C_Vus;
        VCKM(0,2) = C_Vub;
        VCKM(1,0) = C_Vcd;
        VCKM(1,1) = C_Vcs;
        VCKM(1,2) = C_Vcb;
        VCKM(2,0) = C_Vtd;
        VCKM(2,1) = C_Vts;
        VCKM(2,2) = C_Vtb;

        auto norm = (VCKM.adjoint()*VCKM - MatrixXcd::Identity(3,3)).norm();
        if(norm > 1e-5){
            std::cerr << "Your CKM implementation is not unitary!" << std::endl;
            std::cerr << "Your CKM Matrix V is given by \n" << VCKM << std::endl;
            std::cerr << "with adjoint(V)*V = \n" << VCKM.adjoint() * VCKM << std::endl;
            std::cerr << "The norm deviating from 1 is " << norm << std::endl;
            TestResults.at(TestResults.size()-1) = Fail;
        }
    }


	std::cout << std::endl
			<< "This function calculates the masses of the gauge bosons, fermions and Higgs boson and compares them "
			<< "with the parameters defined in SMparam.h."
			<< std::endl;


    TestNames.push_back("Matching gauge boson masses with SMparam.h");
    TestResults.push_back(Pass);
    {
        std::vector<double> gaugeMassesInput;
        gaugeMassesInput.push_back(0);
        gaugeMassesInput.push_back(pow(C_MassW,2));
        gaugeMassesInput.push_back(pow(C_MassW,2));
        gaugeMassesInput.push_back(pow(C_MassZ,2));
        std::sort(gaugeMassesInput.begin(),gaugeMassesInput.end());
        auto GaugeMassCalculated=GaugeMassesSquared(vevTree,0,0);


        std::string prsize_tline1="The SM gauge boson masses squared are : ";
        std::string prsize_tline2="The calculated gauge boson masses squared are : ";
        auto maxlength = std::max(prsize_tline1.size(),prsize_tline2.size());
        auto addline1 = maxlength-prsize_tline1.size();
        auto addline2 = maxlength-prsize_tline2.size();

        std::cout << prsize_tline1 << std::setw(addline1) << " ";
        for(auto x: gaugeMassesInput) std::cout << x << sep;
        std::cout << std::endl;
        std::cout << prsize_tline2 << std::setw(addline2) ;
        for(auto x: GaugeMassCalculated) std::cout << x << sep;
        std::cout << std::endl;

        std::cout << std::endl;

        double sum{0};
        for(std::size_t i{0};i<NGauge;++i)
        {
            sum += std::abs(GaugeMassCalculated.at(i) - gaugeMassesInput.at(i));
        }
        if(sum > ZeroMass)
        {
            TestResults.at(TestResults.size()-1) =Fail;
        }



    }

    TestNames.push_back("Matching lepton masses with SMparam.h");
    TestResults.push_back(Pass);
    TestNames.push_back("Matching quark masses with SMparam.h");
    TestResults.push_back(Pass);
    {
        std::vector<double> leptonMassesInput,quarkMassesInput;
        leptonMassesInput.push_back(0);
        leptonMassesInput.push_back(0);
        leptonMassesInput.push_back(0);
        leptonMassesInput.push_back(pow(C_MassElectron,2));
        leptonMassesInput.push_back(pow(-C_MassElectron,2));
        leptonMassesInput.push_back(pow(C_MassMu,2));
        leptonMassesInput.push_back(pow(-C_MassMu,2));
        leptonMassesInput.push_back(pow(C_MassTau,2));
        leptonMassesInput.push_back(pow(-C_MassTau,2));

        quarkMassesInput.push_back(pow(C_MassUp,2));
        quarkMassesInput.push_back(pow(-C_MassUp,2));
        quarkMassesInput.push_back(pow(C_MassCharm,2));
        quarkMassesInput.push_back(pow(-C_MassCharm,2));
        quarkMassesInput.push_back(pow(C_MassTop,2));
        quarkMassesInput.push_back(pow(-C_MassTop,2));

        quarkMassesInput.push_back(pow(C_MassDown,2));
        quarkMassesInput.push_back(pow(-C_MassDown,2));
        quarkMassesInput.push_back(pow(C_MassStrange,2));
        quarkMassesInput.push_back(pow(-C_MassStrange,2));
        quarkMassesInput.push_back(pow(C_MassBottom,2));
        quarkMassesInput.push_back(pow(-C_MassBottom,2));

        if(NLepton == 0){
            for(std::size_t i=0;i<leptonMassesInput.size();i++) quarkMassesInput.push_back(leptonMassesInput.at(i));
        }


        std::sort(leptonMassesInput.begin(),leptonMassesInput.end());
        std::sort(quarkMassesInput.begin(),quarkMassesInput.end());


        auto leptonMassCalculated=LeptonMassesSquared(vevTree,0);
        auto quarkMassCalculated=QuarkMassesSquared(vevTree,0);

        if(NLepton != 0){
            std::string prsize_tline1="The SM lepton masses squared are : ";
            std::string prsize_tline2="The calculated lepton masses squared are : ";
            auto maxlength = std::max(prsize_tline1.size(),prsize_tline2.size());
            auto addline1 = maxlength-prsize_tline1.size();
            auto addline2 = maxlength-prsize_tline2.size();

            std::cout << prsize_tline1 << std::setw(addline1) << " ";
            for(auto x: leptonMassesInput) std::cout << x << sep;
            std::cout << std::endl;
            std::cout << prsize_tline2 << std::setw(addline2) ;
            for(auto x: leptonMassCalculated) std::cout << x << sep;
            std::cout << std::endl;
            std::cout << std::endl;

            double sum{0};
            for(std::size_t i{0};i<NLepton;++i)
            {
                sum += std::abs(leptonMassCalculated.at(i) - leptonMassesInput.at(i));
            }
            if(sum > ZeroMass)
            {
                TestResults.at(TestResults.size()-2) = Fail;
            }
        }

        std::string prsize_tline1="The SM quark masses squared are : ";
        std::string prsize_tline2="The calculated quark masses squared are : ";
        auto maxlength = std::max(prsize_tline1.size(),prsize_tline2.size());
        auto addline1 = maxlength-prsize_tline1.size();
        auto addline2 = maxlength-prsize_tline2.size();

        std::cout << prsize_tline1 << std::setw(addline1) << " ";
        for(auto x: quarkMassesInput) std::cout << x << sep;
        std::cout << std::endl;
        std::cout << prsize_tline2 << std::setw(addline2) ;
        for(auto x: quarkMassCalculated) std::cout << x << sep;
        std::cout << std::endl << std::endl;

        double sum{0};
        for(std::size_t i{0};i<NQuarks;++i)
        {
            sum += std::abs(quarkMassCalculated.at(i) - quarkMassesInput.at(i));
        }
        if(sum > ZeroMass)
        {
            TestResults.at(TestResults.size()-1) = Fail;
        }

    }


    TestNames.push_back("Correct EW Minimum");
    TestResults.push_back(Pass);
    {
        std::vector<double> CalculatedHiggsVEV, CheckVector, start;
        for(const auto& x: vevTreeMin) start.push_back(0.5*x);
	std::cout << "parStored: " << parStored << std::endl;
	std::cout << "parCTStored: " << parCTStored << std::endl;
	std::cout << "CheckVector: " << CheckVector << std::endl;
	std::cout << "start: " << start << std::endl;
	std::cout << "Minimizer: " << WhichMinimizer << std::endl;

        CalculatedHiggsVEV = Minimizer::Minimize_gen_all_tree_level(Model,parStored,parCTStored,CheckVector,start,WhichMinimizer);

        std::string prsize_tline1="The given VEV configuration at tree-level is : ";
        std::string prsize_tline2="The calculated VEV configuration at tree-level is : ";
        auto maxlength = std::max(prsize_tline1.size(),prsize_tline2.size());
        auto addline1 = maxlength-prsize_tline1.size();
        auto addline2 = maxlength-prsize_tline2.size();

        std::cout << prsize_tline1 << std::setw(addline1) << " ";
        for(auto x: vevTreeMin) std::cout << x << sep;
        std::cout << std::endl;
        std::cout << prsize_tline2 << std::setw(addline2) ;
        for(auto x: CalculatedHiggsVEV) std::cout << x << sep;
        std::cout << std::endl << std::endl;

        double sum{0};
        for(std::size_t i{0};i<nVEV;++i)
        {
            sum += std::abs(std::abs(CalculatedHiggsVEV.at(i)) - std::abs(vevTreeMin.at(i)));
        }
        if(sum > 0.5)
        {
            TestResults.at(TestResults.size()-1) = Fail;
        }
    }


    TestNames.push_back("Tadpole relations fullfilled");
    TestResults.push_back(Pass);
    {
        double SurviveTadpole = 0;
        for(std::size_t i=0;i<NHiggs;i++) SurviveTadpole += std::abs(VTree(vevTree,i+1));

        if(SurviveTadpole > 1e-5){
            std::cout << "The given input parameter does not fulfill the tadpole relations and is not a minimum of the potential."
                      << std::endl
                      << "This may happen if all your parameters are read in from an input file. Try applying the minimum conditions"
                      <<	" in the set_gen function."
                         << std::endl
                         << std::endl;
            TestResults.at(TestResults.size()-1) = Fail;
        }
    }



    TestNames.push_back("Matching Masses between NLO and tree-level");
    TestResults.push_back(Pass);
    {
        std::vector<double> WeinbergNabla,WeinbergHesse;
        WeinbergNabla = WeinbergFirstDerivative();
        WeinbergHesse = WeinbergSecondDerivative();

        VectorXd NablaWeinberg(NHiggs),NablaVCT(NHiggs);
        MatrixXd HesseWeinberg(NHiggs,NHiggs), HesseVCT(NHiggs,NHiggs);
        for(std::size_t i=0;i<NHiggs;i++)
        {
            NablaWeinberg[i] = WeinbergNabla[i];
            for(std::size_t j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);
        }

        for(std::size_t i=0;i<NHiggs;i++){
            NablaVCT(i) = Curvature_Higgs_CT_L1[i];
            for(std::size_t j=0;j<NHiggs;j++){
                NablaVCT(i) += Curvature_Higgs_CT_L2[i][j] * vevTree[j];
                HesseVCT(i,j) = Curvature_Higgs_CT_L2[i][j];
                for(std::size_t k=0;k<NHiggs;k++){
                    NablaVCT(i) += 0.5*Curvature_Higgs_CT_L3[i][j][k] * vevTree[j] *vevTree[k];
                    HesseVCT(i,j) += Curvature_Higgs_CT_L3[i][j][k] * vevTree[k];
                    for(std::size_t l=0;l<NHiggs;l++){
                        NablaVCT(i) += 1.0/6.0 * Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[j] * vevTree[k] * vevTree[l];
                        HesseVCT(i,j) += 0.5* Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k] * vevTree[l];
                    }
                }
            }
        }

        MatrixXd MassMatrix(NHiggs,NHiggs);
        for(std::size_t i=0;i<NHiggs;i++)
        {
            for(std::size_t j=0;j<NHiggs;j++)
            {
                MassMatrix(i,j) = Curvature_Higgs_L2[i][j];
                for(std::size_t k=0;k<NHiggs;k++)
                {
                    MassMatrix(i,j) += Curvature_Higgs_L3[i][j][k]*vevTree[k];
                    for(std::size_t l=0;l<NHiggs;l++) MassMatrix(i,j) += 0.5*Curvature_Higgs_L4[i][j][k][l] * vevTree[k]*vevTree[l];
                }
            }
        }

        SelfAdjointEigenSolver<MatrixXd> esTree(MassMatrix,EigenvaluesOnly);
        SelfAdjointEigenSolver<MatrixXd> esNLO(MassMatrix+HesseVCT+HesseWeinberg,EigenvaluesOnly);

        std::vector<double> TreeMass,NLOMass;

        for(std::size_t i =0;i<NHiggs;i++){
            if(std::abs(esTree.eigenvalues()[i]) < ZeroMass) TreeMass.push_back(0);
            else TreeMass.push_back(esTree.eigenvalues()[i]);
            if(std::abs(esNLO.eigenvalues()[i]) < ZeroMass) NLOMass.push_back(0);
            else NLOMass.push_back(esNLO.eigenvalues()[i]);
        }

        std::cout << "The higgs masses squared at LO | NLO are : " << std::endl;
        for(std::size_t i=0;i<NHiggs;i++){
            std::cout << "m_i^2 = " << TreeMass[i] << " | " << NLOMass[i] << std::endl;
        }

        double sum{0.0};
        for(std::size_t i=0;i<NHiggs;i++){
            double z = std::abs(std::abs(NLOMass[i]) - TreeMass.at(i));
            double n = std::max(std::abs(NLOMass[i]), std::abs(TreeMass[i]));
            if(n != 0) sum += z/n;
        }
        if(sum > 0.5) {
            TestResults.at(TestResults.size()-1) = Fail;
        }
    }


	if(UseVTreeSimplified){
        TestNames.push_back("Checking VTreeSimplified");
        TestResults.push_back(Pass);
		std::default_random_engine randGen(0);
		double ValSimplified = 0, ValTensor = 0, PotentialDifference = 0;
        std::vector<double> VevDummyConfigForMinimiser(nVEV);
		double RNDMax = 500;
		for(int i=0;i<10;i++){
            for(std::size_t j=0;j<nVEV;j++) VevDummyConfigForMinimiser.at(j) = RNDMax*(-1+2*std::generate_canonical<double,std::numeric_limits<double>::digits>(randGen));
            auto VevDummyConfigForPotential=MinimizeOrderVEV(VevDummyConfigForMinimiser);
			ValSimplified = VTree(VevDummyConfigForPotential);
			ValTensor = VTree(VevDummyConfigForPotential,0,true /* = ForceExplicitCalculation */);
			PotentialDifference += std::abs(ValTensor - ValSimplified)/(std::abs(ValTensor) + std::abs(ValSimplified));
		}
		if(PotentialDifference > 1e-5) {
			std::cout << "You provided a simplified version of the tree-level potential but it yields"
					<< " different results for the same input compared to the explicit calculation. "
					<< "Recheck your implementation of the simplified tree-level potential."
					<< std::endl;
            TestResults.at(TestResults.size()-1) = Fail;
		}

	}

	if(UseVCounterSimplified){
        TestNames.push_back("Checking VCounterSimplified");
        TestResults.push_back(Pass);
		std::default_random_engine randGen(0);
		double ValCTSimplified = 0, ValCTTensor = 0, PotentialCTDifference = 0;
        std::vector<double> VevDummyConfigForMinimiser(nVEV);
		double RNDMax = 500;
		for(int i=0;i<10;i++){
            for(std::size_t j=0;j<nVEV;j++) VevDummyConfigForMinimiser.at(j) = RNDMax*(-1+2*std::generate_canonical<double,std::numeric_limits<double>::digits>(randGen));
            auto VevDummyConfigForPotential=MinimizeOrderVEV(VevDummyConfigForMinimiser);
			ValCTSimplified = CounterTerm(VevDummyConfigForPotential);
			ValCTTensor = CounterTerm(VevDummyConfigForPotential,0,	true);
			PotentialCTDifference += std::abs(ValCTTensor - ValCTSimplified)/(std::abs(ValCTSimplified)+std::abs(ValCTTensor));
		}
		if(PotentialCTDifference > 1) {
			std::cout << "You provided a simplified version of the counterterm potential but it yields"
					<< " different results for the same input compared to the explicit calculation. "
					<< "Recheck your implementation of the simplified counterterm potential."
					<< std::endl;
            TestResults.at(TestResults.size()-1) = Fail;
		}

	}



    if(TestNames.size() != TestResults.size())
    {
        std::cout << "TestNames : " << std::endl << TestNames << std::endl;
        std::cout << "TestResults : " << std::endl << TestResults << std::endl;
        std::string errmsg{"Mismatch between the number of labels for the tests and the results."};
        errmsg += "TestNames.size() = " + std::to_string(TestNames.size());
        errmsg +=" TestResults.size() = " + std::to_string(TestResults.size());
        throw std::runtime_error(errmsg);
    }

    std::size_t Passes{0};
    for(const auto& el: TestResults)
    {
        if(el == Pass) Passes++;
    }

    std::cout << "\nTEST | Pass/Fail\n================================\n" << std::endl;
    auto maxsize = TestNames.at(0).size();
    for(const auto& el: TestNames)
    {
        if(el.size() > maxsize) maxsize = el.size();
    }
    maxsize += 5;
    for(std::size_t i{0};i<TestNames.size();++i)
    {
        std::cout << std::setw(maxsize) << std::left << TestNames.at(i)
                  << "| " << TestResults.at(i) << std::endl;
    }

    std::cout << Passes << " tests out of " << TestResults.size() << " passed.\n" << std::endl;
    if(Passes != TestResults.size())
    {
        std::cout << TestResults.size()-Passes << " tests failed. Please check and try again." << std::endl;
    }
    else{
        std::cout << "You're good to go!\n" << std::endl;
    }








}


void Class_Potential_Origin::FindSignSymmetries(){
	SignSymmetries.clear();
	std::vector<double> testvev,testvevPotential;
    for(std::size_t i=0;i<nVEV;i++) testvev.push_back(i+1);
    testvevPotential=MinimizeOrderVEV(testvev);
	double referenceValue = VEff(testvevPotential);

	std::vector<double> vevdummy,vevdummyPotential;
	double VEffDummy;

	// Fill a dummy vector with a certain amout of -1 and look at all possible permutations of it
    for(std::size_t countNegative=1;countNegative<=nVEV;countNegative++){
		std::vector<double> tmpSymmetry;
        for(std::size_t i=0;i<countNegative;i++) tmpSymmetry.push_back(-1);
        for(std::size_t i=countNegative;i<nVEV;i++) tmpSymmetry.push_back(1);

		do{
			vevdummy.clear();
            for(std::size_t i=0;i<nVEV;i++) vevdummy.push_back(tmpSymmetry.at(i)*testvev.at(i));
            vevdummyPotential=MinimizeOrderVEV(vevdummy);
			VEffDummy = VEff(vevdummyPotential);
			if(std::abs(VEffDummy-referenceValue)<=1e-3*std::abs(referenceValue)) SignSymmetries.push_back(tmpSymmetry);
		}while(std::next_permutation(tmpSymmetry.begin(),tmpSymmetry.end()));
	}
}

void Class_Potential_Origin::SetUseTreeLevel(bool val){
	UseTreeLevel=val;
}

std::pair<std::vector<double>, std::vector<double>> Class_Potential_Origin::initModel(std::string linestr){
	std::vector<double> par(nPar), parCT(nParCT);
	resetbools();
	ReadAndSet(linestr,par);
    parCT=initModel(par);

	parStored=par;
	parCTStored=parCT;

	std::pair<std::vector<double>, std::vector<double>> res;
	res.first = par;
	res.second = parCT;

	return res;

}

std::vector<double> Class_Potential_Origin::initModel(const std::vector<double>& par){
	std::vector<double> parCT(nParCT);
	resetbools();
	set_gen(par);
    CalculatePhysicalCouplings();
    parCT=calc_CT();
	set_CT_Pot_Par(parCT);
	CalculateDebye();
	CalculateDebyeGauge();

    parStored = par;
    parCTStored = parCT;

	return parCT;
}

std::vector<double> Class_Potential_Origin::resetScale(const double& newScale)
{
    scale = newScale;
    auto parCT = calc_CT();
    set_CT_Pot_Par(parCT);

    parCTStored = parCT;

    return parCT;
}


std::vector<std::complex<double>> Class_Potential_Origin::QuarkMasses(const std::vector<double>& v) const
{
    std::vector<std::complex<double>> res;
	if(v.size() != nVEV and v.size() != NHiggs){
	  std::string ErrorString = std::string("You have called ") + std::string(__func__)
		+ std::string(" with an invalid vev configuration. Your vev is of dimension ")
	  + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
	  throw std::runtime_error(ErrorString);
	}
	if(v.size() == nVEV and nVEV != NHiggs){
		std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
		    			<< nVEV << " while it should have " << NHiggs
						<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
	  std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      res = QuarkMasses(Transformedv);
      return res ;
	}
    if(!SetCurvatureDone) {
//        SetCurvatureArrays(); std::cout << "Reset of Set Curvature " << std::endl;
        std::string retmes = __func__;
        retmes += " is called before SetCurvatureArrays() is called. \n";
        throw std::runtime_error(retmes);
    }
    MatrixXcd MIJ(NQuarks,NQuarks);
    MIJ = MatrixXcd::Zero(NQuarks,NQuarks);
    double ZeroMass = std::pow(10,-10);
    for(std::size_t i=0;i<NQuarks;i++)
    {
        for(std::size_t j=0;j<NQuarks;j++)
        {

            for(std::size_t k=0;k<NHiggs;k++)
            {
                MIJ(i,j) += Curvature_Quark_F2H1[i][j][k]*v[k];
            }
        }
    }


	ComplexEigenSolver<MatrixXcd> es(MIJ,false);

    for(std::size_t i =0;i<NQuarks;i++)
	{
		auto tmp = es.eigenvalues()[i];
		if(std::abs(tmp) < ZeroMass ) res.push_back(0);
		else res.push_back(tmp);
	}

    return res;
}


std::vector<std::complex<double>> Class_Potential_Origin::LeptonMasses(const std::vector<double>& v) const
{
    std::vector<std::complex<double>> res;
	if(v.size() != nVEV and v.size() != NHiggs){
	  std::string ErrorString = std::string("You have called ") + std::string(__func__)
		+ std::string(" with an invalid vev configuration. Your vev is of dimension ")
	  + std::to_string(v.size()) + std::string(" and it should be ") + std::to_string(NHiggs) + std::string(".");
	  throw std::runtime_error(ErrorString);
	}
	if(v.size() == nVEV and nVEV != NHiggs){
		std::cerr << __func__ << " is being called with a wrong sized vev configuration. It has the dimension of "
		    			<< nVEV << " while it should have " << NHiggs
						<< ". For now this is transformed but please fix this to reduce the runtime." << std::endl;
	  std::vector<double> Transformedv;
      Transformedv=MinimizeOrderVEV(v);
      res = LeptonMasses(Transformedv);
      return res ;
	}
    if(!SetCurvatureDone) {
//        SetCurvatureArrays(); std::cout << "Reset of Set Curvature " << std::endl;
        std::string retmes = __func__;
        retmes += " is called before SetCurvatureArrays();\n";
        throw std::runtime_error(retmes);
    };
    MatrixXcd MIJ(NLepton,NLepton);
    MIJ = MatrixXcd::Zero(NLepton,NLepton);
    double ZeroMass = std::pow(10,-10);
    for(std::size_t i=0;i<NLepton;i++)
    {
        for(std::size_t j=0;j<NLepton;j++)
        {

            for(std::size_t k=0;k<NHiggs;k++)
            {
                MIJ(i,j) += Curvature_Lepton_F2H1[i][j][k]*v[k];
            }
        }
    }


	ComplexEigenSolver<MatrixXcd> es(MIJ,false);

    for(std::size_t i =0;i<NLepton;i++)
	{
		auto tmp = es.eigenvalues()[i];
		if(std::abs(tmp) < ZeroMass ) res.push_back(0);
		else res.push_back(tmp);
	}
    return  res;
}

double Class_Potential_Origin::CalculateRatioAlpha(
		const std::vector<double>& vev_symmetric,
		const std::vector<double>& vev_broken,
        const double& Temp) const{
	double res = 0;
	double PotentialSymmetricPhaseValue = VEff(vev_symmetric,Temp,0);
	double PotentialSymmetricPhaseDeriv = VEff(vev_symmetric,Temp,-1);
	double PotentialBrokenPhaseValue = VEff(vev_broken,Temp,0);
	double PotentialBrokenPhaseDeriv = VEff(vev_broken,Temp,-1);
	res = -(PotentialBrokenPhaseValue-PotentialSymmetricPhaseValue) + Temp*(PotentialBrokenPhaseDeriv-PotentialSymmetricPhaseDeriv);
	//TODO:: Unfinished!
    throw std::runtime_error("The CalculateRatioAlpha function is not finished and you should not be using it!");

	return res;
}


std::vector<double> Class_Potential_Origin::MinimizeOrderVEV(const std::vector<double>& vevMinimizer) const{
     std::vector<double> vevFunction;


     std::size_t count=0;
     for(std::size_t i=0;i<NHiggs;++i)
     {
         if(i==VevOrder[count]) {
             vevFunction.push_back(vevMinimizer.at(count));
             count++;
         }
         else vevFunction.push_back(0);
     }
     return  vevFunction;
}

}

