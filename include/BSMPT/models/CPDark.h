/*
 * ClassTemplate.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete MÃ¼hlleitner

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

/**
  * @file
  */

#pragma once

#include <string>                               // for string
#include <vector>                               // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT{
namespace Models{

/**
 * @brief The Class_Template class
 * Template for implementing a new model
 */
class Class_CPDark : public Class_Potential_Origin
{
public:
  Class_CPDark ();
  virtual
  ~Class_CPDark ();


  // Add here your parameters for the Lagrangian as well as for the counterterm potential
  // Add here your variables in which you will save the Debye correction factors

  double ms11=0., ms22=0., mss=0., L1=0., L2=0., L3=0., L4=0., L5=0., L6=0., L7=0., L8=0., Areal=0., Aimag=0., vh=0.;

  double ms11CT=0., ms22CT=0., mssCT=0., L1CT=0., L2CT=0., L3CT=0., L4CT=0., L5CT=0., L6CT=0., L7CT=0., L8CT=0., ArealCT=0., AimagCT=0.;
  double T1CT =0., T2CT=0., TCBCT=0., TCPCT=0., TSCT=0.;




  void ReadAndSet(const std::string& linestr, std::vector<double>& par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double>& par) override;
  void set_CT_Pot_Par(const std::vector<double>& par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double>& v) const override;
  double VCounterSimplified(const std::vector<double>& v) const override;
  void Debugging(const std::vector<double>& input, std::vector<double>& output) const override;
};

}
}
