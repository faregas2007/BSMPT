// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch.hpp"

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>

#include "GenerateTestCompares/C2HDM.h"
#include <fstream>

const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                              /* lambda_2 = */ 0.274365,
                                              /* lambda_3 = */ 4.71019,
                                              /* lambda_4 = */ -2.23056,
                                              /* Re(lambda_5) = */ -2.43487,
                                              /* Im(lambda_5) = */ 0.124948,
                                              /* Re(m_{12}^2) = */ 2706.86,
                                              /* tan(beta) = */ 4.64487,
                                              /* Yukawa Type = */ 1};

const Compare_C2HDM Expected;

TEST_CASE("Checking NLOVEV for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  std::vector<double> Check;
  auto sol = Minimizer::Minimize_gen_all(modelPointer,
                                         0,
                                         Check,
                                         modelPointer->get_vevTreeMin(),
                                         Minimizer::WhichMinimizerDefault);
  for (std::size_t i{0}; i < sol.size(); ++i)
  {
    auto expected = std::abs(modelPointer->get_vevTreeMin(i));
    auto res      = std::abs(sol.at(i));
    REQUIRE(std::abs(res - expected) <= 1e-4);
  }
}

TEST_CASE("Checking EWPT for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  std::vector<double> Check;
  auto EWPT = Minimizer::PTFinder_gen_all(
      modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
  const double omega_c_expected =
      Expected.EWPTPerSetting.at(Minimizer::WhichMinimizerDefault).vc;
  const double Tc_expected =
      Expected.EWPTPerSetting.at(Minimizer::WhichMinimizerDefault).Tc;
  const std::vector<double> min_expected =
      Expected.EWPTPerSetting.at(Minimizer::WhichMinimizerDefault).EWMinimum;
  REQUIRE(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS);
  REQUIRE(std::abs(omega_c_expected - EWPT.vc) / omega_c_expected <= 1e-4);
  REQUIRE(std::abs(Tc_expected - EWPT.Tc) / Tc_expected <= 1e-4);
  for (std::size_t i{0}; i < EWPT.EWMinimum.size(); ++i)
  {
    auto res      = std::abs(EWPT.EWMinimum.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected != 0)
    {
      UNSCOPED_INFO("Current Option for Minimizer:\t"
                    << Minimizer::WhichMinimizerDefault);
      UNSCOPED_INFO("This ist the position:"
                    << i << "\tFound solution =" << EWPT.EWMinimum.at(i)
                    << "\tExpected solution = " << min_expected.at(i));
      REQUIRE(std::abs(res - expected) / expected <= 1e-4);
    }
    else
    {
      REQUIRE(res <= 1e-4);
    }
  }
}

TEST_CASE("Checking number of CT parameters for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckNumberOfCTParameters(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of VEV labels for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckNumberOfVEVLabels(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE(
    "Checking number of labels for temperature dependend results for C2HDM",
    "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckLegendTemp(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of triple Higgs couplings for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckNumberOfTripleCouplings(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking Gauge Boson masses for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckGaugeBosonMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking fermion and quark masses masses for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckFermionicMasses(*modelPointer);
  REQUIRE(result.first == ModelTests::TestResults::Pass);
  REQUIRE(result.second == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level minimum for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckTreeLevelMin(*modelPointer,
                                              Minimizer::WhichMinimizerDefault);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level tadpoles for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckTadpoleRelations(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking NLO masses matching tree level masses for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckNLOMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking VTreeSimplified for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  if (modelPointer->UseVTreeSimplified)
  {
    modelPointer->initModel(example_point_C2HDM);
    auto result = ModelTests::CheckVTreeSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking VCounterSimplified for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  if (modelPointer->UseVCounterSimplified)
  {
    modelPointer->initModel(example_point_C2HDM);
    auto result = ModelTests::CheckVCounterSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking first derivative of the sum of CT and CW in the C2HDM",
          "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckCTConditionsFirstDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking second derivative of the sum of CT and CW in the C2HDM",
          "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckCTConditionsSecondDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking the identities required to vanish for the CT in the C2HDM",
          "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  auto result = ModelTests::CheckCTIdentities(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking triple higgs NLO couplings in the C2HDM", "[c2hdm]")
{

  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  modelPointer->Prepare_Triple();
  modelPointer->TripleHiggsCouplings();

  auto Check = [](auto result, auto expected) {
    if (expected != 0)
    {
      REQUIRE(std::abs(result - expected) / std::abs(expected) < 1e-4);
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  auto NHiggs = modelPointer->get_NHiggs();
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      for (std::size_t k{0}; k < NHiggs; ++k)
      {
        // INFO("Current Minimizer
        // Set-up:\t"<<Minimizer::WhichMinimizerDefault); INFO("Failed at
        // Tree-Coupling ("
        //      << i << "," << j << "," << k << ")\tFound:\t"
        //      << modelPointer->get_TripleHiggsCorrectionsTreePhysical(i, j, k)
        //      << "\tExpextec:\t" <<
        //      Expected.CheckTripleTree.at(i).at(j).at(k));
        Check(modelPointer->get_TripleHiggsCorrectionsTreePhysical(i, j, k),
              Expected.CheckTripleTree.at(i).at(j).at(k));
        // INFO("Failed at CT-Coupling ("
        //      << i << "," << j << "," << k << ")\tFound:\t"
        //      << modelPointer->get_TripleHiggsCorrectionsCTPhysical(i, j, k)
        //      << "\tExpextec:\t" << Expected.CheckTripleCT.at(i).at(j).at(k));
        Check(modelPointer->get_TripleHiggsCorrectionsCTPhysical(i, j, k),
              Expected.CheckTripleCT.at(i).at(j).at(k));
        // INFO("Failed at CW-Coupling ("
        //      << i << "," << j << "," << k << ")\tFound:\t"
        //      << modelPointer->get_TripleHiggsCorrectionsCWPhysical(i, j, k)
        //      << "\tExpextec:\t" << Expected.CheckTripleCW.at(i).at(j).at(k));
        Check(modelPointer->get_TripleHiggsCorrectionsCWPhysical(i, j, k),
              Expected.CheckTripleCW.at(i).at(j).at(k));
      }
    }
  }
}

TEST_CASE("Check number of calculated CT parameters in the C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckCTNumber(*modelPointer));
}

// TEST_CASE("Create EWPT compares", "[c2hdm-setup]")
//{
//  using namespace BSMPT;
//  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
//      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
//  modelPointer->initModel(example_point_C2HDM);

//  int WhichMin;
//  Minimizer::EWPTReturnType EWPT;

//  std::ofstream outfile("dummy.txt");

//  std::map<int, Minimizer::EWPTReturnType> mdata;
//  for (bool UseGSL : {false, true})
//  {
//    for (bool UseCMAES : {false, true})
//    {
//      for (bool UseNLopt : {false, true})
//      {
//        if (not UseGSL and not UseCMAES and not UseNLopt) continue;

//        WhichMin = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
//        EWPT     = Minimizer::PTFinder_gen_all(modelPointer, 0, 300,
//        WhichMin); mdata[WhichMin] = EWPT; outfile << "EWPTPerSetting[" <<
//        WhichMin
//                << "].Tc = " << mdata[WhichMin].Tc << ";" << std::endl
//                << "EWPTPerSetting[" << WhichMin
//                << "].vc = " << mdata[WhichMin].vc << ";" << std::endl;
//        for (const auto &el : EWPT.EWMinimum)
//        {
//          if (std::abs(el) > 1e-5)
//            outfile << "EWPTPerSetting[" << WhichMin <<
//            "].EWMinimum.push_back("
//                    << el << ");" << std::endl;
//          else
//            outfile << "EWPTPerSetting[" << WhichMin <<
//            "].EWMinimum.push_back("
//                    << 0 << ");" << std::endl;
//        }
//      }
//    }
//  }
//}
