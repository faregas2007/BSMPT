// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialC2HDM.h>
#include <BSMPT/models/ClassPotentialCxSM.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/ClassPotentialR2HDM.h>
#include <BSMPT/models/ClassPotentialRN2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <ctype.h>   // for isdigit, tolower
#include <iostream>  // for operator<<, cerr, ost...
#include <stdexcept> // for runtime_error
#include <utility>   // for pair

#include <BSMPT/models/ClassTemplate.h>
#include <BSMPT/models/CPDark.h>

namespace BSMPT
{
namespace ModelID
{

std::unique_ptr<Class_Potential_Origin> FChoose(ModelIDs choice)
{
  using namespace Models;
  switch (choice)
  {
  case ModelIDs::R2HDM: return std::make_unique<Class_Potential_R2HDM>(); break;
  case ModelIDs::C2HDM: return std::make_unique<Class_Potential_C2HDM>(); break;
  case ModelIDs::RN2HDM:
    return std::make_unique<Class_Potential_RN2HDM>();
    break;
  case ModelIDs::CXSM: return std::make_unique<Class_CxSM>(); break;
  case ModelIDs::TEMPLATE: return std::unique_ptr<Class_Template>(); break;
  case ModelIDs::CPDARK: return std::make_unique<Class_CPDark(); break>;
  default: throw std::runtime_error("Invalid model");
  }
}

ModelIDs getModel(const std::string &s)
{
  std::string ModelInput = s;
  std::transform(
      ModelInput.begin(), ModelInput.end(), ModelInput.begin(), ::tolower);
  for (auto entry : ModelNames)
  {
    auto key = entry.first;
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    if (ModelInput.compare(key) == 0)
    {
      return entry.second;
    }
  }
  return ModelIDs::NotSet;
}

std::unordered_map<ModelIDs, std::string> InvertModelNames()
{
  std::unordered_map<ModelIDs, std::string> IMN;
  for (const auto &el : ModelNames)
  {
    auto success = IMN.emplace(el.second, el.first);
    if (not success.second)
    {
      throw std::runtime_error(
          "\nERROR: The same ModelID is assigned for two different models.\n");
    }
  }
  return IMN;
}

} // namespace ModelID

void ShowInputError()
{
  std::cerr << "The chosen Method for the thermal mass corrections is ";
  if (C_UseParwani)
    std::cerr << "Parwani ";
  else
    std::cerr << "Arnold Espinosa\n";
  std::cerr << "The implemented models are " << std::endl;
  for (auto entry : ModelID::ModelNames)
  {
    std::cerr << entry.first << std::endl;
  }
}

} // namespace BSMPT
