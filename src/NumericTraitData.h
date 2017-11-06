// Copyright 2017 Venelin Mitov
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
//   limitations under the License.

#ifndef NumericTraitData_H_
#define NumericTraitData_H_


#include "ParallelPruningAlgorithm.h"
namespace ppa {
template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of long vectors
  std::vector<NameType> const& names_;
  vec const& z_;
  vec const& se_;
  NumericTraitData(
    std::vector<NameType> const& names,
    vec const& z, vec const& se): names_(names), z_(z), se_(se) {}
};
}

#endif //#ifndef NumericTraitData_H_
