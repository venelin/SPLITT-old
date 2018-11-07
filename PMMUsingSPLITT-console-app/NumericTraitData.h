/**
 *  AbcPMM.h
 *  SPLITT
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of SPLITT: a generic C++ library for Serial and Parallel
 * Lineage Traversal of Trees.
 *
 * SPLITT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * SPLITT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SPLITT.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */
#ifndef NumericTraitData_H_
#define NumericTraitData_H_

#include "./SPLITT.h"
#include <iostream>

namespace PMMUsingSPLITT {

using namespace SPLITT;

template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of long vectors
  std::vector<NameType> const& names_;
  vec const& x_;
  NumericTraitData(
    std::vector<NameType> const& names,
    vec const& x): names_(names), x_(x) {}
};
}
#endif //NumericTraitData_H_
