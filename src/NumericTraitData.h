/*
 *  NumericTraitData.h
 *  SPLiTTree
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of SPLiTTree: a generic C++ library for Serial and Parallel
 * Lineage Traversal of Trees.
 *
 * SPLiTTree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * SPLiTTree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SPLiTTree.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */

#ifndef NumericTraitData_H_
#define NumericTraitData_H_


#include "splittree.h"

namespace splittree {
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
