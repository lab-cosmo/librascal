/**
 * @file   rascal_utility.cc
 *
 * Copyright  2019 COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "rascal_utility.hh"
#include <cstdlib>

#ifdef __GNUG__
#include <cxxabi.h>
#endif

std::string rascal::internal::type_name_demangled(const char * name) {
#ifdef __GNUG__
  int status{0};
  auto demangled = abi::__cxa_demangle(name, nullptr, nullptr, &status);
  auto result = std::string(demangled);
  std::free(demangled);
  return result;
#else
#error "no demangling implementation for this compiler"
#endif
}

void rascal::internal::replace(std::string & string,
                               const std::string & pattern,
                               const std::string & replacement) {
  auto pos = string.find(pattern);
  while (pos != std::string::npos) {
    string.replace(pos, pattern.size(), replacement);
    pos = string.find(pattern, pos + replacement.size());
  }
}
