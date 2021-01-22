/*
 *  borrowr: estimate population average treatment effects with borrowing between data sources.
 *  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
 *  This code is a modified version from the BART R package from April 2019.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
	std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);

#endif
