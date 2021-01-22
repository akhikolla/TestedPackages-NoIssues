// Copyright (c) 2018-2020  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include "spatVector.h"

#ifdef useGDAL
	#include "crs.h"
#endif


SpatHole::SpatHole() {}

SpatHole::SpatHole(std::vector<double> X, std::vector<double> Y) {
	x = X; y = Y;
	extent.xmin = *std::min_element(X.begin(), X.end());
	extent.xmax = *std::max_element(X.begin(), X.end());
	extent.ymin = *std::min_element(Y.begin(), Y.end());
	extent.ymax = *std::max_element(Y.begin(), Y.end());
}

bool SpatPart::addHole(std::vector<double> X, std::vector<double> Y) {
	SpatHole h(X, Y);
	holes.push_back(h);
	// check if inside pol?
	return true;
}


bool SpatPart::addHole(SpatHole h) {
	holes.push_back(h);
	// check if inside pol?
	return true;
}


SpatPart::SpatPart() {}

SpatPart::SpatPart(double X, double Y) {
	x.push_back(X);
	y.push_back(Y);
	extent.xmin = X;
	extent.xmax = X;
	extent.ymin = Y;
	extent.ymax = Y;
}

SpatPart::SpatPart(std::vector<double> X, std::vector<double> Y) {
	x = X; y = Y;
	extent.xmin = *std::min_element(X.begin(), X.end());
	extent.xmax = *std::max_element(X.begin(), X.end());
	extent.ymin = *std::min_element(Y.begin(), Y.end());
	extent.ymax = *std::max_element(Y.begin(), Y.end());
}


SpatGeom::SpatGeom() {}

SpatGeom::SpatGeom(SpatPart p) {
	parts.push_back(p);
	extent = p.extent;
}

SpatGeom::SpatGeom(SpatGeomType g) {
	gtype = g;
}

bool SpatGeom::unite(SpatGeom g) {
	if (parts.size() == 0) {
		parts = g.parts;
		extent = g.extent;
	} else {
		parts.insert(parts.end(), g.parts.begin(), g.parts.end());
		extent.unite(g.extent);
	}
	return true;
}


bool SpatGeom::addPart(SpatPart p) {
	parts.push_back(p);
	if (parts.size() > 1) {
		extent.unite(p.extent);
	} else {
		extent = p.extent;
	}	
	return true;
}

bool SpatGeom::addHole(SpatHole h) {
	long i = parts.size()-1;
	if (i > -1) {
		parts[i].addHole(h);
		return true;
	} else {
		return false;
	}
}


bool SpatGeom::setPart(SpatPart p, unsigned i) {
	parts[i] = p;
	if (parts.size() > 1) {
		extent.unite(p.extent);
	} else {
		extent = p.extent;
	}
	return true;
}

SpatPart SpatGeom::getPart(unsigned i) {
	return parts[i];
}

SpatVector::SpatVector() {
	srs.proj4 = "+proj=longlat +datum=WGS84";
	srs.wkt = "GEOGCS[\"WGS 84\", DATUM[\"WGS_1984\", SPHEROID[\"WGS 84\",6378137,298.257223563]], PRIMEM[\"Greenwich\",0], UNIT[\"degree\",0.0174532925199433]]";
}

SpatVector::SpatVector(SpatGeom g) {
	addGeom(g);
}

/*
SpatVector::SpatVector(const SpatVector &x) {
	srs = x.srs;
	df = SpatDataFrame(x.df);	
}
*/

SpatVector::SpatVector(SpatExtent e, std::string crs) {
	std::vector<double> x = { e.xmin, e.xmin, e.xmax, e.xmax, e.xmin };
	std::vector<double> y = { e.ymin, e.ymax, e.ymax, e.ymin, e.ymin };
	SpatPart p(x, y);
	SpatGeom g(p);
	g.gtype = polygons;
	setGeom(g);
	setSRS( {crs});
}

SpatVector::SpatVector(std::vector<double> x, std::vector<double> y, SpatGeomType g, std::string crs) {
	SpatPart p(x, y);
	SpatGeom geom(p);
	geom.gtype = g;
	setGeom(geom);
	setSRS( {crs} );
}


std::vector<double> SpatVector::getDv(unsigned i) {
	return df.getD(i);
}

std::vector<long> SpatVector::getIv(unsigned i){
	return df.getI(i);
}

std::vector<std::string> SpatVector::getSv(unsigned i){
	return df.getS(i);
}

std::vector<unsigned> SpatVector::getItype(){
	return df.itype;
}

std::vector<unsigned> SpatVector::getIplace(){
	return df.iplace;
}


std::vector<std::string> SpatVector::get_names(){
	return df.get_names();
}

void SpatVector::set_names(std::vector<std::string> s){
	df.set_names(s);
}

unsigned SpatVector::ncol() {
	return df.ncol();
}

unsigned SpatVector::nrow() {
	return geoms.size();
}

size_t SpatVector::size() {
	return geoms.size();
}

bool SpatVector::is_lonlat() {
	return srs.is_lonlat();
}

bool SpatVector::could_be_lonlat() {
	SpatExtent e = getExtent();
	return srs.could_be_lonlat(e);
}


SpatExtent SpatVector::getExtent(){
	return extent;
}


/*
void SpatVector::setPRJ(std::string PRJ){
	crs[0] = PRJ;
}

std::string SpatVector::getPRJ(){
	return crs[0];
}
*/

std::string SpatVector::type(){
	if (size() == 0) {
		return "none";
	} else if (geoms[0].gtype == points) {
		return "points";
	} else if (geoms[0].gtype == lines) {
		return "lines";
	} else if (geoms[0].gtype == polygons) {
		return "polygons";
	} else {
		return("unknown");
	}
}



SpatGeom SpatVector::getGeom(unsigned i) {
	return geoms[i];
}

bool SpatVector::addGeom(SpatGeom p) {
	geoms.push_back(p);
	if (geoms.size() > 1) {
		extent.unite(p.extent);
	} else {
		extent = p.extent;
	}
	return true;
}

bool SpatVector::setGeom(SpatGeom p) {
	geoms.resize(1);
	geoms[0] = p;
	extent = p.extent;
	return true;
}


unsigned SpatVector::nxy() {
	unsigned n = 0;
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		if (g.size() == 0) {
			n++; // empty
		}
		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			n += p.x.size();
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					n += h.x.size();
				}
			}
		}
	}
	return n;
}


std::vector<std::vector<double>> SpatVector::coordinates() {
	std::vector<std::vector<double>> out(2);
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			for (size_t q=0; q < p.x.size(); q++) {
				out[0].push_back( p.x[q] );
				out[1].push_back( p.y[q] );
			}
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					for (size_t q=0; q < h.x.size(); q++) {
						out[0].push_back( h.x[q] );
						out[1].push_back( h.y[q] );
					}
				}
			}
		}
	}
	return out;
}


SpatDataFrame SpatVector::getGeometryDF() {

	SpatDataFrame out;
	out.add_column(1, "geom");
	out.add_column(1, "part");
	out.add_column(0, "x");
	out.add_column(0, "y");
	out.add_column(1, "hole");

	unsigned n = nxy();
	out.resize_rows(n);

	size_t idx = 0;
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		if (g.size() == 0) { // empty
			out.iv[0][idx] = i+1;
			out.iv[1][idx] = 1;
			out.dv[0][idx] = NAN;
			out.dv[1][idx] = NAN;
			out.iv[2][idx] = 0;
			idx++;
		}

		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			for (size_t q=0; q < p.x.size(); q++) {
				out.iv[0][idx] = i+1;
				out.iv[1][idx] = j+1;
				out.dv[0][idx] = p.x[q];
				out.dv[1][idx] = p.y[q];
				out.iv[2][idx] = 0;
				idx++;
			}
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					for (size_t q=0; q < h.x.size(); q++) {
						out.iv[0][idx] = i+1;
						out.iv[1][idx] = j+1;
						out.dv[0][idx] = h.x[q];
						out.dv[1][idx] = h.y[q];
						out.iv[2][idx] = k+1;
						idx++;
					}
				}
			}
		}
	}
	return out;
}

std::string nice_string(const double &x) {
	std::string s = std::to_string(x);
	s.erase(s.find_last_not_of('0') + 1, std::string::npos);
	s.erase(s.find_last_not_of('.') + 1, std::string::npos);
	return s;
}

std::vector<std::string> SpatVector::getGeometryWKT() {

	std::vector<std::string> out(size());
	std::string wkt;
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		size_t n = g.size();
		if (g.gtype == points) {
			if (n > 1) {
				wkt = "MULTIPOINT ";
			} else {
				wkt = "POINT ";				
			}
		} else if (g.gtype == lines) {
			if (n > 1) {
				wkt = "MULTILINESTRING ";
			} else {
				wkt = "LINESTRING ";		
			}
		} else if (g.gtype == polygons) {
			if (n > 1) {
				wkt = "MULTIPOLYGON ";
			} else {
				wkt = "POLYGON ";
			}
		}
		
		if (n == 0) {
			wkt += "EMPTY";
			out[i] = wkt;
			continue;
		}

		if ((g.gtype == polygons) | (n > 1)) { 
			wkt += "(";
		}			
		
		for (size_t j=0; j < n; j++) {
			SpatPart p = g.getPart(j);
			if (j>0) wkt += ",";

			if ((g.gtype == polygons) & (n > 1)) { 
				wkt += "(";
			}			
		
			wkt += "(" + nice_string(p.x[0]) + " " + nice_string(p.y[0]);
			for (size_t q=1; q < p.x.size(); q++) {
				wkt += ", " + nice_string(p.x[q]) + " " + nice_string(p.y[q]);
			}
			wkt += ")";
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					wkt += ",(" + nice_string(h.x[0]) + " " + nice_string(h.y[0]);
					for (size_t q=1; q < h.x.size(); q++) {
						wkt += ", " + nice_string(h.x[q]) + " " + nice_string(h.y[q]);
					}
					wkt += ")";
				}
			}
			if ((g.gtype == polygons) & (n > 1)) { 
				wkt += ")";
			}			
		}
		if ((g.gtype == polygons) | (n > 1)) {
			wkt += ")";
		}
		out[i] = wkt;
	}
	return out;
}


SpatGeomType SpatVector::getGType(std::string &type) {
	if (type == "points") { return points; }
	else if (type == "lines") { return lines; }
	else if (type == "polygons") { return polygons; }
	else { return unknown; }
}


void SpatVector::setGeometry(std::string type, std::vector<unsigned> gid, std::vector<unsigned> part, std::vector<double> x, std::vector<double> y, std::vector<unsigned> hole) {

// it is assumed that values are sorted by gid, part, hole
	unsigned lastgeom = gid[0];
	unsigned lastpart = part[0];
	unsigned lasthole = hole[0];
	bool isHole = lasthole > 0;

	std::vector<double> X, Y;
	SpatGeom g;
	g.gtype = getGType(type);

	for (size_t i=0; i<gid.size(); i++) {
		if ((lastgeom != gid[i]) || (lastpart != part[i]) || (lasthole != hole[i])) {
            if (g.gtype == polygons) {
                if ((X[0] != X[X.size()-1]) || (Y[0] != Y[Y.size()-1])) {
                    X.push_back(X[0]);
                    Y.push_back(Y[0]);
                }
                if (isHole) {
                    SpatHole h(X, Y);
                    g.addHole(h);
                } else {
                    SpatPart p(X, Y);
                    g.addPart(p);
                }
            } else {
                SpatPart p(X, Y);
                g.addPart(p);
            }
			lastpart = part[i];
            lasthole = hole[i];
            isHole = lasthole > 0;
			X.resize(0);
			Y.resize(0);
			if (lastgeom != gid[i]) {
				addGeom(g);
				g.parts.resize(0);
				lastgeom = gid[i];
			}
		}
		X.push_back(x[i]);
		Y.push_back(y[i]);
	}

	if (g.gtype == polygons) {
        if ((X[0] != X[X.size()-1]) || (Y[0] != Y[Y.size()-1])) {
            X.push_back(X[0]);
            Y.push_back(Y[0]);
        }
        if (isHole) {
            SpatHole h(X, Y);
            g.addHole(h);
        } else {
            SpatPart p(X, Y);
            g.addPart(p);
        }
	} else {
		SpatPart p(X, Y);
		g.addPart(p);
	}
	addGeom(g);
}


SpatVector SpatVector::subset_rows(std::vector<int> range) {

	SpatVector out;
	int n = nrow();
	std::vector<unsigned> r;
	for (size_t i=0; i<range.size(); i++) {
	if ((range[i] >= 0) & (range[i] < n)) {
			r.push_back(range[i]);
		}
	}

	for (size_t i=0; i < r.size(); i++) {
		out.addGeom( geoms[r[i]] );
	}
	out.srs = srs;
	out.df = df.subset_rows(r);
	return out;
}


SpatVector SpatVector::subset_rows(int i) {
	std::vector<int> range(1, i);
	SpatVector out = subset_rows(range);
	return out;
}


SpatVector SpatVector::subset_cols(std::vector<int> range) {
	SpatVector out = *this;
	//out.geoms = geoms;
	//out.crs = crs;
	//out.extent = extent;
	int nc = ncol();

	std::vector<unsigned> r;
	for (size_t i=0; i<range.size(); i++) {
	if ((range[i] >= 0) & (range[i] < nc)) {
			r.push_back(range[i]);
		}
	}
	out.df = df.subset_cols(r);
	return out;
}


SpatVector SpatVector::subset_cols(int i) {
	std::vector<int> range(1, i);
	SpatVector out = subset_cols(range);
	return out;
}


SpatVector SpatVector::shift(double x, double y) {

	SpatVector out = *this;

	for (size_t i=0; i < size(); i++) {
		for (size_t j=0; j < geoms[i].size(); j++) {
			for (size_t q=0; q < geoms[i].parts[j].x.size(); q++) {
				out.geoms[i].parts[j].x[q] = geoms[i].parts[j].x[q] + x;
				out.geoms[i].parts[j].y[q] = geoms[i].parts[j].y[q] + y;
			}
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					for (size_t q=0; q < geoms[i].parts[j].holes[k].x.size(); q++) {
						out.geoms[i].parts[j].holes[k].x[q] = geoms[i].parts[j].holes[k].x[q] + x;
						out.geoms[i].parts[j].holes[k].y[q] = geoms[i].parts[j].holes[k].y[q] + y;
					}
					out.geoms[i].parts[j].holes[k].extent.xmin = out.geoms[i].parts[j].holes[k].extent.xmin + x;
					out.geoms[i].parts[j].holes[k].extent.xmax = out.geoms[i].parts[j].holes[k].extent.xmax + x;
					out.geoms[i].parts[j].holes[k].extent.ymin = out.geoms[i].parts[j].holes[k].extent.ymin + y;
					out.geoms[i].parts[j].holes[k].extent.ymax = out.geoms[i].parts[j].holes[k].extent.ymax + y;
				}
			}
			out.geoms[i].parts[j].extent.xmin = out.geoms[i].parts[j].extent.xmin + x;
			out.geoms[i].parts[j].extent.xmax = out.geoms[i].parts[j].extent.xmax + x;
			out.geoms[i].parts[j].extent.ymin = out.geoms[i].parts[j].extent.ymin + y;
			out.geoms[i].parts[j].extent.ymax = out.geoms[i].parts[j].extent.ymax + y;
		}
		out.geoms[i].extent.xmin = out.geoms[i].extent.xmin + x;
		out.geoms[i].extent.xmax = out.geoms[i].extent.xmax + x;
		out.geoms[i].extent.ymin = out.geoms[i].extent.ymin + y;
		out.geoms[i].extent.ymax = out.geoms[i].extent.ymax + y;
	}
	out.extent.xmin = out.extent.xmin + x;
	out.extent.xmax = out.extent.xmax + x;
	out.extent.ymin = out.extent.ymin + y;
	out.extent.ymax = out.extent.ymax + y;
	return out;
}

