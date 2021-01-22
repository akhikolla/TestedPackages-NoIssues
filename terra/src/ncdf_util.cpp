#include "spatRaster.h"
#include "time.h"
#include "recycle.h"

#include "file_utils.h"
#include "string_utils.h"
#include "gdal_info.h"

bool good_ends(std::string const &s) {
	std::vector<std::string> end = {"_bnds", "_bounds", "lat", "lon", "longitude", "latitude"};
	for (size_t i=0; i<end.size(); i++) {
		if (s.length() >= end[i].length()) {
			if (s.compare(s.length() - end[i].length(), s.length(), end[i]) == 0) {
				return false;
			}
		} 
	}
	if (s == "x" || s == "y" || s == "northing" || s == "easting") {
		return false;
	}
	return true;
}

void pickMost(std::vector<std::string> &sd, std::vector<std::string> &name, std::vector<int> &dim1, std::vector<int> &dim2) {
	if (sd.size() < 2) return;
	std::vector<int> ud = dim1;
	std::sort(ud.begin(), ud.end());
	ud.erase(std::unique(ud.begin(), ud.end()), ud.end());
	if (ud.size() > 1) {
		std::vector<std::string> tmpsd, tmpname;
		std::vector<int> tmpdim1, tmpdim2;
		int mx = ud[ud.size()-1];
		for (size_t i=0; i<sd.size(); i++) {
			if (dim1[i] == mx) {
				tmpsd.push_back(sd[i]);
				tmpname.push_back(name[i]);
				tmpdim1.push_back(dim1[i]);
				tmpdim2.push_back(dim2[i]);
			}
		}
		sd = tmpsd;
		name = tmpname;
		dim1 = tmpdim1;
		dim2 = tmpdim2;
	}
}


bool SpatRaster::constructFromSDS(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname, bool ncdf) {

	std::vector<std::vector<std::string>> info = parse_metadata_sds(meta);
	int n = info[0].size();
	std::vector<std::string> sd, varname;

// std::vector<unsigned> varnl;
// for selection based on nlyr
	
	if (info[0].size() == 0) {
		return false;
	}
	// select sds by index
	if (subds[0] >=0) {
		for (size_t i=0; i<subds.size(); i++) {
			if (subds[i] >=0 && subds[i] < n) {
				sd.push_back(info[0][subds[i]]);
				varname.push_back(info[1][i]);
			} else {
				std::string emsg = std::to_string(subds[i]+1) + " is not valid. There are " + std::to_string(info[0].size()) + " subdatasets\n";
				setError(emsg);
				return false;
			}
		}
	// select by name	
	} else if (subdsname[0] != "") {
		for (size_t i=0; i<subdsname.size(); i++) {
			int w = where_in_vector(subdsname[i], info[1]);
			if (w >= 0) {
				sd.push_back(info[0][w]);
				varname.push_back(info[1][w]);
			} else {
				std::string emsg = concatenate(info[1], ", ");
				emsg = subdsname[i] + " not found. Choose one of:\n" + emsg;
				setError(emsg);
				return false;
			}
		}
	// select all
	} else {
		// eliminate sources based on names like "*_bnds" and "lat"
		std::vector<int> rows, cols;
		for (size_t i=0; i<info[1].size(); i++) {
			if (good_ends(info[1][i])) {
				sd.push_back(info[0][i]);
				varname.push_back(info[1][i]);
				rows.push_back(stoi(info[3][i]));
				cols.push_back(stoi(info[4][i]));
			} 
		}
		if (sd.size() == 0) { // all were removed
			std::vector<size_t> nl(n);
			for (size_t i=0; i<nl.size(); i++) {
				nl[i] = stol(info[5][i]);
			}
			size_t mxnl = *max_element(nl.begin(), nl.end());
			for (size_t i=0; i<nl.size(); i++) {
				if (nl[i] == mxnl) {
					sd.push_back(info[0][i]);
					varname.push_back(info[1][i]);
					rows.push_back(stoi(info[3][i]));
					cols.push_back(stoi(info[4][i]));
				}			
			}
		}
		// pick the ones with most rows and then cols
		// to avoid picking the 1 or 2 "row" datasets
		pickMost(sd, varname, rows, cols);
		pickMost(sd, varname, cols, rows);
	}
	
	std::vector<size_t> srcnl;
	size_t cnt=0;
    for (size_t i=0; i < sd.size(); i++) {
		cnt++;
		bool success = constructFromFile(sd[i], {-1}, {""});
		if (success) break;
	}
	std::vector<std::string> skipped, used;
	srcnl.push_back(nlyr());
	used.push_back(varname[0]);				
	SpatRaster out;
    for (size_t i=cnt; i < sd.size(); i++) {
//		printf( "%s\n", sd[i].c_str() );
		bool success = out.constructFromFile(sd[i], {-1}, {""});
		if (success) {
			if (out.compare_geom(*this, false, false)) {
				addSource(out);
				srcnl.push_back(out.nlyr());
				used.push_back(varname[i]);				
			} else {
				skipped.push_back(varname[i]);
			}
		} else {
			skipped.push_back(varname[i]);
		}
	}

	if (skipped.size() > 0) {
		std::string s="skipped sub-datasets (see 'desc(sds=TRUE)'):\n" + skipped[0];
		for (size_t i=1; i<skipped.size(); i++) {
			s += ", " + skipped[i];
			if ((i%3) == 0) s += "\n";
		}
		addWarning(s);
	}

	if (!ncdf) {
		std::vector<std::string> lyrnames;
		for (size_t i=0; i<used.size(); i++) {
			std::vector<std::string> nms = {basename(used[i])};
			recycle(nms, srcnl[i]);
			make_unique_names(nms);
			lyrnames.insert(lyrnames.end(), nms.begin(), nms.end());
			//Rcpp::Rcout << used[i] << std::endl;
			//Rcpp::Rcout << nms.size() << std::endl;		
		}
		if (lyrnames.size() > 0) {
			setNames(lyrnames, false);
		}
	}

	return true;
}



std::vector<int_64> str2int64v(std::string s, std::string delim) {
	std::vector<int_64> out;
	size_t pos = 0;
	while ((pos = s.find(delim)) != std::string::npos) {
		std::string v = s.substr(0, pos);
		s.erase(0, pos + 1);
		out.push_back(std::stoll(v));
	}
	out.push_back(std::stoll(s));
	return out;
}


std::vector<int_64> ncdf_time(const std::vector<std::string> &metadata, std::vector<std::string> vals, std::string &step, std::string &msg) {


	std::vector<int_64> out, bad;
	if (vals.size() < 1) {
		step = "";
		return out;
	}
	
	for (size_t i=0; i<vals.size(); i++) {
		out.push_back(stoll(vals[i]));
	}
	
	bool fu=false;
	bool fc=false;
	std::string origin;
	std::string calendar = "standard";
	for (size_t i=0; i<metadata.size(); i++) {
		if (!fc) {
			std::string pattern = "time#calendar=";
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				calendar = metadata[i];
				calendar.erase(calendar.begin(), calendar.begin()+pattern.size());  
				fc = true;
			}
		}
		if (!fu) {
			std::string pattern = "time#units=";
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				origin = metadata[i];
				origin.erase(origin.begin(), origin.begin()+pattern.size());  
				fu = true;
			}
		}
		if (fc & fu) break;
	}

	bool days = false; 
	bool hours = false;
	bool seconds = false; 
	bool foundorigin = false;

	if (fu) {
		if ((origin.find("hours")) != std::string::npos) {
			hours = true;
		} else if ((origin.find("days")) != std::string::npos) {
			days = true;
		} else if ((origin.find("seconds")) != std::string::npos) {
			seconds = true;		
		} 
		size_t pos;
		if ((pos = origin.find("from")) != std::string::npos) {
			origin.erase(0, pos + 5);
			foundorigin = true;
		} else if ((pos = origin.find("since")) != std::string::npos) {
			origin.erase(0, pos + 6);
			foundorigin = true;	
		}
	}
	
	SpatTime_t offset = 0;
	if (foundorigin) {
		step = "seconds";
		if (days) {
			if (calendar == "noleap" || calendar == "365_day" || calendar == "365 day") { 
				std::vector<int> ymd = getymd(origin);
				for (int_64 &d : out) d = time_from_day_noleap(ymd[0], ymd[1], ymd[2], d);
			} else if (calendar == "360_day" || calendar == "360 day") { 
				std::vector<int> ymd = getymd(origin);
				for (int_64 &d : out) d = time_from_day_360(ymd[0], ymd[1], ymd[2], d);
			} else { 
				if (!(calendar =="gregorian" || calendar =="proleptic_gregorian" || calendar=="standard" || calendar == "julian")) { 
					// julian is perhaps questionable it can mean different things.
					msg = "unknown calendar (assuming standard): " + calendar;
				}
				std::vector<int> ymd = getymd(origin);
				for (int_64 &d : out) d = time_from_day(ymd[0], ymd[1], ymd[2], d);
			}
		} else if (hours) {
			hours_to_time(out, origin);
			//std::vector<int> ymd = getymd(origin);
			//for (int_64 &d : out) d = time_from_hour(ymd[0], ymd[1], ymd[2], d);
		} else if (seconds) {
			offset = get_time_string(origin);
			for (int_64 &d : out) d = d + offset;
		} else {
			step = "raw";
		}
	}

	return out;
}
	
	
//NETCDF_DIM_k=0
//NETCDF_DIM_tile=0
//NETCDF_DIM_time=0
//NETCDF_VARNAME=NVEL

std::vector<std::vector<std::string>> ncdf_names(const std::vector<std::vector<std::string>> &m) {

	std::vector<std::vector<std::string>> out(4);
	out[2] = out[3] = std::vector<std::string>(m.size(), "");

	for (size_t i=0; i<m.size(); i++) {
		std::vector<std::string> b = m[i];
		std::string vname, dim;
		for (size_t j=0; j<b.size(); j++) {
			
			size_t pos = b[j].find("NETCDF_VARNAME");
			if (pos != std::string::npos) {
				vname = b[j].erase(0, pos+15);
				continue;
			} 
			
			pos = b[j].find("NETCDF_DIM_");
			if (pos != std::string::npos) {
				size_t pos = b[j].find("NETCDF_DIM_time");
				if (pos != std::string::npos) {
					out[1].push_back( b[j].erase(0, pos+16) );
				} else {
					dim += b[j].erase(0, pos+11);
				}
				continue;
			}
			pos = b[j].find("units=");
			if (pos != std::string::npos) {
				out[2][i] = b[j].erase(0, pos+6);
				continue;
			}		
			pos = b[j].find("long_name=");
			if (pos != std::string::npos) {
				out[3][i] = b[j].erase(0, pos+10);
				continue;
			}		
			pos = b[j].find("standard_name=");
			if (pos != std::string::npos) {
				if (out[3][i] == "") {
					out[3][i] = b[j].erase(0, pos+14);
				}
			}
		}
		out[0].push_back(vname + dim); 
	}
	return out;
}

void RasterSource::set_names_time_ncdf(std::vector<std::string> metadata, std::vector<std::vector<std::string>> bandmeta, std::string &msg) {

	if (bandmeta.size() == 0) return;

	std::vector<std::vector<std::string>> nms = ncdf_names(bandmeta);

	if (nms[0].size() > 0) {
		make_unique_names(nms[0]);
		names = nms[0];
	}
	unit = nms[2];
	long_names = nms[3];
	if (nms[1].size() > 0) {
		std::string step;
		std::vector<int_64> x = ncdf_time(metadata, nms[1], step, msg);
		if (x.size() == nlyr) {
			time = x;
			timestep = step;
			hasTime = true;
		}
	}

}


/*


std::vector<std::vector<std::string>> metatime(std::vector<std::string> meta) {
	std::vector<std::vector<std::string>> out(meta.size());
	std::string delim = "=";
	for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(delim);
		if (pos != std::string::npos) {
			out[i].push_back(s.erase(pos+1, std::string::npos));
			out[i].push_back(s.erase(0, pos+1));
		} else {
			out[i].push_back(s);
		}
	}
	return out;
}


*/


/*
bool SpatRaster::constructFromSDS(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname) {

	std::vector<std::string> sd; //, nms;
//	std::vector<std::string> dc; //, nms;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";
	for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			sd.push_back(s);
		} //else {
		//	size_t pos = s.find(ddelim);
		//	if (pos != std::string::npos) {
		//		s.erase(0, pos + ddelim.length());
		//		dc.push_back(s);
		//	}
		//}
	}
	if (sd.size() == 0) {
		return false;
	}
	//bool useDC = (dc.size() == sd.size());
	int sdsize = sd.size();
	if (subds[0] >=0) {
		std::vector<std::string> tmp;
		for (size_t i=0; i<subds.size(); i++) {
			if (subds[i] >=0 && subds[i] < sdsize) {
				tmp.push_back(sd[subds[i]]);
			//	if (useDC) {
			//		dc = {dc[subds[0]]};
			//	}
			} else {
				std::string emsg = std::to_string(subds[i]+1) + " is not valid. There are " + std::to_string(sd.size()) + " subdatasets\n";
				setError(emsg);
				return false;
			}
		}
		sd = tmp;		
	} else if (subdsname[0] != "") {
		std::vector<std::string> tmp;
		std::vector<std::string> shortnames = getlastpart(sd, ":");
		for (size_t i=0; i<subdsname.size(); i++) {
			int w = where_in_vector(subdsname[i], shortnames);
			if (w >= 0) {
				tmp.push_back(sd[w]);
			//	if (useDC) {
			//		dc = {dc[w]};
			//	}			
			} else {
				std::string emsg = concatenate(shortnames, ", ");
				emsg = subdsname[i] + " not found. Choose one of:\n" + emsg;
				setError(emsg);
				return false;
			}
		}
		sd = tmp;
	}
	
	bool success = constructFromFile(sd[0], {-1}, {""});
	if (!success) {
		// should continue to the next one  with while
		return false;
	}
	SpatRaster out;
	std::vector<int> skipped;
    for (size_t i=1; i < sd.size(); i++) {
//		printf( "%s\n", sd[i].c_str() );
		success = out.constructFromFile(sd[i], {-1}, {""});
		if (success) {
			if (out.compare_geom(*this, false, false)) {
				addSource(out);
			} else {
				skipped.push_back(i);
			}
		} else {
			if (out.msg.has_error) {
				//setError(out.msg.error);
				//addWarning(out.msg.error);
			}
			//return false;
		}
	}

	for (std::string& s : sd) s = basename_sds(s);
	if (skipped.size() > 0) {
		std::string s="skipped subdatasets (different geometry):";
		for (size_t i=0; i<skipped.size(); i++) {
			s += "\n   " + sd[skipped[i]];
		}
		s += "\nSee 'describe_sds' for more info";
		addWarning(s);
		for (int i=skipped.size()-1; i>0; i--) {
			sd.erase(sd.begin() + skipped[i]);
		}
	}
	success = setNames(sd);
	return true;
}
*/
