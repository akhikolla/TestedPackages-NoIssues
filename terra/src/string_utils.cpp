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

#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <numeric>


std::string double_to_string(double x) { 
	std::string s = std::to_string (x);
	s.erase( s.find_last_not_of('0') + 1, std::string::npos );
	s.erase( s.find_last_not_of('.') + 1, std::string::npos );
	return s;
}


std::vector<char *> string_to_charpnt(std::vector<std::string> s) {
	size_t n = s.size();
	std::vector<char *> out(n + 1);
	for (size_t i = 0; i < n; i++) {
		out[i] = (char *) (s[i].c_str());
	}
	out[n] = NULL;
	return out;
}

//std::vector<char> string_to_char(std::vector<std::string> s) {
//	std::vector<char> charstr(s.c_str(), s.c_str() + s.size() + 1);
//	return charstr;
//}

std::vector<std::string> double_to_string(const std::vector<double> &x, std::string prep) { 
	std::vector<std::string> out(x.size());
	for (size_t i=0; i<x.size(); i++) {
		out[i] = prep + double_to_string (x[i]);
	}
	return out;
}


std::string concatenate(std::vector<std::string> v, std::string delim) {
	for (size_t i=0; i<(v.size()-1); i++) {
		v[i] = v[i] + delim;
	}
	std::string s;
	for (const auto &piece : v) s += piece;
	return s;
}


std::vector<std::string> getlastpart (std::vector<std::string> s, std::string delim) {
	std::vector<std::string> out(s.size());
	for (size_t i=0; i<s.size(); i++) {
		std::size_t pos = s[i].find_last_of(delim);
		out[i] = s[i].substr(pos+1);
	}
	return out;
}


bool in_string(const std::string &x, std::string part) {
	size_t f = x.find(part);
	if (f == std::string::npos) {
		return false;
	} else {
		return true;
	}
}


bool ends_on(std::string const &s, std::string const &end) {
	if (s.length() >= end.length()) {
		if (s.compare(s.length() - end.length(), s.length(), end) == 0) {
			return false;
		} 
	}
	return true;
}


void lowercase(std::string &s) {
	std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

bool is_in_vector(std::string s, std::vector<std::string> ss) {
	//std::set<std::string> sset (ss.begin(), ss.end());
	//return sset.find(s) != sset.end();
	auto it = std::find (ss.begin(), ss.end(), s);
	return (it != ss.end());
}


int where_in_vector(std::string s, std::vector<std::string> ss) {
	int i = -1;
	auto it = std::find (ss.begin(), ss.end(), s);
	if (it != ss.end()) {
		i = std::distance(ss.begin(), it);
	}
	return i;
}


std::string is_in_set_default(std::string s, std::vector<std::string> ss, std::string defvalue, bool tolower) {
	if (tolower) lowercase(s);
	std::set<std::string> sset (ss.begin(), ss.end());
	if (sset.find(s) == sset.end() ) {
		s = defvalue;
	}
	return s;
}


std::vector<std::string> strsplit(std::string s, std::string delimiter){
	std::vector<std::string> out;
	size_t pos = 0;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		out.push_back(s.substr(0, pos));
		s.erase(0, pos + delimiter.length());
	}
	out.push_back(s.substr(0, pos));
	return out;
}


std::vector<double> str2dbl(std::vector<std::string> s) {
	std::vector<double> d (s.size());
	std::transform(s.begin(), s.end(), d.begin(), [](const std::string& val) {
		return std::stod(val);
	});
	return d;
}

std::vector<long> str2long(std::vector<std::string> s) {
	std::vector<long> d (s.size());
	std::transform(s.begin(), s.end(), d.begin(), [](const std::string& val) {
		return std::stol(val);
	});
	return d;
}

std::vector<int> str2int(std::vector<std::string> s) {
	std::vector<int> d (s.size());
	std::transform(s.begin(), s.end(), d.begin(), [](const std::string& val) {
		return std::stoi(val);
	});
	return d;
}

std::vector<std::string> dbl2str(std::vector<double> d) {
	std::vector<std::string> s (d.size());
	std::transform(d.begin(), d.end(), s.begin(),
			[](double i) { return std::to_string(i); }
	);
	return s;
}



// trim from start (in place)
void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
void lrtrim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
std::string lrtrim_copy(std::string s) {
    lrtrim(s);
    return s;
}



void make_valid_names(std::vector<std::string> &s) {
    for (size_t i=0; i<s.size(); i++) {
        lrtrim(s[i]);
        if (s[i] == "") s[i] = "X";
        if (isdigit(s[i][0])) s[i] = "X" + s[i];
//        if ((s[i][0] == ".") & (s[i].size() > 1)) {
//			if (isdigit(s[i][1])) s[i] = "X" + s[i];
//		}
        std::replace(s[i].begin(), s[i].end(), ' ', '.');
    }
}


template <typename T>
std::vector<size_t> order(const std::vector<T> &v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}


void make_unique_names(std::vector<std::string> &s) {
    std::vector<size_t> x = order(s);
    std::sort(s.begin(), s.end());
    std::vector<std::string> ss = s;
    unsigned j = 1;
    for (size_t i=1; i<s.size(); i++) {
        if (s[i] == s[i-1]) {
			if (j==1) {
				ss[i-1] = s[i-1] + "_1";
			}
			j++;
            ss[i] = s[i] + "_" + std::to_string(j);
        } else {
            j = 1;
        }
    }
    for (size_t i=0; i<s.size(); i++) {
        s[x[i]] = ss[i];
    }
}


void str_replace(std::string& s, const std::string& from, const std::string& to) {
    size_t start_pos = s.find(from);
    if (start_pos != std::string::npos) {
		s.replace(start_pos, from.length(), to);
	}
}

size_t str_replace_all(std::string& str, const std::string& from, const std::string& to) {
    size_t count = 0;
    if (from.empty()) return count;
    size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
		count++;
    }
	return count;
}

