std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> openopts);
std::vector<std::vector<std::string>> sdinfo(std::string fname);
std::vector<std::string> get_metadata(std::string filename);
std::vector<std::string> get_metadata_sds(std::string filename);
std::vector<std::vector<std::string>> parse_metadata_sds(std::vector<std::string> meta);

