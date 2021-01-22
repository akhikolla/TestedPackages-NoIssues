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

#include <fstream>
#include <numeric>
#include "spatVector.h"

#ifdef useGDAL
#include "gdal_priv.h"
#endif

#ifdef useRcpp
#include <Rcpp.h>
// Rcpp::depends(RcppProgress)
#include "progress.hpp"
#include "progress_bar.hpp"
#endif

typedef long long int_64;
typedef unsigned long long uint_64;

class SpatCategories {
	public:
		std::vector<double> levels;
		std::vector<std::string> labels;
};


class SpatWindow {
	public:
		SpatExtent full_extent;
		uint_64 full_ncol, full_nrow, off_row, off_col;
		bool expanded = false;
		std::vector<uint_64> expand;
};




class RasterSource {
    private:
//		std::ofstream ofs;
	public:
#ifdef useGDAL
		GDALDataset* gdalconnection;
#endif
		bool open_read=false;
		bool open_write=false;

		RasterSource();

//		void fsopen(std::string filename);
//		bool fswrite(std::vector<double> &v);
//		void fsclose();

		uint_64 ncol, nrow;
		unsigned nlyr;
		unsigned nlyrfile = 0;
		SpatExtent extent;
		bool rotated=false;
		bool flipped=false;
		bool hasWindow=false;
		SpatWindow window;
	
		
		//std::vector<std::string> crs = std::vector<std::string>(2, "");
		SpatSRS srs;
		std::vector<unsigned> layers;
		std::vector<std::string> names;
		std::vector<std::string> long_names;
		std::vector<int_64> time;
		std::string timestep = "seconds";
		bool hasTime = false;
		std::vector<double> depth;
		std::vector<std::string> unit;

		//std::vector< std::vector<double> values;
        std::vector<double> values;
        std::vector<int64_t> ivalues;
        std::vector<bool> bvalues;

		std::vector<bool> hasRange;
		std::vector<double> range_min;
		std::vector<double> range_max;
		std::vector<bool> hasAttributes;
		std::vector<SpatDataFrame> atts;
		std::vector<bool> hasCategories;
		std::vector<SpatCategories> cats;
		std::vector<bool> hasColors;
		std::vector<SpatDataFrame> cols;

		bool memory=true;
		bool hasValues=false;
		std::string filename;
		std::string driver;
		std::string datatype; 
		
		// for native files
		//double NAflag;

		std::vector<bool> has_scale_offset;
		std::vector<double> scale;
		std::vector<double> offset;

//		std::vector<RasterSource> subset(std::vector<unsigned> lyrs);
		RasterSource subset(std::vector<unsigned> lyrs);
		std::vector<double> getValues(unsigned lyr);
		void setRange();
		void resize(unsigned n);
		bool in_order();
		bool combine_sources(const RasterSource &x);

		bool parameters_changed = false;		
		
		void set_names_time_ncdf(std::vector<std::string> metadata, std::vector<std::vector<std::string>> bandmeta, std::string &msg);
		
};


class BlockSize {
	public:
		std::vector<uint_64> row;
		std::vector<uint_64> nrows;
		unsigned n;
};

class SpatRaster {

	protected:
		//SpatExtent extent;
		SpatExtent window;
		//SpatSRS srs;

	public:
	//	std::string prj;

	//	bool GDALregistred = false;

#ifdef useRcpp
		Progress* pbar;
#endif

////////////////////////////////////////////////////
// properties and property-like methods for entire object
////////////////////////////////////////////////////

		std::string name;
		
		std::vector<RasterSource> source;

		BlockSize bs;
		//BlockSize getBlockSize(unsigned n, double frac, unsigned steps=0);
		BlockSize getBlockSize(SpatOptions &opt);
		std::vector<double> mem_needs(SpatOptions &opt);

		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		std::string getError() { return msg.getError(); }
		void addWarning(std::string s) { msg.addWarning(s); }
		void setMessage(std::string s) { msg.setMessage(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }

		//double NA = std::numeric_limits<double>::quiet_NaN();

		uint_64 ncol();
		uint_64 nrow();
		SpatExtent getExtent();
		void setExtent(SpatExtent e);
		void setExtent(SpatExtent ext, bool keepRes=false, std::string snap="");  // also set it for sources?
		SpatVector dense_extent();

		//std::vector<std::string> getCRS();
		//void setCRS(std::vector<std::string> _crs);

		std::string getSRS(std::string x);
		bool setSRS(std::string crs);


/*
#ifdef useGDAL	
		bool setSRS(OGRSpatialReference *poSRS, std::string &msg) {
#endif 
*/

		bool is_lonlat();
		bool could_be_lonlat();
		bool is_global_lonlat();

		std::vector<double> resolution();
		SpatRaster setResolution(double xres, double yres);
		double ncell() { return nrow() * ncol(); }
		double size() { return ncol() * nrow() * nlyr() ; }

		double xres();
		double yres();
		std::vector<double> origin();
		unsigned nlyr();

		// only no values allowed with a single RasterSource
		bool hasValues() { return source[0].hasValues ; };
		std::vector<double> getValues(long lyr= -1);
		
		bool getValuesSource(size_t src, std::vector<double> &out);				
		bool setValues(std::vector<double> _values);
		bool replaceValues(std::vector<double> cells, std::vector<double> _values, int ncols);
		void setRange();


		
////////////////////////////////////////////////////
// property like methods for RasterSources
////////////////////////////////////////////////////
		std::vector<std::string> filenames();
		bool isSource(std::string filename);
		std::vector<bool> inMemory();


////////////////////////////////////////////////////
// property like methods for layers
////////////////////////////////////////////////////

		std::vector<bool> hasRange();
		std::vector<double> range_min();
		std::vector<double> range_max();
		std::vector<std::string> getNames();
		std::vector<std::string> getLongNames();
		bool setLongNames(std::vector<std::string> nms);
		bool setNames(std::vector<std::string> names, bool make_valid=false);
		bool hasTime();
		std::vector<int_64> getTime();
		std::string getTimeStep();
		
		std::vector<std::string> getTimeStr();

		bool setTime(std::vector<int_64> time, std::string step);
		std::vector<double> getDepth();
		bool setDepth(std::vector<double> depths);
		std::vector<std::string> getUnit();
		bool setUnit(std::vector<std::string> units);

////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////

		SpatRaster();
		SpatRaster(unsigned nr, unsigned nc, unsigned nl, SpatExtent ext, std::string crs);
		SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string crs);
		SpatRaster(std::vector<std::string> fname, std::vector<int> subds, std::vector<std::string> subdsname, std::string x);
		SpatRaster(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname);
		SpatRaster(RasterSource s);
		void setSource(RasterSource s);
		void setSources(std::vector<RasterSource> s);
		//SpatRaster(const SpatRaster& x);

        SpatRaster deepCopy();
        SpatRaster geometry(long nlyrs=-1, bool properties=false);

		bool constructFromFile(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname);
		bool constructFromFiles(std::vector<std::string> fnames);
		bool constructFromSDS(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname, bool ncdf);
//		bool constructFromNCDFsds(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname);


		void addSource(SpatRaster x);
		SpatRaster combineSources(SpatRaster x);
		SpatRaster subset(std::vector<unsigned> lyrs, SpatOptions &opt);

////////////////////////////////////////////////////
// helper methods
////////////////////////////////////////////////////

		void gdalogrproj_init(std::string path);

		bool compare_geom(SpatRaster x, bool lyrs, bool crs, bool warncrs=false, bool ext=true, bool rowcol=true, bool res=false);
		bool compare_origin(std::vector<double> x, double tol);
		bool shared_basegeom(SpatRaster &x, double tol, bool test_overlap);

		std::vector<double> cellFromXY (std::vector<double> x, std::vector<double> y);
		double cellFromXY(double x, double y);
		std::vector<double> cellFromRowCol(std::vector<int_64> row, std::vector<int_64> col);
		double cellFromRowCol(int_64 row, int_64 col);
		std::vector<double> cellFromRowColCombine(std::vector<int_64> row, std::vector<int_64> col);
		double cellFromRowColCombine(int_64 row, int_64 col);
		std::vector<double> yFromRow(std::vector<int_64> &row);
		double yFromRow(int_64 row);
		std::vector<double> xFromCol(std::vector<int_64> &col);
		double xFromCol(int_64 col);

		std::vector<int_64> colFromX(std::vector<double> &x);
		int_64 colFromX(double x);
		std::vector<int_64> rowFromY(std::vector<double> &y);
		int_64 rowFromY(double y);
		std::vector<std::vector<double>> xyFromCell( std::vector<double> &cell);
		std::vector<std::vector<double>> xyFromCell( double cell);
		std::vector<std::vector<int_64>> rowColFromCell(std::vector<double> &cell);
		std::vector<int_64> rowColFromY(std::vector<double> &y);
		std::vector<std::vector<int_64>> rowColFromExtent(SpatExtent e);
	
		
        std::vector<unsigned> sourcesFromLyrs(std::vector<unsigned> lyrs);
		int sourceFromLyr(unsigned lyr);
		std::vector<unsigned> findLyr(unsigned lyr);

        std::vector<unsigned> nlyrBySource();
        std::vector<unsigned> lyrsBySource();
        unsigned nsrc();

		void createCategories(unsigned layer);
		std::vector<bool> hasCategories();
		void setCategories(unsigned layer, std::vector<std::string> labs, std::vector<double> levs);
		std::vector<SpatCategories> getCategories();
		SpatCategories getLayerCategories(unsigned layer);

		void createAttributes(unsigned layer);
		std::vector<bool> hasAttributes();
		void setAttributes(unsigned layer, SpatDataFrame df);
		std::vector<SpatDataFrame> getAttributes();
		SpatDataFrame getLayerAttributes(unsigned layer);
		std::vector<bool> hasColors();
		std::vector<SpatDataFrame> getColors();

		double valuesCell(double);
		double valuesCell(int, int);
		std::vector<double> valuesCell(std::vector<double>);
		std::vector<double> valuesRow(int);

		bool isLonLat();
		bool couldBeLonLat();
		bool isGlobalLonLat();


////////////////////////////////////////////////////
// read and write
////////////////////////////////////////////////////

		bool valid_sources(bool files=true, bool rotated=true);
		bool readStart();
		std::vector<double> readValues(uint_64 row, uint_64 nrows, uint_64 col, uint_64 ncols);
		void readChunkMEM(std::vector<double> &out, size_t src, uint_64 row, uint_64 nrows, uint_64 col, uint_64 ncols);

		std::vector<double> readBlock(BlockSize bs, unsigned i);
		std::vector<std::vector<double>> readBlock2(BlockSize bs, unsigned i);
		std::vector<double> readBlockIP(BlockSize bs, unsigned i);		
		std::vector<double> readExtent(SpatExtent e);
		bool readStop();


		bool writeStart(SpatOptions &opt);
		bool writeValues(std::vector<double> &vals, uint_64 startrow, uint_64 nrows, uint_64 startcol, uint_64 ncols);
		bool writeValues2(std::vector<std::vector<double>> &vals, uint_64 startrow, uint_64 nrows, uint_64 startcol, uint_64 ncols);
		bool writeStop();
		bool writeHDR(std::string filename);

		//bool writeStartGDAL(std::string filename, std::string driver, std::string datatype, bool overwrite, SpatOptions &opt);
		bool writeStartGDAL(SpatOptions &opt);		
		bool fillValuesGDAL(double fillvalue);
		bool writeValuesGDAL(std::vector<double> &vals, uint_64 startrow, uint_64 nrows, uint_64 startcol, uint_64 ncols);
		bool writeStopGDAL();

		//bool writeStartBinary(std::string filename, std::string datatype, std::string bandorder, bool overwrite);
		//bool writeValuesBinary(std::vector<double> &vals, unsigned startrow, unsigned nrows, unsigned startcol, unsigned ncols);

		bool writeValuesMem(std::vector<double> &vals, uint_64 startrow, uint_64 nrows, uint_64 startcol, uint_64 ncols);

		// binary (flat) source
		//std::vector<double> readValuesBinary(unsigned src, unsigned row, unsigned nrows, unsigned col, unsigned ncols);
		//std::vector<double> readSampleBinary(unsigned src, unsigned srows, unsigned scols);
		//std::vector<std::vector<double>> readCellsBinary(unsigned src, std::vector<double> cells);

		// gdal source
		std::vector<double> readValuesGDAL(unsigned src, uint_64 row, uint_64 nrows, uint_64 col, uint_64 ncols, int lyr = -1);
		std::vector<double> readGDALsample(unsigned src, int srows, int scols);
		std::vector<std::vector<double>> readRowColGDAL(unsigned src, std::vector<int_64> &rows, const std::vector<int_64> &cols);

		bool readStartGDAL(unsigned src);
		bool readStopGDAL(unsigned src);
		void readChunkGDAL(std::vector<double> &data, unsigned src, uint_64 row, unsigned nrows, uint_64 col, unsigned ncols);

		bool setWindow(SpatExtent x);
		bool removeWindow();
		bool hasWindow();

		void openFS(std::string const &filename);

		SpatRaster writeRaster(SpatOptions &opt);
		//SpatRaster writeRasterGDAL(std::string filename, std::string format, std::string datatype, bool overwrite, SpatOptions &opt);
		//SpatRaster writeRasterBinary(std::string filename, std::string datatype, std::string bandorder, bool overwrite);

		bool canProcessInMemory(SpatOptions &opt);
		uint_64 chunkSize(unsigned n, double frac);

		void fill(double x);

		SpatRaster sources_to_disk(std::vector<std::string> &tmpfs, bool unique, SpatOptions &opt);
		bool sources_from_file();

		bool differentFilenames(std::vector<std::string> outf);


////////////////////////////////////////////////////
// main methods
////////////////////////////////////////////////////

		SpatRaster collapse_sources();
		SpatRaster rectify(std::string method, SpatRaster aoi, unsigned useaoi, bool snap, SpatOptions &opt);
		
        std::vector<std::vector<double>> adjacent(std::vector<double> cells, std::string directions, bool include);
 		SpatRaster aggregate(std::vector<unsigned> fact, std::string fun, bool narm, SpatOptions &opt);
		SpatExtent align(SpatExtent e, std::string snap);
		SpatRaster rst_area(SpatOptions &opt);
		std::vector<double> sum_area(SpatOptions &opt);
		std::vector<std::vector<double>> area_by_value(SpatOptions &opt);

		SpatRaster arith(SpatRaster x, std::string oper, SpatOptions &opt);
		SpatRaster arith(double x, std::string oper, bool reverse, SpatOptions &opt);
		SpatRaster arith(std::vector<double> x, std::string oper, bool reverse, SpatOptions &opt);
		SpatRaster apply(std::vector<unsigned> ind, std::string fun, bool narm, std::vector<std::string> nms, SpatOptions &opt);
		SpatRaster rapply(SpatRaster x, std::string fun, bool narm, SpatOptions &opt);
		std::vector<std::vector<double>> rappvals(SpatRaster x, size_t startrow, size_t nrows);

		SpatVector as_polygons(bool trunc, bool dissolve, bool values, bool narm, SpatOptions &opt);
		SpatVector polygonize(bool trunc, SpatOptions &opt);
		SpatVector as_points(bool values, bool narm, SpatOptions &opt);
		SpatRaster atan_2(SpatRaster x, SpatOptions &opt);

		std::vector<std::vector<double>> bilinearValues(std::vector<double> x, std::vector<double> y);


		SpatRaster buffer(double d, SpatOptions &opt);
		SpatRaster clamp(double low, double high, bool usevalue, SpatOptions &opt);
		SpatRaster cover(SpatRaster x, double value, SpatOptions &opt);
		SpatRaster crop(SpatExtent e, std::string snap, SpatOptions &opt);
		SpatRaster cum(std::string fun, bool narm, SpatOptions &opt);
        SpatRaster disaggregate(std::vector<unsigned> fact, SpatOptions &opt);
		SpatRaster distance(SpatOptions &opt);
		SpatRaster distance(SpatVector p, SpatOptions &opt);

		SpatRaster edges(bool classes, std::string type, unsigned directions, SpatOptions &opt);
		SpatRaster extend(SpatExtent e, SpatOptions &opt);
		std::vector<std::vector<std::vector<double>>> extractVector(SpatVector v, bool touches, std::string method="");
		std::vector<double> getCells(SpatVector v, bool touches, std::string method);

		std::vector<std::vector<double>> extractCell(std::vector<double> &cell);
        std::vector<std::vector<double>> extractXY(std::vector<double> &x, std::vector<double> &y, std::string method);
		
		SpatRaster flip(bool vertical, SpatOptions &opt);
		SpatRaster filler(SpatRaster x, SpatOptions &opt);
		
		SpatRaster focal(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, std::string fun, SpatOptions &opt);
		std::vector<double> focal_values(std::vector<unsigned> w, double fillvalue, int row, int nrows);
		std::vector<std::vector<double>> freq(bool bylayer, bool round, int digits, SpatOptions &opt);
		std::vector<size_t> count(double value, bool bylayer, bool round, int digits, SpatOptions &opt);
		
		bool get_aggregate_dims(std::vector<unsigned> &fact, std::string &message);
		std::vector<unsigned> get_aggregate_dims2(std::vector<unsigned> fact);
		std::vector<std::vector<double> > get_aggregates(std::vector<double> &in, size_t nr, std::vector<unsigned> dim);
//		std::vector<double> compute_aggregates(std::vector<double> &in, size_t nr, std::vector<unsigned> dim, std::function<double(std::vector<double>&, bool)> fun, bool narm);
		SpatDataFrame global(std::string fun, bool narm, SpatOptions &opt);
		SpatDataFrame global_weighted_mean(SpatRaster &weights, std::string fun, bool narm, SpatOptions &opt);

		SpatRaster gridDistance(SpatOptions &opt);
		SpatRaster gridCostDistance(SpatRaster cost, SpatOptions &opt);

		SpatRaster init(std::string value, bool plusone, SpatOptions &opt);
		SpatRaster init(double value, SpatOptions &opt);
		SpatRaster is_in(std::vector<double> m, SpatOptions &opt);
		std::vector<std::vector<double>> is_in_cells(std::vector<double> m, SpatOptions &opt);

		SpatRaster isnot(SpatOptions &opt);
		SpatRaster isnan(SpatOptions &opt);
		SpatRaster isnotnan(SpatOptions &opt);
		SpatRaster isfinite(SpatOptions &opt);
		SpatRaster isinfinite(SpatOptions &opt);

		std::vector<double> line_cells(SpatGeom& g);
		SpatRaster logic(SpatRaster x, std::string oper, SpatOptions &opt);
		SpatRaster logic(bool x, std::string oper, SpatOptions &opt);
		SpatRaster mask(SpatRaster x, bool inverse, double maskvalue, double updatevalue, SpatOptions &opt);
		SpatRaster mask(SpatVector x, bool inverse, double updatevalue, SpatOptions &opt);
		SpatRaster math(std::string fun, SpatOptions &opt);
		SpatRaster math2(std::string fun, unsigned digits, SpatOptions &opt);


		SpatRaster separate(std::vector<double> classes, double keepvalue, double othervalue, SpatOptions &opt);

		SpatRaster modal(std::vector<double> add, std::string ties, bool narm, SpatOptions &opt);

        std::vector<double> polygon_cells(SpatGeom& g);
		SpatRaster quantile(std::vector<double> probs, bool narm, SpatOptions &opt);
		SpatRaster stretch(std::vector<double> minv, std::vector<double> maxv, std::vector<double> minq, std::vector<double> maxq, std::vector<double> smin, std::vector<double> smax, SpatOptions &opt);
		SpatRaster reverse(SpatOptions &opt);

		SpatRaster range(std::vector<double> add, bool narm, SpatOptions &opt);
		//SpatRaster rasterize(SpatVector p, std::vector<double> values, double background, bool update, SpatOptions &opt);
		SpatRaster rasterize(SpatVector x, std::string field, std::vector<double> values, double background, bool update, bool touches, bool inverse, SpatOptions &opt);
		std::vector<double> rasterizeCells(SpatVector &v, bool touches);

		SpatRaster reclassify(std::vector<std::vector<double>> rcl, unsigned right, bool lowest, bool othersNA, SpatOptions &opt);
		SpatRaster reclassify(std::vector<double> rcl, unsigned nc, unsigned right, bool lowest, bool othersNA, SpatOptions &opt);
		SpatRaster classify_layers(std::vector<std::vector<double>> groups, std::vector<double> id, SpatOptions &opt);
		SpatRaster classify_layers(std::vector<double> groups, unsigned nc, std::vector<double> id, SpatOptions &opt);

		std::vector<double> readSample(unsigned src, unsigned srows, unsigned scols);
		SpatRaster rotate(bool left, SpatOptions &opt);

		SpatRaster sampleRegularRaster(unsigned size);
		SpatRaster sampleRandomRaster(unsigned size, bool replace, unsigned seed);
		std::vector<std::vector<double>> sampleRegularValues(unsigned size);
		std::vector<std::vector<double>> sampleRandomValues(unsigned size, bool replace, unsigned seed);

		SpatRaster slope(unsigned neighbors, bool degrees, SpatOptions &opt);

		SpatRaster selRange(SpatRaster x, int z, int recycleby, SpatOptions &opt);

		SpatRaster shift(double x, double y, SpatOptions &opt);
		SpatRaster summary(std::string fun, bool narm, SpatOptions &opt);
		SpatRaster summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt);

		SpatRaster transpose(SpatOptions &opt);
		SpatRaster trig(std::string fun, SpatOptions &opt);
		SpatRaster trim(unsigned padding, SpatOptions &opt);
		std::vector<std::vector<double>> unique(bool bylayer, SpatOptions &opt);
		SpatRaster project1(std::string newcrs, std::string method, SpatOptions &opt);
		SpatRaster project2(SpatRaster &x, std::string method, SpatOptions &opt);
		void project3(SpatRaster &out, std::string method, SpatOptions &opt);

		SpatRaster resample1(SpatRaster &x, const std::string &method, SpatOptions &opt);
		void resample2(SpatRaster &out, const std::string &method, SpatOptions &opt);

#ifdef useGDAL
		bool open_gdal(GDALDatasetH &hDS, int src, SpatOptions &opt);
		bool create_gdalDS(GDALDatasetH &hDS, std::string filename, std::string driver, bool fill, double fillvalue, SpatOptions& opt);
		bool from_gdalMEM(GDALDatasetH hDS, bool set_geometry, bool get_values);
		bool as_gdalvrt(GDALDatasetH &hVRT, SpatOptions &opt);
		//bool as_gdalmem(GDALDatasetH &hVRT);
#endif

		SpatRaster to_memory_copy();
		bool to_memory();

		SpatRaster warp(SpatRaster x, const std::string &method, SpatOptions &opt);
		SpatRaster warpcrs(std::string x, const std::string &method, SpatOptions &opt);

		SpatRaster warper(SpatRaster x, std::string crs, std::string method, bool mask, SpatOptions &opt);
		//SpatRaster tester(bool geom);

		//SpatRaster warp_gdal(SpatRaster x, const std::string &method, SpatOptions &opt);
		//SpatRaster warp_gdal_crs(std::string x, const std::string &method, SpatOptions &opt);
		SpatDataFrame zonal(SpatRaster x, std::string fun, bool narm, SpatOptions &opt);

};

