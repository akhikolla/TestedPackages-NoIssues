#include <Rcpp.h>
//#include "spatRaster.h"
#include "spatRasterMultiple.h"
#include <memory> //std::addressof
//#include "gdal_priv.h"
//# include "gdal_info.h"


Rcpp::List getBlockSizeR(SpatRaster* r, unsigned n, double frac) { 
	SpatOptions opt;
	opt.set_blocksizemp(n);
	opt.set_memfrac(frac);
    BlockSize bs = r->getBlockSize(opt);
	Rcpp::List L = Rcpp::List::create(Rcpp::Named("row") = bs.row, Rcpp::Named("nrows") = bs.nrows, Rcpp::Named("n") = bs.n);
	return(L);
}


bool sameObject(SpatRaster* a, SpatRaster* b) {
//	Rcpp::Rcout << a->source[0].time[0] << std::endl;
//	size_t n = a->source[0].time.size();
//	Rcpp::Rcout << n << std::endl;
//	Rcpp::Rcout << a->nlyr() << std::endl;
//	Rcpp::Rcout << a->source[0].time[n-1] << std::endl;
	return a == b;
}	



Rcpp::List getDataFrame(SpatDataFrame* v) {
	unsigned n = v->ncol();
	Rcpp::List out(n);	
	if (n == 0) {
		return(out);
	} 

	std::vector<std::string> nms = v->names;
	std::vector<unsigned> itype = v->itype;
	for (size_t i=0; i < n; i++) {
		if (itype[i] == 0) {
			out[i] = v->getD(i);
		} else if (itype[i] == 1) {
			out[i] = v->getI(i);
		} else {
			out[i] = v->getS(i);
		}
	}	
	out.names() = nms;
	// todo: deal with NAs in int and str
	return out;
//  Rcpp::df is nice, but no of variables is <= 20, 
//  and no "stringsAsFactors"=false
//	Rcpp::DataFrame result(out);
//	result.attr("names") = v->names();
//	return result;
}	


Rcpp::List getVectorAttributes(SpatVector* v) {
	SpatDataFrame df = v->df;
	Rcpp::List lst = getDataFrame(&df);
	return lst;
}

Rcpp::List getRasterAttributes(SpatRaster* x) {
	Rcpp::List lst;
	if (x->nlyr() > 0) {
		SpatDataFrame df = x->source[0].atts[0];
		lst = getDataFrame(&df);
	}
	return lst;
}



Rcpp::DataFrame getGeometry(SpatVector* v) {
	SpatDataFrame df = v->getGeometryDF();

	Rcpp::DataFrame out = Rcpp::DataFrame::create(
			Rcpp::Named("id") = df.iv[0], 
			Rcpp::Named("part") = df.iv[1], 
			Rcpp::Named("x") = df.dv[0],
			Rcpp::Named("y") = df.dv[1],
			Rcpp::Named("hole") = df.iv[2]
	);
	return out;
}


RCPP_EXPOSED_CLASS(SpatMessages)
RCPP_EXPOSED_CLASS(SpatOptions)
RCPP_EXPOSED_CLASS(SpatExtent)
RCPP_EXPOSED_CLASS(SpatCategories)
RCPP_EXPOSED_CLASS(SpatDataFrame)
RCPP_EXPOSED_CLASS(SpatWindow)
RCPP_EXPOSED_CLASS(RasterSource)
RCPP_EXPOSED_CLASS(SpatRaster)
RCPP_EXPOSED_CLASS(SpatRasterCollection)
RCPP_EXPOSED_CLASS(SpatRasterStack)
RCPP_EXPOSED_CLASS(SpatVector)


RCPP_MODULE(spat){

    using namespace Rcpp;


    class_<SpatExtent>("SpatExtent")
		.constructor()
		.constructor<double, double, double, double>()
		.property("vector", &SpatExtent::asVector)		
		.property("valid", &SpatExtent::valid)		
		.method("align", &SpatExtent::align, "align")		
		.method("intersect", &SpatExtent::intersect, "intersect")		
		.method("unite", &SpatExtent::unite, "unite")		
		.method("as.points", &SpatExtent::asPoints, "as.points")		
		.method("ceil",  &SpatExtent::ceil,  "ceil")		
		.method("compare", &SpatExtent::compare, "compare")		
		.method("floor", &SpatExtent::floor, "floor")		
		.method("round", &SpatExtent::round, "round")		
		.method("union", &SpatExtent::unite, "union")		
	;	

/*
    class_<SpatWindow>("SpatWindow")
		.field_readonly("full_extent", &SpatWindow::full_extent)
		.field_readonly("full_nrow", &SpatWindow::full_nrow)
		.field_readonly("full_ncol", &SpatWindow::full_ncol)
		.field_readonly("off_row", &SpatWindow::off_row)
		.field_readonly("off_col", &SpatWindow::off_col)
		.field_readonly("expand", &SpatWindow::expand)
	;
*/

    class_<SpatMessages>("SpatMessages")
		.constructor()
		//.field("success", &SpatMessages::success)		
		.field("has_error", &SpatMessages::has_error)		
		.field("has_warning", &SpatMessages::has_warning)		
		.method("getError", &SpatMessages::getError)		
		.method("getWarnings", &SpatMessages::getWarnings)
		.method("getMessage", &SpatMessages::getMessage)
	;	
	
    class_<SpatOptions>("SpatOptions")
		.constructor()
		.method("deepcopy", &SpatOptions::deepCopy, "deepCopy")
		.property("tempdir", &SpatOptions::get_tempdir, &SpatOptions::set_tempdir )
		.property("memfrac", &SpatOptions::get_memfrac, &SpatOptions::set_memfrac )
		.property("filenames", &SpatOptions::get_filenames, &SpatOptions::set_filenames )
		.property("filetype", &SpatOptions::get_filetype, &SpatOptions::set_filetype )
		.property("datatype", &SpatOptions::get_datatype, &SpatOptions::set_datatype )
		.property("verbose", &SpatOptions::get_verbose, &SpatOptions::set_verbose )
		.property("NAflag", &SpatOptions::get_NAflag, &SpatOptions::set_NAflag )
		.property("overwrite", &SpatOptions::get_overwrite, &SpatOptions::set_overwrite )
		.property("progress", &SpatOptions::get_progress, &SpatOptions::set_progress)
		.field("ncopies", &SpatOptions::ncopies, "ncopies")

		.property("def_filetype", &SpatOptions::get_def_filetype, &SpatOptions::set_def_filetype )
		.property("def_datatype", &SpatOptions::get_def_datatype, &SpatOptions::set_def_datatype )
		.property("def_verbose", &SpatOptions::get_def_verbose, &SpatOptions::set_def_verbose )

		//.property("def_bandorder", &SpatOptions::get_def_bandorder, &SpatOptions::set_def_bandorder )

		.property("todisk", &SpatOptions::get_todisk, &SpatOptions::set_todisk)
		.field("messages", &SpatOptions::msg, "messages")
		.field("gdal_options", &SpatOptions::gdal_options, "gdal_options")
		.field("names", &SpatOptions::names, "names")
	//	.property("overwrite", &SpatOptions::set_overwrite, &SpatOptions::get_overwrite )
		//.field("gdaloptions", &SpatOptions::gdaloptions)		
	;

    class_<SpatDataFrame>("SpatDataFrame")
		.constructor()

		.field_readonly("itype", &SpatDataFrame::itype)
		.field_readonly("iplace", &SpatDataFrame::iplace)
		
		.property("names", &SpatDataFrame::get_names, &SpatDataFrame::set_names)
		.property("nrow", &SpatDataFrame::nrow, &SpatDataFrame::resize_rows, "nrow")
		.property("ncol", &SpatDataFrame::ncol, &SpatDataFrame::resize_cols, "ncol")
				
		.method("add_column_double", (bool (SpatDataFrame::*)(std::vector<double>, std::string name))( &SpatDataFrame::add_column))
		.method("add_column_long", (bool (SpatDataFrame::*)(std::vector<long>, std::string name))( &SpatDataFrame::add_column))
		.method("add_column_string", (bool (SpatDataFrame::*)(std::vector<std::string>, std::string name))( &SpatDataFrame::add_column))

		.method("remove_column", (bool (SpatDataFrame::*)(std::string field))( &SpatDataFrame::remove_column))
		.method("remove_column", (bool (SpatDataFrame::*)(int i))( &SpatDataFrame::remove_column))
				
		.method("cbind", &SpatDataFrame::cbind)
		.method("rbind", &SpatDataFrame::rbind)
		.method("values", &getDataFrame, "get data.frame")
		.method("unique", &SpatDataFrame::unique)
		.field("messages", &SpatDataFrame::msg, "messages")		
	;

    class_<SpatCategories>("SpatCategories")
		.constructor()
		.field_readonly("levels", &SpatCategories::levels, "levels")
		.field_readonly("labels", &SpatCategories::labels, "labels")		
	;	


    class_<SpatVector>("SpatVector")
		.constructor()	
		.constructor<SpatExtent, std::string>()
		.constructor<std::vector<std::string>>()

		.method("add_column_empty", (void (SpatVector::*)(unsigned dtype, std::string name))( &SpatVector::add_column))
		.method("add_column_double", (bool (SpatVector::*)(std::vector<double>, std::string name))( &SpatVector::add_column))
		.method("add_column_long", (bool (SpatVector::*)(std::vector<long>, std::string name))( &SpatVector::add_column))
		.method("add_column_string", (bool (SpatVector::*)(std::vector<std::string>, std::string name))( &SpatVector::add_column))
		.method("remove_column", (bool (SpatVector::*)(std::string field))( &SpatVector::remove_column))
		.method("remove_column", (bool (SpatVector::*)(int i))( &SpatVector::remove_column))
		.method("remove_holes", &SpatVector::remove_holes, "remove holes")		

		.method("area", &SpatVector::area, "area")		
		.method("as_lines", &SpatVector::as_lines, "as_lines")
		.method("as_points", &SpatVector::as_points, "as_points")
		.method("couldBeLonLat", &SpatVector::could_be_lonlat, "couldBeLonLat") 
		.method("get_crs", &SpatVector::getSRS)
		
		.method("set_crs", (bool (SpatVector::*)(std::string crs))( &SpatVector::setSRS))
		//.method("p	rj", &SpatVector::getPRJ)
		
		.method("distance_self", (SpatDataFrame (SpatVector::*)())( &SpatVector::distance))
		.method("distance_other", (SpatDataFrame (SpatVector::*)(SpatVector, bool))( &SpatVector::distance))
//		.method("distance_other2", (SpatDataFrame (SpatVector::*)(SpatVector))( &SpatVector::distance2))
		.method("extent", &SpatVector::getExtent, "extent")		
		.method("getDF", &getVectorAttributes, "get attributes")
		.method("getGeometry", &getGeometry, "getGeometry")
		.method("getGeometryWKT", &SpatVector::getGeometryWKT, "getGeometryWKT")
		.method("isLonLat", &SpatVector::is_lonlat, "isLonLat")
		.method("length", &SpatVector::length, "length")		

		.field("messages", &SpatVector::msg, "messages")
		.property("names", &SpatVector::get_names, &SpatVector::set_names)
		.method("nrow", &SpatVector::nrow, "nrow")		
		.method("ncol", &SpatVector::ncol, "ncol")		
		.method("project", &SpatVector::project, "project")
		.method("read", &SpatVector::read, "read")		
		.method("setGeometry", &SpatVector::setGeometry, "setGeometry")
		.method("size", &SpatVector::size, "size")		
		.method("subset_cols", ( SpatVector (SpatVector::*)(std::vector<int>))( &SpatVector::subset_cols ))
		.method("subset_rows", ( SpatVector (SpatVector::*)(std::vector<int>))( &SpatVector::subset_rows ))				
		.method("type", &SpatVector::type, "type")		

		.method("write", &SpatVector::write, "write")	
		
		.method("aggregate", &SpatVector::aggregate, "aggregate")	
		.method("disaggregate", &SpatVector::disaggregate, "disaggregate")	
		.method("buffer", &SpatVector::buffer, "buffer")	
		.method("is_valid", &SpatVector::is_valid, "is_valid")	
		.method("make_valid", &SpatVector::make_valid, "make_valid")	
		.method("shift", &SpatVector::shift, "shift")	
#ifdef useGEOS
		.method("buffer2", &SpatVector::buffer2, "buffer2")		
		.method("intersect", &SpatVector::intersect, "intersect")		
#endif
	;


    class_<RasterSource>("RasterSource")	
		.field_readonly("time", &RasterSource::time)
		//.field_readonly("memory", &RasterSource::memory)
	//	.field_readonly("filename", &RasterSource::filename)
		//.field_readonly("driver", &RasterSource::driver)
		//.field_readonly("nrow", &RasterSource::nrow)
		//.field_readonly("ncol", &RasterSource::ncol)
		//.field_readonly("nlyr", &RasterSource::nlyr)
//		.field_readonly("extent", &RasterSource::extent)
//		.field_readonly("hasWindow", &RasterSource::hasWindow)
//		.field_readonly("window", &RasterSource::window)
		//.field_readonly("layers", &RasterSource::layers)
		//.field_readonly("nlyrfile", &RasterSource::nlyrfile)
		//.field_readonly("flipped", &RasterSource::flipped)
		//.field_readonly("rotated", &RasterSource::rotated)
//		.field_readonly("parameters_changed", &RasterSource::parameters_changed)
	;	


    class_<SpatRaster>("SpatRaster")
		.constructor()
	 // .constructor<std::string, int>()
	    .constructor<std::vector<std::string>, std::vector<int>, std::vector<std::string>, std::string>()
		.constructor<std::vector<unsigned>, std::vector<double>, std::string>()


		.field("name", &SpatRaster::name)

		.method("sources_to_disk", &SpatRaster::sources_to_disk, "sources_to_disk")
		.method("mem_needs", &SpatRaster::mem_needs, "mem_needs")

		.method("spatinit", &SpatRaster::gdalogrproj_init, "init")
		
		.method("addSource", &SpatRaster::addSource, "addSource")
		.method("combineSources", &SpatRaster::combineSources, "combineSources")
		.method("compare_geom", &SpatRaster::compare_geom, "compare_geom")
		.method("couldBeLonLat", &SpatRaster::could_be_lonlat, "couldBeLonLat") 
		.method("deepcopy", &SpatRaster::deepCopy, "deepCopy")
		.method("get_crs", &SpatRaster::getSRS)

		.method("set_crs", (bool (SpatRaster::*)(std::string crs))( &SpatRaster::setSRS))
		//.field_readonly("prj", &SpatRaster::prj)
		.property("extent", &SpatRaster::getExtent, &SpatRaster::setExtent )

		.method("setWindow", &SpatRaster::setWindow, "")
		.method("removeWindow", &SpatRaster::removeWindow, "")
		.method("hasWindow", &SpatRaster::hasWindow, "")

		.method("getRasterAtt", &getRasterAttributes, "get attributes")
			
		.property("filenames", &SpatRaster::filenames )

		.method("hasAttributes", &SpatRaster::hasAttributes, "hasAttributes")
		.method("getAttributes", &SpatRaster::getAttributes, "getAttributes")
		.method("setAttributes", &SpatRaster::setAttributes, "setAttributes")
		.method("createAttributes", &SpatRaster::createAttributes, "createAttributes")
		.method("hasCategories", &SpatRaster::hasCategories, "hasCategories")
		.method("getCategories", &SpatRaster::getCategories, "getCategories")
		.method("setCategories", &SpatRaster::setCategories, "setCategories")
		.method("createCategories", &SpatRaster::createCategories, "createCategories")
		.method("hasColors", &SpatRaster::hasColors, "hasColors")
		.method("getColors", &SpatRaster::getColors, "getColors")
		
		.property("hasRange", &SpatRaster::hasRange )
		.property("hasValues", &SpatRaster::hasValues )
		.property("inMemory", &SpatRaster::inMemory )
		.method("isLonLat", &SpatRaster::is_lonlat, "isLonLat")
		.method("isGlobalLonLat", &SpatRaster::is_global_lonlat, "isGlobalLonLat") 

		.property("names", &SpatRaster::getNames)
		.property("long_names", &SpatRaster::getLongNames)
		.method("set_long_names", &SpatRaster::setLongNames)

		.property("hasTime", &SpatRaster::hasTime)
		.property("time", &SpatRaster::getTime)
		.property("timestep", &SpatRaster::getTimeStep)
		.method("settime", &SpatRaster::setTime)
		//.property("timestr", &SpatRaster::getTimeStr)

		.property("depth", &SpatRaster::getDepth)
		.method("set_depth", &SpatRaster::setDepth)
		.property("units", &SpatRaster::getUnit)
		.method("set_units", &SpatRaster::setUnit)

		.method("size", &SpatRaster::size, "size")		
		.method("nrow", &SpatRaster::nrow, "nrow")		
		.method("ncol", &SpatRaster::ncol, "ncol")		
		.method("nsrc", &SpatRaster::nsrc, "nsrc" )	
		.field("messages", &SpatRaster::msg, "messages")
		.method("nlyrBySource", &SpatRaster::nlyrBySource, "nlyrBySource" )		
		.method("nlyr", &SpatRaster::nlyr, "nlyr" )
		.property("origin", &SpatRaster::origin)
		.property("range_min", &SpatRaster::range_min )
		.property("range_max", &SpatRaster::range_max )
		.property("res", &SpatRaster::resolution)
				
// only if RasterSource is exposed
		.field_readonly("source", &SpatRaster::source )

		.method("collapse_sources", &SpatRaster::collapse_sources, "collapse_sources" )
		
		.method("dense_extent", &SpatRaster::dense_extent, "dense_extent")
		.method("setNames", &SpatRaster::setNames, "setNames" )
		.method("setTime", &SpatRaster::setTime, "setTime" )
		.method("setDepth", &SpatRaster::setDepth, "setDepth" )
		.method("setUnit", &SpatRaster::setUnit, "setUnit" )
		.method("set_resolution", &SpatRaster::setResolution, "set resolution")
		.method("subset", &SpatRaster::subset, "subset")
				
		.method("cellFromXY", ( std::vector<double> (SpatRaster::*)(std::vector<double>,std::vector<double>) )( &SpatRaster::cellFromXY ))
		.method("getCells", &SpatRaster::getCells, "getCells")
		
		.method("cellFromRowCol", ( std::vector<double> (SpatRaster::*)(std::vector<int_64>,std::vector<int_64>) )( &SpatRaster::cellFromRowCol ))
		.method("cellFromRowColCombine", ( std::vector<double> (SpatRaster::*)(std::vector<int_64>,std::vector<int_64>) )( &SpatRaster::cellFromRowColCombine ))
		.method("yFromRow", ( std::vector<double> (SpatRaster::*)(std::vector<int_64>&) )( &SpatRaster::yFromRow ))
		.method("xFromCol", ( std::vector<double> (SpatRaster::*)(std::vector<int_64>&) )( &SpatRaster::xFromCol ))
		.method("colFromX", ( std::vector<int_64> (SpatRaster::*)(std::vector<double>&) )( &SpatRaster::colFromX ))
		.method("rowFromY", ( std::vector<int_64> (SpatRaster::*)(std::vector<double>&) )( &SpatRaster::rowFromY ))
		.method("xyFromCell", ( std::vector< std::vector<double> > (SpatRaster::*)(std::vector<double>&) )( &SpatRaster::xyFromCell ))
		.method("rowColFromCell", ( std::vector< std::vector<int_64> > (SpatRaster::*)(std::vector<double>) )( &SpatRaster::rowColFromCell ))
		.method("readStart", &SpatRaster::readStart, "readStart") 
		.method("readStop", &SpatRaster::readStop, "readStop") 
		.method("readValues", &SpatRaster::readValues, "readValues")	
		.method("getValues", &SpatRaster::getValues, "getValues")
		.method("getBlockSize", &getBlockSizeR)
		.method("same", &sameObject)
		.method("setValues", &SpatRaster::setValues)
		//.method("replaceValues", &SpatRaster::replace)
		.method("setRange", &SpatRaster::setRange, "setRange")
		.method("writeStart", &SpatRaster::writeStart, "writeStart") 
		.method("writeStop", &SpatRaster::writeStop, "writeStop") 
		.method("writeValues", &SpatRaster::writeValues, "writeValues") 
		.method("writeRaster", &SpatRaster::writeRaster, "writeRaster")
		.method("canProcessInMemory", &SpatRaster::canProcessInMemory, "canProcessInMemory")
		.method("chunkSize", &SpatRaster::chunkSize, "chunkSize")
		
		.method("adjacent", &SpatRaster::adjacent, "adjacent")
		.method("aggregate", &SpatRaster::aggregate, "aggregate")
		.method("align", &SpatRaster::align, "align")
		.method("apply", &SpatRaster::apply, "apply")
		.method("rapply", &SpatRaster::rapply, "rapply")
		.method("rappvals", &SpatRaster::rappvals, "rappvals")
		.method("arith_rast", ( SpatRaster (SpatRaster::*)(SpatRaster, std::string, SpatOptions&) )( &SpatRaster::arith ))
		.method("arith_numb", ( SpatRaster (SpatRaster::*)(std::vector<double>, std::string, bool, SpatOptions&) )( &SpatRaster::arith ))
		.method("rst_area", &SpatRaster::rst_area, "rst_area")
		.method("sum_area", &SpatRaster::sum_area, "sum_area")
		.method("area_by_value", &SpatRaster::area_by_value, "area_by_value")
		
		.method("as_points", &SpatRaster::as_points, "as_points")
		.method("as_polygons", &SpatRaster::as_polygons, "as_polygons")
		.method("polygonize", &SpatRaster::polygonize, "polygonize")
		
		.method("atan2", &SpatRaster::atan_2, "atan2")

		.method("bilinearValues", &SpatRaster::bilinearValues, "bilin")

		.method("boundaries", &SpatRaster::edges, "edges")
		.method("buffer", &SpatRaster::buffer, "buffer")
		.method("gridDistance", &SpatRaster::gridDistance, "gridDistance")
		.method("rastDistance", ( SpatRaster (SpatRaster::*)(SpatOptions&) )( &SpatRaster::distance), "rastDistance")		
		.method("vectDistance", ( SpatRaster (SpatRaster::*)(SpatVector, SpatOptions&) )( &SpatRaster::distance), "vectDistance")		
		.method("clamp", &SpatRaster::clamp, "clamp")
		.method("classify", ( SpatRaster (SpatRaster::*)(std::vector<double>, unsigned, unsigned, bool, bool, SpatOptions&) )( &SpatRaster::reclassify), "reclassify")		
		//.method("source_collapse", &SpatRaster::collapse, "collapse")
		.method("classify", ( SpatRaster (SpatRaster::*)(std::vector<double>, unsigned, unsigned, bool, bool, SpatOptions&) )( &SpatRaster::reclassify), "reclassify")	
		.method("selRange", &SpatRaster::selRange, "selRange")
		.method("separate", &SpatRaster::separate, "separate")
		
		.method("cover", &SpatRaster::cover, "cover")
		.method("crop", &SpatRaster::crop, "crop")
		.method("cum", &SpatRaster::cum, "cum")
		.method("disaggregate", &SpatRaster::disaggregate, "disaggregate")
		.method("expand", &SpatRaster::extend, "extend")
		.method("extractCell", &SpatRaster::extractCell, "extractCell")
		.method("extractVector", &SpatRaster::extractVector, "extractVector")
		.method("flip", &SpatRaster::flip, "flip")
		.method("focal", &SpatRaster::focal, "focal")
		.method("focalValues", &SpatRaster::focal_values, "focalValues")
		.method("count", &SpatRaster::count, "count")
		.method("freq", &SpatRaster::freq, "freq")
		.method("geometry", &SpatRaster::geometry, "geometry")

		.method("get_aggregates", &SpatRaster::get_aggregates, "get_aggregates")
		.method("get_aggregate_dims", &SpatRaster::get_aggregate_dims2, "get_aggregate_dims")
		.method("global", &SpatRaster::global, "global")
		.method("global_weighted_mean", &SpatRaster::global_weighted_mean, "global weighted mean")
				
		.method("initf", ( SpatRaster (SpatRaster::*)(std::string, bool, SpatOptions&) )( &SpatRaster::init ), "init fun")
		.method("initv", ( SpatRaster (SpatRaster::*)(double, SpatOptions&) )( &SpatRaster::init ), "init value")
		.method("is_in", &SpatRaster::is_in, "isin")
		.method("is_in_cells", &SpatRaster::is_in_cells, "isincells")
		.method("isnan", &SpatRaster::isnan, "isnan")
		.method("isfinite", &SpatRaster::isfinite, "isfinite")
		.method("isinfinite", &SpatRaster::isinfinite, "isinfinite")
		.method("logic_rast", ( SpatRaster (SpatRaster::*)(SpatRaster, std::string, SpatOptions&) )( &SpatRaster::logic ))
		.method("logic_numb", ( SpatRaster (SpatRaster::*)(bool, std::string, SpatOptions&) )( &SpatRaster::logic ))
		.method("mask_raster", ( SpatRaster (SpatRaster::*)(SpatRaster, bool, double, double, SpatOptions&) )( &SpatRaster::mask), "mask raster")
		.method("mask_vector", ( SpatRaster (SpatRaster::*)(SpatVector, bool, double, SpatOptions&) )( &SpatRaster::mask), "mask vector")	
		.method("math", &SpatRaster::math, "math")
		.method("math2", &SpatRaster::math2, "math2")
		.method("modal", &SpatRaster::modal, "modal")
		.method("quantile", &SpatRaster::quantile, "quantile")
		//.method("rasterize", &SpatRaster::rasterize, "rasterize")
		.method("rasterize", &SpatRaster::rasterize, "rasterize")
		.method("reverse", &SpatRaster::reverse, "reverse")
		.method("rotate", &SpatRaster::rotate, "rotate")
		.method("sampleRegularRaster", &SpatRaster::sampleRegularRaster, "sampleRegular")	
		.method("sampleRegularValues", &SpatRaster::sampleRegularValues, "sampleValues")	
		.method("sampleRandomRaster", &SpatRaster::sampleRandomRaster, "sampleRandom")	
		.method("sampleRandomValues", &SpatRaster::sampleRandomValues, "sampleValues")	
		.method("shift", &SpatRaster::shift, "shift")
		.method("slope", &SpatRaster::slope, "slope")
		.method("summary", &SpatRaster::summary, "summary")
		.method("summary_numb", &SpatRaster::summary_numb, "summary_numb")
		.method("transpose", &SpatRaster::transpose, "transpose")
		.method("trig", &SpatRaster::trig, "trig")
		.method("trim", &SpatRaster::trim, "trim")
		.method("unique", &SpatRaster::unique, "unique")

		.method("rectify", &SpatRaster::rectify, "rectify")
		.method("stretch", &SpatRaster::stretch, "stretch")
		.method("warp", &SpatRaster::warper, "warper")
		.method("zonal", &SpatRaster::zonal, "zonal")			
	;

    class_<SpatRasterCollection>("SpatRasterCollection")
		.constructor()
		.field("messages", &SpatRasterCollection::msg, "messages")		
		.field_readonly("x", &SpatRasterCollection::x)
		.method("length", &SpatRasterCollection::size, "size")
		.method("resize", &SpatRasterCollection::resize, "resize")
		.method("add", &SpatRasterCollection::push_back, "push_back")	
		.method("merge", &SpatRasterCollection::merge, "merge")
		//.method("mosaic", &SpatRasterCollection::mosaic, "mosaic")
	;
	
    class_<SpatRasterStack>("SpatDataSet")
		.constructor()
	    .constructor<std::string, std::vector<int>, bool>()
	    .constructor<SpatRaster, std::string>()
		.field("messages", &SpatRasterStack::msg, "messages")

		.method("readStart", &SpatRasterStack::readStart, "readStart")
		.method("readStop", &SpatRasterStack::readStop, "readStop")
		.method("nsds", &SpatRasterStack::nsds, "")
		.method("ncol", &SpatRasterStack::ncol, "")
		.method("nrow", &SpatRasterStack::nrow, "")
		.method("getSRS", &SpatRasterStack::getSRS, "")
		.property("names", &SpatRasterStack::getnames, &SpatRasterStack::setnames)
		.method("add", &SpatRasterStack::push_back, "")
		.method("resize", &SpatRasterStack::resize, "")
		.method("summary", &SpatRasterStack::summary, "summary")
		.method("summary_numb", &SpatRasterStack::summary_numb, "summary_numb")
		.method("getsds", &SpatRasterStack::getsds, "")
		.method("replace", &SpatRasterStack::replace, "")
		.method("subset", &SpatRasterStack::subset, "")
		.method("collapse", &SpatRasterStack::collapse , "")
		.method("extractCell", &SpatRasterStack::extractCell, "extractCell")
		.method("extractVector", &SpatRasterStack::extractVector, "extractVector")
	;
}


/*
SpatRaster SQRT() {
	SpatRaster r = *this;
	std::transform(r.values.begin(), r.values.end(), r.values.begin(), (double(*)(double)) sqrt);
	return r;
}

SpatRaster SQRTfree(SpatRaster* g) {
	SpatRaster r = *g;
	std::transform(r.values.begin(), r.values.end(), r.values.begin(), (double(*)(double)) sqrt);
	return r;
}
*/

