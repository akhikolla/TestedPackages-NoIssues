#include <Rcpp.h>
#include "vcfReader.hpp"

using namespace Rcpp;

class Rvcf : public VcfReader{
 private:
    vector <string> rvcf_chrom;
    vector <int> rvcf_pos;
    Rcpp::DataFrame resultList_;

    void gatherChromPos() {
        for (size_t i = 0; i < this->position_.size(); i++) {
            for (size_t ii = 0; ii < this->position_[i].size(); ii++) {
                rvcf_chrom.push_back(this->chrom_[i]);
                rvcf_pos.push_back(this->position_[i][ii]);
            }
        }

        // VQSLOD to come
        resultList_ = Rcpp::DataFrame::create(_("CHROM") = this->rvcf_chrom,
                                              _("POS") = this->rvcf_pos,
                                              _("refCount") =  this->refCount,
                                              _("altCount") = this->altCount);
    }

 public:
    explicit Rvcf(string fileName) : VcfReader(fileName){
        this->finalize(); // Finalize after remove variantlines
        this->gatherChromPos();
    }
    ~Rvcf() {}
    Rcpp::DataFrame info(){
        return resultList_;
    }
};

// [[Rcpp::plugins(cpp11)]]

//' @title Extract VCF information
//'
//' @description Extract VCF information
//'
//' @param filename VCF file name.
//'
//' @seealso
//' \itemize{
//'   \item \code{extractCoverageFromVcf}
//'   \item \code{extractCoverageFromTxt}
//' }
//'
//' @return A dataframe list with members of haplotypes, proportions and log likelihood of the MCMC chain.
//' \itemize{
//'   \item \code{CHROM} SNP chromosomes.
//'   \item \code{POS} SNP positions.
//'   \item \code{refCount} reference allele count.
//'   \item \code{altCount} alternative allele count.
//' }
//'
//' @export
//'
//' @examples
//' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
//' vcf = extractVcf(vcfFile)
//'
// [[Rcpp::export]]
Rcpp::DataFrame extractVcf(std::string filename) {
    Rvcf vcf(filename);
    return vcf.info();
}
