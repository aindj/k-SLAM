/* Copyright 2012 David Ainsworth
 *
 * This file is part of SLAM
 *
 * SLAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with SLAM.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "main.h"
#include <boost/program_options.hpp>
//#include "Tests.h"
using namespace SLAM;
namespace po = boost::program_options;
int main(int argc, char *argv[]) {
  for(int i=0;i<argc;i++){
    if(i!=0)
      commandLine.append(" ");
    commandLine.append(argv[i]);

  }
  po::options_description allowed("Allowed options");
  po::options_description hidden("Hidden options");
  po::options_description all("All options");
  uint32_t numReads, numReadsAtOnce;
  std::string outFileName, samFileName;
  allowed.add_options()("help", "produce help message")(
      "db", po::value<std::string>(), "SLAM database file which reads "
                                      "will be aligned against")(
      "min-alignment-score",
      po::value<uint32_t>(&scoreThreshold)->default_value(0),
      "alignment score cutoff")(
          "score-fraction-threshold",
          po::value<double>(&scoreFractionThreshold)->default_value(0.95),
          "screen alignments with scores < this*top score")      (
      "match-score", po::value<uint32_t>(&match)->default_value(2),
      "match score")("mismatch-penalty",
                     po::value<uint32_t>(&misMatch)->default_value(3),
                     "mismatch penalty (positive)")(
      "gap-open", po::value<uint32_t>(&gapOpen)->default_value(5),
      "gap opening penalty (positive)")(
      "gap-extend", po::value<uint32_t>(&gapExtend)->default_value(2),
      "gap extend penalty (positive)")(
      "num-reads", po::value<uint32_t>(&numReads)->default_value(UINT32_MAX),
      "Number of reads from R1/R2 File to align")(
      "num-reads-at-once",
      po::value<uint32_t>(&numReadsAtOnce)->default_value(10000000),
      "Reduce RAM usage by only analysing \"arg\" reads at once, "
      "this will increase execution time")(
      "output-file", po::value<std::string>(&outFileName)->default_value(""),
      "write to this file instead of stdout")(
      "sam-file", po::value<std::string>(&samFileName)->default_value(""),
      "write SAM output to this file")(
          "num-alignments", po::value<uint32_t>(&numSAMAlignments)->default_value(10),
          "Number of alignments to report in SAM file")("sam-xa",
          "only output primary alignment lines, use XA field for secondary alignments")
          ("just-align",
                    "only perform alignments, not metagenomics")
//          ("report-cigar","report cigar string in output")
                                       (
      "no-pseudo-assembly", "do not link alignments together");
  ;
  hidden.add_options()("server", "server")(
      "input-file", po::value<std::vector<std::string> >(), "input file")(
      "parse-genbank", "parses genbank files to produce a SLAM database")(
          "parse-fasta", "parses FASTA files to produce a SLAM database")(
      "parse-taxonomy",
      "parses NCBI names.dmp and nodes.dmp file to produce a SLAM "
      "taxonomy database, usage SLAM --parse-taxonomy names.dmp"
      " nodes.dmp --output-file taxonomy")(
      "alignment-only", "Only perform alignment, not metagenomic "
                        "analysis, output in SAM format");
  all.add(allowed).add(hidden);
  po::positional_options_description p;
  p.add("input-file", -1);
  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
  po::notify(vm);
//  if (vm.count("report-cigar")) reportCigar = true;
  if (vm.count("sam-xa")) SAMXA = true;
  if (vm.count("just-align")) justAlign = true;
  if (vm.count("no-pseudo-assembly")) performPseudoAssembly = false;
  if (vm.count("help") || argc == 1) {
    std::cout << "Usage\t"
              << "SLAM [option] --db=DATABASE R1FILE R2FILE\n";
    std::cout << "\tAlign paired reads from R1FILE and R2FILE against DATABASE "
                 "and perform metagenomic analysis\n";
    std::cout << "or\t"
              << "SLAM [option] --db=DATABASE R1FILE\n";
    std::cout << "\tAlign reads from R1FILE against DATABASE "
                 "and perform metagenomic analysis\n";
    std::cout << allowed << "\n";
    return 1;
  }
  if (vm.count("parse-genbank")) {
    log("Parsing Genbank");
    createIndexFromGBFF(vm["input-file"].as<std::vector<std::string> >(),
                        outFileName);
    return 0;
  }
  if (vm.count("parse-fasta")) {
    log("Parsing FASTA");
    createIndexFromFASTA(vm["input-file"].as<std::vector<std::string> >(),
                        outFileName);
    return 0;
  }
  if (vm.count("parse-taxonomy")) {
    log("Parsing taxonomy");
    auto inputs = vm["input-file"].as<std::vector<std::string> >();
    if (inputs.size() != 2) {
      std::cout << "Provide names.dmp and nodes.dmp\n";
      return 1;
    }
    TaxonomyDB taxDB;
    taxDB.writeTaxonomyIndex(outFileName, inputs[0], inputs[1]);
    return 0;
  }
  if (vm.count("input-file")) {
    std::vector<std::string> inputs =
        vm["input-file"].as<std::vector<std::string> >();
    if (inputs.size() == 1) {
      if (numReadsAtOnce != UINT32_MAX) {
        metagenomicAnalysis_Low_Mem(inputs[0], std::string(),
                                    vm["db"].as<std::string>(), outFileName,
                                    samFileName, numReadsAtOnce, numReads);
      } else {
        metagenomicAnalysis(inputs[0], std::string(),
                            vm["db"].as<std::string>(), outFileName,
                            samFileName, numReads);
      }
    }
    if (inputs.size() == 2) {
      if (numReadsAtOnce != UINT32_MAX) {
        metagenomicAnalysis_Low_Mem(inputs[0], inputs[1],
                                    vm["db"].as<std::string>(), outFileName,
                                    samFileName, numReadsAtOnce, numReads);
      } else {
        metagenomicAnalysis(inputs[0], inputs[1], vm["db"].as<std::string>(),
                            outFileName, samFileName, numReads);
      }
    }
  }
  //  return 0;
  //  if (argc > 3 && std::string(argv[1]) == "--server") {
  //    serverMode(argv[2], std::stoi(argv[3]));
  //    return 0;
  //  }
  //  if (std::string(argv[1]) == "--parse-fasta" && argc > 3) {
  //    GenbankIndex index;
  //    log("Parsing fasta");
  //    index.parseFASTA(argv[2], argv[3]);
  //  }
  return 0;
}
