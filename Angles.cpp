#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/util/FileUtil.h>
#include <ausa/util/DynamicBranchVector.h>
#include <TROOT.h>
#include <cxxopts.hpp>
#include <gsl/gsl_matrix.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <filesystem>
#include <cmath>

#ifndef EXEC_VERSION
#define EXEC_VERSION "unstable"
#endif


using namespace std;
using namespace AUSA;
using namespace rapidjson;

int main(int argc, char *argv[]) {
  // parse options
  string setup_path = "setup.json", target_path = "target.json", matcher_path = "matcher.json";
  cxxopts::Options options("Angles", "Angles - Generates distribution of angles for each pixel from AUSAlib setup, target and matcher files");
  options.add_options()
      ("s,setup", "Path to setup file", cxxopts::value<string>()->default_value(setup_path))
      ("t,target", "Path to target file", cxxopts::value<string>()->default_value(target_path))
      ("x,xlocation", "Location of source", cxxopts::value<string>()->default_value("X location of target specified from target file"))
      ("y,ylocation", "Location of source", cxxopts::value<string>()->default_value("Y location of target specified from target file"))
      ("i,implantation", "Implantation depth (nm) - note: target's normal vector assumed parallel to z axis", cxxopts::value<string>()->default_value("Half of target thickness"))
      ("n,iterations", "Max number of iterations to generate angle distribution", cxxopts::value<string>()->default_value("100"))
      ("m,matcher", "Path to matcher file", cxxopts::value<string>()->default_value(matcher_path))
      ("d,detector", "Detector(s) - give multiple detectors with ',' as separator, e.g. '-d U1,U2,U3'", cxxopts::value<vector<string>>())
      ("h,help", "Print usage")
      ("v,version", "Print version")
  ;
  auto result = options.parse(argc, argv);
  if (result.count("help")) {
    cout << options.help() << endl;
    exit(0);
  }
  if (result.count("version")) {
    cout << EXEC_VERSION << endl;
    exit(0);
  }
  if (!result.count("detector")) {
    cout << "Must specificy at least one detector." << endl
         << "Try running " << argv[0] << " -h" << endl;
    exit(-1);
  }
  setup_path = result.count("setup") ? result["setup"].as<string>() : setup_path;
  target_path = result.count("target") ? result["target"].as<string>() : target_path;
  matcher_path = result.count("matcher") ? result["matcher"].as<string>() : matcher_path;

  // load json files
  auto setup = JSON::readSetupFromJSON(setup_path);
  auto target = JSON::readTargetFromJSON(target_path);
  FILE *fp = fopen(matcher_path.c_str(), "rb"); // 2023-02-17: https://rapidjson.org/md_doc_stream.html
  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));
  Document matcher_cfg;
  matcher_cfg.ParseStream(is);

  // get/set implantation depth and define point of particle emission from there
  double implantation_depth = result.count("implantation") ?
      stod(result["implantation"].as<string>()) : 1e6*target.getThickness()/2.; // nm
  implantation_depth /= 1e6;

  double xlocation = result.count("xlocation") ?
      stod(result["xlocation"].as<string>()) : 0.; // nm
  double ylocation = result.count("ylocation") ?
      stod(result["ylocation"].as<string>()) : 0.; // nm

  TVector3 center(xlocation,ylocation,-0.3);

  TVector3 origin = center + (target.getThickness()/2. - implantation_depth)*target.getNormal();
  TVector3 ori = target.getCenter() + (target.getThickness()/2. - implantation_depth)*target.getNormal();
  
  int n = result.count("iterations") ?
      stod(result["iterations"].as<string>()) : 100;

  // logging
  string pwd = filesystem::current_path();
  time_t now = time(nullptr);
  char* datetime = ctime(&now);
  vector<string> args;
  args.reserve(argc);
  for (int i = 0; i < argc; i++) {
    args.emplace_back(argv[i]);
  }
  cout << "# " << datetime
       << "# Output created from within " << pwd << " with the following command" << endl << "# ";
  for (const string &s: args) {
    cout << s << " ";
  }
  cout << endl
       << "# Setup:                         " << setup_path << endl
       << "# Target:                        " << target_path << endl
       << "# Implantation depth:            " << implantation_depth*1e6 << " nm in target of thickness " << 1e6*target.getThickness() << " nm" << endl
       << "# Nr. of max iterations:         " << n << endl
       << "# Matcher:                       " << matcher_path << endl
       << "# Target center=                 (" << target.getCenter().X() << ", " << target.getCenter().Y() << ", " << target.getCenter().Z() << ")" << endl
       << "# Origin: Center of implantation=(" << origin.X() << ", " << origin.Y() << ", " << origin.Z() << ")" << endl
       << "# Origo: Center of target=       (" << ori.X() << ", " << ori.Y() << ", " << ori.Z() << ")" << endl;

  // map disabled strips of various detectors
  map<string, pair<vector<int>, vector<int>>> disabled_strips;
  // this is not pretty, but neither is rapidjson's documentation...
  for (const auto& detector : result["detector"].as<vector<string>>()) {
    auto det_string = detector.c_str();
    vector<int> disabled_front, disabled_back;
    for (auto const &side: {"front", "back"}) {
      if (matcher_cfg.HasMember("DSD")) {
        if (matcher_cfg["DSD"].HasMember(det_string)) {
          if (matcher_cfg["DSD"][det_string].HasMember(side)) {
            if (matcher_cfg["DSD"][det_string][side].HasMember("disable")) {
              if (matcher_cfg["DSD"][det_string][side]["disable"].IsArray()) { // disabled strips can be specified as one or several integers in an array...
                for (auto &el: matcher_cfg["DSD"][det_string][side]["disable"].GetArray()) {
                  if (side == (string) "front") {
                    disabled_front.emplace_back(el.GetInt());
                  } else {
                    disabled_back.emplace_back(el.GetInt());
                  }
                }
              } else { // ...or just as a single integer
                if (side == (string) "front") {
                  disabled_front.emplace_back(matcher_cfg["DSD"][det_string][side]["disable"].GetInt());
                } else {
                  disabled_back.emplace_back(matcher_cfg["DSD"][det_string][side]["disable"].GetInt());
                }
              }
            }
          }
        }
      }
    }
    disabled_strips.insert(make_pair(det_string, make_pair(disabled_front, disabled_back)));
  }

  // find solid angle coverages of each pixel for various detectors and print them
  for (const auto& detector : result["detector"].as<vector<string>>()) {
    auto detptr = setup->getDSSD(detector);
    UInt_t fN = detptr->frontStripCount();
    double max_solid = 0.;

    gsl_matrix *solid_angles = gsl_matrix_alloc(fN, fN);
    auto disabled_front = disabled_strips[detector.c_str()].first;
    auto disabled_back = disabled_strips[detector.c_str()].second;
    for (size_t i = 0; i < solid_angles->size1; i++) {
      for (size_t j = 0; j < solid_angles->size2; j++) {
        double psa;
        double psa_origo;
        if (find(disabled_front.begin(), disabled_front.end(), i + 1) != disabled_front.end()
            || find(disabled_back.begin(), disabled_back.end(), j + 1) != disabled_back.end()) {
          psa = 0.;
        } else {
          psa = detptr->getPixelSolidAngle(i + 1, j + 1, origin);
          psa_origo = detptr->getPixelSolidAngle(i + 1, j + 1, ori); //define max solid with respect to 1000 events from center of target
          if (psa_origo>max_solid){
            max_solid = psa_origo;
          }
        }
        gsl_matrix_set(solid_angles, i, j, psa);
      }
    }

    cout << "# Solid angles of pixels of " << detector.c_str() << endl;
    cout << "# Max solid angle is  " << max_solid << endl;

    cout << "# ";
    for (int j = 0; j < solid_angles->size2; j++) {
      cout << j + 1 << "\t";
    }
    cout << endl;
    for (size_t i = 0; i < solid_angles->size1; i++) {
      cout << "#";
      for (size_t j = 0; j < solid_angles->size2; j++) {
        cout << gsl_matrix_get(solid_angles, i, j);
        if (j + 1 != solid_angles->size2) {
          cout << "\t";
        } else {
          cout << "\t " << i + 1 << endl;
        }
      }
    }

    cout << "#Angle distribution of uniform hits over the detector accounted for solid angle of" << detector.c_str() << endl;
    cout << "# " << "FI" << "\t" << "BI" << "\t" << "Angle" << endl;
    for (size_t i = 0; i < solid_angles->size1; i++) { ///for front strip index
        for (size_t j = 0; j < solid_angles->size2; j++) { //for backstrip index
            double solid = gsl_matrix_get(solid_angles, i, j); //get the solid angle for the given i,j combination
            double relative_solid = solid/max_solid; //gives the relative size of the solid angle 
            double N = floor(n * relative_solid); //number of hits relative to the most prominent pixel
            for (int k=0; k<N; k++){//we loop over each "hit"
                TVector3 pos = detptr->getUniformPixelPosition(i, j);
                TVector3 dir = (pos-ori).Unit();
                auto angle = dir.Angle(-detptr->getNormal());
                cout << i+1 << "\t" << j+1 << "\t" << angle << "\t" << N << "\t" << solid << "\t" << pos.X() << "\t" << pos.Y() << "\t" << pos.Z() << endl;
            }//for k in range(N)
        }//for front index
    }//for back index
  }

  exit(0);
}