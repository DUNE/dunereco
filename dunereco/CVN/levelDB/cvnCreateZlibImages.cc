////////////////////////////////////////////////////////////////////////////////////////////////////
//// Authors:     Saul Alonso-Monsalve, saul.alonso.monsalve@cern.ch
//// Zlib image generator.
////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// std library
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <algorithm>

// Boost, for program options
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/algorithm/string/predicate.hpp"

// ART/fcl stuff
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT stuff
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"

// CVN stuff
#include "dune/CVN/func/CVNImageUtils.h"
//#include "CVN/art/CaffeNetHandler.h"

#include "zlib.h" // compression algorithm

namespace po = boost::program_options;

class Config
{
public:
  Config(const fhicl::ParameterSet& pset):
    fTreeName (pset.get<std::string>("TreeName")),
    fOutputDir (pset.get<std::string>("OutputDir")),
    fSetLog (pset.get<bool>("SetLog")),
    fNEvents (pset.get<unsigned int>("NEvents")),
    fPlaneLimit (pset.get<unsigned int>("PlaneLimit")),
    fTDCLimit (pset.get<unsigned int>("TDCLimit")),
    fReverseViews(pset.get<std::vector<bool> >("ReverseViews"))
  {
  };

  std::string  fTreeName;
  std::string  fTrainingBranchObjectName;
  std::string  fOutputDir;

  /// Number of training examples for each test sample, e.g. 4 for 80/20 split
  /// Number of examples in between progress updates (% complete)

  bool          fSetLog;

  /// Flag to control whether or not we write HDF5 regression features
  /// Limit the number of entries in the tree to consider
  unsigned int  fNEvents;
  /// Limit the number of wires in the output image
  int fPlaneLimit;
  /// Limit the number of TDCs in the output image
  int fTDCLimit;
  /// Views to reverse
  std::vector<bool> fReverseViews;
};

void fill(const Config& config, std::string input)
{

  TChain chain(config.fTreeName.c_str());

  if (boost::ends_with(input,".list")) {
    std::ifstream list_file(input.c_str());
    if (!list_file.is_open()) {
      std::cout << "Could not open " << input << std::endl;
      exit(1);
    }

    std::string ifname;
    while (list_file>>ifname)
      chain.Add(ifname.c_str());
      
  }//end if list file

  else if  (boost::ends_with(input,".root")) {
    chain.Add(input.c_str());
  }//end if root file

  chain.SetMakeClass(1);

  int fInt;
  UInt_t          fPMap_fNWire;
  UInt_t          fPMap_fNTdc;
  std::vector<float>   fPMap_fPEX;
  std::vector<float>   fPMap_fPEY;
  std::vector<float>   fPMap_fPEZ;

  chain.SetBranchAddress("fInt", &fInt);
  chain.SetBranchAddress("fPMap.fNWire", &fPMap_fNWire);
  chain.SetBranchAddress("fPMap.fNTdc", &fPMap_fNTdc);
  chain.SetBranchAddress("fPMap.fPEX", &fPMap_fPEX);
  chain.SetBranchAddress("fPMap.fPEY", &fPMap_fPEY);
  chain.SetBranchAddress("fPMap.fPEZ", &fPMap_fPEZ);

  float fNuEnergy = -1;
  float fLepEnergy = -1;
  float fRecoNueEnergy = -1;
  float fRecoNumuEnergy = -1;
  float fEventWeight = -1;

  int  fNuPDG = -1;
  int  fNProton = -1;
  int  fNPion = -1;
  int  fNPizero = -1;
  int  fNNeutron = -1;

  chain.SetBranchAddress("fNuEnergy", &fNuEnergy);
  chain.SetBranchAddress("fLepEnergy", &fLepEnergy);
  chain.SetBranchAddress("fRecoNueEnergy", &fRecoNueEnergy);
  chain.SetBranchAddress("fRecoNumuEnergy", &fRecoNumuEnergy);
  chain.SetBranchAddress("fEventWeight", &fEventWeight);
 
  chain.SetBranchAddress("fNuPDG", &fNuPDG);
  chain.SetBranchAddress("fNProton", &fNProton);
  chain.SetBranchAddress("fNPion", &fNPion);
  chain.SetBranchAddress("fNPizero", &fNPizero);
  chain.SetBranchAddress("fNNeutron", &fNNeutron);

  unsigned int entries = chain.GetEntries();
  if(config.fNEvents < entries){
    entries = config.fNEvents;
  }
  if(entries <= 0){
    std::cout << "Error: Input tree has no entries." << std::endl;
    exit(4);
  }

  std::cout << "- Will process " << entries << " from the input tree." << std::endl;

  std::srand ( unsigned ( std::time(0) ) );
  std::vector<unsigned int> shuffled;
  for (unsigned int i = 0; i < entries; ++i)
  {
    shuffled.push_back(i);
  }

  if(entries > chain.GetEntries()){
    entries = chain.GetEntries();
  }

  std::string image_path = config.fOutputDir + "/images/";
  std::string info_path = config.fOutputDir + "/info/";

  // create image path
  if (mkdir(image_path.data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
  {
    if( errno == EEXIST ) {
      // alredy exists
    }else{
      // error
      std::cout << "cannot create folder error:" << strerror(errno) << std::endl;
      exit(1);
    }
  }
  // create info path
  if (mkdir(info_path.data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
  {
    if( errno == EEXIST ) {
      // alredy exists
    }else{
      // error
      std::cout << "cannot create folder error:" << strerror(errno) << std::endl;
      exit(1);
    }
  }

  for(unsigned int iEntry = 0; iEntry < entries; ++iEntry)
  {
    unsigned int entry = shuffled[iEntry];
    chain.GetEntry(entry);

    unsigned int nViews = 3;

    // Create a CVNImageUtils object and use it to produce the pixels. The arguments
    // define how large we want the output image to be
    cvn::CVNImageUtils imageUtils(config.fPlaneLimit,config.fTDCLimit,nViews);
    // Since we don't have a PixelMap object, we need to tell it how big it is
    imageUtils.SetPixelMapSize(fPMap_fNWire,fPMap_fNTdc);
    
    std::vector<unsigned char> pixelArray(nViews * config.fPlaneLimit * config.fTDCLimit,0);

    imageUtils.SetLogScale(config.fSetLog);
    imageUtils.SetViewReversal(config.fReverseViews);
    imageUtils.ConvertChargeVectorsToPixelArray(fPMap_fPEX, fPMap_fPEY, fPMap_fPEZ, pixelArray);

    //std::cout << "fNuEnergyi: " << fNuEnergy << std::endl;
    //std::cout << "fRecoNueEnergy: " << fRecoNueEnergy << std::endl;
    //std::cout << "fRecoNumuEnergy: " << fRecoNumuEnergy << std::endl;
    //std::cout << "fEventWeight: " << fEventWeight << std::endl;

    std::cout << "[DEBUG] entry " << entry+1 << " out of " << entries << std::endl;
    std::cout << "[DEBUG] file: " << chain.GetCurrentFile()->GetName() << std::endl;
    //std::cout << "[DEBUG] channels: " << nViews << std::endl;   
    //std::cout << "[DEBUG] planes: " << config.fPlaneLimit << std::endl;   
    //std::cout << "[DEBUG] cells: " << config.fTDCLimit << std::endl;   
    //std::cout << "[DEBUG] channels*planes*cells: " << channels*planes*cells << std::endl; 
    
    std::cout << "[DEBUG] label: " << fInt << std::endl;
    ulong srcLen = nViews * config.fPlaneLimit * config.fTDCLimit; // pixelArray length
    ulong destLen = compressBound(srcLen);     // calculate size of the compressed data               
    char* ostream = (char *) malloc(destLen);  // allocate memory for the compressed data

    int res = compress((Bytef *) ostream, &destLen, (Bytef *) &pixelArray[0], srcLen); // compress pixels 

    // destLen is now the size of actuall buffer needed for compression
    // you don't want to uncompress whole buffer later, just the used part
    //

    // Buffer error

    if(res == Z_BUF_ERROR){
        std::cout << "[DEBUG] Buffer was too small!" << std::endl;
    }

    // Memory error
    
    else if(res ==  Z_MEM_ERROR){
        std::cout << "[DEBUG] Not enough memory for compression!" << std::endl;
    }

    // Compression ok 
   
    else{
        std::cout << "[DEBUG] Compression successful" << std::endl;

        // Create output files 

        std::string image_file_name = image_path + input + std::to_string(iEntry) + ".gz";
        std::string info_file_name = info_path + input + std::to_string(iEntry) + ".info";

        while(1){

            std::ofstream image_file (image_file_name, std::ofstream::binary);
            std::ofstream info_file  (info_file_name);

            if(image_file.is_open() && info_file.is_open()){

                // Write compressed data to file

                image_file.write(ostream, destLen);

                image_file.close(); // close file

                // Write records to file

                // Category

                info_file << fInt << std::endl;

                // Energy

                info_file << fNuEnergy << std::endl;
                info_file << fLepEnergy << std::endl;
                info_file << fRecoNueEnergy << std::endl;
                info_file << fRecoNumuEnergy << std::endl;
                info_file << fEventWeight << std::endl;

                // Topology

                info_file << fNuPDG << std::endl;
                info_file << fNProton << std::endl;
                info_file << fNPion << std::endl;
                info_file << fNPizero << std::endl;         
                info_file << fNNeutron;        

                info_file.close(); // close file

                std::cout << "[DEBUG] Done" << std::endl;
                break;

            }
            else{

                if(image_file.is_open()){
                    image_file.close(); // close file
                }
                else{
                    std::cout << "[DEBUG] Unable to open file: " << image_file_name << std::endl;
                }

                if(info_file.is_open()){
                    info_file.close(); // close file
                }
                else{
                    std::cout << "[DEBUG] Unable to open file: " << info_file_name << std::endl;
                }
            }
       }
        //else std::cout << "[DEBUG] Unable to open files " << image_file_name << " and " << info_file_name << std::endl;
    }
    
    free(ostream);  // free allocated memory

  }
}


po::variables_map getOptions(int argc, char*  argv[], std::string& config,
                                                      std::string& input)
{

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("config,c", po::value<std::string>(&config)->required(),
                                                "configuration file")
    ("input,i", po::value<std::string>(&input)->required(),
                                                "Input data in ROOT file.");
  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

  }
  catch(po::error& e)
  {
    std::cout  << "ERROR: " << e.what() << std::endl;
    exit(1);
  }


  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(1);
  }

  return vm;
}



fhicl::ParameterSet getPSet(std::string configPath)
{

  cet::filepath_first_absolute_or_lookup_with_dot policy("FHICL_FILE_PATH");

  // parse a configuration file; obtain intermediate form
  fhicl::intermediate_table tbl;
  fhicl::parse_document(configPath, policy, tbl);

  // convert to ParameterSet
  fhicl::ParameterSet pset;
  fhicl::make_ParameterSet(tbl, pset);

  return pset;

}

int main(int argc, char* argv[])
{

  std::string configPath, inputPath, outputPath, logPath;
  po::variables_map vm = getOptions(argc, argv, configPath, inputPath);

  Config config(getPSet(configPath));


  fill(config, inputPath);

  return 0;

}



