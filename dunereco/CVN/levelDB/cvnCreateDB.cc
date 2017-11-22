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
#include "art/Utilities/FirstAbsoluteOrLookupWithDotPolicy.h"
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
//#include "CVN/func/TrainingData.h"
//#include "CVN/art/CaffeNetHandler.h"

// Caffe stuff
#include <leveldb/db.h>
#include <leveldb/write_batch.h>
#include <lmdb.h>
#define CPU_ONLY
// Suppress warnings originating in Caffe that we can't do anything about
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "caffe/caffe.hpp"
#pragma GCC diagnostic pop

#include "H5Cpp.h"



namespace po = boost::program_options;


enum LabelingMode
{
    kAll,  ///< Label all interaction types separately
    kNumu, ///< Label numu:1, else 0
    kNue,  ///< Label nue:1, else 0
    kNC,    ///< Label NC:1, else 0
    kEnergy ///< Label is conversion of fNuEnergy to int
};


class Config
{
public:
  Config(const fhicl::ParameterSet& pset):
    fOutputFormat (pset.get<std::string>("OutputFormat")),
    fTreeName (pset.get<std::string>("TreeName")),
    fTrainingBranchObjectName (pset.get<std::string>("TrainingDataBranchName")),
    fTestOutputDir (pset.get<std::string>("TestOutputDir")),
    fTrainOutputDir (pset.get<std::string>("TrainOutputDir")),
    fNTrainPerTest (pset.get<unsigned int>("NTrainPerTest")),
    fProgressInterval (pset.get<unsigned int>("ProgressInterval")),
    fErrorIfExists (pset.get<bool>("ErrorIfExists")),
    fSetLog (pset.get<bool>("SetLog")),
    fCreateIfMissing (pset.get<bool>("CreateIfMissing")),
    fWriteSync (pset.get<bool>("WriteSync")),
    fMaxKeyLength (pset.get<unsigned int>("MaxKeyLength")),
    fWriteBufferSize (pset.get<unsigned int>("WriteBufferSize")),
    fLabeling (pset.get<std::string>("Labeling")),
    fUseGeV (pset.get<bool>("UseGeV")),
    fWriteRegressionHDF5 (pset.get<bool>("WriteRegressionHDF5")),
    fRegressionHDF5NameTrain (pset.get<std::string>("RegressionHDF5NameTrain")),
    fRegressionHDF5NameTest (pset.get<std::string>("RegressionHDF5NameTest")),
    fMaxEnergyForLabel (pset.get<float>("MaxEnergyForLabel")),
    fNEvents (pset.get<unsigned int>("NEvents"))
  {
    if(!fLabeling.compare("all"))      fLabelingMode = kAll;
    if(!fLabeling.compare("numu"))     fLabelingMode = kNumu;
    if(!fLabeling.compare("nue"))      fLabelingMode = kNue;
    if(!fLabeling.compare("nc"))       fLabelingMode = kNC;
    if(!fLabeling.compare("energy"))   fLabelingMode = kEnergy;
  };


  std::string  fOutputFormat;
  std::string  fTreeName;
  std::string  fTrainingBranchObjectName;
  std::string  fTestOutputDir;
  std::string  fTrainOutputDir;

  /// Number of training examples for each test sample, e.g. 4 for 80/20 split
  unsigned int fNTrainPerTest;
  /// Number of examples in between progress updates (% complete)
  unsigned int fProgressInterval;

  bool          fErrorIfExists;
  bool          fSetLog;
  bool          fCreateIfMissing;
  bool          fWriteSync;
  unsigned int  fMaxKeyLength;
  unsigned int  fWriteBufferSize;

  std::string   fLabeling;
  unsigned int  fLabelingMode;

  bool          fUseGeV;

  /// Flag to control whether or not we write HDF5 regression features
  bool          fWriteRegressionHDF5;
  std::string   fRegressionHDF5NameTrain;
  std::string   fRegressionHDF5NameTest;
  float         fMaxEnergyForLabel;
  /// Limit the number of entries in the tree to consider
  unsigned int  fNEvents;
};

class OutputDB {
public:
  OutputDB(std::string sample, const Config& config);
  ~OutputDB();

  void Put(std::string &serializeKey, std::string  &serializeString);

private:
  leveldb::DB* fLevelDB;
  leveldb::WriteOptions fWriteOptions;

  MDB_env *mdb_env;
  MDB_dbi mdb_dbi;
  MDB_val mdb_key, mdb_data;
  MDB_txn *mdb_txn;
  
};

OutputDB::OutputDB(std::string sample, const Config& config) :
  fLevelDB(0),  mdb_env(0), mdb_txn(0) {

  std::string outputDir;
  if (sample=="test") 
    outputDir=config.fTestOutputDir;
  else
    outputDir=config.fTrainOutputDir;

  if (config.fOutputFormat=="LevelDB"){
    leveldb::Options fileOptions;
    fileOptions.error_if_exists   = config.fErrorIfExists;
    fileOptions.create_if_missing = config.fCreateIfMissing;
    fileOptions.write_buffer_size = config.fWriteBufferSize;

    fWriteOptions.sync = config.fWriteSync;

    if(!leveldb::DB::Open(fileOptions, outputDir, &fLevelDB).ok()) {
      std::cout << "Problem opening the database: "
		<< outputDir << std::endl;
      exit(1);
    }
  }

  else if (config.fOutputFormat=="LMDB") {
    mkdir(outputDir.c_str(),0777);
    mdb_env_create(&mdb_env);
    mdb_env_set_mapsize(mdb_env, 10737418240);
    mdb_env_open(mdb_env, outputDir.c_str(), 0, 0777);
    mdb_txn_begin(mdb_env, NULL, 0, &mdb_txn);
    mdb_dbi_open(mdb_txn,NULL, 0, &mdb_dbi);
  }

  else {
    std::cout << "Unrecognized output format " << config.fOutputFormat << std::endl;
    exit(1);
  }

}

OutputDB::~OutputDB() {
  if (mdb_txn)  mdb_txn_commit(mdb_txn);
}

void OutputDB::Put(std::string &serializeKey, std::string  &serializeString) {
  if (fLevelDB) {
    fLevelDB->Put(fWriteOptions, serializeKey, serializeString);
  } //end if LevelDB
  else {//it must be LMDB
    mdb_data.mv_size=serializeString.size();
    mdb_data.mv_data=reinterpret_cast<void*>(&serializeString[0]);
    mdb_key.mv_size=serializeKey.size();
    mdb_key.mv_data=reinterpret_cast<void*>(&serializeKey[0]);
    if ( mdb_put(mdb_txn,mdb_dbi,&mdb_key,&mdb_data,0)!= MDB_SUCCESS){
      std::cout<< "ERROR: Events not loaded correctly" <<std::endl;
    }//end if put fails
  }//end if LMDB
}//end OutputDB::Put

char ConvertToCharDune(float n, bool fSetLog)
{
  float peCorrChunk;
  float truncateCorr;
  float centreScale = 0.7;
  if(fSetLog){
    float scaleFrac=(log(n)/log(1000));
    truncateCorr= ceil(centreScale*scaleFrac*255.0);
      }
  else{
    peCorrChunk = (1000.) / 255.0;
    truncateCorr = ceil((n)/(peCorrChunk));
  }
  if (truncateCorr > 255) return (char)255;
  else return (char)truncateCorr;
}

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

  unsigned int entries = chain.GetEntries();
  if(config.fNEvents < entries){
    entries = config.fNEvents;
  }
  if(entries <= 0){
    std::cout << "Error: Input tree has no entries." << std::endl;
    exit(4);
  }

  std::cout << "- Will process " << entries << " from the input tree." << std::endl;

  OutputDB TrainDB("train",config);
  OutputDB  TestDB( "test",config);

  char* key = new char[config.fMaxKeyLength];
  std::string serializeString;

  //need to shuffle entries...
  std::srand ( unsigned ( std::time(0) ) );
  std::vector<unsigned int> shuffled;
  for (unsigned int i = 0; i < entries; ++i)
  {
    shuffled.push_back(i);
  }

  std::random_shuffle( shuffled.begin(), shuffled.end() );

  // Figure out the size of the train and test samples
  // Call a 'block' a particular set of one test and nTrainPerTest train
  unsigned int blockSize = config.fNTrainPerTest + 1;
  // number of test is the number of blocks, using integer division
  unsigned int nTest = 1 + entries / blockSize;
  // number of training samples is number of blocks times train per test
  unsigned int nTrain    = entries / blockSize * config.fNTrainPerTest;
  // Add on the entries from the last, potentially partial block, minus test
  if (entries % blockSize) nTrain += entries % blockSize - 1;

  // Create an array to hold regression features.
  const unsigned int nRegressionFeatures = 2; // Currently nuEnergy, lepEnergy
//  float regressionDataTest[nTest][nRegressionFeatures];
//  float regressionDataTrain[nTrain][nRegressionFeatures];

  int** regressionDataTest = new int*[nTest];
  for(unsigned int i = 0; i < nTest; ++i) {regressionDataTest[i] = new int[nRegressionFeatures];}
  int** regressionDataTrain = new int*[nTrain];
  for(unsigned int i = 0; i < nTrain; ++i) {regressionDataTrain[i] = new int[nRegressionFeatures];}

  int iTrain = 0;
  int iTest  = 0;

  ////hdf5////

  const char saveFilePath[] = "test.h5";
  const hsize_t ndims = 2;
  const hsize_t ncols = 3;

  hid_t file = H5Fcreate(saveFilePath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  std::cout << "- File created" << std::endl;

  hsize_t dims[ndims] = {0, ncols};
  hsize_t max_dims[ndims] = {H5S_UNLIMITED, ncols};
  hid_t file_space = H5Screate_simple(ndims, dims, max_dims);
  std::cout << "- Dataspace created" << std::endl;

  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_layout(plist, H5D_CHUNKED);
  hsize_t chunk_dims[ndims] = {2, ncols};
  H5Pset_chunk(plist, ndims, chunk_dims);
  std::cout << "- Property list created" << std::endl;

  //hid_t dset = H5Dcreate(file, "dset1", H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(file, "dset1", H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  std::cout << "- Dataset 'dset1' created" << std::endl;

  H5Pclose(plist);
  H5Sclose(file_space);

  for(unsigned int iEntry = 0; iEntry < entries; ++iEntry)
  {
    unsigned int entry = shuffled[iEntry];
    chain.GetEntry(entry);

    // auto maxwireelement_0= std::max_element(fPMap_fPEX.begin(), fPMap_fPEX.end());                                                                                     
    // std::cout<<"maxwire 0: "<<*maxwireelement_0<<std::endl;                                                                                            
    // auto maxwireelement_1= std::max_element(fPMap_fPEY.begin(), fPMap_fPEY.end());                                                                                          
    // std::cout<<"maxwire 1: "<<*maxwireelement_1<<std::endl;                                                                                                                 
    // auto maxwireelement_2= std::max_element(fPMap_fPEZ.begin(), fPMap_fPEZ.end());                         
    // std::cout<<"maxwire 2: "<<*maxwireelement_2<<std::endl;   
        
    caffe::Datum datum; //= PixelMapToDatum(data->fPMap, config.fUseGeV);

    ///////////
    char* pixels = NULL;
    int channels(3), planes(0), cells(0);
    datum.set_channels(channels);
    planes = 500;//fPMap_fNWire;
    cells  = fPMap_fNTdc;
    datum.set_height(planes);
    datum.set_width(cells);
    pixels = new char[channels*planes*cells];
    for (int iChan = 0; iChan < channels; ++iChan)
    {
      for (int iPlane = 0; iPlane < planes; ++iPlane)
      {
        for (int iCell = 0; iCell < cells; ++iCell)
        {
          int i = iCell + cells*(iPlane + planes*iChan);
          float val =0.;
          //std::cout<<iChan<<std::endl;
          if(iChan == 0 ){
            //std::cout<<iChan<<std::endl;
            val=fPMap_fPEX.at(iCell + cells*iPlane);
          }
          if(iChan == 1 ){
            val=0.;//fPMap_fPEY.at(iCell + cells*iPlane);
          }
          if(iChan == 2 ){
            val=fPMap_fPEZ.at(iCell + cells*iPlane);
          }
          // if(val!=0){
          //   std::cout<<val<<std::endl;
          // }
          char pix = ConvertToCharDune(val,config.fSetLog);
          pixels[i] = pix;
        }
      }
    }

    int maxPlaneWithCharge=0;
    for (int iChan = 0; iChan < channels; ++iChan)
    {
      for (int iPlane = 0; iPlane < planes; ++iPlane)
        {
          for (int iCell = 0; iCell < cells; ++iCell)
            {
              float val =0.;
              if(iChan == 1 ){
                val=fPMap_fPEY.at(iCell + cells*iPlane);
                if(val>0){
                  maxPlaneWithCharge=iPlane;
                }
              }
            }
        }
    }
    
    int newMapEdge=500;//1000;
    if(maxPlaneWithCharge>500){//1000){
      newMapEdge=maxPlaneWithCharge+1;
    }
    
    for (int iChan = 0; iChan < channels; ++iChan)
      {
        for (int iPlane = newMapEdge; iPlane > (newMapEdge-planes); --iPlane)
          {
            for (int iCell = 0; iCell < cells; ++iCell)
              {
                int planeCount=planes-iPlane;
                int i = iCell + cells*(planeCount + planes*iChan);
                float val =0.;
                //std::cout<<iChan<<std::endl;
                // if(iChan==1&&iCell==1){
                //   std::cout<<"Real Plane"<<iPlane<<std::endl;
                //   std::cout<<"Fake Plane"<<planeCount<<std::endl;
                // }
                //std::cout<<iCell<<std::endl;
                if(iChan == 1 ){
                  val=fPMap_fPEY.at(iCell + cells*iPlane);
                  char pix = ConvertToCharDune(val,config.fSetLog);
                  pixels[i] = pix;
                }
              }
          }
      }
    
    
    datum.set_data(pixels, channels*planes*cells);
    delete[] pixels;
    //////////
    datum.set_label(fInt);

    datum.SerializeToString(&serializeString);

    if(iEntry % (blockSize))
    {
      snprintf(key, config.fMaxKeyLength, "%08lld", (long long int)iTrain);
      std::string serializeKey(key);

      TrainDB.Put(serializeKey,serializeString);

      regressionDataTrain[iTrain][0] = 1.;
      regressionDataTrain[iTrain][1] = 1.;
      iTrain += 1;

      ////hdf5////
      hsize_t nlines = 1;
      float *buffer = new float[nlines * ncols];
      float **b = new float*[nlines];
      for (hsize_t i = 0; i < nlines; ++i){
        b[i] = &buffer[i * ncols];
      }

      b[0][0] = 0.1;
      b[0][1] = 0.2;
      b[0][2] = 0.3;

    }
    else
    {
      snprintf(key, config.fMaxKeyLength, "%08lld", (long long int)iTest);
      std::string serializeKey(key);

      TestDB.Put(serializeKey,serializeString);

      regressionDataTest[iTest][0] = 1.;
      regressionDataTest[iTest][1] = 1.;
      iTest += 1;
    }
    if(not (iEntry % config.fProgressInterval))
      std::cout << "Fraction complete: "
                << iEntry / (float)entries << std::endl;

  }

  if (config.fWriteRegressionHDF5)
  {

    H5::FloatType type(H5::PredType::IEEE_F32LE);
    std::cout << "Writing HDF5 regression output : "
              << config.fRegressionHDF5NameTest << std::endl;
    H5::H5File h5FileTest(config.fRegressionHDF5NameTest, H5F_ACC_TRUNC );
    hsize_t shape[2];
    shape[0] = nTest;
    shape[1] = nRegressionFeatures;
    H5::DataSpace spaceTest(2, shape);

    H5::DataSet datasetTest = h5FileTest.createDataSet("regression",
                                                    type,
                                                    spaceTest);

    datasetTest.write(regressionDataTest, type);

    std::cout << "Writing HDF5 regression output : "
              << config.fRegressionHDF5NameTrain << std::endl;
    H5::H5File h5FileTrain(config.fRegressionHDF5NameTrain, H5F_ACC_TRUNC );
    shape[0] = nTrain;
    H5::DataSpace spaceTrain(2, shape);

    H5::DataSet datasetTrain = h5FileTrain.createDataSet("regression",
                                                    type,
                                                    spaceTrain);

    datasetTrain.write(regressionDataTrain, type);

  }

  // if(fInt)  delete fInt;
  // if(fPMap_fNWire)  delete fPMap_fNWire;
  // if(fPMap_fNTdc)  delete fPMap_fNTdc;
  // if(fPMap_fPEX)  delete fPMap_fPEX;
  // if(fPMap_fPEY)  delete fPMap_fPEY;
  // if(fPMap_fPEZ)  delete fPMap_fPEZ;

  // Clear up
  delete key;
  for(unsigned int i = 0; i < nTest; ++i) {
    delete [] regressionDataTest[i];
  }
  delete [] regressionDataTest;
  for(unsigned int i = 0; i < nTrain; ++i) {
    delete [] regressionDataTrain[i];
  }
  delete [] regressionDataTrain;

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

  art::FirstAbsoluteOrLookupWithDotPolicy policy("FHICL_FILE_PATH");

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



