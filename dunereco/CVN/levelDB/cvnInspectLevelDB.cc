#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


// Boost, for program options
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"


#include "leveldb/db.h"

#define CPU_ONLY
// Suppress warnings originating in Caffe that we can't do anything about
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "caffe/caffe.hpp"
#pragma GCC diagnostic pop

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;
namespace po = boost::program_options;

po::variables_map getOptions(int argc, char*  argv[], std::string& dir)
{

// Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help", "produce help message")
  ("dir,d", po::value<std::string>(&dir)->required(),
    "leveldb directory");
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



int main(int argc, char* argv[])
{

  std::string directory;
  po::variables_map vm = getOptions(argc, argv, directory);

  leveldb::DB* db;
  leveldb::Options options;
  leveldb::Status status = leveldb::DB::Open(options, directory, &db);

  if (!status.ok()) return 1;

  ofstream outFile("test.txt");

  TH1D* hLabels = new TH1D("hLabels", ";Event Type", 15, 0, 15);
  TH1D* hPE = new TH1D("hPE", "Scaled PE", 256, 0, 256);
  vector<TH2D*> XVec;
  vector<TH2D*> YVec;
  vector<TH2D*> ZVec;

  caffe::Datum datum;

  int total(0);
  leveldb::Iterator* itCount = db->NewIterator(leveldb::ReadOptions());
  for (itCount->SeekToFirst(); itCount->Valid(); itCount->Next())
  {
    total++;
  }
  std::cout << "total: " << total << std::endl;

  int count(0);
  leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next())
  {
    if (count < 10000)
    {
      datum.ParseFromString( it->value().ToString() );
      int channels = datum.channels();
      int height   = datum.height();
      int width    = datum.width();
      int label    = datum.label();
      
      //std::cout<<"channels:"<<channels<<std::endl;
      //std::cout<<"height:"<<height<<std::endl;
      //std::cout<<"width:"<<width<<std::endl;
      //std::cout<<"label:"<<label<<std::endl;

      const char* pixels = datum.data().data();
      for (int iX = 0; iX < channels*height*width; ++iX)
      {
        hPE->Fill( (uint8_t)pixels[iX] );
        hLabels->Fill(label);
      }
      if (count < 100)
      {
        outFile << "key " << it->key().ToString() << endl;
        outFile << "count " << count << endl;
        outFile << "label " << label << endl;

        TString hName = "PixelMap_";
        hName += count;
        hName += "_type";
        hName += label;

        TH2D* xMap = new TH2D(hName + TString("_X"), ";Wire;TDC",
         height, 0, height, width, 0, width);
        XVec.push_back(xMap);
        TH2D* yMap = new TH2D(hName + TString("_Y"), ";Wire;TDC",
         height, 0, height, width, 0, width);
        YVec.push_back(yMap);
        TH2D* zMap = new TH2D(hName + TString("_Z"), ";Wire;TDC",
         height, 0, height, width, 0, width);
        ZVec.push_back(zMap);

        for (int iPix = 0; iPix < channels*height*width; ++iPix)
        {
          int iCell = iPix%width;
          int iPlane = (iPix/width)%height;
          int iChan =  (iPix/width)/height;
          int pixVal = (uint8_t)pixels[iPix];
          if (iChan == 0)
            {
              xMap->SetBinContent(iPlane+1, iCell+1, pixVal);
            }
          if (iChan == 1)
            {
              yMap->SetBinContent(iPlane+1, iCell+1, pixVal);
            }
          if ( iChan == 2 )
            {
              zMap->SetBinContent(iPlane+1, iCell+1, pixVal);
            }
          
          if (pixVal > 0)
            {
              if (iChan == 0)
                {
              outFile << "X " << iCell << "  " << iPlane << endl;
                }
              else if (iChan == 1)
                {
                  outFile << "Y " << iCell << "  " << iPlane << endl;
                }
              else
                {
                  outFile << "Z " << iCell << "  " << iPlane << endl;
                }
            outFile << "[" << pixVal << "]" << endl;
            }
        }
      }
    }
    else
      {
        break;
      }
    ++count;
  }

  TFile* fOut = new TFile("test.root", "recreate");
  fOut->WriteTObject(hPE);
  fOut->WriteTObject(hLabels);
  //gStyle->SetOptStat(0);
  for (size_t i = 0; i < XVec.size(); ++i)
  {
    TString cName = "cPixelMap_";
    cName += i;
    TCanvas* canv = new TCanvas(cName, cName);
    canv->Divide(3,1);
    canv->cd(1);
    std::cout<<XVec[i]->Integral()<<std::endl;
    XVec[i]->Draw("colz");
    canv->cd(2);
    std::cout<<YVec[i]->Integral()<<std::endl;
    YVec[i]->Draw("colz");
    canv->cd(3);
    std::cout<<ZVec[i]->Integral()<<std::endl;
    ZVec[i]->Draw("colz");
    fOut->WriteTObject(canv);
  }


  return 0;
}
