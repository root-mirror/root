#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void TreeReaderSimple() {
	TH1F *myHistogram = new TH1F("h1","ntuple",100,-4,4);

	TFile *myFile = TFile::Open("hsimple.root");
	TTreeReader myHSimpleReader("ntuple", myFile);

	TTreeReaderValue<Float_t> myPx(myHSimpleReader, "px");
	TTreeReaderValue<Float_t> myPy(myHSimpleReader, "py");

	while (myHSimpleReader.Next()){
		myHistogram->Fill(*myPx + *myPy);
	}

	myHistogram->Draw();
}
