/*
 * delaunayTriangulation.cxx
 *
 *  Created on: Jun 30, 2014
 *      Author: dfunke
 *
 *  This test creates a TGraph2D, fills it with 3 points and then performs
 *  the Delaunay triangulation of them.
 *
 *  Because of bug ROOT-XXX the data structures were not properly initialized and no triangle was found
 */

#include "TGraph2D.h"
#include "TGraphDelaunay2D.h"

#include "delaunayTriangulation_bug.h"

#include "TVirtualPad.h"

#include <iostream>

void printDelaunay(const TGraphDelaunay2D & gd){

	auto graph = gd.GetGraph2D();

	for(const auto & triangle : gd){
		printf("[%u](%f,%f) - [%u](%f,%f) - [%u](%f,%f)\n",
				triangle.idx[0], graph->GetX()[triangle.idx[0]], graph->GetY()[triangle.idx[0]],
				triangle.idx[1], graph->GetX()[triangle.idx[1]], graph->GetY()[triangle.idx[1]],
				triangle.idx[2], graph->GetX()[triangle.idx[2]], graph->GetY()[triangle.idx[2]]);
	}

}

int delaunayTriangulation(bool old = false) {

	const int EXP = 4750;
	const bool VERBOSE = false;

	TGraph2D * graph = getGraph();

	//graph->GetHistogram("");

	TGraphDelaunay2D delaunay(graph);

        if (old)
           graph->Draw("tri1 old");
        else
           graph->Draw("tri1");

        if (gPad) gPad->Update(); // to force drawing

	for (int i = 0; i < 100; i++) {
           Double_t pt = 50 + 0.001;
           Double_t eta = 1. + 0.01 * i + 0.001;
           Double_t res = graph->Interpolate(pt, eta);
           Double_t res2 = delaunay.ComputeZ(pt, eta);
           std::cout << eta << " " << res << "  " << res2 << std::endl;
	}

	if(delaunay.GetNdt() == EXP){
           if(VERBOSE) printDelaunay(delaunay);

	return 0;
	}
        else {
           printf("ERROR - Expected: %i\t Gotten: %i\n", EXP, delaunay.GetNdt());
           if(VERBOSE) printDelaunay(delaunay);

           return 4;
	}

}

int main() {
   return delaunayTriangulation();
}
