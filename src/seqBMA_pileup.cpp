#include <stdio.h>
#include <iostream>
#include <stl.h>
#include <get_seqBMA_home.h>
#include <tclap/CmdLine.h>
using namespace TCLAP;

int main(int argc, char** argv)
{
	try{
		// parse command line
		CmdLine cmd("", ' ', "0.1");
		
		UnlabeledValueArg<string> cmpH5fileArg("cmpH5file","aligned data in cmpH5 format",true,"","cmpH5file");
		UnlabeledValueArg<string> outdirArg("outdir","output directory",true,"","outdir");
		ValueArg<double> mapQVthresholdArg("m","mapQVthreshold","minimal mapQV",false,255,"mapQVthreshold");
		ValueArg<string> chemistryArg("r","reagent", "chemistry used for SMRT sequencing, candidates are: \"C2\"", false, "C2", "reagent");
        	
		cmd.add(cmpH5fileArg);
		cmd.add(outdirArg);
		cmd.add(mapQVthresholdArg);
		cmd.add(chemistryArg);
		
		cmd.parse(argc, argv);	
		
		// get parameters
		string cmpH5file = cmpH5fileArg.getValue();
		string outdir = outdirArg.getValue();
		double mapQVthreshold = mapQVthresholdArg.getValue();
		string chemistry = chemistryArg.getValue();
		
		if (mapQVthreshold>255){
			fprintf(stderr, "Error: value of mapQVthreshold should be in 0~255\n");
			return 1;
		}
		if(chemistry!="C2"){
			fprintf(stderr, "Error: reagent should be C2\n");
			return 1;
		}

		// get home path of seqBMA
		string homepath = getHomePath();
		if (homepath==""){
			fprintf(stderr, "Error: can not find SEQBWAHOME\n");
			return 1;
		}
		// excute 	
		char R_cmd[5000];
		sprintf(R_cmd, "Rscript %s/R/R_seqBMA_pileup.R %s %s %s %.0lf",homepath.c_str(), cmpH5file.c_str(), outdir.c_str(), chemistry.c_str(), mapQVthreshold);	
		system(R_cmd);
	}catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}



