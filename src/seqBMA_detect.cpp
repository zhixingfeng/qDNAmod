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
		
		ValueArg<int> n_iterArg("i","interations","number of interations, default is 1",false, 1 ,"number of interations");	
		UnlabeledValueArg<string> nativefileArg("nativefile","native pileup data directory",true,"","native pileup dir");
		UnlabeledValueArg<string> wgafileArg("wgafile","WGA pileup data directory",true,"","WGA pileup dir");
		UnlabeledValueArg<string> priorfileArg("prior","directory containing estimated prior",true,"","prior dir");
		UnlabeledValueArg<string> outdirArg("outdir","output directory",true,"","outdir");
        		
		cmd.add(n_iterArg);
		cmd.add(nativefileArg);
		cmd.add(wgafileArg);
		cmd.add(priorfileArg);
		cmd.add(outdirArg);
		
		cmd.parse(argc, argv);	
		
		// get parameters
		int n_iter = n_iterArg.getValue();
		string nativefile = nativefileArg.getValue();
		string wgafile = wgafileArg.getValue();
		string priorfile = priorfileArg.getValue();
		string outdir = outdirArg.getValue();
		

		// get home path of seqBMA
		string homepath = getHomePath();
		if (homepath==""){
			fprintf(stderr, "Error: can not find SEQBMAHOME\n");
			return 1;
		}
		
		// excute 	
		char R_cmd[5000];
		sprintf(R_cmd, "Rscript %s/R/R_seqBMA_detect.R %s %s %s %s %d",homepath.c_str(), nativefile.c_str(), 
			wgafile.c_str(), priorfile.c_str(), outdir.c_str(), n_iter);	
		system(R_cmd);
	}catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}



