#include <errno.h>
#include <string.h>
#include <sys/stat.h>

#include "Logger.hpp"
#include "RadioRT.hpp"
#include "RadioRT_Parameters.hpp"
#include "ParseLua.hpp"
#include "FileManagement.hpp"

#include "selene/include/selene.h"

void parseParameters(const std::string& filename, RadioRT_Parameters& params);
void showUsage();

int main(int argc, char *argv[]) {
	std::string paramFile = "config/radio-config.lua";

	// Parse commandline arguments.
	if (argc > 2) {
		showUsage();
		exit(0);
	}
	if (argc > 1) {
		for (int iarg = 1; iarg < argc; ++iarg) {
			std::string arg = argv[iarg];
			std::string prefix1("--config=");

			if (!arg.compare(0, prefix1.size(), prefix1))
				paramFile = arg.substr(prefix1.size()).c_str();
			else {
				showUsage();
				exit(0);
			}
		}
	}

	// Run program.
	try {
		RadioRT_Parameters params;
		parseParameters(paramFile, params);
		FileManagement::makeDirectoryPath(params.outputDirectory);

		RadioRT radio(params);
		radio.run();
	}
	catch (std::exception& e) {
		Logger<FileLogPolicy>::Instance().print<SeverityType::FATAL_ERROR>(e.what());
		Logger<ConsoleLogPolicy>::Instance().print<SeverityType::FATAL_ERROR>(e.what());
		std::cout << "See radio.log for more details." << std::endl;
	}

	return 0;
}

void parseParameters(const std::string& filename, RadioRT_Parameters& params) {
	// Create new Lua state and load the lua libraries.
	sel::State luaState{true};
	//luaState.HandleExceptionsWith([](int, std::string msg, std::exception_ptr){ throw std::runtime_error(msg);});

	if (!luaState.Load(filename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + filename);
	}
	else {
		parseLuaVariable(luaState["Parameters"]["torch_params_filename"], params.torchParamFilename);
		parseLuaVariable(luaState["Parameters"]["torch_data_filename"],   params.inputfile);
		parseLuaVariable(luaState["Parameters"]["output_directory"],      params.outputDirectory);

		parseLuaVariable(luaState["Parameters"]["sampling"],              params.sampling);
		parseLuaVariable(luaState["Parameters"]["dopplerShifted"],        params.dopplerShifted);
		parseLuaVariable(luaState["Parameters"]["dopp_shift_phi_inc"],    params.doppShiftPhiIncr);
		parseLuaVariable(luaState["Parameters"]["distance"],              params.dist);
		parseLuaVariable(luaState["Parameters"]["right_ascension"],       params.rightAscension);
		parseLuaVariable(luaState["Parameters"]["declination"],           params.declination);
		parseLuaVariable(luaState["Parameters"]["theta"],                 params.theta);
		parseLuaVariable(luaState["Parameters"]["phi"],                   params.phi);

		parseLuaVariable(luaState["Parameters"]["frequency"],             params.frequency);
		parseLuaVariable(luaState["Parameters"]["bandwidth"],             params.bandwidth);
		parseLuaVariable(luaState["Parameters"]["nchannels"],             params.nchannels);
		parseLuaVariable(luaState["Parameters"]["nlevel"],                params.nlevel);
		parseLuaVariable(luaState["Parameters"]["stokes"],                params.stokes);

		parseLuaVariable(luaState["Parameters"]["turb_broadening"],       params.turb_broad);
		parseLuaVariable(luaState["Parameters"]["vLOS"],                  params.vLOS);

		parseLuaVariable(luaState["Parameters"]["integratingFF"],         params.integratingFF);
		parseLuaVariable(luaState["Parameters"]["integratingRL"],         params.integratingRL);
		parseLuaVariable(luaState["Parameters"]["resolution_scale"],      params.resolutionScale);

		if (params.outputDirectory.back() != '/') {
			params.outputDirectory += "/";
		}
	}

	sel::State luaState2{true};
	//luaState2.HandleExceptionsWith([](int, std::string msg, std::exception_ptr){ throw std::runtime_error(msg);});

	if (!luaState2.Load(params.torchParamFilename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + params.torchParamFilename);
	}
	else {
		parseLuaVariable(luaState2["Parameters"]["Grid"]["no_dimensions"], params.ndims);
		parseLuaVariable(luaState2["Parameters"]["Grid"]["no_cells_x"], params.ncells[0]);
		parseLuaVariable(luaState2["Parameters"]["Grid"]["no_cells_y"], params.ncells[1]);
		parseLuaVariable(luaState2["Parameters"]["Grid"]["no_cells_z"], params.ncells[2]);
		parseLuaVariable(luaState2["Parameters"]["Grid"]["geometry"], params.geometry);
		parseLuaVariable(luaState2["Parameters"]["Grid"]["side_length"], params.sideLength);

		if (params.ncells[0] < 1 || params.ncells[1] < 1 || params.ncells[2] < 1)
			throw std::runtime_error("parseParameters: ncells is invalid.");

		if (params.ndims < 1 || params.ndims > 3)
			throw std::runtime_error("parseParameters: ndims < 1 or ndims > 3.");

		if (params.geometry == "cylindrical" && params.ndims != 2) {
			throw std::runtime_error("parseParameters: cylindrical geometry can only have 2 dimensions.");
		}
		else if (params.geometry == "spherical" && params.ndims != 1) {
			throw std::runtime_error("parseParameters: spherical geometry can only have 1 dimension.");
		}

		if (params.geometry != "cylindrical" && params.geometry != "cartesian")
			throw std::runtime_error("parseParameters: " + params.geometry + " geometry not supported.");
	}
}

void showUsage() {
	std::cout << "Usage: radio [--config=<filename>]" << std::endl;
}
