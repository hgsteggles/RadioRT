#include <errno.h>
#include <string.h>
#include <sys/stat.h>

#include "Logger.hpp"
#include "RadioRT.hpp"
#include "RadioRT_Parameters.hpp"
#include "ParseLua.hpp"

#include "lib/selene/include/selene.h"

//extern void readinfits_3d();
//extern void integ(double slice_ff[][jmax], double slice_gam[][jmax]);
//extern void wfits2d(char filename[], char filetype[], char units[], int imx, int imy, double z[], double pixs, double pixsize, double rb, double rapos, double decpos, double im_max);

void parseParameters(const std::string& filename, RadioRT_Parameters& params);
void showUsage();
int do_mkdir(const char *path, mode_t mode);
int mkpath(const char *path, mode_t mode);

typedef struct stat Stat;
#ifndef lint
/* Prevent over-aggressive optimizers from eliminating ID string */
const char jlss_id_mkpath_c[] = "@(#)$Id: mkpath.c,v 1.13 2012/07/15 00:40:37 jleffler Exp $";
#endif /* lint */

int main(int argc, char *argv[]) {
	std::string paramFile = "config/radioconfig.lua";

	// Parse parameters
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
	luaState.HandleExceptionsWith([](int, std::string msg, std::exception_ptr){ throw std::runtime_error(msg);});

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

		parseLuaVariable(luaState["Parameters"]["rbeam_degrees"],         params.rbeam_degrees);
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

		int rc = mkpath(params.outputDirectory.c_str(), 0777);
		if (rc != 0)
			throw std::runtime_error("Failed to create directory tree: " + params.outputDirectory);

		if (params.outputDirectory.back() != '/') {
			params.outputDirectory += "/";
		}
	}

	sel::State luaState2{true};
	luaState2.HandleExceptionsWith([](int, std::string msg, std::exception_ptr){ throw std::runtime_error(msg);});

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

		if (params.geometry != "cylindrical")
			throw std::runtime_error("parseParameters: only cylindrical geometry supported.");
	}
}

void showUsage() {
	std::cout << "Usage: radio [--config=<filename>]" << std::endl;
}

int do_mkdir(const char *path, mode_t mode)
{
	Stat            st;
	int             status = 0;

	if (stat(path, &st) != 0)
	{
		/* Directory does not exist. EEXIST for race condition */
		if (mkdir(path, mode) != 0 && errno != EEXIST)
			status = -1;
	}
	else if (!S_ISDIR(st.st_mode))
	{
		errno = ENOTDIR;
		status = -1;
	}

	return(status);
}

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
int mkpath(const char *path, mode_t mode)
{
	char           *pp;
	char           *sp;
	int             status;
	char           *copypath = strdup(path);

	status = 0;
	pp = copypath;
	while (status == 0 && (sp = strchr(pp, '/')) != 0)
	{
		if (sp != pp)
		{
			/* Neither root nor double slash in path */
			*sp = '\0';
			status = do_mkdir(copypath, mode);
			*sp = '/';
		}
		pp = sp + 1;
	}
	if (status == 0)
		status = do_mkdir(path, mode);
	free(copypath);
	return (status);
}
