-- RadioRT Parameters
Parameters = {
	torch_params_filename =      "config/torch-config.lua",
	torch_data_filename =        "config/example-data.txt",
	output_directory =           "tmp/",

	sampling =                   4.0,
	dopplerShifted =             false,
	dopp_shift_phi_inc =         0.01,

	distance =                   1.5,
	right_ascension =            308.94617775,
	declination =                41.379058775,
	theta =                      45.0,
	phi =                        0.0,

	frequency =                  5.0e9,
	bandwidth =                  0,
	nchannels =                  1,
	nlevel =                     0,
	stokes =                     1,

	turb_broadening =            1.0e6,
	vLOS =                       0,

	integratingFF =              true,
	integratingRL =              false,
	resolution_scale =           2,
}
