solution "Examples"
	configurations{"debug", "release"}
	project "DynMeansExample"
		kind "ConsoleApp"
		language "C++"
		location "build"
		files {"maindm.cpp"}
		links {"lpsolve55"}
		includedirs{"/usr/local/include/eigen3", "/usr/local/include/dynmeans"}
		configuration "debug"
			flags{"Symbols", "ExtraWarnings"}
			buildoptions{"-std=c++0x"}
		configuration "release"
			flags{"Optimize"}
			buildoptions{"-std=c++0x"}
	project "SpecDynMeansExample"
		kind "ConsoleApp"
		language "C++"
		location "build"
		files {"mainsdm.cpp"}
		links {"lpsolve55", "gurobi_c++", "gurobi55"}
		includedirs{"/usr/local/include/eigen3", "/opt/gurobi550/linux64/include", "/usr/local/include/specdynmeans"}
		libdirs{"/opt/gurobi550/linux64/lib"}
		configuration "debug"
			flags{"Symbols", "ExtraWarnings"}
			buildoptions{"-std=c++0x"}
		configuration "release"
			flags{"Optimize"}
			buildoptions{"-std=c++0x"}
