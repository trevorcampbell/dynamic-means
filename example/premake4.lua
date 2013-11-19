solution "DynMeansExample"
	configurations{"debug", "release"}
	project "DynMeansExample"
		kind "ConsoleApp"
		language "C++"
		location "build"
		files {"main.cpp"}
		links {"lpsolve55"}
		includedirs{"/usr/local/include/eigen3", "/usr/local/include/dynmeans"}
		configuration "debug"
			flags{"Symbols", "ExtraWarnings"}
			buildoptions{"-std=c++0x"}
		configuration "release"
			flags{"Optimize"}
			buildoptions{"-std=c++0x"}

