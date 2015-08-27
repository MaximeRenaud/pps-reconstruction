# pps-reconstruction
Reconstruction code for the test beam performed at CERN's SPS on behalf of the CT-PPS project

## Building the code
First of all, you need to ensure CMake is properly installed on your machine.

Once done, create a new directory (for instance <path_to_your_pps-reconstruction_folder>/build) and 'cd' into it.
Then, type 'cmake <path_to_your_pps-reconstruction_folder>' (or 'cmake ..' if you followed the convention above) to generate all dependencies/makefiles/...
Finally, type 'make' to build the library.

In order to build any script in the test/ subdirectory, you have to 
1) insert it (if brand new) into the test/CMakeLists.txt file by adding:
> add_test(\<name_of_your_file_without_cpp_extension\>)

2) build it (type "make <name_of_your_file_without_cpp_extension>")
