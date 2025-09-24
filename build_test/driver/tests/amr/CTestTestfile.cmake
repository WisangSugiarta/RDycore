# CMake generated Testfile for 
# Source directory: /Users/wsugiarta/Documents/github/RDycore/driver/tests/amr
# Build directory: /Users/wsugiarta/Documents/github/RDycore/build_test/driver/tests/amr
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(amr_c_np_1_basic "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/bin/mpiexec" "--oversubscribe" "-n" "1" "/Users/wsugiarta/Documents/github/RDycore/build_test/driver/tests/../rdycore_amr" "amr_dx1.yaml" "-dm_plex_transform_type" "refine_sbr" "-dm_plex_transform_active" "adapt" "-dm_view" "-dm_fine_view")
set_tests_properties(amr_c_np_1_basic PROPERTIES  _BACKTRACE_TRIPLES "/Users/wsugiarta/Documents/github/RDycore/driver/tests/amr/CMakeLists.txt;17;add_test;/Users/wsugiarta/Documents/github/RDycore/driver/tests/amr/CMakeLists.txt;0;")
