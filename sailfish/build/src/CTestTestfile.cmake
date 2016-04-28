# CMake generated Testfile for 
# Source directory: /home/hirak/quark/sailfish/src
# Build directory: /home/hirak/quark/sailfish/build/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(unit_tests "/usr/bin/cmake" "-DTOPLEVEL_DIR=/home/hirak/quark/sailfish" "-P" "/home/hirak/quark/sailfish/cmake/UnitTests.cmake")
ADD_TEST(simple_test "/usr/bin/cmake" "-DTOPLEVEL_DIR=/home/hirak/quark/sailfish" "-P" "/home/hirak/quark/sailfish/cmake/SimpleTest.cmake")
