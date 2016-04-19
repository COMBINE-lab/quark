# CMake generated Testfile for 
# Source directory: /home/hirak/RapCompress/sailfish/src
# Build directory: /home/hirak/RapCompress/sailfish/build/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(unit_tests "/usr/bin/cmake" "-DTOPLEVEL_DIR=/home/hirak/RapCompress/sailfish" "-P" "/home/hirak/RapCompress/sailfish/cmake/UnitTests.cmake")
ADD_TEST(simple_test "/usr/bin/cmake" "-DTOPLEVEL_DIR=/home/hirak/RapCompress/sailfish" "-P" "/home/hirak/RapCompress/sailfish/cmake/SimpleTest.cmake")
