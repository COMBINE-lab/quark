FILE(REMOVE_RECURSE
  "CMakeFiles/libtbb"
  "CMakeFiles/libtbb-complete"
  "libtbb-prefix/src/libtbb-stamp/libtbb-install"
  "libtbb-prefix/src/libtbb-stamp/libtbb-mkdir"
  "libtbb-prefix/src/libtbb-stamp/libtbb-download"
  "libtbb-prefix/src/libtbb-stamp/libtbb-update"
  "libtbb-prefix/src/libtbb-stamp/libtbb-patch"
  "libtbb-prefix/src/libtbb-stamp/libtbb-configure"
  "libtbb-prefix/src/libtbb-stamp/libtbb-build"
  "libtbb-prefix/src/libtbb-stamp/libtbb-reconfigure"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/libtbb.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
