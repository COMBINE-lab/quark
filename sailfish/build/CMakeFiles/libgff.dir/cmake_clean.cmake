FILE(REMOVE_RECURSE
  "CMakeFiles/libgff"
  "CMakeFiles/libgff-complete"
  "libgff-prefix/src/libgff-stamp/libgff-install"
  "libgff-prefix/src/libgff-stamp/libgff-mkdir"
  "libgff-prefix/src/libgff-stamp/libgff-download"
  "libgff-prefix/src/libgff-stamp/libgff-update"
  "libgff-prefix/src/libgff-stamp/libgff-patch"
  "libgff-prefix/src/libgff-stamp/libgff-configure"
  "libgff-prefix/src/libgff-stamp/libgff-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/libgff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
