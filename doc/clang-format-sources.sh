#!/bin/bash

# doc: https://clang.llvm.org/docs/ClangFormatStyleOptions.html
projectrootdir="$(dirname pwd)"
cd $projectrootdir
git ls-files -- '*.cpp' '*.h' | grep "^\(include\|src\|python\|test\)" |  xargs clang-format -i --style=file
