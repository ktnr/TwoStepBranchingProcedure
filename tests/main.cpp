#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// How to test private methods?
////#define private public
// Use instead: https://stackoverflow.com/questions/27778908/define-private-to-public-in-c. Doesn't work.
// Try https://www.codeproject.com/Tips/5249547/How-to-Unit-Test-a-Private-Function-in-Cplusplus.
// Current solution: just declare the respective methods public.