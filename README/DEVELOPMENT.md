# ROOT Development Practice

## Overview

The development of ROOT almost exclusively happens using the [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests)
model of github. A pull request (PR) should contain a set focused changes
organized in one or more [atomic commits](https://en.wikipedia.org/wiki/Atomic_commit#Revision_control).
PRs should be well-documented and well-tested in order to allow other community
members to use, maintain and modify. If the PR contains performance-critical
code consider writing a benchmark against the [rootbench repository](https://github.com/root-project/rootbench).


## Quality Assurance

Each contribution should contain developer documentation in the form of code
comments and sufficient amount of tests in the form of unit and/or integration
tests. Unit tests are relatively small and quick programs focused to check if
small pieces of code and API work as expected. Integration tests are checks
which ensure the synergy between different (unit tested) components. Put in
practice, unit tests verify (member) function behavior whereas integration tests
check classes and their cooperation. The boundary between both kinds of testing
is blurred.

ROOT has support for both kinds of tests in the [roottest repository](https://github.com/root-project/roottest)
and supports "inline" unit tests in each component's `test` folder. Unit testing
uses the [GTest and GMock](https://github.com/google/googletest) infrastructure
along with small ROOT-specific extensions located in
[ROOTUnitTestSupport](../test/unit_testing_support). The documentation of GTest
and GMock is rather extensive and we will describe some of the features of
ROOTUnitTestSupport. In order to write an inline unit test, add a new file in the
nearest to the tested component's `test` folder and call `ROOT_ADD_GTEST` in the
`CMakeLists.txt` file.

In many cases using standard GTest facility is sufficient to write a good test.
However, sometimes we want to test the error conditions of an interface and
its diagnostics. For example,

```cpp

void MyROOTFunc(unsigned x)
{
  if (x == 0) {
    Error("MyROOTFunc", "x should be greater than 0!");
    return;
  }
  /* some work */
}

```

In order to check if the error is printed on the right spot we can write the
following test:

```cpp

#include "ROOTUnitTestSupport.h"

#include "gtest/gtest.h"

TEST(MyROOTFunc, ErrorCases)
{
  ROOT_EXPECT_ERROR(MyROOTFunc(0), "MyROOTFunc", "x should be greater than 0!");
  // Also ROOT_EXPECT_WARNING, ROOT_EXPECT_INFO, ROOT_EXPECT_NODIAG and ROOT_EXPECT_SYSERROR available.
}

```

