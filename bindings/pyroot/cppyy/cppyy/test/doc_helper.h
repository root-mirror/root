#ifndef CPPYY_TEST_DOCHELPER_H
#define CPPYY_TEST_DOCHELPER_H

#include <exception>
#include <string>


class SomeError : public std::exception {
public:
    explicit SomeError(const std::string& msg) : fMsg(msg) {}
    SomeError(const SomeError& s) : fMsg(s.fMsg) {}

    const char* what() const throw() override { return fMsg.c_str(); }

private:
    std::string fMsg;
};

class SomeOtherError : public SomeError {
public:
    explicit SomeOtherError(const std::string& msg) : SomeError(msg) {}
    SomeOtherError(const SomeOtherError& s) : SomeError(s) {}
};

namespace DocHelper {

void throw_an_error(int i);

} // namespace DocHelper

#endif // !CPPYY_TEST_DOCHELPER_H
