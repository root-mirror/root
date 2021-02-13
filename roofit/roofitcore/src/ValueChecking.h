/*
 * ValueChecking.h
 *
 *  Created on: 03.08.2020
 *      Author: shageboeck
 */

#ifndef ROOFIT_ROOFITCORE_INC_VALUECHECKING_H_
#define ROOFIT_ROOFITCORE_INC_VALUECHECKING_H_

#include <TSystem.h>

#include <exception>
#include <vector>
#include <string>
#include <sstream>

class CachingError : public std::exception {
  public:
    CachingError(const std::string& newMessage) :
      std::exception(),
      _messages()
    {
      _messages.push_back(newMessage);
    }

    CachingError(CachingError&& previous, const std::string& newMessage) :
    std::exception(),
    _messages{std::move(previous._messages)}
    {
      _messages.push_back(newMessage);
    }

    const char* what() const noexcept override {
      std::stringstream out;
      out << "**Computation/caching error** in\n";

      std::string indent;
      for (auto it = _messages.rbegin(); it != _messages.rend(); ++it) {
        std::string message = *it;
        auto pos = message.find('\n', 0);
        while (pos != std::string::npos) {
          message.insert(pos+1, indent);
          pos = (message.find('\n', pos+1));
        }

        out << indent << message << "\n";
        indent += " ";
      }

      out << std::endl;

      std::string* ret = new std::string(out.str()); //Make it survive this method

      return ret->c_str();
    }


  private:
    std::vector<std::string> _messages;
};


class FormatPdfTree {
  public:
    template <class T,
    typename std::enable_if<std::is_base_of<RooAbsArg, T>::value>::type* = nullptr >
    FormatPdfTree& operator<<(const T& arg) {
      _stream << arg.ClassName() << "::" << arg.GetName() << " " << &arg << " ";
      arg.printArgs(_stream);
      return *this;
    }

    template <class T,
    typename std::enable_if< ! std::is_base_of<RooAbsArg, T>::value>::type* = nullptr >
    FormatPdfTree& operator<<(const T& arg) {
      _stream << arg;
      return *this;
    }

    operator std::string() const {
      return _stream.str();
    }

    std::ostream& stream() {
      return _stream;
    }

  private:
    std::ostringstream _stream;
};



#endif /* ROOFIT_ROOFITCORE_INC_VALUECHECKING_H_ */
