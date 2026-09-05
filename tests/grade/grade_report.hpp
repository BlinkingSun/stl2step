#ifndef GRADE_REPORT_HPP
#define GRADE_REPORT_HPP

#include "grade_compare.hpp"

#include <string>

namespace grade {

std::string writeJson(const GradeDocument& d);
std::string writeMd(const GradeDocument& d);
bool writeFile(const std::string& path, const std::string& text);

}  // namespace grade

#endif
