// POSIX process/env shims for MSVC. Include is a no-op on POSIX.
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_POSIX_COMPAT_HPP
#define STL2STEP_POSIX_COMPAT_HPP

#ifdef _WIN32
#include <cstdio>
#include <cstdlib>

inline FILE* popen(const char* command, const char* mode) {
    return ::_popen(command, mode);
}

inline int pclose(FILE* stream) {
    return ::_pclose(stream);
}

inline int setenv(const char* name, const char* value, int overwrite) {
    if (!name) return -1;
    if (!overwrite && std::getenv(name) != nullptr) return 0;
    return ::_putenv_s(name, value ? value : "");
}

inline int unsetenv(const char* name) {
    if (!name) return -1;
    return ::_putenv_s(name, "");
}
#endif

#endif
