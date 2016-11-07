/**
 * This header's purpose is mainly to replace boost's filesystem library. The dir_iterator gives access to files within
 * a directory using tinydir, the other methods are self-explanatory.
 *
 * @file filesystem.h
 * @brief Definition of the dir_iterator and some file system utilities.
 * @author clonker
 * @date 14.10.16
 */

#ifndef READDY_MAIN_FILESYSTEM_H
#define READDY_MAIN_FILESYSTEM_H

#include <string>
#include <memory>

namespace readdy {
namespace util {
namespace fs {

constexpr char separator =
#if READDY_WINDOWS
        '\\';
#else
        '/';
#endif

struct dir_iterator {

    dir_iterator(const std::string& path);
    bool has_next() const;
    std::string base_name() const;
    std::string next();
    virtual ~dir_iterator();

    dir_iterator(const dir_iterator&) = delete;
    dir_iterator& operator=(const dir_iterator&) = delete;
    dir_iterator(dir_iterator&&);
    dir_iterator& operator=(dir_iterator&&);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

};

std::string current_path();

bool exists(const std::string& path);

bool is_file(const std::string &path);

bool is_directory(const std::string &path);

}
}
}
#endif //READDY_MAIN_FILESYSTEM_H