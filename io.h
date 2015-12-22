#ifndef IO_H
#define IO_H

#include <string>
#include <vector>

int write_to_file(const std::string, const std::vector <double> &,
                  const std::vector <double> &, size_t = 1);
int write_to_file(const std::string, const std::vector <double> &,
                  const std::vector <double> &, const std::vector<double> &,
                  size_t = 1);
int write_to_file_slvecs(const std::string, const std::vector <double> &,
                         const std::vector <double> &, const std::vector<double> &,
                         size_t = 1);

#endif /*OI _H */
