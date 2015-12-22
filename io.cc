// io.cc
//
// Copyright (C) 2015 - Atri Bhattacharya
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include "io.h"

using namespace std;
using boost::format;

int write_to_file(const string fname, const vector<double> & x,
                  const vector<double> & y, size_t skips)
{
    ofstream fout(fname.c_str());

    if(fout.bad())
    {
        return 1;
    }

    for (unsigned int i = 0; i < x.size(); i += skips)
    {
        fout.precision(7);
        fout.width(10);
        fout << fixed << x.at(i);
        fout.width(30);
        fout << fixed << y.at(i)
             << endl;
    }

    fout.close();
    return 0;

}

int write_to_file(const string fname, const vector<double> & x,
                  const vector<double> & y,  const vector<double> & z,
                  size_t skips)
{
	string fbase    = fname.substr(0, fname.find_last_of('.'));
	vector<string> ax_file = { fbase + "-xax.dat", fbase + "-yax.dat",
	                           fbase + "-zax.dat" };

	size_t numfiles  = ax_file.size();

	vector<ofstream> ax_fout(3);

	for (size_t f = 0; f < numfiles; ++f)
	{
		ax_fout[f].open(ax_file[f].c_str());
	}

	for (double xx : x)
		ax_fout[0] << format("%|.5E|\n") % xx;

	for (double yy : y)
		ax_fout[1] << format("%|.5E|\n") % yy;

	for (unsigned int k = 0; k < z.size(); ++k)
	{
		if (! ((k) % (z.size()/x.size()/y.size())))
			ax_fout[2] << endl;
		ax_fout[2] << format("%|15.5f|") % z[k];
	}
	
	ofstream fout(fname.c_str());

    if(fout.bad())
    {
        return 1;
    }

//  Write data after every skips iterations to reduce file size
    for (unsigned int i = 0; i < x.size(); i += skips)
    {
		for (unsigned int j = 0; j < y.size(); j += 1)
		{
/*			fout.precision(7);
			fout.width(10);
			fout << fixed << x.at(i);
			fout.width(30);
			fout << fixed << y.at(j);
*/
			for (unsigned int k = 0; k < z.size() / (x.size() * y.size()); ++k)
			{
				fout << format("%|15.5E|") 
					    % z.at(k + j * y.size() + i * x.size());
			}
			fout << endl;
		}
		fout << endl;
    }

    fout.close();
    return 0;

}

int write_to_file_slvecs(const string fname, const vector<double> & x,
                         const vector<double> & y,  const vector<double> & z,
                         size_t skips)
{
	ofstream fout(fname.c_str());

    if(fout.bad())
    {
        return 1;
    }

//  Write data after every skips iterations to reduce file size
    for (unsigned int i = 0; i < x.size(); i += skips)
    {
		fout.precision(7);
		fout.width(10);
		fout << fixed << x.at(i);
		fout.width(30);
		fout << fixed << y.at(i);
		fout.width(30);
		fout << fixed << z.at(i) << endl;
    }

    fout.close();
    return 0;

}
