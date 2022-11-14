//------------------------------------------------------------------------------
// Mongoose/Tests/Mongoose_Test_IO_exe.cpp
//------------------------------------------------------------------------------

// Mongoose Graph Partitioning Library, Copyright (C) 2017-2018,
// Scott P. Kolodziej, Nuri S. Yeralan, Timothy A. Davis, William W. Hager
// Mongoose is licensed under Version 3 of the GNU General Public License.
// Mongoose is also available under other licenses; contact authors for details.
// SPDX-License-Identifier: GPL-3.0-only

//------------------------------------------------------------------------------

#include "Mongoose_Test.hpp"

using namespace Mongoose;

#undef LOG_ERROR
#undef LOG_WARN
#undef LOG_INFO
#undef LOG_TEST
#define LOG_ERROR 1
#define LOG_WARN 1
#define LOG_INFO 0
#define LOG_TEST 1

int main(int argn, char** argv)
{
    SuiteSparse_start();

    // Set Logger to report all messages
    Logger::setDebugLevel(All);
    
    if (argn != 3)
    {
        // Wrong number of arguments - return error
        LogError("Usage: mongoose_test_io <MM-input-file.mtx> <1 for valid graph, 0 for invalid>");
        SuiteSparse_finish();
        return EXIT_FAILURE;
    }

    // Read in input file name
    std::string inputFile = std::string(argv[1]);

    // Read in whether this file should produce a valid graph
    bool validGraph = static_cast<bool>(atoi(argv[2]));

    // Run the I/O test
    int status = runIOTest(inputFile, validGraph);

    SuiteSparse_finish();

    return status;
}
