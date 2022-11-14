//------------------------------------------------------------------------------
// Mongoose/Tests/Mongoose_Test_IO.cpp
//------------------------------------------------------------------------------

// Mongoose Graph Partitioning Library, Copyright (C) 2017-2018,
// Scott P. Kolodziej, Nuri S. Yeralan, Timothy A. Davis, William W. Hager
// Mongoose is licensed under Version 3 of the GNU General Public License.
// Mongoose is also available under other licenses; contact authors for details.
// SPDX-License-Identifier: GPL-3.0-only

//------------------------------------------------------------------------------

#include <string>
#include "Mongoose_IO.hpp"
#include "Mongoose_Test.hpp"

using namespace Mongoose;

int runIOTest(const std::string &inputFile, bool validGraph)
{
    LogTest("Running I/O Test on " << inputFile);

    Graph *G = read_graph(inputFile);

    if (validGraph)
    {
        assert(G != NULL);    // A valid graph should not be null
        assert(G->n > 0);     // A valid graph should have 
        assert(G->nz >= 0);   // At least 1 edge
        assert(G->p != NULL); // Column pointers should not be null
        assert(G->i != NULL); // Row numbers should not be null
        G->~Graph();
    }
    else
    {
        assert(G == NULL);
    }

    // Also try with C-style string
    Graph *G2 = read_graph(inputFile.c_str());

    if (validGraph)
    {
        assert(G2 != NULL);    // A valid graph should not be null
        assert(G2->n > 0);     // A valid graph should have 
        assert(G2->nz >= 0);   // At least 1 edge
        assert(G2->p != NULL); // Column pointers should not be null
        assert(G2->i != NULL); // Row numbers should not be null
        G2->~Graph();
    }
    else
    {
        assert(G2 == NULL);
    }
    
    LogTest("I/O Test Completed Successfully");

    return EXIT_SUCCESS;
}
