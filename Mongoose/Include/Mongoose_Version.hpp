/* ========================================================================== */
/* === Include/Mongoose_Version.hpp ========================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * Mongoose Graph Partitioning Library, Copyright (C) 2017-2018,
 * Scott P. Kolodziej, Nuri S. Yeralan, Timothy A. Davis, William W. Hager
 * Mongoose is licensed under Version 3 of the GNU General Public License.
 * Mongoose is also available under other licenses; contact authors for details.
 * SPDX-License-Identifier: GPL-3.0-only
 * -------------------------------------------------------------------------- */

#pragma once

#include <string>

// Configuration information from CMake
#define Mongoose_VERSION_MAJOR 3
#define Mongoose_VERSION_MINOR 0
#define Mongoose_VERSION_PATCH 0
#define Mongoose_DATE "Nov 12, 2022"

namespace Mongoose
{

int major_version();
int minor_version();
int patch_version();
std::string mongoose_version();

} // end namespace Mongoose
