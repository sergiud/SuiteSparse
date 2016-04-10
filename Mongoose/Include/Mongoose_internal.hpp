#ifndef MONGOOSE_INTERNAL_HPP
#define MONGOOSE_INTERNAL_HPP

#ifndef MONGOOSE_ONE
#define MONGOOSE_ONE 1.0
#endif

#ifndef MONGOOSE_ZERO
#define MONGOOSE_ZERO 0.0
#endif

#ifndef MONGOOSE_MIN2
#define MONGOOSE_MIN2(x,y) ((x)<(y) ? (x) : (y))
#endif

#ifndef MONGOOSE_MAX2
#define MONGOOSE_MAX2(x,y) ((x)>(y) ? (x) : (y))
#endif

#include "Mongoose.hpp"

#endif