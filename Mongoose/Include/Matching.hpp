#ifndef MATCHING_HPP_
#define MATCHING_HPP_

#include "Mongoose_internal.hpp"

namespace SuiteSparse_Mongoose
{

enum MatchType
{
    MatchType_Orphan = 0,
    MatchType_Standard = 1,
    MatchType_Brotherly = 2,
    MatchType_Community = 3
};

void match(Graph *G, Options *O);

void matching_Random(Graph *G, Options *O);
void matching_HEM(Graph *G, Options *O);
void matching_PA(Graph *G, Options *O);
void matching_DavisPA(Graph *G, Options *O);
void matching_Cleanup(Graph *G, Options *O);

}

/* Mongoose Matching Macros */
#ifndef MONGOOSE_MATCHING_MACROS
#define MONGOOSE_MATCHING_MACROS

#define MONGOOSE_IS_MATCHED(a) (matching[(a)] > 0)

#define MONGOOSE_GETMATCH(a) (matching[(a)]-1)

#define MONGOOSE_MATCH(a,b,t)                           \
{                                                       \
    matching[(a)] = (b)+1;                              \
    matching[(b)] = (a)+1;                              \
    invmatchmap[cn] = (a);                              \
    matchtype[(a)] = matchtype[(b)] = (t);              \
    matchmap[(a)] = matchmap[(b)] = cn;                 \
    cn++;                                               \
}                                                       \

#define MONGOOSE_COMMUNITY_MATCH(a,b,t)                 \
{                                                       \
    Int vm[4] = {-1,-1,-1,-1};                          \
    vm[0] = a;                                          \
    vm[1] = MONGOOSE_GETMATCH(vm[0]);                   \
    vm[2] = MONGOOSE_GETMATCH(vm[1]);                   \
    vm[3] = MONGOOSE_GETMATCH(vm[2]);                   \
                                                        \
    bool is3Way = (vm[0] == vm[3]);                     \
    if(is3Way)                                          \
    {                                                   \
        matching[vm[1]] = a+1;                          \
        MONGOOSE_MATCH(vm[2], b, t);                    \
    }                                                   \
    else                                                \
    {                                                   \
        matching[b] = matching[a];                      \
        matching[a] = b+1;                              \
        matchmap[b] = matchmap[a];                      \
        matchtype[b] = t;                               \
    }                                                   \
}                                                       \

#endif

#endif
