
#include "Interop.hpp"
#include "Matching.hpp"
#include "BoundaryHeap.hpp"
#include "cs.hpp"
#include "stdio.h"
#include "string.h"

namespace SuiteSparse_Mongoose
{

/* Sanitizes a matrix market input file into a CSparse3 input file. */
void sanitize
(
    const char *mmFilename,
    const char *csFilename
);

//-----------------------------------------------------------------------------
// load a 1-indexed triplet matrix from a flat text file
//-----------------------------------------------------------------------------
cs *cs_load2 (FILE *f)
{
    if(!f) return NULL;

    double i, j ;   /* use double for integers to avoid csi conflicts */
    double x ;
    cs *T ;
    T = cs_spalloc (0, 0, 1, 1, 1) ;                    /* allocate result */

    double m, n, nz ;
    fscanf (f, "%lg %lg %lg\n", &m, &n, &nz) ;

    while (fscanf (f, "%lg %lg %lg\n", &i, &j, &x) == 3)
    {
        if(i == j) continue;
        if (!cs_entry (T, (csi) i-1, (csi) j-1, x)) return (cs_spfree (T)) ;
    }

    /* Override cs_entry with what the proclaimed values are,
     * since we strip the diagonal from the input. */
    T->m = (csi) m;
    T->n = (csi) n;

    return (T) ;
}

//-----------------------------------------------------------------------------
// Take a matrix market format file, sanitize it, and use the sanitized file
// to load a CSparse3 compressed column sparse matrix and returns it.
//-----------------------------------------------------------------------------
Graph *readGraphFromMM
(
    const char *mmFilename
)
{
    if(!mmFilename) return NULL;

    /* Build the name of the output file. */
    char outputFile[256];
    strcpy(outputFile, mmFilename);
    strcat(outputFile, ".s");

    /* Sanitize the MM file and save the result in the output file. */
    sanitize(mmFilename, outputFile);

    /* Open the output file. */
    FILE *csFile = fopen(outputFile, "r");
    if(!csFile) return NULL;
    cs *A = cs_load2(csFile);
    fclose(csFile);
    if(!A) return NULL;

    /* Compress to column pointer form. */
    cs *G = cs_compress(A);
    cs_spfree(A);

    /* Brain transplant CSparse3 matrix into an empty Mongoose Graph. */
    Graph *returner = CSparse3ToGraph(G);

    /* Return the Mongoose Graph. */
    return returner;
}

void sanitize
(
    const char *mmFilename,
    const char *csFilename
)
{
    FILE *graph = fopen(mmFilename, "r");
    if(!graph) return;
    FILE *output = fopen(csFilename, "w");
    if(!output) return;

    bool isFirst = true;
    char line[256];
    while(fgets(line, 256, graph))
    {
        if(line[0] == '%') continue;
        else if(isFirst)
        {
            fprintf(output, "%s", line);
            isFirst = false;
            continue;
        }
        else
        {
            fprintf(output, "%s", line);
        }
    }

    fclose(output);
    fclose(graph);
}

/* Configure a CSparse3 matrix from an existing Mongoose Graph. */
cs *GraphToCSparse3(Graph *G, bool copy)
{
    cs *A = (cs*) SuiteSparse_malloc(1, sizeof(cs));
    if(!A) return NULL;
    A->n = G->cs_n;
    A->m = G->cs_m;
    A->nz = G->cs_nz;
    A->nzmax = G->cs_nzmax;
    if(!copy)
    {
        A->p = (ptrdiff_t*) G->p;
        A->i = (ptrdiff_t*) G->i;
        A->x = G->x;
    }
    else
    {
        Int n = G->n;
        Int nz = G->nz;
        A->p = (ptrdiff_t*) SuiteSparse_malloc((n+1), sizeof(ptrdiff_t));
        A->i = (ptrdiff_t*) SuiteSparse_malloc(nz, sizeof(ptrdiff_t));
        A->x = (Double*) SuiteSparse_malloc(nz, sizeof(Double));
        if(!A->p || !A->i || !A->x)
        {
            cs_spfree(A);
            return NULL;
        }

        for(Int k=0; k<=n; k++)
        {
            A->p[k] = (ptrdiff_t) G->p[k];
        }
        for(Int p=0; p<nz; p++)
        {
            A->i[p] = (ptrdiff_t) G->i[p];
            A->x[p] = G->x[p];
        }
    }

    return A;
}

/* Create a new Mongoose Graph from an existing CSparse3 matrix. */
Graph *CSparse3ToGraph(cs *G, bool resetEW, bool resetNW)
{
    Graph *returner = (Graph*) SuiteSparse_calloc(1, sizeof(Graph));
    if(!returner) return NULL;
    new (returner) Graph();

    /* Brain-transplant the graph to the new representation. */
    returner->cs_n = G->n;
    returner->cs_m = G->m;
    returner->cs_nz = G->nz;
    returner->cs_nzmax = G->nzmax;
    returner->n = MONGOOSE_MAX2(G->n, G->m);
    returner->nz = G->p[G->n];
    returner->p = (Int*) G->p;
    returner->i = (Int*) G->i;
    returner->x = G->x;

    /* Allocate edge weights if necessary. */
    bool attachEdgeWeights = false;
    if(!returner->x || resetEW)
    {
        Int nz = returner->nz;
        returner->x = (Weight*) SuiteSparse_malloc(nz, sizeof(Weight));
        attachEdgeWeights = true;        
    }

    /* Allocate node weights if necessary. */
    bool attachNodeWeights = false;
    if(!returner->w || resetNW)
    {
        Int n = returner->n;
        returner->w = (Weight*) SuiteSparse_malloc(n, sizeof(Weight));
        attachNodeWeights = true;        
    }

    /* If we failed to attach weights, free the graph and return. */
    if((attachEdgeWeights && !returner->x) || (attachNodeWeights && !returner->w))
    {
        /* Undo the brain transplant, free the Graph skeleton, and return. */
        returner->p = NULL;
        returner->i = NULL;
        returner->x = NULL;
        returner->~Graph();
        SuiteSparse_free(returner);
        return NULL;
    }

    if(attachEdgeWeights)
    {
        Int nz = returner->nz;
        for(Int p=0; p<nz; p++) returner->x[p] = 1.0;
    }

    if(attachNodeWeights)
    {
        int n = returner->n;
        for(Int k=0; k<n; k++) returner->w[k] = 1.0;
    }

    return returner;
}

#if 0
void writeDot
(
    Graph *G,
    Options *O,
    const char *comment,
    DOTRenderSetting in_setting
)
{
    if(!O->writeDot) return;

    Int n = G->n;
    Int *Gp = G->p;
    Int *Gi = G->i;
    Int *bhIndex = G->bhIndex;
    bool *partition = G->partition;
    Int *externalDegree = G->externalDegree;
    Int *matching = G->matching;
    Int *matchtype = G->matchtype;

    DOTRenderSetting setting = in_setting;
    if((!partition || !bhIndex) && setting == CutSet) setting = GraphOnly;
    if(!matching && setting == Matching) setting = GraphOnly;

    char filename[256];
    sprintf(filename, "%s-%03ld-%s.dot", O->problemName, O->graphSID++, comment);
    FILE *dot = fopen(filename, "w");

    fprintf(dot, "graph G {\n");
    fprintf(dot, "dim=2;\n");
    fprintf(dot, "bgcolor=black;\n");
    fprintf(dot, "node [shape=point,fontcolor=white,color=white];\n");
    fprintf(dot, "edge [color=%s];\n", setting == Matching ? "\"#333333\"" : "white");

    /* Write out the vertex specifiers. */
    for(Int k=0; k<n; k++)
    {
        fprintf(dot, "%ld ", k);
        switch(setting)
        {
            case CutSet:
            {
                bool kPartition = partition[k];
                const char *borderColor = externalDegree[k] ? "red" : partition[k] ? "green" : "blue";
                const char *fillColor = kPartition ? "green" : "blue";
                fprintf(dot, "[color=%s,fillcolor=%s]", borderColor, fillColor);
                break;
            }

            case Matching:
            {
                const char *colors[4] = {"red", "blue", "green", "purple"};
                fprintf(dot, "[color=%s]", colors[matchtype[k]]);
                break;
            }

            case GraphOnly: default:
                break;
        }
        fprintf(dot, ";\n");
    }

    /* Write out the edge specifiers. */
    for(Int k=0; k<n; k++)
    {
        for(Int p=Gp[k]; p<Gp[k+1]; p++)
        {
            Int neighbor = Gi[p];

            /* After the first coarsening, row indices are no longer sorted in ascending order.
             * However, the underlying matrix representation is still symmetric.
             * We can get by with only traversing the lower triangular part. */
            if(neighbor > k) continue;

            fprintf(dot, "%ld -- %ld ", k, neighbor);
            switch(setting)
            {
                case CutSet:
                {
                    bool kPartition = partition[k];
                    bool onSameSide = (kPartition == partition[neighbor]);
                    fprintf(dot, "[color=%s]", onSameSide ? partition[k] ? "green" : "blue" : "red");
                    break;
                }

                case Matching:
                    if(GET_MATCH(k) == neighbor) fprintf(dot, "[color=white]");
                    break;

                case GraphOnly: default:
                    break;
            }
            fprintf(dot, ";\n");
        }
    }

    fprintf(dot, "}\n");

    fclose(dot);
}

void heapWriteDot
(
    Graph *G,
    Options *O,
    Int *heap,
    Int size,
    Weight *gains,
    Int vHighlight
)
{
    if(!O->writeDot) return;

    char filename[256];
    sprintf(filename, "%s-%03ld-%s.dot", O->problemName, O->graphSID++, "heap");
    FILE *dot = fopen(filename, "w");

    fprintf(dot, "graph G{\n");
    fprintf(dot, "bgcolor=black;\n");
    fprintf(dot, "node [shape=circle,fontcolor=white,color=white,fillcolor=black];\n");
    fprintf(dot, "edge [color=white];\n");

    /* Nodes */
    for(Int k=0; k<size; k++)
    {
        fprintf(dot, "%ld [label=\"%ld (%f)\"", k, heap[k], gains[heap[k]]);
        if(k == vHighlight) fprintf(dot, ", color=red");
        fprintf(dot, "];\n");
    }

    /* Edges */
    for(Int k=0; k<size; k++)
    {
        if(LEFT_CHILD(k) < size)  fprintf(dot, "%ld -- %ld;\n", k, LEFT_CHILD(k));
        if(RIGHT_CHILD(k) < size) fprintf(dot, "%ld -- %ld;\n", k, RIGHT_CHILD(k));
    }

    fprintf(dot, "}\n");

    fclose(dot);
}
#endif

}