#include "precompile.h"
#include "Bool.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include "MPMesh.h"
#include <ctime>
#include "configure.h"

namespace CSG
{

extern HANDLE _output;
extern clock_t t0;


void DebugInfo(char* str, clock_t& t0);
void StdOutput(const char* str);


extern "C" CSG_API MPMesh2* BooleanOperationFeito(MPMesh2* mesh1, MPMesh2* mesh2, int op, HANDLE hd)
{
    _output = hd;
    StdOutput("Hi");

    

    return NULL;
}

extern "C" CSG_API MPMesh2* Convert2MPMesh2(const GS::BaseMesh* mesh)
{
    return NULL;
}

extern "C" CSG_API GS::BaseMesh* Convert2BaseMesh2(const MPMesh2* mesh)
{
    return NULL;
}


}
