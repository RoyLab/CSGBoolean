#include "NewMeshBoolOp.h"
#include "Bool.h"
#include <Windows.h>
extern bool bDebug;

namespace GS
{
NewMeshBoolOp::NewMeshBoolOp(void)
{
}


NewMeshBoolOp::~NewMeshBoolOp(void)
{
}

BoolOp* NewMeshBoolOp::GetInstance()
{
    static NewMeshBoolOp boolOp;
    return &boolOp;
}

BaseMesh* NewMeshBoolOp::DoCompute(CSGExprNode* node)
{
    return CSG::BooleanOperation(node, GetStdHandle(STD_OUTPUT_HANDLE), bDebug);
}

}
