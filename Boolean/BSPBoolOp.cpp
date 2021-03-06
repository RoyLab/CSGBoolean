#include "BSPBoolOp.h"
#include "PlaneMesh.h"
#include "BSPTree.h"
#include "OctTree.h"
#include "FixedPlaneMesh.h"
#include "FixedBSPTree.h"
#include <iostream>
#include <ctime>
#include "BSPOctTree.h"
#include <stack>
#include <cctype>
#include <Bool.h>

namespace GS{

static void GSOutputTimeLog(const wchar_t* ch)
{
#ifdef _DEBUG
	std::wstring debug;
	debug += ch;
	debug += std::to_wstring(clock());
	debug += L"\n";
	OutputDebugString(debug.c_str());
#endif 
}


BSPBoolOp:: ~BSPBoolOp()
{

}


BoolOp* BSPBoolOp::GetInstance()
{
    static BSPBoolOp boolOp;
    return &boolOp;
}
//#define SWITCH_WR
#ifdef SWITCH_WR
BaseMesh* BSPBoolOp::ComputeBoolean(BaseMesh* mesh1, BaseMesh* mesh2, BOOL_OP op)
{
    PlaneMesh* pPlaneMesh1 = mesh1->ToPlaneMesh();
    PlaneMesh* pPlaneMesh2 = mesh2->ToPlaneMesh();
    BaseMesh* pMesh = ComputeBoolean(pPlaneMesh1, pPlaneMesh2, op);
    delete pPlaneMesh1;
    delete pPlaneMesh2;
    return pMesh;
}
#else
BaseMesh* BSPBoolOp::ComputeBoolean(BaseMesh* mesh1, BaseMesh* mesh2, BOOL_OP op)
{
	Box3 bbox(mesh1->AABB());
	bbox.IncludeBox(mesh2->AABB());
    
    mesh1->NormalizeCoord(&bbox);
    mesh2->NormalizeCoord(&bbox);

    FixedPlaneMesh* pPlaneMesh1 = ToFixedPlaneMesh(mesh1);
    FixedPlaneMesh* pPlaneMesh2 = ToFixedPlaneMesh(mesh2);
    BaseMesh* pMesh = ComputeBoolean(pPlaneMesh1, pPlaneMesh2, op);
    delete pPlaneMesh1;
    delete pPlaneMesh2;

    mesh1->DenormalizeCoord();
    mesh2->DenormalizeCoord();
    return pMesh;
}
#endif
BaseMesh* BSPBoolOp::ComputeBoolean(FixedPlaneMesh* mesh1, FixedPlaneMesh* mesh2, BOOL_OP op)
{
	GSOutputTimeLog(L"start: ");
    FixedBSPTree::SET_OP BSPOp;
    switch (op)
    {
        case eUnion:
            BSPOp = FixedBSPTree::OP_UNION;
            break;
        case eIntersect:
            BSPOp = FixedBSPTree::OP_INTERSECT;
            break;
        case eDiff:
            BSPOp = FixedBSPTree::OP_DIFFERENCE;
            break;
        default :
            break;
    }

	Box3 bbox(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5);
    FixedBSPTree* pTree1 =  mesh1->ToBSPTree();
	GSOutputTimeLog(L"tree1, build: ");

    pTree1->FormSubHyperPlane(bbox);
	//pTree1->OutputDebugInfo("D:\\x.txt");
     
	GSOutputTimeLog(L"tree1, sp: ");

    FixedBSPTree* pTree2 = mesh2->ToBSPTree();
	GSOutputTimeLog(L"tree2, build: ");
    if (BSPOp == FixedBSPTree::OP_DIFFERENCE)
    {
        pTree2->Negate();
        BSPOp = FixedBSPTree::OP_INTERSECT;
    }
    pTree2->FormSubHyperPlane(bbox);
	GSOutputTimeLog(L"tree2, sp: ");
	//pTree2->OutputDebugInfo("D:\\y.txt");

    FixedBSPTree* pResultTree = pTree1->Merge(pTree2, BSPOp);
	GSOutputTimeLog(L"merged: ");
	//pResultTree->OutputDebugInfo("D:\\z.txt");
    delete pTree2;
    delete pTree1;
    //FixedPlaneMesh* pPlaneMesh1 = new FixedPlaneMesh(pTree1, mesh1->Color());//debug2

    // delete pTree1;//debug2
    FixedPlaneMesh* pPlaneMesh1 = new FixedPlaneMesh(pResultTree, mesh1->Color());
    delete pResultTree;
	GSOutputTimeLog(L"toPlaneMesh: ");

	pPlaneMesh1->SetAABB(mesh1->AABB());
    BaseMesh* pMesh = pPlaneMesh1->ToBaseMesh();
	GSOutputTimeLog(L"toBaseMesh: ");

    pMesh->DenormalizeCoord();
    pMesh->GenAABB(true);
    //BaseMesh* pMesh = mesh1->ToBaseMesh();//debug1
    if (pMesh && pMesh->PrimitiveCount() ==0 )
	{
        delete pMesh;
        pMesh = NULL;
    }
    delete pPlaneMesh1;
    return pMesh;
    return nullptr;
}

BaseMesh* BSPBoolOp::ComputeBoolean(PlaneMesh* mesh1, PlaneMesh* mesh2, BOOL_OP op)
{
    BSPTree::SET_OP BSPOp;
    switch (op)
    {
        case eUnion:
            BSPOp = BSPTree::OP_UNION;
            break;
        case eIntersect:
            BSPOp = BSPTree::OP_INTERSECT;
            break;
        case eDiff:
            BSPOp = BSPTree::OP_DIFFERENCE;
            break;
        default :
            break;
    }
	Box3 bbox(mesh1->AABB());
	bbox.IncludeBox(mesh2->AABB());

    BSPTree* pTree1 =  mesh1->ToBSPTree();
    pTree1->FormSubHyperPlane(bbox);
	pTree1->OutputDebugInfo("D:\\x1.txt");

    BSPTree* pTree2 = mesh2->ToBSPTree();
    if (BSPOp == BSPTree::OP_DIFFERENCE)
    {
        pTree2->Negate();
        BSPOp = BSPTree::OP_INTERSECT;
    }
    pTree2->FormSubHyperPlane(bbox);
	pTree2->OutputDebugInfo("D:\\y1.txt");

    BSPTree* pResultTree = pTree1->Merge(pTree2, BSPOp);
	pResultTree->OutputDebugInfo("D:\\z1.txt");
    delete pTree2;
    delete pTree1;

   PlaneMesh* pPlaneMesh1 = new PlaneMesh(pResultTree);
   //PlaneMesh* pPlaneMesh1 = new PlaneMesh(pTree1);
    //delete pTree1;
    BaseMesh* pMesh = pPlaneMesh1->ToBaseMesh();
    if (pMesh && pMesh->PrimitiveCount() ==0 )
	{
        delete pMesh;
        pMesh = NULL;
    }
    delete pPlaneMesh1;

    return pMesh;
}

BaseMesh*  BSPBoolOp::Evalute(std::vector<BaseMesh*>& meshList, std::string& postfix) 
{
    return NULL;
}

///////////////////////////////////////////////////////////////////////////////////

LBSPBoolOp:: ~LBSPBoolOp()
{
}


BoolOp* LBSPBoolOp::GetInstance()
{
    static LBSPBoolOp boolOp;
    return &boolOp;
}

BaseMesh* LBSPBoolOp::ComputeBoolean(BaseMesh* mesh1,  BaseMesh* mesh2, BOOL_OP op)
{
    FixedBSPTree::SET_OP BSPOp;
    switch (op)
    {
        case eUnion:
            if (!mesh1) return mesh2;
            if (!mesh2) return mesh1;
            BSPOp = FixedBSPTree::OP_UNION;
            break;
        case eIntersect:
            if (!mesh1 || !mesh2) return NULL;
            BSPOp = FixedBSPTree::OP_INTERSECT;
            break;
        case eDiff:
            if (!mesh1) return NULL;
            if (!mesh2) return mesh1;
            BSPOp = FixedBSPTree::OP_DIFFERENCE;
            break;
        default :
            break;
    }
    BSPOctree* bspOctree = new BSPOctree(BSPOp);
    BaseMesh* result = new BaseMesh;
    bspOctree->BSPOperation(mesh1, mesh2, &result);
    delete bspOctree;
    return result;
}

CSG::MPMesh2* LBSPBoolOp::ComputeBoolean2(CSG::MPMesh2* mesh1,  CSG::MPMesh2* mesh2, BOOL_OP op)
{
    int opnum = -1;
    switch (op)
    {
    case GS::eUnion:
        opnum = 0;
        break;
    case GS::eIntersect:
        opnum = 1;
        break;
    case GS::eDiff:
        opnum = 2;
        break;
    default:
        break;
    }
    return CSG::BooleanOperationFeito(mesh1, mesh2, opnum, GetStdHandle(STD_OUTPUT_HANDLE));
}

BaseMesh*  LBSPBoolOp::Evalute(std::vector<BaseMesh*>& meshList, std::string& postfix)
{
    return Evalute2(meshList, postfix);
}

BaseMesh*  LBSPBoolOp::Evalute2(std::vector<BaseMesh*>& meshList, std::string& postfix)
{
    std::stack<CSG::MPMesh2*> temp;

    auto c0 = clock();
	int k = 0;
    for (unsigned i = 0; i < postfix.size(); i ++)
    {
        if (postfix[i] == ' ') continue;
        if (isdigit(postfix[i]))
        {
            std::string index;
            index += postfix[i];
            while (isdigit(postfix[++i]))
                index += postfix[i];

            temp.push(CSG::Convert2MPMesh2(meshList[atoi(index.c_str())]));
            i--;
        }
        else 
        {
			wchar_t s[32];
			wsprintf(s, L"%d...", k++);
			WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), s, wcslen(s), 0, 0);
            CSG::MPMesh2 *B = temp.top();
            temp.pop();
            CSG::MPMesh2 *A = temp.top(), *res;
            temp.pop();
            switch (postfix[i])
            {
            case '+':
                res = ComputeBoolean2(A, B, eUnion);
                break;
            case '-':
                res = ComputeBoolean2(A, B, eDiff);
                break;
            case '*':
                res = ComputeBoolean2(A, B, eIntersect);
                break;
            default:
                assert(0);
                break;
            }
            ReleaseMPMesh2(A);
            ReleaseMPMesh2(B);
            temp.push(res);
        }
    }
    auto t = clock()-c0;
    wchar_t ch[32];
    wsprintf(ch, L"\n...total time...%d\n", t);
    WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), ch, wcslen(ch), 0, 0);
    BaseMesh* finalres = CSG::Convert2BaseMesh2(temp.top());
    ReleaseMPMesh2(temp.top());

    return finalres;
}

BaseMesh*  LBSPBoolOp::Evalute1(std::vector<BaseMesh*>& meshList, std::string& postfix) 
{
    std::stack<BaseMesh*> temp;
    std::list<BaseMesh*> garbage;

    auto c0 = clock();
	int k = 0;
    for (unsigned i = 0; i < postfix.size(); i ++)
    {
        if (postfix[i] == ' ') continue;
        if (isdigit(postfix[i]))
        {
            std::string index;
            index += postfix[i];
            while (isdigit(postfix[++i]))
                index += postfix[i];

            temp.push(meshList[atoi(index.c_str())]);
            i--;
        }
        else 
        {
			wchar_t s[32];
			wsprintf(s, L"%d...", k++);
			WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), s, wcslen(s), 0, 0);

            BaseMesh *B = temp.top();
            temp.pop();
            BaseMesh *A = temp.top(), *res;
            temp.pop();
            switch (postfix[i])
            {
            case '+':
                res = ComputeBoolean(A, B, eUnion);
                break;
            case '-':
                res = ComputeBoolean(A, B, eDiff);
                break;
            case '*':
                res = ComputeBoolean(A, B, eIntersect);
                break;
            default:
                assert(0);
                break;
            }

			auto itrA = std::find(garbage.begin(), garbage.end(), A);
			if (itrA != garbage.end())
			{
				delete *itrA;
				garbage.erase(itrA);
			}

			auto itrB = std::find(garbage.begin(), garbage.end(), B);
			if (itrB != garbage.end())
			{
				delete *itrB;
				garbage.erase(itrB);
			}

            temp.push(res);
            garbage.push_back(res);
        }
    }
    auto t = clock()-c0;
    wchar_t ch[32];
    wsprintf(ch, L"\n...total time...%d\n", t);
    WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), ch, wcslen(ch), 0, 0);
    BaseMesh* finalres = garbage.back();
    garbage.pop_back();

    for (auto itr: garbage)
    {
        if (itr) delete itr;
    }

    return finalres;
}


}

