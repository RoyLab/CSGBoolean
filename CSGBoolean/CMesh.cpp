#include "precompile.h"
#include "CMesh.h"
#include "BaseMesh.h"

namespace CSG
{
    CSGMesh::CSGMesh(GS::BaseMesh* pMesh):
        ID(-1), pOrigin(pMesh), bInverse(false)
    {
        int n = (int)pMesh->VertexCount();
        for (int i = 0; i < n; i++)
            mVertex.push_back(pMesh->Vertex(i).pos);

        n = (int)pMesh->PrimitiveCount();
        for (int i = 0; i < n; i++)
        {
            auto &index = pMesh->TriangleInfo(i);
            mTriangle.emplace_back(index.VertexId, index.Normal);
        }

        unsigned vn = mVertex.size();
        mAABB.Clear();
        for (unsigned i = 0; i < vn; i++)
            mAABB.IncludePoint(mVertex[i]);
    }

    CSGMesh::~CSGMesh(void)
    {
    }

    GS::BaseMesh* ConverteToBaseMesh(GS::ListOfvertices &vertices)
    {
        GS::BaseMesh* pResult = new GS::BaseMesh(false);
        const unsigned n = vertices.size();
        GS::double3 v[3];
        for (unsigned i = 0; i < n; i+=3)
        {
            v[0] = vertices[i].pos;
            v[1] = vertices[i+1].pos;
            v[2] = vertices[i+2].pos;
            pResult->AddTriangle(v);
        }
        
        return pResult;
    }

}
