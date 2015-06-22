#include "precompile.h"
#include "Bool.h"
#include "BaseMesh.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include "MPMesh.h"
#include <ctime>
#include "configure.h"
#include "COctree.h"
#include "IsectTriangle.h"
#include "isect.h"

namespace CSG
{

extern HANDLE _output;
extern clock_t t0;


void DebugInfo(char* str, clock_t& t0);
void StdOutput(const char* str);


class FeitoISectZone
{
public:
    struct ISVertexInfo;
	typedef std::list<ISVertexInfo>::iterator	ISVertexItr;

	struct ISVertexInfo
	{
        ISVertexInfo(){reset();}
        bool is_valid() {return pos > -1;}
        void reset() {pos = -1;}
		int	 pos; // -1 is invalid
		ISVertexItr	next;
	};

    struct SegData
    {
        MPMesh2::FaceHandle face;
        ISVertexItr points[2];
    };

    struct VertexData
    {
        Vec3d coord;
        bool used;
        VertexData():used(false){}
        VertexData(const Vec3d& v):used(false), coord(v){}
    };

    FeitoISectZone(MPMesh2* mesh, MPMesh2::FaceHandle fh):
    face(fh), thisMesh(mesh)
    {
		// 三个角点在vertices的最后三个
        int id[3];
		auto fvItr = thisMesh->fv_begin(face);
        id[0] = (fvItr++)->idx();
        id[1] = (fvItr++)->idx();
        id[2] = (fvItr)->idx();

        ISVertexInfo info;

        vertices.push_back(info);
		corner[0] = --vertices.end();
        vertices.push_back(info);
		corner[1] = --vertices.end();
        vertices.push_back(info);
		corner[2] = --vertices.end();

		// add to result mesh.  #WR#
		MPMesh2::VertexHandle tmp;
		for (int i = 0; i < 3; i++)
		{
            auto &prop = thisMesh->property(thisMesh->VertexIndexPropHandle, thisMesh->vertex_handle(id[i]));
			if (prop == -1) 
				prop = add_point(thisMesh->point(thisMesh->vertex_handle(id[i])));
			corner[i]->pos = prop.value;
		}
		zone_record.push_back(this);
    }

    ~FeitoISectZone(){}

    	void InsertSegment(ISVertexItr v0, ISVertexItr v1, FeitoISectZone* tri2)
	{
		SegData seg;
        seg.points[0] = v0;
        seg.points[1] = v1;
		seg.face = tri2->face;
		segs.push_back(seg);
	}

	std::list<MPMesh2::FaceHandle> coplanarTris;
	std::list<SegData> segs;
    std::list<ISVertexInfo> vertices;
    MPMesh2::FaceHandle face;

    ISVertexItr corner[3];
    MPMesh2 *thisMesh;

private:
    int Id;

public:
    static void clearList(){verticesList.clear();}

    static void clear()
    {
        for (auto itr:zone_record) SAFE_RELEASE(itr);
        zone_record.clear();
        clearList();
    }

    static const VertexData& point(unsigned id){return verticesList[id];}

    static unsigned add_point(Vec3d& vec)
    {
        verticesList.emplace_back(vec);
        return verticesList.size()-1;
        return 0;
    }


    static std::vector<FeitoISectZone*> zone_record;

private:
    static std::vector<VertexData> verticesList;

};

std::vector<FeitoISectZone::VertexData> FeitoISectZone::verticesList;
std::vector<FeitoISectZone*> FeitoISectZone::zone_record;

inline void GetCorners(MPMesh2* pMesh, MPMesh2::FaceHandle fhandle, Vec3d*& v0, Vec3d *&v1, Vec3d *&v2)
{
	auto fvItr = pMesh->fv_begin(fhandle);
	v0 = &(pMesh->verticesList[fvItr->idx()]);fvItr++;
	v1 = &(pMesh->verticesList[fvItr->idx()]);fvItr++;
	v2 = &(pMesh->verticesList[fvItr->idx()]);
}

void GetLeafNodes(OctreeNode*, std::list<OctreeNode*>&, int);

inline bool CompareVertex(FeitoISectZone::ISVertexItr& ref, const Vec3d& vec)
{
	auto refcpy = ref;
    while (!refcpy->is_valid()) refcpy = refcpy->next;
	return IsEqual(FeitoISectZone::point(refcpy->pos).coord, vec);
}

bool IsVertexExisted(FeitoISectZone* tri, Vec3d& vec, FeitoISectZone::ISVertexItr& ref)
{
	// we do not check if it is a corner point #WR#
	// 但是，如果这个ref本身就在这个里面，那么问题就来了
	// 因为所有的插入点都是在前面的，所以不用检查到最后一个：a, b, c, d, v0, v1, v2.
	const auto end = tri->corner[0];
	for (auto itr = tri->vertices.begin(); itr != end; itr++)
	{
		if (CompareVertex(itr, vec))
		{
			ref = itr;
			return true;
		}
	}
	return false;
}

FeitoISectZone::ISVertexItr InsertPoint(FeitoISectZone* tri, VertexPos pos, FeitoISectZone::ISVertexItr ref)
{
	// 追溯到最源头的迭代器
	while (!ref->is_valid()) ref = ref->next;

    Vec3d vec = FeitoISectZone::point(ref->pos).coord;
	FeitoISectZone::ISVertexItr output = ref;
	FeitoISectZone::ISVertexInfo info;
	info.next = ref;
	if (pos == INNER)
	{
		if (!IsVertexExisted(tri, vec, output))
		{
			tri->vertices.push_front(info);
			output = tri->vertices.begin();
		}
		else
		{
			if (!output->is_valid() || output->pos != ref->pos)
			{
				output->reset();
				output->next = ref;
			}
		}
	}
	else if (pos <= EDGE_2)
	{
		auto mesh = tri->thisMesh;
		auto ffItr = mesh->ff_begin(tri->face);
		switch (pos)
		{
		case CSG::EDGE_0:
			ffItr ++;
			ffItr ++;
			break;
		case CSG::EDGE_1:
			break;
		case CSG::EDGE_2:
			ffItr ++;
			break;
		default:
			assert(0);
		}

		FeitoISectZone*& other = mesh->property(mesh->SurfacePropHandle, *ffItr);
		if (!other) other = new FeitoISectZone(mesh, *ffItr);
		InsertPoint(other, INNER, ref);
		return InsertPoint(tri, INNER, ref);
	}
	else
	{
		// get the iterator of vertex #WR#
		switch (pos)
		{
		case CSG::VER_0:
			if (tri->corner[0]->pos == ref->pos) break;
			tri->corner[0]->reset();
			tri->corner[0]->next = ref;
			output = tri->corner[0];
			break;
		case CSG::VER_1:
			if (tri->corner[1]->pos == ref->pos) break;
			tri->corner[1]->reset();
			tri->corner[1]->next = ref;
			output = tri->corner[1];
			break;
		case CSG::VER_2:
			if (tri->corner[2]->pos == ref->pos) break;
			tri->corner[2]->reset();
			tri->corner[2]->next = ref;
			output = tri->corner[2];
			break;
		}
	}
		
	return output;
}

FeitoISectZone::ISVertexItr InsertPoint(FeitoISectZone* tri, VertexPos pos, Vec3d& vec)
{
	FeitoISectZone::ISVertexItr output;
	if (pos == INNER)
	{
		if (!IsVertexExisted(tri, vec, output))
		{
			FeitoISectZone::ISVertexInfo info;
            info.pos = FeitoISectZone::add_point(vec);
			tri->vertices.push_front(info);
			output = tri->vertices.begin();
		}
	}
	else if (pos <= EDGE_2)
	{
		auto mesh = tri->thisMesh;
		auto ffItr = mesh->ff_begin(tri->face);
		switch (pos)
		{
		case CSG::EDGE_0:
			ffItr ++;
			ffItr ++;
			break;
		case CSG::EDGE_1:
			break;
		case CSG::EDGE_2:
			ffItr ++;
			break;
		default: 
			assert(0);
		}

		FeitoISectZone*& other = mesh->property(mesh->SurfacePropHandle, *ffItr);
#ifdef _DEBUG_
		auto fvItr = mesh->fv_begin(tri->face);
		if (pos == EDGE_1) fvItr++;
		if (pos == EDGE_2) {fvItr++; fvItr++;}
			
		auto f0 = *fvItr;
		fvItr = mesh->fv_begin(*ffItr);
		assert(*fvItr != f0); fvItr++;
		assert(*fvItr != f0); fvItr++;
		assert(*fvItr != f0);
#endif
		if (!other) other = new FeitoISectZone(mesh, *ffItr);
		output = InsertPoint(other, INNER, vec);
		return InsertPoint(tri, INNER, output);
	}
	else
	{
		// get the iterator of vertex #WR#
		switch (pos)
		{
		case CSG::VER_0:
			output = tri->corner[0];
			break;
		case CSG::VER_1:
			output = tri->corner[1];
			break;
		case CSG::VER_2:
			output = tri->corner[2];
			break;
		}
	}
	return output;
}

void ISectTest2(Octree<MPMesh2>* pOctree, bool bInverse)
{
	assert(pOctree);
	std::list<OctreeNode*> leaves;
	GetLeafNodes(pOctree->Root, leaves, NODE_COMPOUND);

	std::set<GS::IndexPair> antiOverlapMap;
	for (auto leaf: leaves)
	{
		unsigned i, j, ni, nj;
		MPMesh2 *mesh1 = pOctree->pMesh[0];
		MPMesh2 *mesh2 = pOctree->pMesh[1];
		MPMesh2::FaceHandle tri1, tri2;
		Vec3d *v0,*v1,*v2, nv,*u0,*u1,*u2,nu,start,end;
		MPMesh::FVIter fvItr;
		int isISect;
		int startT(0), endT(0);
		VertexPos startiT, startjT, endiT, endjT;
		FeitoISectZone::ISVertexItr vP1, vP2;

        int triId[2];
        GS::IndexPair iPair;
        auto &tt1 = leaf->TriangleTable[0];
        auto &tt2 = leaf->TriangleTable[1];
        ni = tt1.size();
        nj = tt2.size();

		for (i = 0; i < ni; i++)
		{
			for (j = 0; j < nj; j++)
			{
				tri1 = tt1[i];
				tri2 = tt1[j];
                triId[0] = tri1.idx();
                triId[1] = tri2.idx();
				GS::MakeIndex(triId, iPair);
				if (antiOverlapMap.find(iPair) != antiOverlapMap.end()) continue;
				else antiOverlapMap.insert(iPair);

				// intersection test main body
				GetCorners(mesh1, tri1, v0, v1, v2);
				GetCorners(mesh2, tri2, u0, u1, u2);

				nv = mesh1->normal(tri1);
				nu = mesh2->normal(tri2);
							
				startT = 0; endT = 0; // return to Zero.

				isISect = TriTriIntersectTest(*v0, *v1, *v2, nv,
					*u0, *u1, *u2, nu, startT, endT, start, end);

				if (isISect < 0) continue;

				FeitoISectZone **si = &mesh1->property(mesh1->SurfacePropHandle, tri1);
				FeitoISectZone **sj = &mesh2->property(mesh2->SurfacePropHandle, tri2);

				if (!*si) *si = new FeitoISectZone(mesh1, tri1);
				if (!*sj) *sj = new FeitoISectZone(mesh2, tri2);

				if (isISect == 0)
				{
					(*si)->coplanarTris.emplace_back(tri2);
					(*sj)->coplanarTris.emplace_back(tri1);
					continue;
				}

				startiT = VertexPos(startT & 0xffff);
				startjT = VertexPos(startT >> 16);

				endiT = VertexPos(endT & 0xffff);
				endjT = VertexPos(endT >> 16);

				if (IsEqual(start, end))
				{
					// 点相交
					Vec3d point = (start+end)/2;

					// 最后一个参数表示，可能存在两个以上的插入点
					vP1 = InsertPoint(*si, startiT, point);
					InsertPoint(*si, endiT, vP1);
					InsertPoint(*sj, startjT, vP1);
					InsertPoint(*sj, endjT, vP1);
				}
				else
				{
					// 线相交
					double d = OpenMesh::dot(OpenMesh::cross(nv, nu), end-start);
					if (bInverse) d = -d;
					if (d > 0)
					{
						vP1 = InsertPoint(*si, startiT, start);
						vP2 = InsertPoint(*si, endiT, end);
						(*si)->InsertSegment(vP1, vP2, *sj);
						vP1 = InsertPoint(*sj, startjT, vP1);
						vP2 = InsertPoint(*sj, endjT, vP2);
						(*sj)->InsertSegment(vP2, vP1, *si);
					}
					else
					{
						vP1 = InsertPoint(*si, startiT, start);
						vP2 = InsertPoint(*si, endiT, end);
						(*si)->InsertSegment(vP2, vP1, *sj);
						vP1 = InsertPoint(*sj, startjT, vP1);
						vP2 = InsertPoint(*sj, endjT, vP2);
						(*sj)->InsertSegment(vP1, vP2, *si);
					}
				}
			}
		}
	}
}


// 0 Union:
// 1 Intersect
// 2 Diff 



extern "C" CSG_API MPMesh2* BooleanOperationFeito(MPMesh2* mesh1, MPMesh2* mesh2, int op, HANDLE hd)
{
    _output = hd;
    StdOutput("Hi");

    MPMesh2* res = nullptr;
    MPMesh2 *arrMesh[2] = {mesh1, mesh2};
    Octree<MPMesh2>* pOctree = BuildOctree2(arrMesh, 2);
    ISectTest2(pOctree, (op==2)?true:false);
    //Tesselation();
    //Classification();

    FeitoISectZone::clear();
    delete pOctree;
    return res;
}

inline OpenMesh::Vec3d convert_double3(GS::double3& vec)
{
	return OpenMesh::Vec3d(vec.x, vec.y, vec.z);
}

extern "C" CSG_API MPMesh2* Convert2MPMesh2(const GS::BaseMesh* pMesh)
{
    if (!pMesh) return nullptr;

	MPMesh2 *res = new MPMesh2(pMesh);

    int n = (int)pMesh->VertexCount();
    res->BBox.Clear();
	MPMesh::VertexHandle vhandle;
	res->verticesList = new Vec3d[n];
    for (int i = 0; i < n; i++)
	{
		res->verticesList[i] = convert_double3(pMesh->Vertex(i).pos);
        vhandle = res->add_vertex(res->verticesList[i]);
		//res->property(res->PointInOutTestPropHandle, vhandle) = 0;
		res->BBox.IncludePoint(res->point(vhandle));
	}

	std::vector<MPMesh::VertexHandle>  face_vhandles(3);
    n = (int)pMesh->PrimitiveCount();
    for (int i = 0; i < n; i++)
    {
        auto &info = pMesh->TriangleInfo(i);
		face_vhandles[0] = res->vertex_handle(info.VertexId[0]);
		face_vhandles[1] = res->vertex_handle(info.VertexId[1]);
		face_vhandles[2] = res->vertex_handle(info.VertexId[2]);
        res->add_face(face_vhandles);
    }

	res->update_normals();

	return res;
}

extern "C" CSG_API GS::BaseMesh* Convert2BaseMesh2(MPMesh2* mesh)
{
    if (!mesh) return nullptr;
    GS::BaseMesh* res = new GS::BaseMesh;

    for (auto face = mesh->faces_begin(); face != mesh->faces_end(); face++)
    {
		GS::double3 v[3];
		Vec3d *v0, *v1, *v2;
        MPMesh2::FaceHandle face2 = face;
		GetCorners(mesh, face2, v0, v1, v2);
		v[0] = Vec3dToDouble3(*v0);
		v[1] = Vec3dToDouble3(*v1);
		v[2] = Vec3dToDouble3(*v2);
		res->AddTriangle(v);
    }

    return NULL;
}


}
