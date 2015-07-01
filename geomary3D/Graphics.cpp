#include "stdafx.h"
#include "Graphics.h"
#include "SceneParser.h"
#include "typedefs.h"
#include "Global.h"
#include "BoolOp.h"
#include <fstream>

std::vector<std::string> expression;

bool GetLoop(std::string& s)
{
    bool flag = false;
    char c[20];
    int end;
    std::ifstream script("D:\\bool\\boolnum.txt");
	if (!script) return false;
    script >> end;
    script.close();
	if (end == 0) return false;

	end /= 20;
	s = "0";
    for (int i = 1 ; i < end+1; i++)
    {
		for (int j = 0; j < 20; j++)
		{
			if (flag) s += "+";
			else s += "-";
			sprintf_s(c, 20, "%d", i  + j* 40);
			s += c;
		}
		flag = !flag;
	}
	return true;
}

CGraphics::CGraphics()
	:mViewPortSize(0, 0)
	, mCursorPos(0, 0)
{
	m_pD3D =0;
    std::ifstream script("D:\\bool\\boolconfig.txt");
    if (!script) return;
    char eval[1000];
    while (!script.eof())
    {
        script.getline(eval, 1000);
        if (strlen(eval))
            expression.emplace_back(eval);
    }
    script.close();
	expression.emplace_back();
    if (!GetLoop(expression.back())) expression.pop_back();
}

CGraphics::~CGraphics()
{
	if (m_pD3D)
	{
		delete m_pD3D;
		m_pD3D = 0 ; 
	}
}

void CGraphics::Shutdown()
{
	mShader.Release();
	if (m_pD3D)
	{
		m_pD3D->Shutdown();
		delete m_pD3D;
		m_pD3D =0;
	}
	GS::Global::Shutdown();
}

bool CGraphics::Initialize(int screenWidth, int screenHeight, HWND hwnd)
{
	bool result;
	GS::Global::Initialize();
	mViewPortSize = GS::int2(screenWidth,screenHeight );
	// Create the Direct3D object.
	m_pD3D = new GS::VirtualDeviceDx11();
	if(!m_pD3D)
	{
		return false;
	}
    if (screenWidth == 0)
		screenWidth =1;
	if (screenHeight == 0)
		screenHeight =1;
	// Initialize the Direct3D object.
	result = m_pD3D->Initialize(screenWidth, screenHeight, VSYNC_ENABLED, hwnd, false, SCREEN_DEPTH, SCREEN_NEAR);
	if(!result)
	{
		MessageBox(hwnd, L"Could not initialize Direct3D.", L"Error", MB_OK);
		return false;
	}
	result = mShader.Initialize(m_pD3D->GetDevice(), hwnd);
	if (!result)
	{
		MessageBox(hwnd, L"Could not initialize Shader.", L"Error", MB_OK);
		return false;
	}
    return true;
}

void  CGraphics::ResizeWindow(int width, int height)
{

	if (width == 0 || height ==0 )
		return ; 
	mViewPortSize = GS::int2(width, height);
	if (m_pD3D)
		m_pD3D->Resize(width, height);

    //mCamera.SetProjectionMode(true);
    //mCamera.SetPosition(GS::float3(0,0,-60), GS::float3(0,0,0), GS::float3(0,1,0));
    mCamera.ComputeFovProjectMatrix(PI/4.0f, (float)width/height, 0.01f, 1000.0f);
}

inline void recordVec4(std::ofstream& f, GS::float4& vec)
{
    f << vec.x << '\t';
    f << vec.y << '\t';
    f << vec.z << '\t';
    f << vec.w << '\n';
}

inline void recoverVec4(std::ifstream& f, GS::float4& vec)
{
    f >> vec.x;
    f >> vec.y;
    f >> vec.z;
    f >> vec.w;
}

void CGraphics::SnapCam(int i)
{
    char file[32];
    sprintf_s(file, "D:\\bool\\cam%d.cam", i);
    std::ofstream f(file);
    if (!f) return;
    auto vec = mCamera.Eye();
    f << vec.x << '\t';
    f << vec.y << '\t';
    f << vec.z << '\n';

    vec = mCamera.Target();
    f << vec.x << '\t';
    f << vec.y << '\t';
    f << vec.z << '\n';

    f << mCamera.mNearPlane << '\t';
    f << mCamera.mFarPlane << '\t';
    f << mCamera.mOrthoRange << '\t';
    f << mCamera.mbPerspective << '\n';

    auto vec2 = mCamera.mViewMatrix[0];
    recordVec4(f, vec2);
    vec2 = mCamera.mViewMatrix[1];
    recordVec4(f, vec2);
    vec2 = mCamera.mViewMatrix[2];
    recordVec4(f, vec2);
    vec2 = mCamera.mViewMatrix[3];
    recordVec4(f, vec2);

    vec2 = mCamera.mProjectionMatrix[0];
    recordVec4(f, vec2);
    vec2 = mCamera.mProjectionMatrix[1];
    recordVec4(f, vec2);
    vec2 = mCamera.mProjectionMatrix[2];
    recordVec4(f, vec2);
    vec2 = mCamera.mProjectionMatrix[3];
    recordVec4(f, vec2);
    
    wchar_t str[] = L"Write Cam Succeeded!\n";
    WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), str, wcslen(str), 0, 0);

    f.close();
}

void CGraphics::RecoverCam(int i)
{
    char file[32];
    sprintf_s(file, "D:\\bool\\cam%d.cam", i);
    std::ifstream f(file);
    if (!f) return;
    GS::float3 eye;
    f >> eye.x;
    f >> eye.y;
    f >> eye.z;
    mCamera.mEye = eye;

    GS::float3 target;
    f >> target.x;
    f >> target.y;
    f >> target.z;
    mCamera.mTarget = target;

    f >> mCamera.mNearPlane;
    f >> mCamera.mFarPlane;
    f >> mCamera.mOrthoRange;
    f >> mCamera.mbPerspective;

    recoverVec4(f, mCamera.mViewMatrix[0]);
    recoverVec4(f, mCamera.mViewMatrix[1]);
    recoverVec4(f, mCamera.mViewMatrix[2]);
    recoverVec4(f, mCamera.mViewMatrix[3]);

    recoverVec4(f, mCamera.mProjectionMatrix[0]);
    recoverVec4(f, mCamera.mProjectionMatrix[1]);
    recoverVec4(f, mCamera.mProjectionMatrix[2]);
    recoverVec4(f, mCamera.mProjectionMatrix[3]);

    wchar_t str[] = L"Import Cam Succeeded!\n";
    WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), str, wcslen(str), 0, 0);

    f.close();
}

bool CGraphics::Frame()
{
	bool result;
	
	if (!m_pD3D)
		return false; 

	// Render the graphics scene.
	result = Render();
	if(!result)
	{
		return false;
	}

	return true;
}

extern int _renderState;

bool CGraphics::Flip(float x, float y)
{
    auto dir = mCamera.Direction();
    auto up = mCamera.Up();
    auto right = GS::normalize(GS::cross(dir, up));
    auto up2 = GS::normalize(GS::cross(right, dir));

    auto dright = right*0.01*x;
    auto dup = up2*0.01*y;
    auto diff = dright+dup;

    mCenterPos += diff;
    mCamera.SetPosition(mCamera.Eye()+diff, mCamera.Target()+diff);
    return true;
}
extern bool bHintLine;

bool CGraphics::Render()
{
	if(!m_pD3D)
		return false;
	GS::float4x4 worldMatrix (GS::id4x4());
	// Clear the buffers to begin the scene.
    GS::float4x4 worldMat = GS::inverse(mCamera.ViewMatrix());

	const GS::float4& diffuse1 = mLightMgr.DefaultLight(0).DiffuseColor();
	const GS::float4& specular1 = mLightMgr.DefaultLight(0).SpecularColor();
	GS::float3 lightDir1 = normalize(mLightMgr.DefaultLight(0).Direction());

	const GS::float4& diffuse2 = mLightMgr.DefaultLight(1).DiffuseColor();
	GS::float3 lightDir2 = normalize(mLightMgr.DefaultLight(1).Direction());
	ComputeLightDirection(mLightMgr.DefaultLight(1), worldMat, lightDir2);
	float SpecularFactor = mLightMgr.DefaultLight(0).SpecularFactor();
	m_pD3D->BeginScene(1.0f, 1.0f, 1.0f, 1.0f);
	mShader.SetShaderMatrix(m_pD3D->GetDeviceContext(), worldMatrix, mCamera.ViewMatrix(), mCamera.ProjectionMatrix(), mCamera.Eye());
	mShader.SetShaderLights(m_pD3D->GetDeviceContext(),mLightMgr.AmbientColor(), diffuse1, specular1, lightDir1, diffuse2, lightDir2, SpecularFactor);

    bool Pass = (_renderState == 1)?false:true;
	mShader.Render( m_pD3D->GetDeviceContext(), Pass);
    m_pD3D->SetRenderState(Pass);
	mModelMgr.Render(m_pD3D->GetDevice(), m_pD3D->GetDeviceContext());

    if (_renderState == 2)
    {
        Pass = false;
	    mShader.Render( m_pD3D->GetDeviceContext(), Pass);
        m_pD3D->SetRenderState(Pass);

        if (bHintLine)
        {
            mModelMgr.RenderOrigList(m_pD3D->GetDevice(), m_pD3D->GetDeviceContext());
        }
        else mModelMgr.Render(m_pD3D->GetDevice(), m_pD3D->GetDeviceContext());
    }

	// Present the rendered scene to the screen.  
    m_pD3D->EndScene();  
  
    return true; 
}

inline
void CGraphics::ComputeLightDirection(const GS::DirectionLight& light, const GS::float4x4& worldMat, GS::float3& Lightdir)
{

	Lightdir = GS::mul(worldMat,GS::point(light.Direction())).xyz;
   // Lightdir = light.Direction();
	Lightdir = GS::normalize(Lightdir);
	
	
}

void CGraphics::ShowOctree(bool bShow)
{
    if (!bShow)
        mModelMgr.ClearResults();
    else 
        mModelMgr.CreateOctree();


}


void CGraphics::Union()
{
    mModelMgr.BoolOperation(eUnion, GS::eMeshBool);
}

void CGraphics::Difference()
{
    mModelMgr.BoolOperation(eDifference, GS::eMeshBool);
}

void CGraphics::Intersect()
{
    mModelMgr.BoolOperation(eIntersect, GS::eMeshBool);
}

 void CGraphics::EvaluateBoolExpression()
 {
     bool flag = false;
     static unsigned count = 0;
     mModelMgr.EvaluateBoolExpression(expression[count]);
     count = (count+1)%expression.size();
 }


void  CGraphics::BSPUnion()
{
    mModelMgr.BoolOperation(eUnion, GS::eBSPBool);
}

void  CGraphics::BSPIntersect()
{
    mModelMgr.BoolOperation(eIntersect, GS::eBSPBool);
}


void  CGraphics::BSPDifference()
{
    mModelMgr.BoolOperation(eDifference, GS::eBSPBool);
}


void CGraphics::LocalizedBSPUnion()
{
    mModelMgr.BoolOperation(eUnion, GS::eLBSPBool);
}
    
void CGraphics::LocalizedBSPIntersect()
{
    mModelMgr.BoolOperation(eIntersect, GS::eLBSPBool);
}
    
void CGraphics::LocalizedBSPDifference()
{
    mModelMgr.BoolOperation(eDifference, GS::eLBSPBool);
}


bool CGraphics::LoadSceneFromStream(char* fbxFileName)
{

	GS::float3 pos(0, 0, -10);
	GS::float3 target(0, 0, 1);
	GS::float3 up(0,1, 0);
	//GS::float3 pos(31.696, -85.373, 105.02);
	//GS::float3 target(20.960, 14.08, 3.013);
	//GS::float3 up(-0.07, 0.709, 0.700);
	mCamera.SetPosition(pos, target, up);
	SceneParser* parser = SceneParser::GetParser(fbxFileName);
	if (parser == NULL)
		return false ;
	if (!parser->Initialize())
		return false;
	parser->LoadScene(fbxFileName);
	std::vector<GS::BaseMesh*>  meshList; 
	std::vector<GS::Light*>     lights; 
	parser->ProcessScene(meshList, lights, mCamera);
	mModelMgr.Clear();
	mBoundingBox.Clear();
	for (int i = 0; i < meshList.size(); i++)
	{

		mModelMgr.Add(meshList[i]);
		mBoundingBox.IncludeBox(meshList[i]->AABB());
	}
    //mCamera.SetOrthoRange(mBoundingBox);
	mModelMgr.Invalidate();
	meshList.clear();
	lights.clear();
	delete parser;
	return true;
}

void CGraphics::ZoomView(short zDelta)
{
	 
		 float fPercent = 60/ 100.0;
        if (zDelta > 0) 
            fPercent = -fPercent / 4.0;
        fPercent += 1.0;  
		float fieldWidth = mCamera.HorizontalExtent();
		float fieldHeight = mCamera.VerticalExtent();
        fieldWidth  *= fPercent;
        fieldHeight *= fPercent;
		mCamera.SetExtends(fieldWidth, fieldHeight);
	
}

void CGraphics::StartOrbitView(int x, int y)
{
	mCursorPos.x  = x;
	mCursorPos.y = FlipYAxis(y);
	mCenterPos = mCamera.Target();

}

void CGraphics::UpdateOrbitView(int x, int y)
{
	int  dx = x - mCursorPos.x;
	int  dy = FlipYAxis(y)- mCursorPos.y;
	
	int sx = GS::min(mViewPortSize.x,  mViewPortSize.y);
	if (sx < 1) {
		 sx = 1;
     }
    double pixelAngleScale = PI_2 / ((double)(sx) * 0.4) ;
    double ax = (double)(dx) * pixelAngleScale;
    double ay = -(double)(dy) * pixelAngleScale;
	mCamera.Orbit(mCenterPos, ax, ay);
    //wchar_t str[32];
    //swprintf_s(str, L"%f and %f\n", ax, ay);
    //OutputDebugString(str);
	mCursorPos.x  = x;
    mCursorPos.y = FlipYAxis(y);
      
}

inline 
int CGraphics::FlipYAxis(int y)
{
	return mViewPortSize.y - y -1;
}
