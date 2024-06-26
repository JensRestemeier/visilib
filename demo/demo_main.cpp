/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2021 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/

#ifdef _WIN32
#include <windows.h>
#endif
#ifdef USE_GLUT
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include <imgui.h>
#include <examples/imgui_impl_glut.h>
#include <examples/imgui_impl_opengl2.h>

#include <string>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "xmmintrin.h"
#include "pmmintrin.h"

#include "helper_triangle_mesh_container.h"
#include "silhouette_container_embree.h"
#include "helper_synthetic_mesh_builder.h"
#include "helper_geometry_scene_reader.h"
#include "demo_viewer_glut.h"
#include "demo_debug_visualisation_gl.h"
#include "demo_helper.h"

using namespace visilib;
using namespace std;
using namespace visilibDemo;

visilib::SimpleStats visilib::stats;

namespace visilibDemo
{
    class VisilibDemoMain;

    class VisilibDemoMain
    {
        bool show_demo_window;
        int wc = 10;
        int rc = 0;
        std::string output;
    public:
        VisilibDemoMain()
        {
            show_demo_window = false;
        }
        ~VisilibDemoMain()
        {
            printf("done\n");
        }

        bool init()
        {
            _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
            _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
            readConfig("config.txt");
            forceDisplay = true;
#ifdef USE_GLUT
            gluLookAt(1, 1, 1,
                0, 0, 0,
                1, 0, 0);
#endif
            return initScene(mDemoConfiguration.sceneIndex);
        }
        bool initScene(int s)
        {
            srand(0);

            delete debugger;
            debugger = new HelperVisualDebugger();

            delete meshContainer;
            meshContainer = DemoHelper::createScene(s, mDemoConfiguration.globalScaling);
            
            delete occluderSet;
            occluderSet = DemoHelper::createOccluderSet(meshContainer);

            return true;
        }

        void resolveVisibility()
        {
            VisibilityExactQueryConfiguration config;
            config.silhouetteOptimization = mDemoConfiguration.silhouetteOptimisation;
            config.hyperSphereNormalization = mDemoConfiguration.normalization;
            config.precision = mDemoConfiguration.precisionType;
            config.representativeLineSampling = mDemoConfiguration.representativeLineSampling;
            config.detectApertureOnly = mDemoConfiguration.detectApertureOnly;
            config.resolveMode = mDemoConfiguration.resolveMode;
#if EMBREE 
            config.useEmbree = mDemoConfiguration.embree;
#endif
            result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, debugger);
        }

        void animate()
        {
            if (animated || step)
            {
                mDemoConfiguration.phi += 0.005f;
                mDemoConfiguration.eta += 0.001f;
                forceDisplay = true;
                step = false;
            }

            if (forceDisplay)
            {
                stats.Reset();
#if 1
                std::stringstream out;

                std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
                std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
#endif

                DemoHelper::generatePolygon(v0, mDemoConfiguration.vertexCount0, mDemoConfiguration.scaling, mDemoConfiguration.phi - (float)M_PI, mDemoConfiguration.globalScaling);
                DemoHelper::generatePolygon(v1, mDemoConfiguration.vertexCount1, mDemoConfiguration.scaling, mDemoConfiguration.phi, mDemoConfiguration.globalScaling);

                resolveVisibility();

#if 1
                std::cout.rdbuf(coutbuf);
                output = out.str();
#endif

                forceDisplay = false;
            }
#ifdef USE_GLUT
            glutPostRedisplay();
#endif
        }

        void writeConfig(const std::string & filename)
        {
            mDemoConfiguration.writeConfig(filename);
        }

        void readConfig(const std::string & filename)
        {
            mDemoConfiguration.readConfig(filename);
        }

        void nextPrecisionType()
        {
            switch (mDemoConfiguration.precisionType)
            {
            case VisibilityExactQueryConfiguration::FLOAT:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::DOUBLE;
                break;
#define NEXT_ITEM VisibilityExactQueryConfiguration::DOUBLE
#ifdef EXACT_ARITHMETIC
            case NEXT_ITEM:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::EXACT;
                break;
#undef NEXT_ITEM
#define NEXT_ITEM VisibilityExactQueryConfiguration::EXACT
#endif
#ifdef ENABLE_GMP
            case NEXT_ITEM:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::GMP_FLOAT;
                break;
            case VisibilityExactQueryConfiguration::GMP_FLOAT:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::GMP_RATIONAL;
                break;
#undef NEXT_ITEM
#define NEXT_ITEM VisibilityExactQueryConfiguration::GMP_RATIONAL
#endif
#ifdef ENABLE_REALEXPR
            case NEXT_ITEM:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::REAL_EXPR;
                break;
#undef NEXT_ITEM
#define NEXT_ITEM VisibilityExactQueryConfiguration::REAL_EXPR
#endif
#ifdef ENABLE_MPFR
            case NEXT_ITEM:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::MPFR;
                break;
#undef NEXT_ITEM
#define NEXT_ITEM VisibilityExactQueryConfiguration::MPFR
#endif
            default:
                mDemoConfiguration.precisionType = VisibilityExactQueryConfiguration::FLOAT;
                break;
#undef NEXT_ITEM
            }
        }

#if 1
        void showGUI()
        {
            if (ImGui::BeginMainMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    if (ImGui::MenuItem("Quit", NULL, false, true)) {
                        exit(0);
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMainMenuBar();
            }

            if (show_demo_window)
            {
                ImGui::ShowDemoWindow(&show_demo_window);
            }

            {
                ImGui::Begin("Visilib 1.0. Demo application");

                forceDisplay |= ImGui::Checkbox("silhouette optimisation", &mDemoConfiguration.silhouetteOptimisation);
                forceDisplay |= ImGui::Checkbox("normalization", &mDemoConfiguration.normalization);

                const char* resolve_mode_items[] = { "Recursive", "NonRecursive", "Compare" };
                int resolveMode = mDemoConfiguration.resolveMode;
                forceDisplay |= ImGui::Combo("Resolve", &resolveMode, resolve_mode_items, IM_ARRAYSIZE(resolve_mode_items));
                mDemoConfiguration.resolveMode = (VisibilityExactQueryConfiguration::ResolveMode)resolveMode;

#if EMBREE
                if (ImGui::Checkbox("embree ray tracing", mDemoConfiguration.embree)) {
                    initScene(mDemoConfiguration.sceneIndex);
                    forceDisplay = true;
                }
#endif
                forceDisplay |= ImGui::Checkbox("representative line sampling strategy", &mDemoConfiguration.representativeLineSampling);
                forceDisplay |= ImGui::Checkbox("detect aperture only", &mDemoConfiguration.detectApertureOnly);

                const char* precision_items[] = { 
                    "FLOAT"
                    , "DOUBLE"
#ifdef EXACT_ARITHMETIC  
                    , "EXACT"
#endif
#ifdef ENABLE_GMP
                    , "GMP_FLOAT"
                    , "GMP_RATIONAL"
#endif
#ifdef ENABLE_REALEXPR
                    , "REAL_EXPR"
#endif
#ifdef ENABLE_MPFR
                    , "MPFR"
#endif
                };
                int precisionType = mDemoConfiguration.precisionType;
                forceDisplay |= ImGui::Combo("Precision", &precisionType, precision_items, IM_ARRAYSIZE(precision_items));
                mDemoConfiguration.precisionType = (VisibilityExactQueryConfiguration::PrecisionType)precisionType;

                ImGui::SameLine();
                ImGui::Text(DemoConfiguration::toStr(mDemoConfiguration.precisionType).c_str()); // validate that the enum and the string list match up...

                ImGui::Text("select scene"); ImGui::SameLine();
                if (ImGui::Button("prev ")) {
                    if (mDemoConfiguration.sceneIndex > 0)
                    {
                        mDemoConfiguration.sceneIndex --;
                    }
                    else
                    {
                        mDemoConfiguration.sceneIndex = 9;
                    }
                    initScene(mDemoConfiguration.sceneIndex);
                    forceDisplay = true;
                }
                ImGui::SameLine();
                ImGui::Text("%i", mDemoConfiguration.sceneIndex);
                ImGui::SameLine();
                if (ImGui::Button("next ")) {
                    if (mDemoConfiguration.sceneIndex < 9) 
                    {
                        mDemoConfiguration.sceneIndex ++;
                    }
                    else
                    {
                        mDemoConfiguration.sceneIndex = 0;
                    }
                    initScene(mDemoConfiguration.sceneIndex);
                    forceDisplay = true;
                }

                int tmp = mDemoConfiguration.vertexCount1;
                forceDisplay |= ImGui::SliderInt("number of vertices of query polygon", &tmp, 1, 12);
                mDemoConfiguration.vertexCount1 = tmp;

                ImGui::Text("adjust global scaling %i%%", (int)round(mDemoConfiguration.globalScaling * 100.0f)); ImGui::SameLine();
                if (ImGui::Button("+")) {
                    mDemoConfiguration.globalScaling *= 2;
                    forceDisplay = true;
                    setViewPortScaling(mDemoConfiguration.globalScaling);
                    initScene(mDemoConfiguration.sceneIndex);
                }
                ImGui::SameLine();
                if (ImGui::Button("-")) {
                    mDemoConfiguration.globalScaling /= 2;
                    forceDisplay = true;
                    setViewPortScaling(mDemoConfiguration.globalScaling);
                    initScene(mDemoConfiguration.sceneIndex);
                }

                if (ImGui::Button("change geometry mode")) {
                    drawGeometryType++;
                    drawGeometryType = drawGeometryType % 4;
                    forceDisplay = true;
                }
                ImGui::SameLine();
                ImGui::Text("%i", drawGeometryType);

                if (ImGui::Button("step")) {
                    step = true;
                    animated = false;
                    forceDisplay = true;
                }
                ImGui::SameLine();
                forceDisplay |= ImGui::Checkbox("animation", &animated);

                if (ImGui::Button("write config")) {
                    std::stringstream ss;

                    std::cout << "Save config.txt" << std::endl;
                    ss << "config_" << wc++ << ".txt";
                    writeConfig(ss.str());
                }
                ImGui::SameLine();
                if (ImGui::Button("open config")) {
                    std::stringstream ss;

                    std::cout << "Read config.txt" << std::endl;
                    ss << "config_" << rc++ << ".txt";
                    if (rc > wc)
                        rc = 0;
                    readConfig(ss.str());
                    initScene(mDemoConfiguration.sceneIndex);
                }

                ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
                ImGui::End();
            }

            {
                ImGui::Begin("Output");
                ImGui::Text(output.c_str());
                ImGui::End();
            }

            { 
                ImGui::Begin("Stats");
                ImGui::Text("resolveInternal %i", stats.resolveInternal);
                ImGui::Text("findNextEdge %i", stats.findNextEdge);
                ImGui::Text("maxRecursionLevel %i", stats.maxRecursionLevel);
                ImGui::Text("indexCount %i", stats.indexCount);
                ImGui::Text("polytopes %i", stats.polytopes);
                ImGui::End();
            }
        }
#endif
        void display()
        {
            stats.indexCount = 0;
            DemoDebugVisualisationGl::display(debugger, *meshContainer, v0, v1, result, drawGeometryType);
        }

        void displaySettings()
        {
            mDemoConfiguration.displaySettings();
        }

        void writeHelp()
        {
            std::cout << "Visilib 1.0. Demo application " << std::endl;

            std::cout << "  s: enable/disable silhouette optimisation" << std::endl;
            std::cout << "  n: enable/disable nomalization" << std::endl;
            std::cout << "  e: change precision type" << std::endl;
#if EMBREE
            std::cout << "  g: enable/disable embree ray tracing" << std::endl;
#endif
            std::cout << "  r: enable/disable representative line sampling strategy" << std::endl;
            std::cout << "  f: enable/disable detect aperture only" << std::endl;

            // std::cout << "  n: enable/disable Plucker normalization" << std::endl;
            // std::cout << "  f: enable/disable fast silhouette rejection test" << std::endl;

            std::cout << "  x: change scene " << std::endl;
            std::cout << "  +/-: increase/decrease scaling of query polygons" << std::endl;
            std::cout << "  1/2: increase/decrease number of vertices of query polygons" << std::endl;

            std::cout << "  [*] / [/] adjust global scaling" << std::endl;

            std::cout << "  w: write config" << std::endl;
            std::cout << "  o: open config" << std::endl;

            std::cout << "  space: start/pause animation" << std::endl;
            std::cout << "  Enter: show/hide geometry" << std::endl;

            std::cout << "  h: write this help" << std::endl;
        }

        void keyboard(unsigned char key, int, int)
        {
            std::stringstream ss;

            switch (key)
            {
            case 27:  // The escape key
            case 'Q':
            case 'q':
                exit(0);   // Simply exit
                break;

            case '2':
                if (mDemoConfiguration.vertexCount1 < 12)
                    mDemoConfiguration.vertexCount1++;

                forceDisplay = true;
                break;
            case '1':
                if (mDemoConfiguration.vertexCount1 > 1)
                    mDemoConfiguration.vertexCount1--;

                forceDisplay = true;
                break;

            case '+':
                if (mDemoConfiguration.scaling < 1.00f)
                    mDemoConfiguration.scaling += 0.01f;
                forceDisplay = true;
                break;

            case '-':
                if (mDemoConfiguration.scaling > 0.02f)
                    mDemoConfiguration.scaling -= 0.01f;

                forceDisplay = true;
                break;

            case '*':
                mDemoConfiguration.globalScaling *= 2;
                forceDisplay = true;
                setViewPortScaling(mDemoConfiguration.globalScaling);
                initScene(mDemoConfiguration.sceneIndex);

                break;

            case '/':
                mDemoConfiguration.globalScaling /= 2;

                forceDisplay = true;
                setViewPortScaling(mDemoConfiguration.globalScaling);
                initScene(mDemoConfiguration.sceneIndex);

                break;

            case 'h':

                writeHelp();

                displaySettings();
                break;
            case 's':
                mDemoConfiguration.silhouetteOptimisation = !mDemoConfiguration.silhouetteOptimisation;
                forceDisplay = true;
                break;
            case 'f':
                mDemoConfiguration.detectApertureOnly = !mDemoConfiguration.detectApertureOnly;
                forceDisplay = true;

                break;

            case 'x':
                mDemoConfiguration.sceneIndex++;
                if (mDemoConfiguration.sceneIndex > 9)
                    mDemoConfiguration.sceneIndex = 0;
                initScene(mDemoConfiguration.sceneIndex);
                forceDisplay = true;
                break;
            case 'r':
                mDemoConfiguration.representativeLineSampling = !mDemoConfiguration.representativeLineSampling;
                forceDisplay = true;
                break;


            case 'e':
                nextPrecisionType();
                std::cout << "  [Arithmetic: " << DemoConfiguration::toStr(mDemoConfiguration.precisionType) << "]" << std::endl;
                forceDisplay = true;

                break;
#if EMBREE
            case 'g':
                mDemoConfiguration.embree = !mDemoConfiguration.embree;

                initScene(mDemoConfiguration.sceneIndex);
                forceDisplay = true;

                break;
#endif

            case 'w':
                std::cout << "Save config.txt" << std::endl;
                ss << "config_" << wc++ << ".txt";
                writeConfig(ss.str());
                break;

            case 'o':
                std::cout << "Read config.txt" << std::endl;
                ss << "config_" << rc++ << ".txt";
                if (rc > wc)
                    rc = 0;
                readConfig(ss.str());
                initScene(mDemoConfiguration.sceneIndex);
                forceDisplay = true;
                break;

            case 'n':
                mDemoConfiguration.normalization = !mDemoConfiguration.normalization;

                forceDisplay = true;
                break;

            case 32:
                animated = !animated;
                break;

            case 13:
                drawGeometryType++;
                drawGeometryType = drawGeometryType % 4;
                break;
            }
        }
    private:
        std::vector<float> v0;
        std::vector<float> v1;

        GeometryOccluderSet* occluderSet = nullptr;
        HelperTriangleMeshContainer* meshContainer = nullptr;
        HelperVisualDebugger* debugger = nullptr;
        VisibilityResult result = UNKNOWN;

        DemoConfiguration mDemoConfiguration;
        bool forceDisplay = true;
        bool animated = false;
        bool step = false;
        int drawGeometryType = 0;
    };
}

static VisilibDemoMain* demo = nullptr;

void display()
{
#if 1
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGLUT_NewFrame();
#endif
    demo->display();
#if 1
    demo->showGUI();

    ImGui::Render();
    ImGuiIO& io = ImGui::GetIO();
    glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
#endif
#ifdef USE_GLUT
    glutSwapBuffers();
#endif
}

void keyboard(unsigned char key, int x, int y)
{
#if 1
    ImGui_ImplGLUT_KeyboardFunc(key, x, y);

    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureKeyboard)
#endif
    {
        demo->keyboard(key, x, y);
    }
}

void animate()
{
    demo->animate();
}

#if 1
const DWORD MS_VC_EXCEPTION = 0x406D1388;
#pragma pack(push,8)
typedef struct tagTHREADNAME_INFO
{
    DWORD dwType; // Must be 0x1000.
    LPCSTR szName; // Pointer to name (in user addr space).
    DWORD dwThreadID; // Thread ID (-1=caller thread).
    DWORD dwFlags; // Reserved for future use, must be zero.
} THREADNAME_INFO;
#pragma pack(pop)
void SetThreadName(DWORD dwThreadID, const char* threadName) {
    THREADNAME_INFO info;
    info.dwType = 0x1000;
    info.szName = threadName;
    info.dwThreadID = dwThreadID;
    info.dwFlags = 0;
#pragma warning(push)
#pragma warning(disable: 6320 6322)
    __try {
        RaiseException(MS_VC_EXCEPTION, 0, sizeof(info) / sizeof(ULONG_PTR), (ULONG_PTR*)&info);
    }
    __except (EXCEPTION_EXECUTE_HANDLER) {
    }
#pragma warning(pop)
}
#endif

int main(int argc, char** argv)
{
#ifdef WIN32
    HRESULT r = SetThreadDescription( GetCurrentThread(), L"Main Thread" );
#endif

    demo = new VisilibDemoMain();

    //demo->writeHelp();

#ifdef USE_GLUT
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowSize(800, 800);

    glutInitWindowPosition(100, 100);

    glutCreateWindow("Visilib demo");
#if 1
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer bindings
    ImGui_ImplGLUT_Init();
    ImGui_ImplGLUT_InstallFuncs();
    ImGui_ImplOpenGL2_Init();
#endif

    glutDisplayFunc(display);

    glutKeyboardFunc(keyboard);

    glutIdleFunc(animate);

    glutReshapeFunc(zprReshape);

    glutMouseFunc(zprMouse);

    glutMotionFunc(zprMotion);
#endif
    if (!demo->init())
    {
        std::cout << "Error reading geometry files. Exit" << std::endl;
        return 1;
    }
#ifdef USE_GLUT
    glutMainLoop();
#else
    animate();
#endif
    delete demo;
#if 1
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();
#endif
    return 0;
}

#ifdef _WIN32
// WinMain declaration to run without console subsystem on Windows
static int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, char*, int nShowCmd)
{
    return main(__argc, __argv);
}
#endif
