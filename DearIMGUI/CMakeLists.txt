# adds a static library for DearIMGUI

set( IMGUI_DIR "E:/External/imgui" )

include_directories(${IMGUI_DIR})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories( E:/External/freeglut-2.8.0/include/ )

add_library( DearIMGUI STATIC
	${IMGUI_DIR}/examples/imgui_impl_opengl2.cpp
	${IMGUI_DIR}/examples/imgui_impl_opengl2.h
	${IMGUI_DIR}/examples/imgui_impl_glut.cpp
	${IMGUI_DIR}/examples/imgui_impl_glut.h
	${IMGUI_DIR}/imgui.cpp
	${IMGUI_DIR}/imgui_demo.cpp
	${IMGUI_DIR}/imgui_draw.cpp
	${IMGUI_DIR}/imgui_node_graph_test.cpp
	${IMGUI_DIR}/imgui_widgets.cpp
	${IMGUI_DIR}/imconfig.h
	${IMGUI_DIR}/imgui.h
	${IMGUI_DIR}/imgui_internal.h
	${IMGUI_DIR}/imstb_rectpack.h
	${IMGUI_DIR}/imstb_textedit.h
	${IMGUI_DIR}/imstb_truetype.h
)

set_target_properties(DearIMGUI PROPERTIES LINKER_LANGUAGE CXX)
