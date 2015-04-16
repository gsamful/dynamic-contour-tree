CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project (LINK_CUT_TREE)
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
set(SRC_DIR ${CMAKE_SOURCE_DIR})
set(DYNAMIC_TREE_SRC_DIR ${CMAKE_SOURCE_DIR}/dynamic-tree)

set(SPLAY_TREE
	${DYNAMIC_TREE_SRC_DIR}/splay_tree.h
	${DYNAMIC_TREE_SRC_DIR}/statistics.h
	${DYNAMIC_TREE_SRC_DIR}/navigator.h
)

set(LINK_CUT_TREE
	${SPLAY_TREE}
	${DYNAMIC_TREE_SRC_DIR}/link_cut_tree.h
)

set(EULER_TREE
	${SPLAY_TREE}
	${DYNAMIC_TREE_SRC_DIR}/euler_tree.h
)

set(CONTOUR_TREE
	${SRC_DIR}/terrain.h
	${SRC_DIR}/vertex.h
	${SRC_DIR}/contour_tree.h
)

set(DYNAMIC_CONTOUR_TREE
	${LINK_CUT_TREE}
	${EULER_TREE}
	${SRC_DIR}/dynamic_contour_tree.h
	${SRC_DIR}/dynamic_contour_tree_interfaces.h
)

add_executable(contour_tree_test contour_tree_test_unit.cpp ${CONTOUR_TREE})
add_executable(dynamic_contour_tree_test dynamic_contour_tree_test_unit.cpp ${DYNAMIC_CONTOUR_TREE})
add_executable(kinetic_test kinetic_test.cpp ${DYNAMIC_CONTOUR_TREE})