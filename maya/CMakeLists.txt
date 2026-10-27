if(NOT WIN32 OR NOT MSVC)
  message(FATAL_ERROR "The Maya adapter currently targets Windows x64 / MSVC")
endif()
if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  message(FATAL_ERROR "Use a 64-bit compiler for Maya 2023")
endif()
if(MSVC_VERSION LESS 1920 OR MSVC_VERSION GREATER_EQUAL 1930)
  message(FATAL_ERROR "Maya 2023 requires Visual Studio 2019 / v142. Configure with -T v142, or use the VS2019 developer command prompt.")
endif()
set(MAYA_ROOT "" CACHE PATH "Maya 2023 SDK root containing include/ and lib/")
if(NOT MAYA_ROOT)
  message(FATAL_ERROR "Set MAYA_ROOT to the Maya 2023 devkitBase or installation directory")
endif()
find_path(CURVENET_MAYA_INCLUDE maya/MPxDeformerNode.h PATHS "${MAYA_ROOT}/include" NO_DEFAULT_PATH REQUIRED)
file(STRINGS "${CURVENET_MAYA_INCLUDE}/maya/MTypes.h" maya_api_line REGEX "^#define MAYA_API_VERSION ")
if(NOT maya_api_line MATCHES "20230000")
  message(FATAL_ERROR "Use the Maya 2023 SDK (API 20230000)")
endif()
foreach(library Foundation OpenMaya OpenMayaAnim)
  find_library(CURVENET_MAYA_${library} ${library} PATHS "${MAYA_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
endforeach()

add_library(curvenet MODULE src/maya_deformer.cpp)
target_include_directories(curvenet PRIVATE "${CURVENET_MAYA_INCLUDE}")
target_link_libraries(curvenet PRIVATE CurveNet::core ${CURVENET_MAYA_Foundation} ${CURVENET_MAYA_OpenMaya} ${CURVENET_MAYA_OpenMayaAnim})
target_compile_definitions(curvenet PRIVATE NT_PLUGIN REQUIRE_IOSTREAM _BOOL)
target_link_options(curvenet PRIVATE /export:initializePlugin /export:uninitializePlugin)
set_target_properties(curvenet PROPERTIES PREFIX "" SUFFIX ".mll")
foreach(configuration "" _DEBUG _RELEASE _RELWITHDEBINFO _MINSIZEREL)
  set_target_properties(curvenet PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY${configuration} "${CMAKE_CURRENT_BINARY_DIR}/plug-ins"
    LIBRARY_OUTPUT_DIRECTORY${configuration} "${CMAKE_CURRENT_BINARY_DIR}/plug-ins")
endforeach()

# Produce an unpacked module next to the binary. No SDK or asset is redistributed.
add_custom_command(TARGET curvenet POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/scripts" "${CMAKE_CURRENT_BINARY_DIR}/tests"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/scripts/curvenet_maya.py" "${CMAKE_CURRENT_BINARY_DIR}/scripts/curvenet_maya.py"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/tests/maya2023_smoke.py" "${CMAKE_CURRENT_BINARY_DIR}/tests/maya2023_smoke.py"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/docs/maya2023.md" "${CMAKE_CURRENT_BINARY_DIR}/README.md"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_BINARY_DIR}/LICENSE"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/docs/EIGEN-LICENSE.txt" "${CMAKE_CURRENT_BINARY_DIR}/EIGEN-LICENSE.txt")
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/CurveNet.mod"
  CONTENT "+ MAYAVERSION:2023 PLATFORM:win64 CurveNet 0.2.1 ${CMAKE_CURRENT_BINARY_DIR}\nPYTHONPATH +:= scripts\n")
