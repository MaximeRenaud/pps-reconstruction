# Tests needing a ROOT instance
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /usr/share/root/cmake)
find_package(ROOT)

if (ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIRS})   # Sometimes [...]_DIR
#   message("Notice me senpai !!!! : ${ROOT_INCLUDE_DIRS}\n")
  link_directories(${ROOT_LIBRARY_DIR})
#   message("Me too senpai !!! : ${ROOT_LIBRARY_DIR}\n")
#   message("And us senpai !!! : 1) ${ROOT_LIBRARIES}\n              2) ${ROOT_COMPONENT_LIBRARIES}\n")  #faulty !!
  
  if(${ROOT_VERSION} LESS 6.0)
    set(GCC_COMPILE_FLAGS "-Wno-shadow -fPIC")
  else()
    set(GCC_COMPILE_FLAGS "-Wno-shadow -fPIC -std=c++11")
  endif()
  add_definitions(${GCC_COMPILE_FLAGS})
endif()

function(ADD_TEST exec)
  add_executable(${exec} ${PROJECT_SOURCE_DIR}/test/${exec}.cpp $<TARGET_OBJECTS:reader_lib>)
  target_link_libraries(${exec} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES}) # component_libraries may be empty ?
  set_property(TARGET ${exec} PROPERTY EXCLUDE_FROM_ALL true)
endfunction()

if (ROOT_FOUND)
  add_test(write_tree_sorted)
  add_test(NoiseIssue)
  add_test(write_tree_gastof)
  add_test(gastof_simple)
endif()

#add_test(testdb)
#set_property(TARGET testdb PROPERTY LINK_FLAGS "-lsqlite3")

