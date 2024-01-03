include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# interface function for setting common library properties
#
function(setLibProperties targetname outputname)
    set_target_properties(${targetname} PROPERTIES
        OUTPUT_NAME ${outputname})
endfunction(setLibProperties)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)

set(scipsdpsources
    scipsdp/SdpVarmapper.c
    scipsdp/SdpVarfixer.c
    scipsdp/cons_sdp.c
    scipsdp/cons_savedsdpsettings.c
    scipsdp/cons_savesdpsol.c
    scipsdp/relax_sdp.c
    scipsdp/disp_sdpiterations.c
    scipsdp/disp_sdpavgiterations.c
    scipsdp/disp_sdpfastsettings.c
    scipsdp/disp_sdppenalty.c
    scipsdp/disp_sdpunsolved.c
    scipsdp/prop_sdpredcost.c
    scipsdp/branch_sdpmostfrac.c
    scipsdp/branch_sdpmostinf.c
    scipsdp/branch_sdpobjective.c
    scipsdp/branch_sdpinfobjective.c
    scipsdp/heur_sdpfracround.c
    scipsdp/heur_sdpinnerlp.c
    scipsdp/heur_sdpfracdiving.c
    scipsdp/heur_sdprand.c
    scipsdp/prop_sdpsymmetry.c
    scipsdp/prop_sdpobbt.c
    scipsdp/prop_companalcent.c
    scipsdp/reader_cbf.c
    scipsdp/reader_sdpa.c
    scipsdp/sdpsymmetry.c
    scipsdp/table_relaxsdp.c
    scipsdp/table_slater.c
    scipsdp/scipsdpdefplugins.c
    sdpi/sdpi.c
    sdpi/sdpsolchecker.c
    sdpi/solveonevarsdp.c
    sdpi/lapack_interface.c
    sdpi/sdpiclock.c
    scipsdpgithash.c
)

set(objscipsdpsources
)

set(sdpisources
   ${sdpi}
   )

set(symsources
   ${sym}
   )

set(sdpiheaders
    sdpi/configf77.h
    sdpi/sdpisolver.h
    sdpi/lapack_interface.h
    sdpi/sdpiclock.h
    sdpi/type_sdpiclock.h
    sdpi/struct_sdpiclock.h
    sdpi/type_sdpi.h
    sdpi/sdpi.h
    sdpi/sdpsolchecker.h
    sdpi/solveonevarsdp.h
)

set(objscipsdpheaders
)

set(scipsdpheaders
    scipsdp/SdpVarmapper.h
    scipsdp/SdpVarfixer.h
    scipsdp/cons_sdp.h
    scipsdp/cons_savedsdpsettings.h
    scipsdp/cons_savesdpsol.h
    scipsdp/relax_sdp.h
    scipsdp/disp_sdpiterations.h
    scipsdp/disp_sdpavgiterations.h
    scipsdp/disp_sdpfastsettings.h
    scipsdp/disp_sdppenalty.h
    scipsdp/disp_sdpunsolved.h
    scipsdp/prop_sdpredcost.h
    scipsdp/branch_sdpmostfrac.h
    scipsdp/branch_sdpmostinf.h
    scipsdp/branch_sdpobjective.h
    scipsdp/branch_sdpinfobjective.h
    scipsdp/heur_sdpfracround.h
    scipsdp/heur_sdpinnerlp.h
    scipsdp/heur_sdpfracdiving.h
    scipsdp/heur_sdprand.h
    scipsdp/prop_sdpobbt.h
    scipsdp/prop_companalcent.h
    scipsdp/prop_sdpsymmetry.h
    scipsdp/reader_cbf.h
    scipsdp/reader_sdpa.h
    scipsdp/sdpsymmetry.h
    scipsdp/table_relaxsdp.h
    scipsdp/table_slater.h
    scipsdp/scipsdpdefplugins.h
)

# all source files should be compiled with a c++ compiler
if(CXXONLY)
    set_source_files_properties(scipsdp/main.c ${scipsdpsources} ${objscipsdpsources} ${sdpisources} ${sym} PROPERTIES LANGUAGE CXX)

    # for the clang compiler this suppresses the warnings about treating 'c' input as 'c++' when CXXONLY is enabled
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
       add_compile_options(-x c++)
    endif()
endif()

add_library(libscipsdp ${scipsdpsources} ${objscipsdpsources} ${sdpisources} ${sym})
setLibProperties(libscipsdp "scipsdp")

#put binary in bin directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
if(SHARED)
    target_link_libraries(libscipsdp PRIVATE ${ZIMPL_PIC_LIBRARIES} ${SDPS_PIC_LIBRARIES} ${SCIP_LIBRARIES} ${SYM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} m)
    add_executable(scipsdp scipsdp/main.c ${scipsdpsources} ${objscipsdpsources} ${sdpisources} ${sym})
    target_link_libraries(scipsdp ${ZLIB_LIBRARIES} ${Readline_LIBRARY} ${GMP_LIBRARIES} ${ZIMPL_LIBRARIES} ${SDPS_LIBRARIES} ${SCIP_LIBRARIES} ${SYM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} m)
else()
    target_link_libraries(libscipsdp PRIVATE ${ZIMPL_LIBRARIES} ${SDPS_LIBRARIES} ${SYM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} m)
    add_executable(scipsdp scipsdp/main.c ${scipsdpsources} ${objscipsdpsources} ${sdpisources} ${sym})
    target_link_libraries(scipsdp libscipsdp ${SCIP_LIBRARIES})
endif()

target_compile_definitions(scipsdp PRIVATE EXTERN=extern)

#add_dependencies(libscipsdp scipsdp_update_githash)
#add_dependencies(scipsdp scipsdp_update_githash)

set_target_properties(libscipsdp PROPERTIES
    VERSION ${SCIPSDP_VERSION_MAJOR}.${SCIPSDP_VERSION_MINOR}.${SCIPSDP_VERSION_PATCH}
    SOVERSION ${SCIPSDP_VERSION_MAJOR}.${SCIPSDP_VERSION_MINOR})

# install the header files of scip
install(FILES ${sdpiheaders} DESTINATION include/sdpi)
install(FILES ${objscipsdpheaders} DESTINATION include/objscip)
install(FILES ${scipsdpheaders} DESTINATION include/scip)

# install the binary and the library to appropriate lcoations and add them to an export group
install(TARGETS scipsdp libscipsdp EXPORT scipsdp-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

# Add all targets to the build-tree export set
export(TARGETS scipsdp libscipsdp
  FILE "${CMAKE_BINARY_DIR}/scipsdp-targets.cmake")
