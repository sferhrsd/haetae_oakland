# shake only

add_executable(${LIB_NAME2}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME2}-main ${LIB_NAME2})

add_executable(${LIB_NAME3}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME3}-main ${LIB_NAME3})

add_executable(${LIB_NAME5}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME5}-main ${LIB_NAME5})

# shake + aes

add_executable(${LIB_NAME2_AES}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME2_AES}-main ${LIB_NAME2_AES})

add_executable(${LIB_NAME3_AES}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME3_AES}-main ${LIB_NAME3_AES})

add_executable(${LIB_NAME5_AES}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME5_AES}-main ${LIB_NAME5_AES})

# blake only:

add_executable(${LIB_NAME2_BLAKE}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME2_BLAKE}-main ${LIB_NAME2_BLAKE})

add_executable(${LIB_NAME3_BLAKE}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME3_BLAKE}-main ${LIB_NAME3_BLAKE})

add_executable(${LIB_NAME5_BLAKE}-main main.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME5_BLAKE}-main ${LIB_NAME5_BLAKE})



add_executable(${LIB_NAME2}-test test.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
target_link_libraries(${LIB_NAME2}-test ${LIB_NAME2})

#add_executable(${LIB_NAME2}-fft-test fft_test.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
#target_link_libraries(${LIB_NAME2}-fft-test ${LIB_NAME2} fftw3 m)

#add_executable(${LIB_NAME3}-fft-test fft_test.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
#target_link_libraries(${LIB_NAME3}-fft-test ${LIB_NAME3} fftw3 m)

#add_executable(${LIB_NAME5}-fft-test fft_test.c ${CMAKE_SOURCE_DIR}/src/randombytes.c)
#target_link_libraries(${LIB_NAME5}-fft-test ${LIB_NAME5} fftw3 m)