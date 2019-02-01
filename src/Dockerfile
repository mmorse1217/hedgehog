FROM mobo

#TEST try compiling ebi #yolo
VOLUME ${MOBO_DIR}/src/ebi
#ADD . ${MOBO_DIR}/src/ebi
ENV DOCKER 1
WORKDIR ${MOBO_DIR}/src/ebi/

ENV PVFMM_INC ${PVFMM_INC}/include
ENV PVFMM_LIB /usr/local/share/pvfmm ###????
CMD ["make", "BOPT=g"]


