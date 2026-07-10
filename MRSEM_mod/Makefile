#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

USE_HDF5=1
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5



HEADERS = OPS_MRSEM_kernels.h OPS_MRSEM_eddy_functions.h OPS_MRSEM_constants.h TBL_data.h

OTHER_FILES = 
OPS_GENERATED = OPS_MRSEM_ops.cpp 
OPS_FILES = OPS_MRSEM.cpp 


APP=OPS_MRSEM
MAIN_SRC=OPS_MRSEM

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app_legacy
