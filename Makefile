##=======================================================================
###############################################
# COPYRIGHT: See COPYRIGHT.txt                #
# 2015 by MCS, Argonne National Laboratory.   #
###############################################
##=======================================================================

##=======================================================================
##   PLEASE SET THESE VARIABLES BEFORE COMPILING
##=======================================================================

HDF5PATH	= /home/sdi/Install/hdf5-1.10.3-install
#HDF5PATH	= /home/sdi/Install/hdf5-1.12.1-install

##=======================================================================
##   DIRECTORY TREE
##=======================================================================

LIB		= lib
OBJ		= obj
SRC		= src
INC		= include

##=======================================================================
##   COMPILERS
##=======================================================================

CC 		= gcc
MPICC 		= mpicc

##=======================================================================
##   FLAGS
##=======================================================================

#ROCCIFLAGS         = -I$(ROCCIPATH)/include -L$(ROCCIPATH)/lib -Wl,-rpath,$(ROCCIPATH)/lib

HDF5FLAGS	= -I$(HDF5PATH)/include #$(HDF5PATH)/lib/libhdf5.a

##=======================================================================
##   TARGETS
##=======================================================================

OBJS		= $(OBJ)/H5Z_ROCCI.o

SHARED		= libhdf5ROCCI.so
STATIC		= libhdf5ROCCI.a

all: 		$(LIB)/$(SHARED) $(LIB)/$(STATIC)

$(OBJ)/%.o:	$(SRC)/%.c
		@mkdir -p $(OBJ)
		$(CC) -c -fPIC -g -O3 -I./include $(HDF5FLAGS) $(ROCCIFLAGS) $< -o $@
		
$(LIB)/$(SHARED):	$(OBJS)
		@mkdir -p $(LIB)
		$(CC) -O3 -g -shared -o $(LIB)/$(SHARED) $(OBJS) $(ROCCIFLAGS) -L$(HDF5PATH)/lib -lc -lhdf5

$(LIB)/$(STATIC):	$(OBJS)
		@mkdir -p $(LIB)
		$(RM) $@
		$(AR) -cvq $@ $(OBJS)
		
install: $(LIB)/$(SHARED) $(LIB)/$(STATIC)
		install -d $(ROCCIPATH)/lib  $(ROCCIPATH)/include
		install $(INC)/* $(ROCCIPATH)/include/
		install $(LIB)/* $(ROCCIPATH)/lib/

uninstall:
		$(RM) $(ROCCIPATH)/$(LIB)/libhdf5ROCCI.* $(ROCCIPATH)/$(INC)/H5Z_ROCCI.h

clean:
		rm -f $(OBJ)/*.o $(LIB)/*

.PHONY:		$(SHARED) $(STATIC) install uninstall clean


