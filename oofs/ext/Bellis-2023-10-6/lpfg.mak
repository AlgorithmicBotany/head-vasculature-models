CC=g++

SRC_DIR=./include
LIB_DIR=./lib
OBJ_DIR=.

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXX_FLAGS +=-O3 -std=c++17 -Wno-return-type-c-linkage -fPIC
    CPPFLAGS +=
    LIBS +=
    SHARED_LIB_FLAGS +=-shared -Xlinker -Bsymbolic
    INC +=-I$(SRC_DIR) -I$(LIB_DIR)
endif
ifeq ($(UNAME_S),Darwin)
    CXX_FLAGS +=-arch x86_64 -O3 -std=c++17 -stdlib=libc++ -Wno-return-type-c-linkage
    CPPFLAGS +=
    LIBS +=-framework OpenGL
    SHARED_LIB_FLAGS +=-dynamiclib
    INC +=-I$(SRC_DIR) -I$(LIB_DIR)
endif

SRC_FILES=vec2f.cpp\
          vec3f.cpp\
		  vec2i.cpp\
		  vec3i.cpp\
		  vec4f.cpp\
		  mat3f.cpp\
		  mat4f.cpp\
		  barycentric_coordinates.cpp\
		  vec_utility.cpp\
          math_utility.cpp\
		  phyllotaxis.cpp\
		  timestamp_primordia_file_io.cpp\
		  minimal_resistance_rodkaew.cpp\
		  cage_deformation_2D.cpp\
		  cage_deformation_sequence_file_io.cpp\
		  cage_weight_utility.cpp\
		  volume_of_revolution.cpp\
		  packed_contiguous_utility.cpp\
		  edge.cpp\
		  curve.cpp\
		  triangle.cpp\
		  triangle_mesh.cpp\
		  triangle_mesh_2D.cpp\
		  random_utility.cpp\
		  receptacle_utility.cpp\
		  solid_angle_clustering.cpp
		  
SOURCES=$(addprefix $(SRC_DIR)/, $(SRC_FILES))

OBJECTS=$(patsubst %.cpp,%.o,$(notdir $(SOURCES)))

%.o: $(SRC_DIR)/%.cpp
	$(CC) $(INC) $(CPPFLAGS) $(CXX_FLAGS) -c -o $@ $<

lsys.so: lsys.o $(OBJECTS)
	$(CC) $(INC) $(CPPFLAGS) $(CXX_FLAGS) $(SHARED_LIB_FLAGS) -o lsys.so lsys.o $(OBJECTS) $(LIBS)

lsys.o: lsys.ii
	$(CC) $(INC) $(CPPFLAGS) $(CXX_FLAGS) -c lsys.ii -Wno-return-type-c-linkage
	
clean:
	rm *.o *.ii *.so
