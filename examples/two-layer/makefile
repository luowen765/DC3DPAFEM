
MFEM_DIR = ../../contrib/mfem_package/mfem-4.5
include $(MFEM_DIR)/config/config.mk

#gnu compiler
CXX		= $(MFEM_CXX)
#Ubuntu
CXXFLAGS	= $(MFEM_FLAGS)

INCLUDES	= $(wildcard ../../src/*.h)
SRCS		= $(wildcard ../../src/*.cpp)
OBJS		= $(patsubst %.cpp, %.o, $(SRCS))

DC3D_INCLUDE	= -I../../src/

%.o : %.cpp
	@echo "DC3D is compiling "$<"..."
	@$(CXX) $(CXXFLAGS) $(DC3D_INCLUDE) $(MFEM_INCLUDE)  -c $< -o $@ 

DC3D: $(OBJS)
	@$(CXX) $(CXXFLAGS) -o DC3D $(OBJS) $(MFEM_LIBS) 


clean:
	@rm -rf $(OBJS)
	
cleanall: 
	@rm -rf $(OBJS) DC3D
	
echo:	
	@echo $(MFEM_CXX)
