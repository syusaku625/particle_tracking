COMPILER=icpc

#local environment
DEST = ~/bin
LIBDIR = /mnt/c/share/lib
TARDIR = ./bin
OBJDIR = ./obj
SRCDIR = ./

CPPFLAGS = -Wall -O3 -g -std=c++0x -lstdc++fs -DHDF5
INCLUDE = -I./include

#GLOG
#GLOG      = $(LIBDIR)/glog
#LDFLAGS += -L$(GLOG)/lib -lglog -lpthread
#INCLUDE += -I$(GLOG)/include

#HDF5
HDF5 = $(LIBDIR)/HDF5
LDFLAGS += -L$(HDF5)/lib -hdf5_hl_cpp -lhdf5_cpp -hdf5_hl -lhdf5 -lz
INCLUDE += -I$(HDF5)/include

#boost
#BOOST     = $(LIBDIR)/boost
#LDFLAGS += -L$(BOOST)/lib -lboost_graph 
#INCLUDE += -I$(BOOST)/include

EIGEN      = $(LIBDIR)/eigen
INCLUDE += -I$(EIGEN)/include/eigen3

TARGET = $(TARDIR)/run

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))
DEPENDS = $(OBJECTS:.o=d)

$(TARGET): $(OBJECTS) $(LIBS)
		@[ -d $(TARDIR) ] || mkdir -p $(TARDIR)
		$(COMPILER) -o $@ $^ $(LDFLAGS) 

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
		@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
		$(COMPILER) $(CPPFLAGS) $(INCLUDE) -o $@ -c $<

all: clean $(TARGET)

clean:
	rm -f $(OBJECTS) $(DEPENDS) $(TARGET)

install:$(TARGET)
	install -s $(TARGET) $(DEST)

-include $(DEPENDS)