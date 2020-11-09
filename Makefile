CXX= clang++
LINKER= clang++
LDFLAGS= `root-config --libs`
CPPFLAGS= `root-config --cflags`
CXXFLAGS= -Wno-deprecated -g
LNKFLAGS= -g
XBDIR=info
#DIR_INC=-I.
DIR_INC=-I$(XBDIR)
EXEC=qfs

all: $(EXEC)

OBJ=  qfs.o

SRC_XBTOOLS=  headers.hh	\
	info.hh	\
	tree.hh	\
	libs.hh

DEPS = $(patsubst %,$(XBDIR)/%,$(SRC_XBTOOLS))

# $(warning $(DEPS))

%.o : %.cc $(DEPS)
	@$(MAKEDEPEND)
	@${CXX} ${CPPFLAGS} ${CXXFLAGS} $(DIR_INC) -c $< -o $@
	@echo "	CXX $@"

%.o : %.C $(DEPS)
	@$(MAKEDEPEND)
	@${CXX} ${CPPFLAGS} ${CXXFLAGS} $(DIR_INC) -c $< -o $@
	@echo "	CXX $@"

$(EXEC): $(OBJ)
	@${LINKER} ${LDFLAGS} ${LNKFLAGS} -o $(EXEC) $(DIR_INC) $(OBJ)
	@echo "	COMP $(EXEC)"

clean:
	rm -f $(OBJ)
	rm -f $(EXEC)
