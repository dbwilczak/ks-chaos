CAPD = ../CAPD/build/bin/
INCLUDE = `$(CAPD)capd-config --cflags` -I.
LIBS = `$(CAPD)capd-config --libs` -lpthread
FLAGS = -O2 -s $(INCLUDE)
CXX = g++

# files to be compiled 
SRC = $(wildcard progs/*.cpp)
OBJS = $(patsubst progs/%.cpp,obj/%.o,$(SRC)) 
DEPS = $(patsubst progs/%.cpp,dep/%.d,$(SRC))
PROGS = $(patsubst progs/%.cpp,%,$(SRC))

all: $(PROGS)

include $(DEPS)

%: obj/%.o 
	@echo $(CXX) $< -o $@ 
	@$(CXX) $< -o $@ $(LIBS)

obj/%.o: progs/%.cpp dep/%.d
	@echo $(CXX) -c $< -o $@
	@$(CXX) $(FLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f obj/*.o dep/*.d $(PROGS)

tar:
	make clean
	rm -f supplement.tgz
	tar cvfz supplement.tgz *

dep/%.d: */%.cpp 
	@$(CXX) $(FLAGS) -MM -MT obj/$*.o $< > $@; $(include $@)
