# ----------------------------------------
# inputs 

FASTJETLOCATION=/Users/ntran/Documents/Research/CMS/PhysicsAnalysis/fastjet/fastjet-install/
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags`

# ----------------------------------------

FILES = Qjets.o QjetsPlugin.o

all: $(FILES) lib/libQjets.a anaSubstructure

lib/libQjets.a: $(FILES) $(FILES:.cc=.o)
	ar cq lib/libQjets.a $(FILES)

lib/libQjets.so: $(FILES) $(FILES:.cc=.o)
	gcc $(SOFLAGS) -o $@ $(FILES)

anaSubstructure: anaSubstructure.cpp $(FILES)
	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict -lEnergyCorrelator -lQjets -lSoftDrop -L./lib anaSubstructure.cpp
    
%.o: %.C %.h
	gcc -fPIC -O3 -c $(FASTJETINC) $< -o $@ 

clean:
	rm $(FILES) anaSubstructure lib/*