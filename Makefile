
FC = gfortran -fcoarray=lib -O3
BUILD = build
moddir = $(BUILD)/mod
objdir = $(BUILD)/obj
app = $(BUILD)/ising_parallel
srcdir = src
srcfiles = mod_parallel.f90 parameters.f90
mainfile = main.f90
objfiles = $(patsusbt %,$(objdir)/%.o,$(srcfiles))

$(app): build/obj/mod_parallel.o build/obj/parameters.o build/obj/main.o
	gfortran $^ -o $@ -lcaf_mpi 

$(objdir)/mod_parallel.o: $(srcdir)/mod_parallel.f90 $(objdir) $(moddir)
	$(FC) -c -J $(moddir) $< -o $@


$(objdir)/parameters.o: $(srcdir)/parameters.f90 $(objdir) $(moddir)
	$(FC) -c -J $(moddir) $< -o $@

build/obj/main.o: src/main.f90 $(objdir)
	$(FC) -c -I $(moddir) $< -o $@

$(BUILD):
	mkdir -p $@


$(objdir): $(BUILD)
	mkdir -p $@

$(moddir): $(BUILD)
	mkdir -p $@

clean :
	rm -r $(BUILD)/*

run:
	echo input/parameters.nml | LD_LIBRARY_PATH=/home/jose/OpenCoarrays/build/lib cafrun -n 4 $(app)
