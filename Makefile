
FC = caf
BUILD = build
moddir = $(BUILD)/mod
objdir = $(BUILD)/obj
app = $(BUILD)/ising_parallel


$(app): build/obj/mod_parallel.o build/obj/main.o
	$(FC) $^ -o $@

build/obj/mod_parallel.o: src/mod_parallel.f90 $(objdir) $(moddir)
	$(FC) -c -J $(moddir) $< -o $@


build/obj/main.o: src/main.f90 $(objdir)
	$(FC) -c -I $(moddir) $< -o $@


$(BUILD):
	mkdir -p $@


$(objdir): $(BUILD)
	mkdir -p $@

$(moddir): $(BUILD)
	mkdir -p $@


run:
	cafrun -n 1 $(app)
