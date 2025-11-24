.PHONY: all release binary clean

binary: build/minimera

all: build/minimera build/minimera.1 build/minimera.sif

# Config ----------------------------------------------------------------------
lisps := $(shell ffind '\.(asd|lisp)$$')

# Build -----------------------------------------------------------------------
build/asdf-manifest: Makefile minimera.asd
	mkdir -p build/
	sbcl --disable-debugger --quit --eval '(ql:write-asdf-manifest-file "build/asdf-manifest")'

build/minimera: $(lisps) Makefile build/asdf-manifest build-binary.sh
	mkdir -p build/
	./build-binary.sh

build/minimera.1: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "build-manual.lisp" --quit

build/minimera.sif: build/minimera build/minimera.1 contrib/minimera.def
	mkdir -p build/
	singularity build --force --fakeroot build/minimera.sif contrib/minimera.def


# Releases --------------------------------------------------------------------
release: build/minimera build/minimera.1 build/minimera.sif build-release.sh Makefile
	mkdir -p build/
	./build-release.sh

# Clean -----------------------------------------------------------------------
clean:
	rm -r ./build/
