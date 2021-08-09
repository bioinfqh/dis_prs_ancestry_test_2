
bin_PROGRAMS = ancestry relfinder relfinder_vcf
DISTCHECK_CONFIGURE_FLAGS=LDFLAGS=-L/opt/local/lib CXXFLAGS=-I/opt/local/include
ancestry_SOURCES = ancestry.cpp Data23.cpp Ancestry_params.cpp CmdLine.cpp gzstream.cpp FreqData.cpp
relfinder_SOURCES = relfinder.cpp Data23.cpp Ancestry_params.cpp CmdLine.cpp gzstream.cpp FreqData.cpp
relfinder_vcf_SOURCES = relfinder_vcf.cpp DataVCF.cpp VCF_params.cpp CmdLine.cpp gzstream.cpp

noinst_HEADERS = gzstream.h FreqData.h DataVCF.h VCF_params.h CmdLine.h Data23.h Ancestry_params.h
