# which tests to run
TESTS = test-qktrace


# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs
test_qktrace_SOURCES = test-qktrace.cc

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
