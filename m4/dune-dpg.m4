dnl -*- autoconf -*-
# Macros needed to find dune-dpg and dependent libraries.  They are called by
# the macros in ${top_src_dir}/dependencies.m4, which is generated by
# "dunecontrol autogen"

# Additional checks needed to build dune-dpg
# This macro should be invoked by every module which depends on dune-dpg, as
# well as by dune-dpg itself
AC_DEFUN([DUNE_DPG_CHECKS])

# Additional checks needed to find dune-dpg
# This macro should be invoked by every module which depends on dune-dpg, but
# not by dune-dpg itself
AC_DEFUN([DUNE_DPG_CHECK_MODULE],
[
  DUNE_CHECK_MODULES([dune-dpg],[dpg/dpg.hh])
])
