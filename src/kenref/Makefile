USE=core tools bias

#generic makefile
include ../maketools/make.module

# KEnRef's scalar type. This module and libkenref_core must agree on it.
CPPFLAGS += -DKENREF_DOUBLE

# The KEnRefBias constructor is hosted in the KEnRef repository and compiled here (see the header comment
# in KEnRefBias_setup.cpp), so this module needs KEnRef's include paths. They are resolved ONCE, at
# configure time, into Makefile.conf as KENREF_PLUMED_CPPFLAGS -- which make.module above has already
# -included. Do NOT turn this back into $(shell pkg-config ...): `sudo make install` re-parses this file
# as root with PKG_CONFIG_PATH stripped, so a parse-time lookup returns empty and the build dies. Baking
# mirrors how the LINK flags already flow (configure -> LIBS -> Makefile.conf's DYNAMIC_LIBS).
#
# Only the goals that actually COMPILE this module may demand the flags. Everything else has to work on a
# tree configured WITHOUT kenref, because PLUMED's own tooling runs those goals across every module dir
# regardless of what is enabled -- `maketools/plumedcheck` shells out to `make show_used_modules` (with
# stderr discarded) to learn the module's USE= list, and `make clean`/`links` run on unconfigured trees.
# When an $(error) fires there, make prints nothing on stdout, plumedcheck silently records NO used
# modules for kenref, and then reports every `#include "dir/file.h"` in this module as
# `include_non_used_module` -- with cppcheck losing its include configuration on top of that. The
# diagnostic points at the includes, but the cause is here.
#
# Whitelist rather than blacklist, so a query/analysis goal added to make.module later is safe by default.
kenref_building_goals := all obj lib install
ifeq ($(strip $(MAKECMDGOALS)),)
  kenref_needs_flags := yes                 # no explicit goal => make.module's default goal, `all`
else ifeq ($(strip $(filter $(kenref_building_goals),$(MAKECMDGOALS))),)
  kenref_needs_flags := no
else
  kenref_needs_flags := yes
endif

ifneq ($(strip $(KENREF_PLUMED_CPPFLAGS)),)
  CPPFLAGS := $(KENREF_PLUMED_CPPFLAGS) $(CPPFLAGS)
else ifneq ($(strip $(KENREF_SRC_REPO)),)
  # No configure-resolved value: fall back to a KEnRef SOURCE checkout. Mirrors kenref_plumed.pc's layout.
  CPPFLAGS := -I$(KENREF_SRC_REPO)/include_core -I$(KENREF_SRC_REPO)/include_plumedinterface -I$(KENREF_SRC_REPO)/src $(CPPFLAGS)
else ifeq ($(kenref_needs_flags),yes)
  $(error KEnRef: KENREF_PLUMED_CPPFLAGS is empty -- re-run ./configure with the installed kenref prefix on \
    PKG_CONFIG_PATH (it bakes the flags into Makefile.conf), or set KENREF_SRC_REPO=/path/to/KEnRef)
endif

# The forwarded KEnRef source file's own quoted PLUMED includes ("tools/Communicator.h" etc.) resolve
# relative to that file's directory, not this one: `-I.` picks up the USE= header symlinks created here,
# `-I..` picks up PLUMED's src/ root.
CPPFLAGS += -I. -I..
