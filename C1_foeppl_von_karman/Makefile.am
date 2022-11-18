# A few file definitions
#-----------------------

# Define the sources
sourcesFvK =  \
foeppl_von_karman_elements.cc \
foeppl_von_karman_curved_bell_elements.cc

# Define the headers
headersFvK =  \
foeppl_von_karman.h \
foeppl_von_karman_elements.h \
foeppl_von_karman_curved_bell_elements.h

# Define name of library
libnameFvK = C1_foeppl_von_karman

# Combine headers and sources
headersFvK_and_sourcesFvK = $(headersFvK) $(sourcesFvK)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libC1_foeppl_von_karman.la

# Sources that the library depends on: 
#-------------------------------------
libC1_foeppl_von_karman_la_SOURCES = $(headersFvK_and_sourcesFvK)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
C1_foeppl_von_karmaninclude_HEADERS = $(libnameFvK).h 
# C1_basisincludedir = $(includedir) 
C1_foeppl_von_karmanincludedir = $(includedir)

#Here's the subdirectory where the actual header files are placed
library_C1_foeppl_von_karmanincludedir=$(includedir)/C1_foeppl_von_karman

#These are the header files that are to be placed in subdirectory
library_C1_foeppl_von_karmaninclude_HEADERS=$(headersFvK) 


# Required libraries -- [assuming that we want to link against stuff in generic
#-------------------     add other oomph-lib libraries if you need them....]
# Include path for library headers -- need to refer to include files
# in their respective source directories as they will not have been 
# installed yet! 
AM_CPPFLAGS += -I$(top_builddir)/src/generic


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libnameFvK).h: dummy_$(libnameFvK).h

dummy_$(libnameFvK).h:  $(headersFvK)
	echo  $(libnameFvK) $(headersFvK) > all_$(libnameFvK).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libnameFvK).aux > $(libnameFvK).h
	rm all_$(libnameFvK).aux

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_C1_foeppl_von_karmanincludedir) && rm -f $(headersFvK) )
	(echo "$(headersFvK)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_C1_foeppl_von_karmanincludedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Tidy up
#--------
clean-local:
	rm -f $(libnameFvK).h


