if did_filetype()	" filetype already set..
 finish		" ..don't do these checks
endif
if getline(1) =~ '^#! FIELDS .*$' || getline(1) =~ '^#! SET .*$'
  setfiletype plumedf
endif
