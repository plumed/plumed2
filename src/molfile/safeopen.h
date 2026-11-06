/*
 * Local addition to the VMD plugins (not part of the VMD distribution):
 * open a file for writing without making it writable by everybody.
 * fopen(path,"w") asks for the mode 0666 and relies on the umask of the
 * process; here the mode is 0644, whatever the umask is.
 */
#ifndef __PLUMED_molfile_safeopen_h
#define __PLUMED_molfile_safeopen_h

#include <stdio.h>
#if !defined(_WIN32)
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

/* mode is the fopen mode of a write open ("w", "wt", "wb") */
static inline FILE *molfile_fopen_write(const char *path, const char *mode) {
#if defined(_WIN32)
  return fopen(path, mode);
#else
  const int fdesc = open(path, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fdesc < 0) return NULL;
  FILE *f = fdopen(fdesc, mode);
  if (!f) close(fdesc);
  return f;
#endif
}

#endif
