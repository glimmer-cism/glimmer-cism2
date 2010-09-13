/** This file uses libdl to redirect the glide interface subroutines to the
 * desired implementations */
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<dlfcn.h>
#include<config.inc>

#define NAMELEN 512
void *glide_if_lib = NULL;
char glide_if_lib_name[NAMELEN] = "libglide.so";
int glide_if_lib_is_open = 0;

#define glide_if_set_lib FC_FUNC (glide_if_set_lib, GLIDE_IF_SET_LIB)
/* this function is called from glide_if_config in the glide_if_setup module
 */
void
glide_if_set_lib (const char * libname, int namelength)
{
  int i;
  if (namelength >= NAMELEN) {
    abort();
  }

  glide_if_lib_name[namelength] = '\0';
  strncpy (glide_if_lib_name, libname, namelength);
  for (i = 0; i < NAMELEN; i++) {
    if (glide_if_lib_name[i] == ' ') {
      glide_if_lib_name[i] = '\0';
      break;
    }
  }
}

/* only open the library once */
void *
glide_if_dlopen()
{
  const char *err;
  if (glide_if_lib_is_open) {
    if (glide_if_lib == NULL) {
      abort();
    }
  }
  else {
    dlerror();
    glide_if_lib=dlopen(glide_if_lib_name, RTLD_LAZY);
    if (glide_if_lib == NULL) {
      printf ("%s\n", dlerror());
      abort();
    }
    glide_if_lib_is_open = 1;
  }
  return glide_if_lib;
}

void
glide_if_dlclose ()
{
  if (glide_if_lib_is_open) {
    if (glide_if_lib == NULL) {
      abort();
    }
    dlclose(glide_if_lib);
    glide_if_lib = NULL;
    glide_if_lib_is_open = 0;
  }
  else {
    if (glide_if_lib != NULL) {
      abort();
    }
  }
}

/* These are the implementations of the subroutines whose interfaces are given
 * in glide_if.inc */
#define glide_if_config FC_FUNC_ (glide_config, GLIDE_CONFIG)
#define glide_if_initialise FC_FUNC_ (glide_initialise, GLIDE_INITIALISE)
#define glide_if_tstep_p1 FC_FUNC_ (glide_tstep_p1, GLIDE_TSTEP_P1)
#define glide_if_tstep_p2 FC_FUNC_ (glide_tstep_p2, GLIDE_TSTEP_P2)
#define glide_if_tstep_p3 FC_FUNC_ (glide_tstep_p3, GLIDE_TSTEP_P3)
#define glide_if_finalise FC_FUNC_ (glide_finalise, GLIDE_FINALISE)
#define glide_if_nc_fillall FC_FUNC_ (glide_nc_fillall, GLIDE_NC_FILLALL)
#define glide_if_lithot FC_FUNC_ (setup_lithot, SETUP_LITHOT)
#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

typedef void (*onearg_t) (void *);
typedef void (*twoarg_t) (void *, void *);

void
glide_if_config (void * model, void * config)
{
  void *lib;
  twoarg_t config_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  config_fn = (twoarg_t) dlsym(lib, QUOTE(glide_if_config));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  config_fn (model, config);
}


void
glide_if_initialise (void * model)
{
  void *lib;
  onearg_t initialise_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  initialise_fn = (onearg_t) dlsym(lib, QUOTE(glide_if_initialise));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  initialise_fn (model);
}

void
glide_if_tstep_p1 (void * model, void * time)
{
  void *lib;
  twoarg_t tstep_p1_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  tstep_p1_fn = (twoarg_t) dlsym(lib, QUOTE(glide_if_tstep_p1));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  tstep_p1_fn (model, time);
}

void
glide_if_tstep_p2 (void * model)
{
  void *lib;
  onearg_t tstep_p2_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  tstep_p2_fn = (onearg_t) dlsym(lib, QUOTE(glide_if_tstep_p2));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  tstep_p2_fn (model);
}

void
glide_if_tstep_p3 (void * model)
{
  void *lib;
  onearg_t tstep_p3_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  tstep_p3_fn = (onearg_t) dlsym(lib, QUOTE(glide_if_tstep_p3));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  tstep_p3_fn (model);
}

void
glide_if_finalise (void * model, void * crash)
{
  void *lib;
  twoarg_t finalise_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  finalise_fn = (twoarg_t) dlsym(lib, QUOTE(glide_if_finalise));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  finalise_fn (model, crash);
  glide_if_dlclose();
}

void
glide_if_nc_fillall (void * model)
{
  void *lib;
  onearg_t nc_fillall_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  nc_fillall_fn = (onearg_t) dlsym(lib, QUOTE(glide_if_nc_fillall));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  nc_fillall_fn (model);
}

void
glide_if_lithot (void * model)
{
  void *lib;
  onearg_t lithot_fn;
  const char *err;

  lib=glide_if_dlopen();
  dlerror();
  lithot_fn = (onearg_t) dlsym(lib, QUOTE(glide_if_lithot));
  err=dlerror();
  if (err) {
    printf("%s\n", err);
    abort();
  }
  lithot_fn (model);
}
