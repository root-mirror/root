#if !defined(GLEW_BUILD)
#ifndef __building_module
# define __building_module(X) 0
#endif
#  if !__building_module(OpenGL)
#  error "You shouldn't #include gl.h directly. Please use TGLIncludes.h instead."
#  endif
#endif
