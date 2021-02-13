#ifdef ROOT_TX11GL
// FIXME: TX11GL.h should not include glx.h. It should use glew.h. However if
// this becomes the case all sort of other implementation deficiencies occur.
// They are fixable but require a lot of efforts and intricate understanding of
// our GL implementation.
#include_next <GL/gl.h>
#else
#if !defined(GLEW_BUILD)
#ifndef __has_feature
# define __has_feature(X) 0
#endif
#  if !__has_feature(modules)
#  error "You shouldn't #include gl.h directly. Please use TGLIncludes.h instead."
#  endif
#endif
#endif
