/*
 * Copyright (c) 2001,2000,1999 Sasha Vasko <sasha at aftercode.net>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef _WIN32
#include "win32/config.h"
#else
#include "config.h"
#endif

#define LOCAL_DEBUG
#undef DEBUG_SL2XIMAGE
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
 
#ifdef _WIN32
# include "win32/afterbase.h"
#else
# include "afterbase.h"
#endif
#include "asvisual.h"
#include "scanline.h"

#if defined(XSHMIMAGE) && !defined(X_DISPLAY_MISSING) 
# include <sys/ipc.h>
# include <sys/shm.h>
# include <X11/extensions/XShm.h>
#else
# undef XSHMIMAGE
#endif

#if defined(HAVE_GLX) && !defined(X_DISPLAY_MISSING) 
# include <GL/gl.h>
# include <GL/glx.h>
#else
# undef HAVE_GLX
#endif



#ifndef X_DISPLAY_MISSING
static int  get_shifts (unsigned long mask);
static int  get_bits (unsigned long mask);

void _XInitImageFuncPtrs(XImage*);

int
asvisual_empty_XErrorHandler (Display * dpy, XErrorEvent * event)
{
    return 0;
}
/***************************************************************************/
#if defined(LOCAL_DEBUG) && !defined(NO_DEBUG_OUTPUT)
Status
debug_AllocColor( const char *file, int line, ASVisual *asv, Colormap cmap, XColor *pxcol )
{
	Status sret ;
	sret = XAllocColor( asv->dpy, cmap, pxcol );
	show_progress( " XAllocColor in %s:%d has %s -> cmap = %lX, pixel = %lu(%8.8lX), color = 0x%4.4X, 0x%4.4X, 0x%4.4X",
				   file, line, (sret==0)?"failed":"succeeded", (long)cmap, (unsigned long)(pxcol->pixel), (unsigned long)(pxcol->pixel), pxcol->red, pxcol->green, pxcol->blue );
	return sret;
}
#define ASV_ALLOC_COLOR(asv,cmap,pxcol)  debug_AllocColor(__FILE__, __LINE__, (asv),(cmap),(pxcol))
#else
#define ASV_ALLOC_COLOR(asv,cmap,pxcol)  XAllocColor((asv)->dpy,(cmap),(pxcol))
#endif

#else
#define ASV_ALLOC_COLOR(asv,cmap,pxcol)   0
#endif   /* ndef X_DISPLAY_MISSING */

/**********************************************************************/
/* returns the maximum number of true colors between a and b          */
long
ARGB32_manhattan_distance (long a, long b)
{
	register int d = (int)ARGB32_RED8(a)   - (int)ARGB32_RED8(b);
	register int t = (d < 0 )? -d : d ;

	d = (int)ARGB32_GREEN8(a) - (int)ARGB32_GREEN8(b);
	t += (d < 0)? -d : d ;
	d = (int)ARGB32_BLUE8(a)  - (int)ARGB32_BLUE8(b);
	return (t+((d < 0)? -d : d)) ;
}



/***************************************************************************
 * ASVisual :
 * encoding/decoding/querying/setup
 ***************************************************************************/
int get_bits_per_pixel(Display *dpy, int depth)
{
#if 0
#ifndef X_DISPLAY_MISSING
 	register ScreenFormat *fmt = dpy->pixmap_format;
 	register int i;

 	for (i = dpy->nformats + 1; --i; ++fmt)
 		if (fmt->depth == depth)
 			return(fmt->bits_per_pixel);
#endif
#endif
	if (depth <= 4)
	    return 4;
	if (depth <= 8)
	    return 8;
	if (depth <= 16)
	    return 16;
	return 32;
 }

/* ********************* ASVisual ************************************/
ASVisual *_set_default_asvisual( ASVisual *new_v );

#ifndef X_DISPLAY_MISSING
static XColor black_xcol = { 0, 0x0000, 0x0000, 0x0000, DoRed|DoGreen|DoBlue };
static XColor white_xcol = { 0, 0xFFFF, 0xFFFF, 0xFFFF, DoRed|DoGreen|DoBlue };

static void find_useable_visual( ASVisual *asv, Display *dpy, int screen,
	                             Window root, XVisualInfo *list, int nitems,
								 XSetWindowAttributes *attr )
{
	int k ;
	int (*oldXErrorHandler) (Display *, XErrorEvent *) =
						XSetErrorHandler (asvisual_empty_XErrorHandler);
	Colormap orig_cmap = attr->colormap ;

	for( k = 0  ; k < nitems ; k++ )
	{
		Window       w = None, wjunk;
		unsigned int width, height, ujunk ;
		int          junk;
		/* try and use default colormap when possible : */
		if( orig_cmap == None )
		{
  			if( list[k].visual == DefaultVisual( dpy, (screen) ) )
			{
				attr->colormap = DefaultColormap( dpy, screen );
				LOCAL_DEBUG_OUT( "Using Default colormap %lX", attr->colormap );
			}else
			{
				attr->colormap = XCreateColormap( dpy, root, list[k].visual, AllocNone);
				LOCAL_DEBUG_OUT( "DefaultVisual is %p, while ours is %p, so Created new colormap %lX", DefaultVisual( dpy, (screen) ), list[k].visual, attr->colormap );
			}
		}
		ASV_ALLOC_COLOR( asv, attr->colormap, &black_xcol );
		ASV_ALLOC_COLOR( asv, attr->colormap, &white_xcol );
		attr->border_pixel = black_xcol.pixel ;

/*fprintf( stderr, "checking out visual ID %d, class %d, depth = %d mask = %X,%X,%X\n", list[k].visualid, list[k].class, list[k].depth, list[k].red_mask, list[k].green_mask, list[k].blue_mask 	);*/
		w = XCreateWindow (dpy, root, -10, -10, 10, 10, 0, list[k].depth, CopyFromParent, list[k].visual, CWColormap|CWBorderPixel, attr );
		if( w != None && XGetGeometry (dpy, w, &wjunk, &junk, &junk, &width, &height, &ujunk, &ujunk))
		{
			/* don't really care what's in it since we do not use it anyways : */
			asv->visual_info = list[k] ;
			XDestroyWindow( dpy, w );
			asv->colormap = attr->colormap ;
			asv->own_colormap = (attr->colormap != DefaultColormap( dpy, screen ));
			asv->white_pixel = white_xcol.pixel ;
			asv->black_pixel = black_xcol.pixel ;
			break;
		}
		if( orig_cmap == None )
		{
			if( attr->colormap != DefaultColormap( dpy, screen ))
				XFreeColormap( dpy, attr->colormap );
			attr->colormap = None ;
		}
	}
	XSetErrorHandler(oldXErrorHandler);
}
#endif

/* Main procedure finding and querying the best visual */
Bool
query_screen_visual_id( ASVisual *asv, Display *dpy, int screen, Window root, int default_depth, VisualID visual_id, Colormap cmap )
{
#ifndef X_DISPLAY_MISSING
	int nitems = 0 ;
	/* first  - attempt locating 24bpp TrueColor or DirectColor RGB or BGR visuals as the best cases : */
	/* second - attempt locating 32bpp TrueColor or DirectColor RGB or BGR visuals as the next best cases : */
	/* third  - lesser but still capable 16bpp 565 RGB or BGR modes : */
	/* forth  - even more lesser 15bpp 555 RGB or BGR modes : */
	/* nothing nice has been found - use whatever X has to offer us as a default :( */
	int i ;

	XVisualInfo *list = NULL;
	XSetWindowAttributes attr ;
	static XVisualInfo templates[] =
		/* Visual, visualid, screen, depth, class      , red_mask, green_mask, blue_mask, colormap_size, bits_per_rgb */
		{{ NULL  , 0       , 0     , 24   , TrueColor  , 0xFF0000, 0x00FF00  , 0x0000FF , 0            , 0 },
		 { NULL  , 0       , 0     , 24   , TrueColor  , 0x0000FF, 0x00FF00  , 0xFF0000 , 0            , 0 },
		 { NULL  , 0       , 0     , 24   , TrueColor  , 0x0     , 0x0       , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 32   , TrueColor  , 0xFF0000, 0x00FF00  , 0x0000FF , 0            , 0 },
		 { NULL  , 0       , 0     , 32   , TrueColor  , 0x0000FF, 0x00FF00  , 0xFF0000 , 0            , 0 },
		 { NULL  , 0       , 0     , 32   , TrueColor  , 0x0     , 0x0       , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , TrueColor  , 0xF800  , 0x07E0    , 0x001F   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , TrueColor  , 0x001F  , 0x07E0    , 0xF800   , 0            , 0 },
		 /* big endian or MBR_First modes : */
		 { NULL  , 0       , 0     , 16   , TrueColor  , 0x0     , 0xE007    , 0x0      , 0            , 0 },
		 /* some misrepresented modes that really are 15bpp : */
		 { NULL  , 0       , 0     , 16   , TrueColor  , 0x7C00  , 0x03E0    , 0x001F   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , TrueColor  , 0x001F  , 0x03E0    , 0x7C00   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , TrueColor  , 0x0     , 0xE003    , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 15   , TrueColor  , 0x7C00  , 0x03E0    , 0x001F   , 0            , 0 },
		 { NULL  , 0       , 0     , 15   , TrueColor  , 0x001F  , 0x03E0    , 0x7C00   , 0            , 0 },
		 { NULL  , 0       , 0     , 15   , TrueColor  , 0x0     , 0xE003    , 0x0      , 0            , 0 },
/* no suitable TrueColorMode found - now do the same thing to DirectColor :*/
		 { NULL  , 0       , 0     , 24   , DirectColor, 0xFF0000, 0x00FF00  , 0x0000FF , 0            , 0 },
		 { NULL  , 0       , 0     , 24   , DirectColor, 0x0000FF, 0x00FF00  , 0xFF0000 , 0            , 0 },
		 { NULL  , 0       , 0     , 24   , DirectColor, 0x0     , 0x0       , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 32   , DirectColor, 0xFF0000, 0x00FF00  , 0x0000FF , 0            , 0 },
		 { NULL  , 0       , 0     , 32   , DirectColor, 0x0000FF, 0x00FF00  , 0xFF0000 , 0            , 0 },
		 { NULL  , 0       , 0     , 32   , DirectColor, 0x0     , 0x0       , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , DirectColor, 0xF800  , 0x07E0    , 0x001F   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , DirectColor, 0x001F  , 0x07E0    , 0xF800   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , DirectColor, 0x0     , 0xE007    , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , DirectColor, 0x7C00  , 0x03E0    , 0x001F   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , DirectColor, 0x001F  , 0x03E0    , 0x7C00   , 0            , 0 },
		 { NULL  , 0       , 0     , 16   , DirectColor, 0x0     , 0xE003    , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 15   , DirectColor, 0x7C00  , 0x03E0    , 0x001F   , 0            , 0 },
		 { NULL  , 0       , 0     , 15   , DirectColor, 0x001F  , 0x03E0    , 0x7C00   , 0            , 0 },
		 { NULL  , 0       , 0     , 15   , DirectColor, 0x0     , 0xE003    , 0x0      , 0            , 0 },
		 { NULL  , 0       , 0     , 0    , 0          , 0       , 0         , 0        , 0            , 0 },
		} ;
#endif /*ifndef X_DISPLAY_MISSING */
	if( asv == NULL )
		return False ;
	memset( asv, 0x00, sizeof(ASVisual));

	asv->dpy = dpy ;

#ifndef X_DISPLAY_MISSING
	memset( &attr, 0x00, sizeof( attr ));
	attr.colormap = cmap ;

	if( visual_id == 0 )
	{
		for( i = 0 ; templates[i].depth != 0 ; i++ )
		{
			int mask = VisualScreenMask|VisualDepthMask|VisualClassMask ;

			templates[i].screen = screen ;
			if( templates[i].red_mask != 0 )
				mask |= VisualRedMaskMask;
			if( templates[i].green_mask != 0 )
				mask |= VisualGreenMaskMask ;
			if( templates[i].blue_mask != 0 )
				mask |= VisualBlueMaskMask ;

			if( (list = XGetVisualInfo( dpy, mask, &(templates[i]), &nitems ))!= NULL )
			{
				find_useable_visual( asv, dpy, screen, root, list, nitems, &attr );
				XFree( list );
				list = NULL ;
				if( asv->visual_info.visual != NULL )
					break;
			}
		}
	}else
	{
		templates[0].visualid = visual_id ;
		if( (list = XGetVisualInfo( dpy, VisualIDMask, &(templates[0]), &nitems )) != NULL )
		{
			find_useable_visual( asv, dpy, screen, root, list, nitems, &attr );
		 	XFree( list );
			list = NULL ;
		}
		if( asv->visual_info.visual == NULL )
			show_error( "Visual with requested ID of 0x%X is unusable - will try default instead.", visual_id );
	}

	if( asv->visual_info.visual == NULL )
	{  /* we ain't found any decent Visuals - that's some crappy screen you have! */
		register int vclass = 6 ;
		while( --vclass >= 0 )
			if( XMatchVisualInfo( dpy, screen, default_depth, vclass, &(asv->visual_info) ) )
				break;
		if( vclass < 0 )
			return False;
		/* try and use default colormap when possible : */
		if( asv->visual_info.visual == DefaultVisual( dpy, screen ) )
			attr.colormap = DefaultColormap( dpy, screen );
		else
			attr.colormap = XCreateColormap( dpy, root, asv->visual_info.visual, AllocNone);
		ASV_ALLOC_COLOR( asv, attr.colormap, &black_xcol );
		ASV_ALLOC_COLOR( asv, attr.colormap, &white_xcol );
		asv->colormap = attr.colormap ;
		asv->own_colormap = (attr.colormap != DefaultColormap( dpy, screen ));
		asv->white_pixel = white_xcol.pixel ;
		asv->black_pixel = black_xcol.pixel ;
	}
	if( get_output_threshold() >= OUTPUT_VERBOSE_THRESHOLD )
	{
		fprintf( stderr, "Selected visual 0x%lx: depth %d, class %d\n RGB masks: 0x%lX, 0x%lX, 0x%lX, Byte Ordering: %s\n",
				 (unsigned long)asv->visual_info.visualid,
				 asv->visual_info.depth,
				 asv->visual_info.class,
				 (unsigned long)asv->visual_info.red_mask,
				 (unsigned long)asv->visual_info.green_mask,
				 (unsigned long)asv->visual_info.blue_mask,
				 (ImageByteOrder(asv->dpy)==MSBFirst)?"MSBFirst":"LSBFirst" );
	}
#else
	asv->white_pixel = ARGB32_White ;
	asv->black_pixel = ARGB32_Black ;
#endif /*ifndef X_DISPLAY_MISSING */
	return True;
}

ASVisual *
create_asvisual_for_id( Display *dpy, int screen, int default_depth, VisualID visual_id, Colormap cmap, ASVisual *reusable_memory )
{
	ASVisual *asv = reusable_memory ;
#ifndef X_DISPLAY_MISSING
    Window root = dpy?RootWindow(dpy,screen):None;
#endif /*ifndef X_DISPLAY_MISSING */

	if( asv == NULL )
        asv = safecalloc( 1, sizeof(ASVisual) );
#ifndef X_DISPLAY_MISSING
    if( dpy )
    {
        if( query_screen_visual_id( asv, dpy, screen, root, default_depth, visual_id, cmap ) )
        {   /* found visual - now off to decide about color handling on it : */
            if( !setup_truecolor_visual( asv ) )
            {  /* well, we don't - lets try and preallocate as many colors as we can but up to
                * 1/4 of the colorspace or 12bpp colors, whichever is smaller */
                setup_pseudo_visual( asv );
                if( asv->as_colormap == NULL )
                    setup_as_colormap( asv );
            }
        }else
        {
            if( reusable_memory != asv )
                free( asv );
            asv = NULL ;
        }
    }
#endif /*ifndef X_DISPLAY_MISSING */
	_set_default_asvisual( asv );
	return asv;
}

ASVisual *
create_asvisual( Display *dpy, int screen, int default_depth, ASVisual *reusable_memory )
{
	VisualID visual_id = 0;
	char *id_env_var ;

	if( (id_env_var = getenv( ASVISUAL_ID_ENVVAR )) != NULL )
		visual_id = strtol(id_env_var,NULL,16);

	return create_asvisual_for_id( dpy, screen, default_depth, visual_id, None, reusable_memory );
}


void
destroy_asvisual( ASVisual *asv, Bool reusable )
{
	if( asv )
	{
		if( get_default_asvisual() == asv )
			_set_default_asvisual( NULL );
#ifndef X_DISPLAY_MISSING
	 	if( asv->own_colormap )
	 	{
	 		if( asv->colormap )
	 			XFreeColormap( asv->dpy, asv->colormap );
	 	}
	 	if( asv->as_colormap )
		{
	 		free( asv->as_colormap );
			if( asv->as_colormap_reverse.xref != NULL )
			{
				if( asv->as_colormap_type == ACM_12BPP )
					destroy_ashash( &(asv->as_colormap_reverse.hash) );
				else
					free( asv->as_colormap_reverse.xref );
			}
		}
#ifdef HAVE_GLX
		if( asv->glx_scratch_gc_direct )
			glXDestroyContext(asv->dpy, asv->glx_scratch_gc_direct );
		if( asv->glx_scratch_gc_indirect )
			glXDestroyContext(asv->dpy, asv->glx_scratch_gc_indirect );
#endif
		if( asv->scratch_window ) 
			XDestroyWindow( asv->dpy, asv->scratch_window );

#endif /*ifndef X_DISPLAY_MISSING */
		if( !reusable )
			free( asv );
	}
}

int as_colormap_type2size( int type )
{
	switch( type )
	{
		case ACM_3BPP :
			return 8 ;
		case ACM_6BPP :
			return 64 ;
		case ACM_12BPP :
			return 4096 ;
		default:
			return 0 ;
	}
}

Bool
visual2visual_prop( ASVisual *asv, size_t *size_ret,
								   unsigned long *version_ret,
								   unsigned long **data_ret )
{
	int cmap_size = 0 ;
	unsigned long *prop ;
	size_t size;

	if( asv == NULL || data_ret == NULL)
		return False;

	cmap_size = as_colormap_type2size( asv->as_colormap_type );

	if( cmap_size > 0 && asv->as_colormap == NULL )
		return False ;
	size = (1+1+2+1+cmap_size)*sizeof(unsigned long);
	prop = safemalloc( size ) ;
#ifndef X_DISPLAY_MISSING
	prop[0] = asv->visual_info.visualid ;
	prop[1] = asv->colormap ;
	prop[2] = asv->black_pixel ;
	prop[3] = asv->white_pixel ;
	prop[4] = asv->as_colormap_type ;
	if( cmap_size > 0 )
	{
		register int i;
		for( i = 0 ; i < cmap_size ; i++ )
			prop[i+5] = asv->as_colormap[i] ;
	}
	if( size_ret )
		*size_ret = size;
#endif /*ifndef X_DISPLAY_MISSING */
	if( version_ret )
		*version_ret = (1<<16)+0;                        /* version is 1.0 */
	*data_ret = prop ;
	return True;
}

Bool
visual_prop2visual( ASVisual *asv, Display *dpy, int screen,
								   size_t size,
								   unsigned long version,
								   unsigned long *data )
{
#ifndef X_DISPLAY_MISSING
	XVisualInfo templ, *list ;
	int nitems = 0 ;
	int cmap_size = 0 ;
#endif /*ifndef X_DISPLAY_MISSING */

	if( asv == NULL )
		return False;

	asv->dpy = dpy ;

	if( size < (1+1+2+1)*sizeof(unsigned long) ||
		(version&0x00FFFF) != 0 || (version>>16) != 1 || data == NULL )
		return False;

	if( data[0] == None || data[1] == None ) /* we MUST have valid colormap and visualID !!!*/
		return False;

#ifndef X_DISPLAY_MISSING
	templ.screen = screen ;
	templ.visualid = data[0] ;

	list = XGetVisualInfo( dpy, VisualScreenMask|VisualIDMask, &templ, &nitems );
	if( list == NULL || nitems == 0 )
		return False;   /* some very bad visual ID has been requested :( */

	asv->visual_info = *list ;
	XFree( list );

	if( asv->own_colormap && asv->colormap )
		XFreeColormap( dpy, asv->colormap );

	asv->colormap = data[1] ;
	asv->own_colormap = False ;
	asv->black_pixel = data[2] ;
	asv->white_pixel = data[3] ;
	asv->as_colormap_type = data[4];

	cmap_size = as_colormap_type2size( asv->as_colormap_type );

	if( cmap_size > 0 )
	{
		register int i ;
		if( asv->as_colormap )
			free( asv->as_colormap );
		asv->as_colormap = safemalloc( cmap_size );
		for( i = 0 ; i < cmap_size ; i++ )
			asv->as_colormap[i] = data[i+5];
	}else
		asv->as_colormap_type = ACM_None ;     /* just in case */
#else

#endif /*ifndef X_DISPLAY_MISSING */
	return True;
}

Bool
setup_truecolor_visual( ASVisual *asv )
{
#ifndef X_DISPLAY_MISSING
	XVisualInfo *vi = &(asv->visual_info) ;

	if( vi->class != TrueColor )
		return False;

#ifdef HAVE_GLX
	if( glXQueryExtension (asv->dpy, NULL, NULL))
	{
		int val = False;
		glXGetConfig(asv->dpy, vi, GLX_USE_GL, &val);		
		if( val ) 
		{
			asv->glx_scratch_gc_indirect = glXCreateContext (asv->dpy, &(asv->visual_info), NULL, False);
			if( asv->glx_scratch_gc_indirect ) 
			{	
				set_flags( asv->glx_support, ASGLX_Available );
				if( glXGetConfig(asv->dpy, vi, GLX_RGBA, &val) == 0 )
					if( val ) set_flags( asv->glx_support, ASGLX_RGBA );
				if( glXGetConfig(asv->dpy, vi, GLX_DOUBLEBUFFER, &val) == 0 )
					if( val ) set_flags( asv->glx_support, ASGLX_DoubleBuffer );
				if( glXGetConfig(asv->dpy, vi, GLX_DOUBLEBUFFER, &val) == 0 )
					if( val ) set_flags( asv->glx_support, ASGLX_DoubleBuffer );
				
				if( (asv->glx_scratch_gc_direct = glXCreateContext (asv->dpy, &(asv->visual_info), NULL, True)) != NULL ) 
					if( !glXIsDirect( asv->dpy, asv->glx_scratch_gc_direct ) )
					{	
						glXDestroyContext(asv->dpy, asv->glx_scratch_gc_direct );
						asv->glx_scratch_gc_direct = NULL ;
					}
#if 0                          /* that needs some more research :  */
				/* under Cygwin that seems to be 40% faster then regular XImage for some reason */
				set_flags( asv->glx_support, ASGLX_UseForImageTx );
#endif
			}
		}	 
	}	 
#endif

	asv->BGR_mode = ((vi->red_mask&0x0010)!=0) ;
	asv->rshift = get_shifts (vi->red_mask);
	asv->gshift = get_shifts (vi->green_mask);
	asv->bshift = get_shifts (vi->blue_mask);
	asv->rbits = get_bits (vi->red_mask);
	asv->gbits = get_bits (vi->green_mask);
	asv->bbits = get_bits (vi->blue_mask);
	asv->true_depth = vi->depth ;
	asv->msb_first = (ImageByteOrder(asv->dpy)==MSBFirst);

	if( asv->true_depth == 16 && ((vi->red_mask|vi->blue_mask)&0x8000) == 0 )
		asv->true_depth = 15;
	/* setting up conversion handlers : */
	switch( asv->true_depth )
	{
		case 24 :
		case 32 :
			asv->color2pixel_func     = (asv->BGR_mode)?color2pixel32bgr:color2pixel32rgb ;
			asv->pixel2color_func     = (asv->BGR_mode)?pixel2color32bgr:pixel2color32rgb ;
			asv->ximage2scanline_func = ximage2scanline32 ;
			asv->scanline2ximage_func = scanline2ximage32 ;
		    break ;
/*		case 24 :
			scr->color2pixel_func     = (bgr_mode)?color2pixel24bgr:color2pixel24rgb ;
			scr->pixel2color_func     = (bgr_mode)?pixel2color24bgr:pixel2color24rgb ;
			scr->ximage2scanline_func = ximage2scanline24 ;
			scr->scanline2ximage_func = scanline2ximage24 ;
		    break ;
  */	case 16 :
			asv->color2pixel_func     = (asv->BGR_mode)?color2pixel16bgr:color2pixel16rgb ;
			asv->pixel2color_func     = (asv->BGR_mode)?pixel2color16bgr:pixel2color16rgb ;
			asv->ximage2scanline_func = ximage2scanline16 ;
			asv->scanline2ximage_func = scanline2ximage16 ;
		    break ;
		case 15 :
			asv->color2pixel_func     = (asv->BGR_mode)?color2pixel15bgr:color2pixel15rgb ;
			asv->pixel2color_func     = (asv->BGR_mode)?pixel2color15bgr:pixel2color15rgb ;
			asv->ximage2scanline_func = ximage2scanline15 ;
			asv->scanline2ximage_func = scanline2ximage15 ;
		    break ;
	}
#endif /*ifndef X_DISPLAY_MISSING */
	return (asv->ximage2scanline_func != NULL) ;
}

ARGB32 *
make_reverse_colormap( unsigned long *cmap, size_t size, int depth, unsigned short mask, unsigned short shift )
{
	unsigned int max_pixel = 0x01<<depth ;
	ARGB32 *rcmap = safecalloc( max_pixel, sizeof( ARGB32 ) );
	register int i ;

	for( i = 0 ; i < (int)size ; i++ )
		if( cmap[i] < max_pixel )
			rcmap[cmap[i]] = MAKE_ARGB32( 0xFF, (i>>(shift<<1))& mask, (i>>(shift))&mask, i&mask);
	return rcmap;
}

ASHashTable *
make_reverse_colorhash( unsigned long *cmap, size_t size, int depth, unsigned short mask, unsigned short shift )
{
	ASHashTable *hash = create_ashash( 0, NULL, NULL, NULL );
	register unsigned int i ;

	if( hash )
	{
		for( i = 0 ; i < size ; i++ )
			add_hash_item( hash, (ASHashableValue)cmap[i], (void*)((intptr_t)MAKE_ARGB32( 0xFF, (i>>(shift<<1))& mask, (i>>(shift))&mask, i&mask)) );
	}
	return hash;
}

void
setup_pseudo_visual( ASVisual *asv  )
{
#ifndef X_DISPLAY_MISSING
	XVisualInfo *vi = &(asv->visual_info) ;

	/* we need to allocate new usable list of colors based on available bpp */
	asv->true_depth = vi->depth ;
	if( asv->as_colormap == NULL )
	{
		if( asv->true_depth < 8 )
			asv->as_colormap_type = ACM_3BPP ;
		else if( asv->true_depth < 12 )
			asv->as_colormap_type = ACM_6BPP ;
		else
			asv->as_colormap_type = ACM_12BPP ;
	}
	/* then we need to set up hooks : */
	switch( asv->as_colormap_type )
	{
		case ACM_3BPP:
			asv->ximage2scanline_func = ximage2scanline_pseudo3bpp ;
			asv->scanline2ximage_func = scanline2ximage_pseudo3bpp ;
			asv->color2pixel_func = color2pixel_pseudo3bpp ;
		    break ;
		case ACM_6BPP:
			asv->ximage2scanline_func = ximage2scanline_pseudo6bpp ;
			asv->scanline2ximage_func = scanline2ximage_pseudo6bpp ;
			asv->color2pixel_func = color2pixel_pseudo6bpp ;
		    break ;
		default:
			asv->as_colormap_type = ACM_12BPP ;
		case ACM_12BPP:
			asv->ximage2scanline_func = ximage2scanline_pseudo12bpp ;
			asv->scanline2ximage_func = scanline2ximage_pseudo12bpp ;
			asv->color2pixel_func = color2pixel_pseudo12bpp ;
		    break ;
	}
	if( asv->as_colormap != NULL )
	{
		if( asv->as_colormap_type == ACM_3BPP || asv->as_colormap_type == ACM_6BPP )
		{
			unsigned short mask = 0x0003, shift = 2 ;
			if( asv->as_colormap_type==ACM_3BPP )
			{
				mask = 0x0001 ;
				shift = 1 ;
			}
			asv->as_colormap_reverse.xref = make_reverse_colormap( asv->as_colormap,
															  as_colormap_type2size( asv->as_colormap_type ),
															  asv->true_depth, mask, shift );
		}else if( asv->as_colormap_type == ACM_12BPP )
		{
			asv->as_colormap_reverse.hash = make_reverse_colorhash( asv->as_colormap,
															  as_colormap_type2size( asv->as_colormap_type ),
															  asv->true_depth, 0x000F, 4 );
		}
	}
#endif /*ifndef X_DISPLAY_MISSING */
}

#ifndef X_DISPLAY_MISSING
static unsigned long*
make_3bpp_colormap( ASVisual *asv )
{
	XColor colors_3bpp[8] =
	/* list of non-white, non-black colors in order of decreasing importance: */
	{   { 0, 0, 0xFFFF, 0, DoRed|DoGreen|DoBlue, 0},
		{ 0, 0xFFFF, 0, 0, DoRed|DoGreen|DoBlue, 0},
		{ 0, 0, 0, 0xFFFF, DoRed|DoGreen|DoBlue, 0},
	 	{ 0, 0xFFFF, 0xFFFF, 0, DoRed|DoGreen|DoBlue, 0},
	    { 0, 0, 0xFFFF, 0xFFFF, DoRed|DoGreen|DoBlue, 0},
	    { 0, 0xFFFF, 0, 0xFFFF, DoRed|DoGreen|DoBlue, 0}} ;
	unsigned long *cmap ;

	cmap = safemalloc( 8 * sizeof(unsigned long) );
	/* fail safe code - if any of the alloc fails - colormap entry will still have
	 * most suitable valid value ( black or white in 1bpp mode for example ) : */
	cmap[0] = cmap[1] = cmap[2] = cmap[3] = asv->black_pixel ;
	cmap[7] = cmap[6] = cmap[5] = cmap[4] = asv->white_pixel ;
	if( ASV_ALLOC_COLOR( asv, asv->colormap, &colors_3bpp[0] ))  /* pure green */
		cmap[0x02] = cmap[0x03] = cmap[0x06] = colors_3bpp[0].pixel ;
	if( ASV_ALLOC_COLOR( asv, asv->colormap, &colors_3bpp[1] ))  /* pure red */
		cmap[0x04] = cmap[0x05] = colors_3bpp[1].pixel ;
	if( ASV_ALLOC_COLOR( asv, asv->colormap, &colors_3bpp[2] ))  /* pure blue */
		cmap[0x01] = colors_3bpp[2].pixel ;
	if( ASV_ALLOC_COLOR( asv, asv->colormap, &colors_3bpp[3] ))  /* yellow */
		cmap[0x06] = colors_3bpp[3].pixel ;
	if( ASV_ALLOC_COLOR( asv, asv->colormap, &colors_3bpp[4] ))  /* cyan */
		cmap[0x03] = colors_3bpp[4].pixel ;
	if( ASV_ALLOC_COLOR( asv, asv->colormap, &colors_3bpp[5] ))  /* magenta */
		cmap[0x05] = colors_3bpp[5].pixel ;
	return cmap;
}

static unsigned long*
make_6bpp_colormap( ASVisual *asv, unsigned long *cmap_3bpp )
{
	unsigned short red, green, blue ;
	unsigned long *cmap = safemalloc( 0x0040*sizeof( unsigned long) );
	XColor xcol ;

	cmap[0] = asv->black_pixel ;

	xcol.flags = DoRed|DoGreen|DoBlue ;
	for( blue = 0 ; blue <= 0x0003 ; blue++ )
	{
		xcol.blue = (0xFFFF*blue)/3 ;
		for( red = 0 ; red <= 0x0003 ; red++ )
		{	                                /* red has highier priority then blue */
			xcol.red = (0xFFFF*red)/3 ;
/*			green = ( blue == 0 && red == 0 )?1:0 ; */
			for( green = 0 ; green <= 0x0003 ; green++ )
			{                                  /* green has highier priority then red */
				unsigned short index_3bpp = ((red&0x0002)<<1)|(green&0x0002)|((blue&0x0002)>>1);
				unsigned short index_6bpp = (red<<4)|(green<<2)|blue;
				xcol.green = (0xFFFF*green)/3 ;

				if( (red&0x0001) == ((red&0x0002)>>1) &&
					(green&0x0001) == ((green&0x0002)>>1) &&
					(blue&0x0001) == ((blue&0x0002)>>1) )
					cmap[index_6bpp] = cmap_3bpp[index_3bpp];
				else
				{
					if( ASV_ALLOC_COLOR( asv, asv->colormap, &xcol) != 0 )
						cmap[index_6bpp] = xcol.pixel ;
					else
						cmap[index_6bpp] = cmap_3bpp[index_3bpp] ;
				}
			}
		}
	}
	return cmap;
}

static unsigned long*
make_9bpp_colormap( ASVisual *asv, unsigned long *cmap_6bpp )
{
	unsigned long *cmap = safemalloc( 512*sizeof( unsigned long) );
	unsigned short red, green, blue ;
	XColor xcol ;

	cmap[0] = asv->black_pixel ;               /* just in case  */

	xcol.flags = DoRed|DoGreen|DoBlue ;
	for( blue = 0 ; blue <= 0x0007 ; blue++ )
	{
		xcol.blue = (0xFFFF*blue)/7 ;
		for( red = 0 ; red <= 0x0007 ; red++ )
		{	                                /* red has highier priority then blue */
			xcol.red = (0xFFFF*red)/7 ;
			for( green = 0 ; green <= 0x0007 ; green++ )
			{                                  /* green has highier priority then red */
				unsigned short index_6bpp = ((red&0x0006)<<3)|((green&0x0006)<<1)|((blue&0x0006)>>1);
				unsigned short index_9bpp = (red<<6)|(green<<3)|blue;
				xcol.green = (0xFFFF*green)/7 ;

				if( (red&0x0001) == ((red&0x0002)>>1) &&
					(green&0x0001) == ((green&0x0002)>>1) &&
					(blue&0x0001) == ((blue&0x0002)>>1) )
					cmap[index_9bpp] = cmap_6bpp[index_6bpp];
				else
				{
					if( ASV_ALLOC_COLOR( asv, asv->colormap, &xcol) != 0 )
						cmap[index_9bpp] = xcol.pixel ;
					else
						cmap[index_9bpp] = cmap_6bpp[index_6bpp] ;
				}
			}
		}
	}
	return cmap;
}

static unsigned long*
make_12bpp_colormap( ASVisual *asv, unsigned long *cmap_9bpp )
{
	unsigned long *cmap = safemalloc( 4096*sizeof( unsigned long) );
	unsigned short red, green, blue ;
	XColor xcol ;

	cmap[0] = asv->black_pixel ;               /* just in case  */

	xcol.flags = DoRed|DoGreen|DoBlue ;
	for( blue = 0 ; blue <= 0x000F ; blue++ )
	{
		xcol.blue = (0xFFFF*blue)/15 ;
		for( red = 0 ; red <= 0x000F ; red++ )
		{	                                /* red has highier priority then blue */
			xcol.red = (0xFFFF*red)/15 ;
			for( green = 0 ; green <= 0x000F ; green++ )
			{                                  /* green has highier priority then red */
				unsigned short index_9bpp = ((red&0x000E)<<5)|((green&0x000E)<<2)|((blue&0x000E)>>1);
				unsigned short index_12bpp = (red<<8)|(green<<4)|blue;
				xcol.green = (0xFFFF*green)/15 ;

				if( (red&0x0001) == ((red&0x0002)>>1) &&
					(green&0x0001) == ((green&0x0002)>>1) &&
					(blue&0x0001) == ((blue&0x0002)>>1) )
					cmap[index_12bpp] = cmap_9bpp[index_9bpp];
				else
				{
					if( ASV_ALLOC_COLOR( asv, asv->colormap, &xcol) != 0 )
						cmap[index_12bpp] = xcol.pixel ;
					else
						cmap[index_12bpp] = cmap_9bpp[index_9bpp] ;
				}
			}
		}
	}
	return cmap;
}
#endif /*ifndef X_DISPLAY_MISSING */

void
setup_as_colormap( ASVisual *asv )
{
#ifndef X_DISPLAY_MISSING
	unsigned long *cmap_lower, *cmap ;

	if( asv == NULL || asv->as_colormap != NULL )
		return ;

	cmap = make_3bpp_colormap( asv );
	if( asv->as_colormap_type == ACM_3BPP )
	{
		asv->as_colormap = cmap ;
		asv->as_colormap_reverse.xref = make_reverse_colormap( cmap, 8, asv->true_depth, 0x0001, 1 );
		return ;
	}
	cmap_lower = cmap ;
	cmap = make_6bpp_colormap( asv, cmap_lower );
	free( cmap_lower );
	if( asv->as_colormap_type == ACM_6BPP )
	{
		asv->as_colormap = cmap ;
		asv->as_colormap_reverse.xref = make_reverse_colormap( cmap, 64, asv->true_depth, 0x0003, 2 );
	}else
	{
		cmap_lower = cmap ;
		cmap = make_9bpp_colormap( asv, cmap_lower );
		free( cmap_lower );
		cmap_lower = cmap ;
		cmap = make_12bpp_colormap( asv, cmap_lower );
		free( cmap_lower );

		asv->as_colormap = cmap ;
		asv->as_colormap_reverse.hash = make_reverse_colorhash( cmap, 4096, asv->true_depth, 0x000F, 4 );
	}
#endif /*ifndef X_DISPLAY_MISSING */
}

/*********************************************************************/
/* handy utility functions for creation of windows/pixmaps/XImages : */
/*********************************************************************/
Window
create_visual_window( ASVisual *asv, Window parent,
					  int x, int y, unsigned int width, unsigned int height,
					  unsigned int border_width, unsigned int wclass,
 					  unsigned long mask, XSetWindowAttributes *attributes )
{
#ifndef X_DISPLAY_MISSING
	XSetWindowAttributes my_attr ;
	int depth = 0;

	if( asv == NULL || parent == None )
		return None ;
LOCAL_DEBUG_OUT( "Colormap %lX, parent %lX, %ux%u%+d%+d, bw = %d, class %d",
				  asv->colormap, parent, width, height, x, y, border_width,
				  wclass );
	if( attributes == NULL )
	{
		attributes = &my_attr ;
		memset( attributes, 0x00, sizeof(XSetWindowAttributes));
		mask = 0;
	}

	if( width < 1 )
		width = 1 ;
	if( height < 1 )
		height = 1 ;

	if( wclass == InputOnly )
	{
		border_width = 0 ;
		if( (mask&INPUTONLY_LEGAL_MASK) != mask )
			show_warning( " software BUG detected : illegal InputOnly window's mask 0x%lX - overriding", mask );
		mask &= INPUTONLY_LEGAL_MASK ;
	}else
	{
		depth = asv->visual_info.depth ;
		if( !get_flags(mask, CWColormap ) )
		{
			attributes->colormap = asv->colormap ;
			set_flags(mask, CWColormap );
		}
		if( !get_flags(mask, CWBorderPixmap ) )
		{
			attributes->border_pixmap = None ;
			set_flags(mask, CWBorderPixmap );
		}

		clear_flags(mask, CWBorderPixmap );
		if( !get_flags(mask, CWBorderPixel ) )
		{
			attributes->border_pixel = asv->black_pixel ;
			set_flags(mask, CWBorderPixel );
		}
		/* If the parent window and the new window have different bit
		** depths (such as on a Solaris box with 8bpp root window and
		** 24bpp child windows), ParentRelative will not work. */
		if ( get_flags(mask, CWBackPixmap) && attributes->background_pixmap == ParentRelative &&
			 asv->visual_info.visual != DefaultVisual( asv->dpy, DefaultScreen(asv->dpy) ))
		{
			clear_flags(mask, CWBackPixmap);
		}
	}
	LOCAL_DEBUG_OUT( "parent = %lX, mask = 0x%lX, VisualID = 0x%lX, Border Pixel = %ld, colormap = %lX",
					  parent, mask, asv->visual_info.visual->visualid, attributes->border_pixel, attributes->colormap );
	return XCreateWindow (asv->dpy, parent, x, y, width, height, border_width, depth,
						  wclass, asv->visual_info.visual,
	                      mask, attributes);
#else
	return None ;
#endif /*ifndef X_DISPLAY_MISSING */

}


GC
create_visual_gc( ASVisual *asv, Window root, unsigned long mask, XGCValues *gcvalues )
{
   	GC gc = NULL ;

#ifndef X_DISPLAY_MISSING
	if( asv )
	{
		XGCValues scratch_gcv ;
		if( asv->scratch_window == None ) 
			asv->scratch_window = create_visual_window( asv, root, -20, -20, 10, 10, 0, InputOutput, 0, NULL );
		if( asv->scratch_window != None )
			gc = XCreateGC( asv->dpy, asv->scratch_window, gcvalues?mask:0, gcvalues?gcvalues:&scratch_gcv );
	}
#endif
	return gc;
}

Pixmap
create_visual_pixmap( ASVisual *asv, Window root, unsigned int width, unsigned int height, unsigned int depth )
{
#ifndef X_DISPLAY_MISSING
	Pixmap p = None ;
	if( asv != NULL )
	{	
		if( root == None ) 
			root = RootWindow(asv->dpy,DefaultScreen(asv->dpy));
		if( depth==0 )
			depth = asv->true_depth ;
		p = XCreatePixmap( asv->dpy, root, MAX(width,(unsigned)1), MAX(height,(unsigned)1), depth );
	}
	return p;
#else
	return None ;
#endif /*ifndef X_DISPLAY_MISSING */
}

void
destroy_visual_pixmap( ASVisual *asv, Pixmap *ppmap )
{
	if( asv && ppmap )
		if( *ppmap )
		{
#ifndef X_DISPLAY_MISSING
			XFreePixmap( asv->dpy, *ppmap );
			*ppmap = None ;
#endif
		}
}

#ifndef X_DISPLAY_MISSING
static int
quiet_xerror_handler (Display * dpy, XErrorEvent * error)
{
    return 0;
}

#endif

int
get_dpy_drawable_size (Display *drawable_dpy, Drawable d, unsigned int *ret_w, unsigned int *ret_h)
{
	int result = 0 ;
#ifndef X_DISPLAY_MISSING
	if( d != None && drawable_dpy != NULL ) 
	{
		Window        root;
		unsigned int  ujunk;
		int           junk;
		int           (*oldXErrorHandler) (Display *, XErrorEvent *) = XSetErrorHandler (quiet_xerror_handler);
		result = XGetGeometry (drawable_dpy, d, &root, &junk, &junk, ret_w, ret_h, &ujunk, &ujunk);
		XSetErrorHandler (oldXErrorHandler);
	}
#endif
	if ( result == 0)
	{
		*ret_w = 0;
		*ret_h = 0;
		return 0;
	}
	return 1;
}

Bool
get_dpy_window_position (Display *window_dpy, Window root, Window w, int *px, int *py, int *transparency_x, int *transparency_y)
{
	Bool result = False ;
	int x = 0, y = 0, transp_x = 0, transp_y = 0 ;
#ifndef X_DISPLAY_MISSING
	if( window_dpy != NULL && w != None ) 
	{
		Window wdumm;
		int rootHeight = XDisplayHeight(window_dpy, DefaultScreen(window_dpy) );
		int rootWidth = XDisplayWidth(window_dpy, DefaultScreen(window_dpy) );

		if( root == None ) 
			root = RootWindow(window_dpy,DefaultScreen(window_dpy));
			
		result = XTranslateCoordinates (window_dpy, w, root, 0, 0, &x, &y, &wdumm);
		if( result ) 
		{
			/* taking in to consideration virtual desktopping */
			result = (x < rootWidth && y < rootHeight );
			if( result )
			{
				unsigned int width = 0, height = 0;
				get_dpy_drawable_size (window_dpy, w, &width, &height);				
				result = (x + width > 0 && y+height > 0) ; 
			}

			for( transp_x = x ; transp_x < 0 ; transp_x += rootWidth ); 			
			for( transp_y = y ; transp_y < 0 ; transp_y += rootHeight ); 			
			while( transp_x > rootWidth ) transp_x -= rootWidth ; 
			while( transp_y > rootHeight ) transp_y -= rootHeight ; 
		}
	}
#endif
	if( px ) 
		*px = x;
	if( py ) 
		*py = y;
	if( transparency_x ) 
		*transparency_x = transp_x ; 
	if( transparency_y ) 
		*transparency_y = transp_y ; 
	return result;
}


#ifndef X_DISPLAY_MISSING
static unsigned char *scratch_ximage_data = NULL ;
static int scratch_use_count = 0 ;
static size_t scratch_ximage_allocated_size = 0;  
#endif 
static size_t scratch_ximage_max_size = ASSHM_SAVED_MAX*2;  /* maximum of 512 KBytes is default  */  
static size_t scratch_ximage_normal_size = ASSHM_SAVED_MAX;  /* normal usage of scratch pool is 256 KBytes is default  */  

int
set_scratch_ximage_max_size( int new_max_size )
{
	int tmp = scratch_ximage_max_size ;
	scratch_ximage_max_size = new_max_size ;
	return tmp;
}

int
set_scratch_ximage_normal_size( int new_normal_size )
{
	int tmp = scratch_ximage_normal_size ;
	scratch_ximage_normal_size = new_normal_size ;
	return tmp;
}

#ifndef X_DISPLAY_MISSING
static void*
get_scratch_data(size_t size)
{
	if( scratch_ximage_max_size < size || scratch_use_count > 0) 
		return NULL;
	if( scratch_ximage_allocated_size < size ) 
	{
		scratch_ximage_allocated_size = size ;
		scratch_ximage_data = realloc( scratch_ximage_data, size );
	}
	
	++scratch_use_count;
	return scratch_ximage_data ; 
}

static Bool
release_scratch_data( void *data )
{
	if( scratch_use_count == 0 || data != scratch_ximage_data )
		return False;
	--scratch_use_count ;
	if( scratch_use_count == 0 )
	{
		/* want to deallocate if too much is allocated ? */
		
	}	 
	return True;
}	 
#endif 

#ifdef XSHMIMAGE

int	(*orig_XShmImage_destroy_image)(XImage *ximage) = NULL ;

typedef struct ASXShmImage
{
	XImage 			*ximage ;
	XShmSegmentInfo *segment ;
	int 			 ref_count ;
	Bool			 wait_completion_event ;
	unsigned int 	 size ;
  ASVisual *asv ;
}ASXShmImage;

typedef struct ASShmArea
{
	unsigned int 	 size ;
	char *shmaddr ;
	int shmid ;
	struct ASShmArea *next, *prev ;
}ASShmArea;

static ASHashTable	*xshmimage_segments = NULL ;
static ASHashTable	*xshmimage_images = NULL ;
/* attempt to reuse 256 Kb of shmem - no reason to reuse more than that,
 * since most XImages will be in range of 20K-100K */
static ASShmArea  *shm_available_mem_head = NULL ;
static int shm_available_mem_used = 0 ;

static Bool _as_use_shm_images = False ;

void really_destroy_shm_area( char *shmaddr, int shmid )
{
	shmdt (shmaddr);
	shmctl (shmid, IPC_RMID, 0);
	LOCAL_DEBUG_OUT("XSHMIMAGE> DESTROY_SHM : freeing shmid = %d, remaining in cache = %d bytes ", shmid, shm_available_mem_used );
}

void remove_shm_area( ASShmArea *area, Bool free_resources )
{
	if( area )
	{
		if( area == shm_available_mem_head )
			shm_available_mem_head = area->next ;
		if( area->next )
			area->next->prev = area->prev ;
		if( area->prev )
			area->prev->next = area->next ;
		shm_available_mem_used -= area->size ;
		if( free_resources )
			really_destroy_shm_area( area->shmaddr, area->shmid );
		else
		{
			LOCAL_DEBUG_OUT("XSHMIMAGE> REMOVE_SHM : reusing shmid = %d, size %d, remaining in cache = %d bytes ", area->shmid, area->size, shm_available_mem_used );
		}
		free( area );
	}

}

void flush_shm_cache( )
{
	if( xshmimage_images ) 
		destroy_ashash( &xshmimage_images );
	if( xshmimage_segments )
		destroy_ashash( &xshmimage_segments );
	while( shm_available_mem_head != NULL )
		remove_shm_area( shm_available_mem_head, True );
}

void save_shm_area( char *shmaddr, int shmid, int size )
{
	ASShmArea *area;

	if( shm_available_mem_used+size >= ASSHM_SAVED_MAX )
	{
	  	really_destroy_shm_area( shmaddr, shmid );
		return ;
	}

	shm_available_mem_used+=size ;
	area = safecalloc( 1, sizeof(ASShmArea) );

	area->shmaddr = shmaddr ;
	area->shmid = shmid ;
	area->size = size ;
	LOCAL_DEBUG_OUT("XSHMIMAGE> SAVE_SHM : saving shmid = %d, size %d, remaining in cache = %d bytes ", area->shmid, area->size, shm_available_mem_used );

	area->next = shm_available_mem_head ;
	if( shm_available_mem_head )
		shm_available_mem_head->prev = area ;
	shm_available_mem_head = area ;
}

char *get_shm_area( int size, int *shmid )
{
	ASShmArea *selected = NULL, *curr = shm_available_mem_head;

	while( curr != NULL )
	{
		if( curr->size >= size && curr->size < (size * 4)/3 )
		{
			if( selected == NULL )
				selected = curr ;
			else if( selected->size > curr->size )
				selected = curr ;
		}
		curr = curr->next ;
	}
	if( selected != NULL )
	{
		char *tmp = selected->shmaddr ;
		*shmid = selected->shmid ;
		remove_shm_area( selected, False );
		return tmp ;
	}

	*shmid = shmget (IPC_PRIVATE, size, IPC_CREAT|0666);
	return shmat (*shmid, 0, 0);
}

void
destroy_xshmimage_segment(ASHashableValue value, void *data)
{
	ASXShmImage *img_data = (ASXShmImage*)data ;
	if( img_data->segment != NULL )
	{
		LOCAL_DEBUG_OUT( "XSHMIMAGE> FREE_SEG : img_data = %p : segent to be freed: shminfo = %p ", img_data, img_data->segment );
		XShmDetach (img_data->asv->dpy, img_data->segment);
		save_shm_area( img_data->segment->shmaddr, img_data->segment->shmid, img_data->size );
		free( img_data->segment );
		img_data->segment = NULL ;
		if( img_data->ximage == NULL )
			free( img_data );
	}else
	{
		LOCAL_DEBUG_OUT( "XSHMIMAGE> FREE_SEG : img_data = %p : segment data is NULL already value = %ld!!", img_data, value );
	}
}

Bool destroy_xshm_segment( ShmSeg shmseg )
{
	if( xshmimage_segments )
	{
		if(remove_hash_item( xshmimage_segments, AS_HASHABLE(shmseg), NULL, True ) == ASH_Success)
		{
			LOCAL_DEBUG_OUT( "XSHMIMAGE> REMOVE_SEG : segment %ld removed from the hash successfully!", shmseg );
			return True ;
		}
		LOCAL_DEBUG_OUT( "XSHMIMAGE> ERROR : could not find segment %ld(0x%lX) in the hash!", shmseg, shmseg );
	}else
	{
		LOCAL_DEBUG_OUT( "XSHMIMAGE> ERROR : segments hash is %p!!", xshmimage_segments );
	}

	return False ;
}


void
destroy_xshmimage_image(ASHashableValue value, void *data)
{
	ASXShmImage *img_data = (ASXShmImage*)data ;
	if( img_data->ximage != NULL )
	{
		if( orig_XShmImage_destroy_image )
			orig_XShmImage_destroy_image( img_data->ximage );
		else
			XFree ((char *)img_data->ximage);
		LOCAL_DEBUG_OUT( "XSHMIMAGE> FREE_XIM : ximage freed: img_data = %p, xim = %p", img_data, img_data->ximage);
		img_data->ximage = NULL ;
		if( img_data->segment != NULL && !img_data->wait_completion_event )
		{
			if( destroy_xshm_segment( img_data->segment->shmseg ) )
				return ;
			img_data->segment = NULL ;
		}
		if( img_data->segment == NULL )
			free( img_data );
	}
}

Bool enable_shmem_images_for_visual (ASVisual *asv)
{
#ifndef DEBUG_ALLOCS
	if( asv && asv->dpy && XShmQueryExtension (asv->dpy) )
	{
		_as_use_shm_images = True ;
		if( xshmimage_segments == NULL )
			xshmimage_segments = create_ashash( 0, NULL, NULL, destroy_xshmimage_segment );
		if( xshmimage_images == NULL )
			xshmimage_images = create_ashash( 0, pointer_hash_value, NULL, destroy_xshmimage_image );
	}else
#endif
		_as_use_shm_images = False ;
	return _as_use_shm_images;
}

Bool enable_shmem_images ()
{
	return enable_shmem_images_for_visual (get_default_asvisual());
}


void disable_shmem_images()
{
	_as_use_shm_images = False ;
}

Bool 
check_shmem_images_enabled()
{
	return _as_use_shm_images ;
}


int destroy_xshm_image( XImage *ximage )
{
	if( xshmimage_images )
	{
		if( remove_hash_item( xshmimage_images, AS_HASHABLE(ximage), NULL, True ) != ASH_Success )
		{
			if (ximage->data != NULL)
				free ((char *)ximage->data);
			if (ximage->obdata != NULL)
				free ((char *)ximage->obdata);
			XFree ((char *)ximage);
			LOCAL_DEBUG_OUT( "XSHMIMAGE> FREE_XIM : ximage freed: xim = %p", ximage);
		}
	}
	return 1;
}

unsigned long 
ximage2shmseg( XImage *xim )
{
	void *vptr = NULL ;
	if( get_hash_item( xshmimage_images, AS_HASHABLE(xim), &vptr ) == ASH_Success )		
	{
		ASXShmImage *data = (ASXShmImage *)vptr ;
		if( data->segment ) 
			return data->segment->shmseg;	
	}	
	return 0; 
}	 

void registerXShmImage( ASVisual *asv, XImage *ximage, XShmSegmentInfo* shminfo )
{
	ASXShmImage *data = safecalloc( 1, sizeof(ASXShmImage));
	LOCAL_DEBUG_OUT( "XSHMIMAGE> CREATE_XIM : img_data = %p : image created: xiom = %p, shminfo = %p, segment = %d, data = %p", data, ximage, shminfo, shminfo->shmid, ximage->data );
  data->asv = asv ;
	data->ximage = ximage ;
	data->segment = shminfo ;
	data->size = ximage->bytes_per_line * ximage->height ;

	orig_XShmImage_destroy_image = ximage->f.destroy_image ;
	ximage->f.destroy_image = destroy_xshm_image ;

	add_hash_item( xshmimage_images, AS_HASHABLE(ximage), data );
	add_hash_item( xshmimage_segments, AS_HASHABLE(shminfo->shmseg), data );
}

void *
check_XImage_shared( XImage *xim )
{
	ASXShmImage *img_data = NULL ;
	if( _as_use_shm_images )
	{
		ASHashData hdata ;
		if(get_hash_item( xshmimage_images, AS_HASHABLE(xim), &hdata.vptr ) != ASH_Success)
			img_data = NULL ;
		else
			img_data = hdata.vptr ;
	}
	return img_data ;
}

Bool ASPutXImage( ASVisual *asv, Drawable d, GC gc, XImage *xim,
                  int src_x, int src_y, int dest_x, int dest_y,
				  unsigned int width, unsigned int height )
{
	ASXShmImage *img_data = NULL ;
	if( xim == NULL || asv == NULL )
		return False ;

	if( ( img_data = check_XImage_shared( xim )) != NULL )
	{
/*		LOCAL_DEBUG_OUT( "XSHMIMAGE> PUT_XIM : using shared memory Put = %p", xim ); */
		if( XShmPutImage( asv->dpy, d, gc, xim, src_x, src_y, dest_x, dest_y,width, height, True ) )
		{
			img_data->wait_completion_event = True ;
			return True ;
		}
	}
/*	LOCAL_DEBUG_OUT( "XSHMIMAGE> PUT_XIM : using normal Put = %p", xim ); */
	return XPutImage( asv->dpy, d, gc, xim, src_x, src_y, dest_x, dest_y,width, height );
}

XImage *ASGetXImage( ASVisual *asv, Drawable d,
                  int x, int y, unsigned int width, unsigned int height,
				  unsigned long plane_mask )
{
	XImage *xim = NULL ;

	if( asv == NULL || d == None )
		return NULL ;
	if( _as_use_shm_images && width*height > 4000)
	{
		unsigned int depth ;
		Window        root;
		unsigned int  ujunk;
		int           junk;
		if(XGetGeometry (asv->dpy, d, &root, &junk, &junk, &ujunk, &ujunk, &ujunk, &depth) == 0)
			return NULL ;

		xim = create_visual_ximage(asv,width,height,depth);
		XShmGetImage( asv->dpy, d, xim, x, y, plane_mask );

	}else
		xim = XGetImage( asv->dpy, d, x, y, width, height, plane_mask, ZPixmap );
	return xim ;
}
#else

Bool enable_shmem_images (){return False; }
void disable_shmem_images(){}
void *check_XImage_shared( XImage *xim ) {return NULL ; }

Bool ASPutXImage( ASVisual *asv, Drawable d, GC gc, XImage *xim,
                  int src_x, int src_y, int dest_x, int dest_y,
				  unsigned int width, unsigned int height )
{
#ifndef X_DISPLAY_MISSING
	if( xim == NULL || asv == NULL )
		return False ;
	return XPutImage( asv->dpy, d, gc, xim, src_x, src_y, dest_x, dest_y,width, height );
#else
	return False;
#endif
}

XImage * ASGetXImage( ASVisual *asv, Drawable d,
                  int x, int y, unsigned int width, unsigned int height,
				  unsigned long plane_mask )
{
#ifndef X_DISPLAY_MISSING
	if( asv == NULL || d == None )
		return NULL ;
	return XGetImage( asv->dpy, d, x, y, width, height, plane_mask, ZPixmap );
#else
	return NULL ;
#endif	
}

#endif                                         /* XSHMIMAGE */

#ifndef X_DISPLAY_MISSING
int
My_XDestroyImage (XImage *ximage)
{
	if( !release_scratch_data(ximage->data) )
		if (ximage->data != NULL)
			free (ximage->data);
	if (ximage->obdata != NULL)
		free (ximage->obdata);
	XFree (ximage);
	return 1;
}
#endif /*ifndef X_DISPLAY_MISSING */


XImage*
create_visual_ximage( ASVisual *asv, unsigned int width, unsigned int height, unsigned int depth )
{
#ifndef X_DISPLAY_MISSING
	register XImage *ximage = NULL;
	unsigned long dsize;
	char         *data;
	int unit ;

	if( asv == NULL )
		return NULL;

#if 0
	unit = asv->dpy->bitmap_unit;
#else
	if( depth == 0 ) 
		unit = (asv->true_depth+7)&0x0038;
	else
		unit = (depth+7)&0x0038;
	if( unit == 24 )
		unit = 32 ;
#endif
#ifdef XSHMIMAGE
	if( _as_use_shm_images && width*height > 4000 )
	{
		XShmSegmentInfo *shminfo = safecalloc( 1, sizeof(XShmSegmentInfo));

		ximage = XShmCreateImage (asv->dpy, asv->visual_info.visual,
			                      (depth==0)?asv->visual_info.depth/*true_depth*/:depth,
								  ZPixmap, NULL, shminfo,
								  MAX(width,(unsigned int)1), MAX(height,(unsigned int)1));
		if( ximage == NULL )
			free( shminfo );
		else
		{
			shminfo->shmaddr = ximage->data = get_shm_area( ximage->bytes_per_line * ximage->height, &(shminfo->shmid) );
			if( shminfo->shmid == -1 )
			{
				static int shmem_failure_count = 0 ;
			    show_warning( "unable to allocate %d bytes of shared image memory", ximage->bytes_per_line * ximage->height ) ;
				if( ximage->bytes_per_line * ximage->height < 100000 || ++shmem_failure_count > 10 )
				{
					show_error( "too many shared memory failures - disabling" ) ;
					_as_use_shm_images = False ;
				}
				free( shminfo );
				shminfo = NULL ;
				XFree( ximage );
				ximage = NULL ;
			}else
			{
				shminfo->readOnly = False;
				XShmAttach (asv->dpy, shminfo);
				registerXShmImage( asv, ximage, shminfo );
			}
		}
	}
#endif
	if( ximage == NULL )
	{
		ximage = XCreateImage (asv->dpy, asv->visual_info.visual, (depth==0)?asv->visual_info.depth/*true_depth*/:depth, ZPixmap,
                                       0, NULL, MAX(width,(unsigned int)1), MAX(height,(unsigned int)1),
                                       unit, 0);
		if (ximage != NULL)
		{
			_XInitImageFuncPtrs (ximage);
			ximage->obdata = NULL;
			ximage->f.destroy_image = My_XDestroyImage;
			dsize = ximage->bytes_per_line*ximage->height;
	    	if (((data = (char *)safemalloc (dsize)) == NULL) && (dsize > 0))
			{
				XFree ((char *)ximage);
				return (XImage *) NULL;
			}
			ximage->data = data;
		}
	}
	return ximage;
#else
	return NULL ;
#endif /*ifndef X_DISPLAY_MISSING */
}
/* this is the vehicle to use static allocated buffer for temporary XImages 
 * in order to reduce XImage meory allocation overhead */
XImage*
create_visual_scratch_ximage( ASVisual *asv, unsigned int width, unsigned int height, unsigned int depth )
{
#ifndef X_DISPLAY_MISSING
	register XImage *ximage = NULL;
	char         *data;
	int unit ;

	if( asv == NULL )
		return NULL;

#if 0
	unit = asv->dpy->bitmap_unit;
#else
	if( depth == 0 )
		unit = (asv->true_depth+7)&0x0038;
	else
		unit = (depth+7)&0x0038;
	if( unit == 24 )
		unit = 32 ;
#endif

	/* for shared memory XImage we already do caching - no need for scratch ximage */	   
#ifdef XSHMIMAGE
	if( _as_use_shm_images )
		return create_visual_ximage( asv, width, height, depth );
#endif
		   
	if( ximage == NULL )
	{
		ximage = XCreateImage (asv->dpy, asv->visual_info.visual, 
                                       (depth==0)?asv->visual_info.depth/*true_depth*/:depth, ZPixmap, 
                                       0, NULL, MAX(width,(unsigned int)1), MAX(height,(unsigned int)1),
                                       unit, 0);
		if (ximage != NULL)
		{
			data = get_scratch_data(ximage->bytes_per_line * ximage->height);
			if( data == NULL ) 
			{
				XFree ((char *)ximage);	
				return create_visual_ximage( asv, width, height, depth );/* fall back */
			}	 
			_XInitImageFuncPtrs (ximage);
			ximage->obdata = NULL;
			ximage->f.destroy_image = My_XDestroyImage;
			ximage->data = data;
		}
	}
	return ximage;
#else
	return NULL ;
#endif /*ifndef X_DISPLAY_MISSING */
}



/****************************************************************************/
/* Color manipulation functions :                                           */
/****************************************************************************/


/* misc function to calculate number of bits/shifts */
#ifndef X_DISPLAY_MISSING
static int
get_shifts (unsigned long mask)
{
	register int  i = 1;

	while (mask >> i)
		i++;

	return i - 1;							   /* can't be negative */
}

static int
get_bits (unsigned long mask)
{
	register int  i;

	for (i = 0; mask; mask >>= 1)
		if (mask & 1)
			i++;

	return i;								   /* can't be negative */
}
#endif
/***************************************************************************/
/* Screen color format -> AS color format conversion ; 					   */
/***************************************************************************/
/***************************************************************************/
/* Screen color format -> AS color format conversion ; 					   */
/***************************************************************************/
/* this functions convert encoded color values into real pixel values, and
 * return half of the quantization error so we can do error diffusion : */
/* encoding scheme : 0RRRrrrr rrrrGGgg ggggggBB bbbbbbbb
 * where RRR, GG and BB are overflow bits so we can do all kinds of funky
 * combined adding, note that we don't use 32'd bit as it is a sign bit */

#ifndef X_DISPLAY_MISSING
static inline void
query_pixel_color( ASVisual *asv, unsigned long pixel, CARD32 *r, CARD32 *g, CARD32 *b )
{
	XColor xcol ;
	xcol.flags = DoRed|DoGreen|DoBlue ;
	xcol.pixel = pixel ;
	if( XQueryColor( asv->dpy, asv->colormap, &xcol ) != 0 )
	{
		*r = xcol.red>>8 ;
		*g = xcol.green>>8 ;
		*b = xcol.blue>>8 ;
	}
}
#endif /*ifndef X_DISPLAY_MISSING */


CARD32 color2pixel32bgr(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	*pixel = ARGB32_RED8(encoded_color)|(ARGB32_GREEN8(encoded_color)<<8)|(ARGB32_BLUE8(encoded_color)<<16);
	return 0;
}
CARD32 color2pixel32rgb(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	*pixel = encoded_color&0x00FFFFFF;
	return 0;
}
CARD32 color2pixel24bgr(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	*pixel = encoded_color&0x00FFFFFF;
	return 0;
}
CARD32 color2pixel24rgb(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	*pixel = encoded_color&0x00FFFFFF;
	return 0;
}
CARD32 color2pixel16bgr(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	register CARD32 c = encoded_color ;
    *pixel = ((c&0x000000F8)<<8)|((c&0x0000FC00)>>5)|((c&0x00F80000)>>19);
	return (c>>1)&0x00030103;
}
CARD32 color2pixel16rgb(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	register CARD32 c = encoded_color ;
    *pixel = ((c&0x00F80000)>>8)|((c&0x0000FC00)>>5)|((c&0x000000F8)>>3);
	return (c>>1)&0x00030103;
}
CARD32 color2pixel15bgr(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	register CARD32 c = encoded_color ;
    *pixel = ((c&0x000000F8)<<7)|((c&0x0000F800)>>6)|((c&0x00F80000)>>19);
	return (c>>1)&0x00030303;
}
CARD32 color2pixel15rgb(ASVisual *asv, CARD32 encoded_color, unsigned long *pixel)
{
	register CARD32 c = encoded_color ;
    *pixel = ((c&0x00F80000)>>9)|((c&0x0000F800)>>6)|((c&0x000000F8)>>3);
	return (c>>1)&0x00030303;
}

CARD32 color2pixel_pseudo3bpp( ASVisual *asv, CARD32 encoded_color, unsigned long *pixel )
{
	register CARD32 c = encoded_color ;
	*pixel = asv->as_colormap[((c>>25)&0x0008)|((c>>16)&0x0002)|((c>>7)&0x0001)];
	return (c>>1)&0x003F3F3F;
}

CARD32 color2pixel_pseudo6bpp( ASVisual *asv, CARD32 encoded_color, unsigned long *pixel )
{
	register CARD32 c = encoded_color ;
	*pixel = asv->as_colormap[((c>>22)&0x0030)|((c>>14)&0x000C)|((c>>6)&0x0003)];
	return (c>>1)&0x001F1F1F;
}

CARD32 color2pixel_pseudo12bpp( ASVisual *asv, CARD32 encoded_color, unsigned long *pixel )
{
	register CARD32 c = encoded_color ;
	*pixel = asv->as_colormap[((c>>16)&0x0F00)|((c>>10)&0x00F0)|((c>>4)&0x000F)];
	return (c>>1)&0x00070707;
}

void pixel2color32rgb(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color32bgr(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color24rgb(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color24bgr(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color16rgb(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color16bgr(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color15rgb(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}
void pixel2color15bgr(ASVisual *asv, unsigned long pixel, CARD32 *red, CARD32 *green, CARD32 *blue)
{}

void ximage2scanline32(ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register CARD32 *a = sl->alpha+sl->offset_x;
	int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x);
	register CARD32 *src = (CARD32*)xim_data ;
/*	src += sl->offset_x; */
/*fprintf( stderr, "%d: ", y);*/

#ifdef WORDS_BIGENDIAN
	if( !asv->msb_first )
#else
	if( asv->msb_first )
#endif
	{
		while (--i >= 0)
		{
			b[i] = (src[i]>>24)&0x0ff;
			g[i] = (src[i]>>16)&0x0ff;
			r[i] = (src[i]>>8)&0x0ff;
			a[i] = src[i]&0x0ff;
/*			fprintf( stderr, "[%d->%8.8X %8.8X %8.8X %8.8X = %8.8X]", i, r[i], g[i], b[i], a[i], src[i]);*/
		}
	}else
	{
		while (--i >= 0)
		{
			a[i] = (src[i]>>24)&0x0ff;
			r[i] = (src[i]>>16)&0x0ff;
			g[i] = (src[i]>>8)&0x0ff;
			b[i] =  src[i]&0x0ff;
		}
	}
/*fprintf( stderr, "\n");*/
}

void ximage2scanline16( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD16 *src = (CARD16*)xim_data ;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
#ifdef WORDS_BIGENDIAN
	if( !asv->msb_first )
#else
	if( asv->msb_first )
#endif
		do
		{
#define ENCODE_MSBF_565(r,gh3,gl3,b)	(((gh3)&0x0007)|((gl3)&0xE000)|((r)&0x00F8)|((b)&0x1F00))
			r[i] =  (src[i]&0x00F8);
			g[i] = ((src[i]&0x0007)<<5)|((src[i]&0xE000)>>11);
			b[i] =  (src[i]&0x1F00)>>5;
		}while( --i >= 0);
	else
		do
		{
#define ENCODE_LSBF_565(r,g,b) (((g)&0x07E0)|((r)&0xF800)|((b)&0x001F))
			r[i] =  (src[i]&0xF800)>>8;
			g[i] =  (src[i]&0x07E0)>>3;
			b[i] =  (src[i]&0x001F)<<3;
		}while( --i >= 0);

}
void ximage2scanline15( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD16 *src = (CARD16*)xim_data ;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
#ifdef WORDS_BIGENDIAN
	if( !asv->msb_first )
#else
	if( asv->msb_first )
#endif
		do
		{
#define ENCODE_MSBF_555(r,gh2,gl3,b)	(((gh2)&0x0003)|((gl3)&0xE000)|((r)&0x007C)|((b)&0x1F00))
			r[i] =  (src[i]&0x007C)<<1;
			g[i] = ((src[i]&0x0003)<<6)|((src[i]&0xE000)>>10);
			b[i] =  (src[i]&0x1F00)>>5;
		}while( --i >= 0);
	else
		do
		{
#define ENCODE_LSBF_555(r,g,b) (((g)&0x03E0)|((r)&0x7C00)|((b)&0x001F))
			r[i] =  (src[i]&0x7C00)>>7;
			g[i] =  (src[i]&0x03E0)>>2;
			b[i] =  (src[i]&0x001F)<<3;
		}while( --i >= 0);
}

#ifndef X_DISPLAY_MISSING

void
ximage2scanline_pseudo3bpp( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;

	do
	{
		unsigned long pixel = XGetPixel( xim, i, y );
		ARGB32 c = asv->as_colormap_reverse.xref[pixel] ;
		if( c == 0 )
			query_pixel_color( asv, pixel, r+i, g+i, b+i );
		else
		{
			r[i] =  ARGB32_RED8(c);
			g[i] =  ARGB32_GREEN8(c);
			b[i] =  ARGB32_BLUE8(c);
		}
	}while( --i >= 0);

}

void
ximage2scanline_pseudo6bpp( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;

	if( xim->bits_per_pixel == 8 )
	{
		register CARD8 *src = (CARD8*)xim_data ;
		do
		{
			ARGB32 c = asv->as_colormap_reverse.xref[src[i]] ;
			if( c == 0 )
				query_pixel_color( asv, src[i], r+i, g+i, b+i );
			else
			{
				r[i] =  ARGB32_RED8(c);
				g[i] =  ARGB32_GREEN8(c);
				b[i] =  ARGB32_BLUE8(c);
			}
		}while( --i >= 0);

	}else
		do
		{
			unsigned long pixel = XGetPixel( xim, i, y );
			ARGB32 c = asv->as_colormap_reverse.xref[pixel] ;
			if( c == 0 )
				query_pixel_color( asv, pixel, r+i, g+i, b+i );
			else
			{
				r[i] =  ARGB32_RED8(c);
				g[i] =  ARGB32_GREEN8(c);
				b[i] =  ARGB32_BLUE8(c);
			}
		}while( --i >= 0);
}

void
ximage2scanline_pseudo12bpp( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;

	if( xim->bits_per_pixel == 16 )
	{
		register CARD16 *src = (CARD16*)xim_data ;
		do
		{
            ASHashData hdata ;
            ARGB32 c ;
            if( get_hash_item( asv->as_colormap_reverse.hash, AS_HASHABLE((unsigned long)src[i]), &hdata.vptr ) != ASH_Success )
				query_pixel_color( asv, src[i], r+i, g+i, b+i );
			else
			{
                c = hdata.c32;
				r[i] =  ARGB32_RED8(c);
				g[i] =  ARGB32_GREEN8(c);
				b[i] =  ARGB32_BLUE8(c);
			}
		}while( --i >= 0);

	}else
		do
		{
			unsigned long pixel = XGetPixel( xim, i, y );
            ASHashData hdata ;
			ARGB32 c ;
            if( get_hash_item( asv->as_colormap_reverse.hash, (ASHashableValue)pixel, &hdata.vptr ) != ASH_Success )
				query_pixel_color( asv, pixel, r+i, g+i, b+i );
			else
			{
                c = hdata.c32;
				r[i] =  ARGB32_RED8(c);
				g[i] =  ARGB32_GREEN8(c);
				b[i] =  ARGB32_BLUE8(c);
			}
		}while( --i >= 0);
}
#endif

void scanline2ximage32( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register CARD32 *a = sl->alpha+sl->offset_x;
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x);
	register CARD32 *src = (CARD32*)xim_data;
/*	src += sl->offset_x ; */
/*fprintf( stderr, "%d: ", y);*/
#ifdef WORDS_BIGENDIAN
	if( !asv->msb_first )
#else
	if( asv->msb_first )
#endif
		while( --i >= 0)
		{ 
			src[i] = (b[i]<<24)|(g[i]<<16)|(r[i]<<8)|a[i];
/*			fprintf( stderr, "[%d->%8.8X %8.8X %8.8X %8.8X = %8.8X]", i, r[i], g[i], b[i], a[i], src[i]);  */
		}
	else
		while( --i >= 0) src[i] = (a[i]<<24)|(r[i]<<16)|(g[i]<<8)|b[i];
/*fprintf( stderr, "\n");*/
#ifdef DEBUG_SL2XIMAGE
	i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x);
	src = (CARD32*)xim_data;
	src += sl->offset_x;
	printf( "%d: xim->width = %d, sl->width = %d, sl->offset = %d: ", y, xim->width, sl->width, sl->offset_x );
	while(--i>=0 )	printf( "%8.8lX ", src[i] );
	printf( "\n" );
#endif
}

void scanline2ximage16( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD16 *src = (CARD16*)xim_data ;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register CARD32 c = (r[i]<<20) | (g[i]<<10) | (b[i]);
#ifdef WORDS_BIGENDIAN
	if( !asv->msb_first )
#else
	if( asv->msb_first )
#endif
		do
		{
			src[i] = ENCODE_MSBF_565((c>>20),(c>>15),(c<<1),(c<<5));
			if( --i < 0 )
				break;
			/* carry over quantization error allow for error diffusion:*/
			c = ((c>>1)&0x00300403)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
/*fprintf( stderr, "c = 0x%X, d = 0x%X, c^d = 0x%X\n", c, d, c^d );*/
			}
		}while(1);
	else
		do
		{
			src[i] = ENCODE_LSBF_565((c>>12),(c>>7),(c>>3));
			if( --i < 0 )
				break;
			/* carry over quantization error allow for error diffusion:*/
			c = ((c>>1)&0x00300403)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
/*fprintf( stderr, "c = 0x%X, d = 0x%X, c^d = 0x%X\n", c, d, c^d );*/
			}
		}while(1);
}

void scanline2ximage15( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD16 *src = (CARD16*)xim_data ;
    register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register CARD32 c = (r[i]<<20) | (g[i]<<10) | (b[i]);
#ifdef WORDS_BIGENDIAN
	if( !asv->msb_first )
#else
	if( asv->msb_first )
#endif
		do
		{
			src[i] = ENCODE_MSBF_555((c>>21),(c>>16),c/*(c>>2)*/,(c<<5));
			if( --i < 0 )
				break;
			/* carry over quantization error allow for error diffusion:*/
			c = ((c>>1)&0x00300C03)+((r[i]<<20) | (g[i]<<10) | (b[i]));
/*fprintf( stderr, "%s:%d src[%d] = 0x%4.4X, c = 0x%X, color[%d] = #%2.2X%2.2X%2.2X\n", __FUNCTION__, __LINE__, i+1, src[i+1], c, i, r[i], g[i], b[i]);*/
			{
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
/*fprintf( stderr, "%s:%d c = 0x%X, d = 0x%X, c^d = 0x%X\n", __FUNCTION__, __LINE__, c, d, c^d );*/
			}
		}while(1);
	else
	{
		do
		{
			src[i] = ENCODE_LSBF_555((c>>13),(c>>8),(c>>3));
			if( --i < 0 )
				break;
			/* carry over quantization error allow for error diffusion:*/
			c = ((c>>1)&0x00300C03)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
			}
		}while(1);
	}
}

#ifndef X_DISPLAY_MISSING
void
scanline2ximage_pseudo3bpp( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD32 c = (r[i]<<20) | (g[i]<<10) | (b[i]);

	do
	{
		XPutPixel( xim, i, y, asv->as_colormap[((c>>25)&0x0008)|((c>>16)&0x0002)|((c>>7)&0x0001)] );
		if( --i < 0 )
			break;
		c = ((c>>1)&0x03F0FC3F)+((r[i]<<20) | (g[i]<<10) | (b[i]));
		{/* handling possible overflow : */
			register CARD32 d = c&0x300C0300 ;
			if( d )
			{
				if( c&0x30000000 )
					d |= 0x0FF00000;
				if( c&0x000C0000 )
					d |= 0x0003FC00 ;
				if( c&0x00000300 )
					d |= 0x000000FF ;
				c ^= d;
			}
		}
	}while(i);
}

void
scanline2ximage_pseudo6bpp( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD32 c = (r[i]<<20) | (g[i]<<10) | (b[i]);

	if( xim->bits_per_pixel == 8 )
	{
		register CARD8 *dst = (CARD8*)xim_data ;
		do
		{
			dst[i] = asv->as_colormap[((c>>22)&0x0030)|((c>>14)&0x000C)|((c>>6)&0x0003)];
			if( --i < 0 )
				break;
			c = ((c>>1)&0x01F07C1F)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{/* handling possible overflow : */
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
			}
		}while(i);
	}else
	{
		do
		{
			XPutPixel( xim, i, y, asv->as_colormap[((c>>22)&0x0030)|((c>>14)&0x000C)|((c>>6)&0x0003)] );
			if( --i < 0 )
				break;
			c = ((c>>1)&0x01F07C1F)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{/* handling possible overflow : */
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
			}
		}while(i);
	}
}

void
scanline2ximage_pseudo12bpp( ASVisual *asv, XImage *xim, ASScanline *sl, int y,  register unsigned char *xim_data )
{
	register CARD32 *r = sl->xc1+sl->offset_x, *g = sl->xc2+sl->offset_x, *b = sl->xc3+sl->offset_x;
	register int i = MIN((unsigned int)(xim->width),sl->width-sl->offset_x)-1;
	register CARD32 c = (r[i]<<20) | (g[i]<<10) | (b[i]);

	if( xim->bits_per_pixel == 16 )
	{
		register CARD16 *dst = (CARD16*)xim_data ;
		do
		{
			dst[i] = asv->as_colormap[((c>>16)&0x0F00)|((c>>10)&0x00F0)|((c>>4)&0x000F)];
			if( --i < 0 )
				break;
			c = ((c>>1)&0x00701C07)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{/* handling possible overflow : */
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
			}
		}while(i);
	}else
	{
		do
		{
			XPutPixel( xim, i, y, asv->as_colormap[((c>>16)&0x0F00)|((c>>10)&0x00F0)|((c>>4)&0x000F)] );
			if( --i < 0 )
				break;
			c = ((c>>1)&0x00701C07)+((r[i]<<20) | (g[i]<<10) | (b[i]));
			{/* handling possible overflow : */
				register CARD32 d = c&0x300C0300 ;
				if( d )
				{
					if( c&0x30000000 )
						d |= 0x0FF00000;
					if( c&0x000C0000 )
						d |= 0x0003FC00 ;
					if( c&0x00000300 )
						d |= 0x000000FF ;
					c ^= d;
				}
			}
		}while(i);
	}
}

#endif /* ifndef X_DISPLAY_MISSING */

