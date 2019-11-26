#ifndef ROOT_RVersion
#define ROOT_RVersion

/* Version information automatically generated by installer. */

/*
 * These macros can be used in the following way:
 *
 *    #if ROOT_VERSION_CODE >= ROOT_VERSION(2,23,4)
 *       #include <newheader.h>
 *    #else
 *       #include <oldheader.h>
 *    #endif
 *
*/

#define ROOT_RELEASE "6.19/01"
#define ROOT_RELEASE_DATE "May 29 2019"
#define ROOT_RELEASE_TIME "10:59:40"
#define ROOT_VERSION(a,b,c) (((a) << 16) + ((b) << 8) + (c))
#define ROOT_VERSION_CODE ROOT_VERSION(6,19,1) /* 398081 */

#endif
